// DiskMetrology.swift — sub-pixel limb fit → the TRUE (full-circle) disk centre and
// radius, ported from imu_fusion/disk_metrology.py (gradient-peak sub-pixel edge on
// radial rays + least-squares circle).  CPU version for the MVP; the production app
// moves the ray sweep to a Metal kernel.
//
// Fits the FULL circle from the bright-limb arc, so the centre is correct even for a
// partly-lit Moon — NOT the centroid of the lit blob (which is pulled toward the
// bright limb, ~9′ ≈ 17 km of error at first quarter).  © 2026 MIT.

import Foundation
import simd

public struct DiskFit {
    public let cx: Double
    public let cy: Double
    public let radiusPx: Double
    public let rmsePx: Double
    public let nPoints: Int
}

/// A single-channel luminance image.
public struct GrayImage {
    public let width: Int
    public let height: Int
    public let pixels: [Float]            // row-major, length width*height
    public init(width: Int, height: Int, pixels: [Float]) {
        self.width = width; self.height = height; self.pixels = pixels
    }
    @inline(__always) func at(_ x: Int, _ y: Int) -> Float {
        pixels[y * width + x]
    }
    /// Bilinear sample (clamped).
    func sample(_ x: Double, _ y: Double) -> Double {
        let x0 = min(max(Int(x.rounded(.down)), 0), width - 2)
        let y0 = min(max(Int(y.rounded(.down)), 0), height - 2)
        let fx = x - Double(x0), fy = y - Double(y0)
        let a = Double(at(x0, y0)),   b = Double(at(x0 + 1, y0))
        let c = Double(at(x0, y0 + 1)), d = Double(at(x0 + 1, y0 + 1))
        return a*(1-fx)*(1-fy) + b*fx*(1-fy) + c*(1-fx)*fy + d*fx*fy
    }
}

public enum DiskMetrology {

    private static func circleFit(_ pts: [SIMD2<Double>]) -> (cx: Double, cy: Double, r: Double)? {
        // Kasa algebraic fit: minimise |x²+y² - (2a x + 2b y + c)|.
        var Sxx = 0.0, Sxy = 0.0, Syy = 0.0, Sx = 0.0, Sy = 0.0
        var Sxz = 0.0, Syz = 0.0, Sz = 0.0
        let n = Double(pts.count)
        for p in pts {
            let z = p.x*p.x + p.y*p.y
            Sxx += p.x*p.x; Sxy += p.x*p.y; Syy += p.y*p.y
            Sx += p.x; Sy += p.y; Sxz += p.x*z; Syz += p.y*z; Sz += z
        }
        // Solve 3×3 [ [Sxx Sxy Sx],[Sxy Syy Sy],[Sx Sy n] ] [2a,2b,c]^T = [Sxz,Syz,Sz]
        // by Cramer's rule.
        func det3(_ a:Double,_ b:Double,_ c:Double,_ d:Double,_ e:Double,_ f:Double,_ g:Double,_ h:Double,_ i:Double)->Double{
            a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g)
        }
        let D  = det3(Sxx,Sxy,Sx, Sxy,Syy,Sy, Sx,Sy,n)
        let Da = det3(Sxz,Sxy,Sx, Syz,Syy,Sy, Sz,Sy,n)
        let Db = det3(Sxx,Sxz,Sx, Sxy,Syz,Sy, Sx,Sz,n)
        let Dc = det3(Sxx,Sxy,Sxz, Sxy,Syy,Syz, Sx,Sy,Sz)
        if abs(D) < 1e-9 { return nil }
        let a = (Da/D)/2, b = (Db/D)/2, c = Dc/D
        let r = (c + a*a + b*b).squareRoot()
        return (a, b, r)
    }

    /// Rough centre/radius from a brightness threshold (night Moon = high contrast).
    /// The production path uses the gradient sky-limb RANSAC seed (disk_metrology.py)
    /// for dim / partial phases.
    private static func seed(_ img: GrayImage) -> (cx: Double, cy: Double, r: Double) {
        var sorted = img.pixels; sorted.sort()
        let sky = Double(sorted[sorted.count / 5])                 // 20th pct
        let hi  = Double(sorted[sorted.count * 99 / 100])          // 99th pct
        let thr = sky + 0.4 * (hi - sky)
        var sx = 0.0, sy = 0.0, n = 0.0
        for y in 0..<img.height { for x in 0..<img.width {
            if Double(img.at(x, y)) > thr { sx += Double(x); sy += Double(y); n += 1 }
        }}
        if n < 1 { return (Double(img.width)/2, Double(img.height)/2, Double(img.width)/6) }
        return (sx/n, sy/n, (n / .pi).squareRoot())
    }

    /// Sub-pixel limb fit.  Returns the full-circle centre & radius.
    public static func subpixelLimb(_ img: GrayImage, nRays: Int = 720,
                                    halfWindowPx: Double = 12) -> DiskFit? {
        var sorted = img.pixels; sorted.sort()
        let sky = Double(sorted[sorted.count / 5])
        let bright = Double(sorted[sorted.count * 99 / 100])
        var (cx, cy, r) = seed(img)
        var pts: [SIMD2<Double>] = []

        for _ in 0..<3 {
            pts.removeAll(keepingCapacity: true)
            for k in 0..<nRays {
                let th = 2 * .pi * Double(k) / Double(nRays)
                let cth = cos(th), sth = sin(th)
                // sample the profile across the expected limb radius
                var prof: [Double] = []
                var rr = -halfWindowPx
                while rr < halfWindowPx {
                    prof.append(img.sample(cx + (r+rr)*cth, cy + (r+rr)*sth))
                    rr += 0.5
                }
                let m = prof.count
                let inMean = prof[0..<max(1, m/3)].reduce(0,+) / Double(max(1, m/3))
                let outMean = prof[(m - max(1, m/3))...].reduce(0,+) / Double(max(1, m/3))
                if inMean < bright*0.5 || outMean > sky + 0.25*(bright - sky) { continue }
                // Falling-edge gradient g[i] = prof[i-1] − prof[i+1] (positive at the
                // bright→dark limb); take its peak and parabola-refine to sub-pixel.
                var g = [Double](repeating: 0, count: m)
                for i in 1..<(m-1) { g[i] = prof[i-1] - prof[i+1] }
                var kBest = 1; var gBest = -Double.greatestFiniteMagnitude
                for i in 2..<(m-2) { if g[i] > gBest { gBest = g[i]; kBest = i } }
                if kBest <= 1 || kBest >= m-2 { continue }
                let denom = g[kBest-1] - 2*g[kBest] + g[kBest+1]
                let dsub = abs(denom) > 1e-9 ? 0.5 * (g[kBest-1] - g[kBest+1]) / denom : 0.0
                let rEdge = r - halfWindowPx + (Double(kBest) + dsub) * 0.5
                pts.append(SIMD2<Double>(cx + rEdge*cth, cy + rEdge*sth))
            }
            if pts.count < 12 { break }
            if let fit = circleFit(pts) { cx = fit.cx; cy = fit.cy; r = fit.r }
        }
        if pts.count < 12 { return nil }
        var se = 0.0
        for p in pts { let d = hypot(p.x - cx, p.y - cy) - r; se += d*d }
        return DiskFit(cx: cx, cy: cy, radiusPx: r,
                       rmsePx: (se / Double(pts.count)).squareRoot(), nPoints: pts.count)
    }

    /// Plate scale from the fitted radius and the Moon's known angular size.
    public static func arcsecPerPixel(radiusPx: Double, moonAngularRadiusArcsec: Double) -> Double {
        moonAngularRadiusArcsec / radiusPx
    }
}
