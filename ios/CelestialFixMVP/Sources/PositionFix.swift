// PositionFix.swift — a dependency-free 2-DOF Gauss-Newton position fix from ONE
// altitude line of position plus a position prior (the dead-reckoned / last-known
// seed).  One altitude sight is a line; the prior pins the along-line direction, so
// the pair is a well-posed 2×2 solve with a real covariance.
//
// This is the MVP stand-in for the GTSAM factor graph: the exact same intercept-
// method gradient (d alt/d position = (sinAz, cosAz)/R) that celestial_factor_graph
// uses analytically.  Scaling to many shots + IMU is where GTSAM replaces this.
// © 2026 MIT.

import Foundation

public struct Fix {
    public let latDeg: Double
    public let lonDeg: Double
    /// Covariance in the local ENU tangent (metres²): [[σ_ee, σ_en],[σ_en, σ_nn]].
    public let covEN: (ee: Double, en: Double, nn: Double)
    public var oneSigmaMajorKm: Double {
        let a = covEN.ee, b = covEN.en, c = covEN.nn
        let tr = a + c, det = a * c - b * b
        let l = tr / 2 + ((tr * tr) / 4 - det).squareRoot()
        return l.squareRoot() / 1000.0
    }
    public var iterations: Int
}

public enum PositionFix {

    private static let R = 6_371_000.0                   // Earth radius, m
    private static func rad(_ d: Double) -> Double { d * .pi / 180 }

    /// Solve for (lat, lon) from one altitude sight + a position prior.
    /// - measAltDeg: measured altitude of the Moon centre.
    /// - dec/gha: Moon geographic position at the shot time.
    /// - sigmaAltDeg: 1σ of the sight.
    /// - prior: seed lat/lon and its 1σ (km); the horizontal disambiguator.
    public static func solve(measAltDeg: Double, decDeg: Double, ghaDeg: Double,
                             sigmaAltDeg: Double,
                             priorLatDeg: Double, priorLonDeg: Double,
                             priorSigmaKm: Double,
                             maxIter: Int = 12) -> Fix {
        var lat = priorLatDeg, lon = priorLonDeg
        let mPerDegLat = R * .pi / 180.0
        let wAlt = 1.0 / rad(sigmaAltDeg)                 // weight = 1/σ (rad)
        let wPos = 1.0 / (priorSigmaKm * 1000.0)          // weight = 1/σ (m)

        var iter = 0
        // Info matrix accumulators (in ENU metres), reused for covariance.
        var Iee = 0.0, Ien = 0.0, Inn = 0.0
        for _ in 0..<maxIter {
            iter += 1
            let mPerDegLon = mPerDegLat * cos(rad(lat))
            // Altitude residual + gradient (intercept method).
            let altPred = CelestialMath.predictedAltitude(latDeg: lat, lonDeg: lon,
                                                          decDeg: decDeg, ghaDeg: ghaDeg)
            let az = rad(CelestialMath.predictedAzimuth(latDeg: lat, lonDeg: lon,
                                                        decDeg: decDeg, ghaDeg: ghaDeg))
            // d(alt)/d(east,north) = (sin Az, cos Az)/R  [rad per metre]
            let gE = sin(az) / R, gN = cos(az) / R
            let rAlt = rad(altPred - measAltDeg)          // residual, rad

            // Prior residual + gradient in ENU metres.
            let dLatM = (lat - priorLatDeg) * mPerDegLat
            let dLonM = (lon - priorLonDeg) * mPerDegLon

            // Normal equations J^T W J  x = -J^T W r, with x = [dE, dN] (metres).
            Iee = wAlt*wAlt*gE*gE + wPos*wPos
            Inn = wAlt*wAlt*gN*gN + wPos*wPos
            Ien = wAlt*wAlt*gE*gN
            var bE = -(wAlt*wAlt*gE*rAlt + wPos*wPos*dLonM)
            var bN = -(wAlt*wAlt*gN*rAlt + wPos*wPos*dLatM)

            let det = Iee*Inn - Ien*Ien
            if abs(det) < 1e-18 { break }
            let dE = (Inn*bE - Ien*bN) / det
            let dN = (-Ien*bE + Iee*bN) / det

            lon += (dE / mPerDegLon)
            lat += (dN / mPerDegLat)
            _ = (bE, bN)
            if hypot(dE, dN) < 0.5 { break }              // <0.5 m step → converged
        }
        // Covariance = inverse of the final information matrix (ENU metres²).
        let det = Iee*Inn - Ien*Ien
        let cee = Inn/det, cnn = Iee/det, cen = -Ien/det
        return Fix(latDeg: lat, lonDeg: lon, covEN: (cee, cen, cnn), iterations: iter)
    }
}
