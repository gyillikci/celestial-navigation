// PipelineTests.swift — device-independent checks of the MVP math.  Run in Xcode
// (⌘U).  These mirror the Python reference behaviour in imu_fusion.  © 2026 MIT.

import XCTest
import simd
@testable import CelestialFixMVP

final class PipelineTests: XCTestCase {

    func testMoonEphemerisSane() {
        // 2026-03-24 12:00 UTC (the study's canonical epoch).
        var c = DateComponents(); c.year = 2026; c.month = 3; c.day = 24
        c.hour = 12; c.timeZone = TimeZone(identifier: "UTC")
        let d = Calendar(identifier: .gregorian).date(from: c)!
        let m = MoonEphemeris.position(d)
        XCTAssertLessThan(abs(m.declinationDeg), 29.0)          // |Dec| ≤ ~28.6°
        XCTAssertTrue(m.ghaDeg >= 0 && m.ghaDeg < 360)
        XCTAssertTrue(m.distanceKm > 350_000 && m.distanceKm < 410_000)
        let rad = MoonEphemeris.angularRadiusArcsec(distanceKm: m.distanceKm)
        XCTAssertTrue(rad > 880 && rad < 1010)                 // ~15–17′
    }

    func testAltitudeAtSubpointIs90() {
        // A body is overhead at (lat=dec, lon=-gha).
        let dec = 20.0, gha = 45.0
        let alt = CelestialMath.predictedAltitude(latDeg: dec, lonDeg: -gha,
                                                  decDeg: dec, ghaDeg: gha)
        XCTAssertEqual(alt, 90.0, accuracy: 1e-6)
    }

    func testPositionFixSatisfiesTheSightAndMovesOffThePrior() {
        let dec = 15.0, gha = 30.0
        let trueLat = 40.0, trueLon = -5.0
        let measAlt = CelestialMath.predictedAltitude(latDeg: trueLat, lonDeg: trueLon,
                                                      decDeg: dec, ghaDeg: gha)
        // Prior displaced ~0.4° from truth.
        let fix = PositionFix.solve(measAltDeg: measAlt, decDeg: dec, ghaDeg: gha,
                                    sigmaAltDeg: 0.03,
                                    priorLatDeg: trueLat + 0.4, priorLonDeg: trueLon - 0.3,
                                    priorSigmaKm: 50)
        // The altitude sight is satisfied at the fix.
        let altAtFix = CelestialMath.predictedAltitude(latDeg: fix.latDeg, lonDeg: fix.lonDeg,
                                                       decDeg: dec, ghaDeg: gha)
        XCTAssertEqual(altAtFix, measAlt, accuracy: 0.01)
        // The fix moved off the prior (the sight pulled it onto the LOP).
        XCTAssertGreaterThan(abs(fix.latDeg - (trueLat + 0.4))
                             + abs(fix.lonDeg - (trueLon - 0.3)), 0.01)
        XCTAssertGreaterThan(fix.covEN.ee, 0)
        XCTAssertGreaterThan(fix.covEN.nn, 0)
        XCTAssertLessThan(fix.oneSigmaMajorKm, 60)              // bounded by the prior
    }

    // MARK: sub-pixel limb on a synthetic disk

    private func syntheticDisk(w: Int, h: Int, cx: Double, cy: Double, r: Double,
                               leftShadow: Bool) -> GrayImage {
        var px = [Float](repeating: 8, count: w * h)           // dark sky
        for y in 0..<h { for x in 0..<w {
            let d = hypot(Double(x) - cx, Double(y) - cy)
            // erf-like soft edge over ~2 px
            let t = max(-1.0, min(1.0, (d - r) / 2.0))
            var v = 8.0 + 200.0 * (0.5 * (1 - t))
            if leftShadow && Double(x) < cx { v = 8.0 + (v - 8.0) * 0.05 }  // unlit half
            px[y * w + x] = Float(v)
        }}
        return GrayImage(width: w, height: h, pixels: px)
    }

    func testSubpixelLimbRecoversFullDisk() {
        let img = syntheticDisk(w: 700, h: 560, cx: 355, cy: 300, r: 160, leftShadow: false)
        let fit = DiskMetrology.subpixelLimb(img)
        XCTAssertNotNil(fit)
        guard let fit else { return }
        XCTAssertEqual(fit.cx, 355, accuracy: 0.6)
        XCTAssertEqual(fit.cy, 300, accuracy: 0.6)
        XCTAssertEqual(fit.radiusPx, 160, accuracy: 0.8)
        XCTAssertLessThan(fit.rmsePx, 0.6)
    }

    func testFullCircleCentreNotLitCentroidOnHalfPhase() {
        // Half the disk in shadow: the full-circle centre must still be the true
        // centre, NOT the (biased) lit-blob centroid.
        let img = syntheticDisk(w: 700, h: 560, cx: 355, cy: 300, r: 160, leftShadow: true)
        let fit = DiskMetrology.subpixelLimb(img)
        XCTAssertNotNil(fit)
        guard let fit else { return }
        XCTAssertEqual(fit.cx, 355, accuracy: 1.5)             // true centre
        // Lit-blob centroid would sit well to the right of centre.
        var sx = 0.0, n = 0.0
        for y in 0..<img.height { for x in 0..<img.width {
            if img.at(x, y) > 100 { sx += Double(x); n += 1 } } }
        let litCx = sx / n
        XCTAssertGreaterThan(litCx - 355, 30)                  // centroid is biased
        XCTAssertLessThan(abs(fit.cx - 355), abs(litCx - 355)) // fit is far better
    }
}
