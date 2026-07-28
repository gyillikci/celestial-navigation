// FixViewModel.swift — glue: on a shutter, run limb → altitude sight → position fix.
// © 2026 MIT.

import Foundation
import CoreLocation
import simd

@MainActor
public final class FixViewModel: ObservableObject {

    public struct Result {
        public let fix: Fix
        public let measAltDeg: Double
        public let sigmaAltDeg: Double
        public let radiusPx: Double
        public let arcsecPerPixel: Double
        public let rmsePx: Double
        public let time: Date
        public let moon: GeographicPosition
    }

    @Published public var status: String = "Point the tele camera at the Moon, then Capture."
    @Published public var result: Result?

    private let capture: CaptureController
    private let location = CLLocationManager()

    // Seed position (dead-reckoned / last-known); the horizontal disambiguator.
    public var priorLat = 51.5
    public var priorLon = 0.0
    public var priorSigmaKm = 30.0

    public init(capture: CaptureController) {
        self.capture = capture
        location.requestWhenInUseAuthorization()
        location.desiredAccuracy = kCLLocationAccuracyKilometer
        location.startUpdatingLocation()
    }

    public func onCapture() {
        // Refresh the prior from GPS if we have it (still just a seed — the sight
        // is what earns the accuracy).
        if let loc = location.location {
            priorLat = loc.coordinate.latitude
            priorLon = loc.coordinate.longitude
        }
        let now = Date()
        capture.captureLuminance { [weak self] img, gravity in
            Task { @MainActor in self?.process(img, gravity, now) }
        }
    }

    private func process(_ img: GrayImage?, _ gravity: SIMD3<Double>?, _ now: Date) {
        guard let img else { status = "No camera frame yet — try again."; return }
        guard let gravity else { status = "Waiting for motion sensors…"; return }
        guard let fit = DiskMetrology.subpixelLimb(img) else {
            status = "Couldn't fit the Moon's limb. Fill more of the frame, steady the shot."
            return
        }
        let moon = MoonEphemeris.position(now)
        let angRadius = MoonEphemeris.angularRadiusArcsec(distanceKm: moon.distanceKm)
        let arcsecPerPx = DiskMetrology.arcsecPerPixel(radiusPx: fit.radiusPx,
                                                       moonAngularRadiusArcsec: angRadius)
        let principal = SIMD2<Double>(Double(img.width) / 2, Double(img.height) / 2)
        let measAlt = CelestialMath.altitudeSightDeg(
            gravityDevice: gravity,
            moonCentrePx: SIMD2<Double>(fit.cx, fit.cy),
            principalPx: principal, arcsecPerPixel: arcsecPerPx)
        let sigmaAlt = CelestialMath.altitudeSigmaDeg(arcsecPerPixel: arcsecPerPx)

        let fix = PositionFix.solve(
            measAltDeg: measAlt, decDeg: moon.declinationDeg, ghaDeg: moon.ghaDeg,
            sigmaAltDeg: sigmaAlt,
            priorLatDeg: priorLat, priorLonDeg: priorLon, priorSigmaKm: priorSigmaKm)

        result = Result(fix: fix, measAltDeg: measAlt, sigmaAltDeg: sigmaAlt,
                        radiusPx: fit.radiusPx, arcsecPerPixel: arcsecPerPx,
                        rmsePx: fit.rmsePx, time: now, moon: moon)
        status = String(format: "Fix: %.4f°, %.4f°  (±%.1f km, 1σ)",
                        fix.latDeg, fix.lonDeg, fix.oneSigmaMajorKm)
    }
}
