// CelestialMath.swift — altitude prediction + turning one Moon photo into an
// altitude sight, using the IMU gravity vector for the horizontal and the Moon's
// KNOWN angular size for the plate scale (no camera intrinsics needed).
//
// Mirrors imu_fusion/astro.predicted_altitude and the synthetic-horizon idea in
// iphone_model.py.  © 2026 MIT.

import Foundation
import simd

public enum CelestialMath {

    static func rad(_ d: Double) -> Double { d * .pi / 180 }
    static func deg(_ r: Double) -> Double { r * 180 / .pi }

    /// Predicted geometric altitude (deg) of a body at geographic position
    /// (dec, gha) seen from (lat, lon).  East longitude positive; H = GHA + lon.
    public static func predictedAltitude(latDeg: Double, lonDeg: Double,
                                         decDeg: Double, ghaDeg: Double) -> Double {
        let H = rad(ghaDeg + lonDeg)
        let lat = rad(latDeg), dec = rad(decDeg)
        let sinAlt = sin(lat) * sin(dec) + cos(lat) * cos(dec) * cos(H)
        return deg(asin(max(-1, min(1, sinAlt))))
    }

    /// Predicted azimuth (deg, from North through East) — used for the analytic
    /// LOP gradient in PositionFix.
    public static func predictedAzimuth(latDeg: Double, lonDeg: Double,
                                        decDeg: Double, ghaDeg: Double) -> Double {
        let H = rad(ghaDeg + lonDeg)
        let lat = rad(latDeg), dec = rad(decDeg)
        let y = -cos(dec) * sin(H)
        let x = cos(lat) * sin(dec) - sin(lat) * cos(dec) * cos(H)
        var az = deg(atan2(y, x))
        if az < 0 { az += 360 }
        return az
    }

    /// One altitude SIGHT from a single frame.
    ///
    /// The measured altitude of the Moon's centre above the gravity horizontal =
    /// (boresight elevation from the IMU) + (the Moon centre's offset from the
    /// image principal point, projected onto the image "up" = the direction gravity
    /// projects onto the sensor, scaled by the plate scale).
    ///
    /// - gravityDevice: `CMDeviceMotion.gravity` as a unit vector in the DEVICE
    ///   frame (points toward the ground). CoreMotion device axes: +x right,
    ///   +y up, +z out of the screen; the rear camera looks along −z.
    /// - moonCentrePx / principalPx: pixel coordinates (origin top-left, +y DOWN).
    /// - arcsecPerPixel: from DiskMetrology (Moon known size / fitted radius).
    public static func altitudeSightDeg(gravityDevice g0: SIMD3<Double>,
                                        moonCentrePx: SIMD2<Double>,
                                        principalPx: SIMD2<Double>,
                                        arcsecPerPixel: Double) -> Double {
        let g = simd_normalize(g0)
        let up = -g                                   // local up in device frame
        let boresight = SIMD3<Double>(0, 0, -1)       // rear camera, device frame
        // Elevation of the boresight above the horizontal plane.
        let boreElevDeg = deg(asin(max(-1, min(1, simd_dot(boresight, up)))))

        // Image "up": project device up onto the sensor plane (device x,y).  The
        // sensor y axis points DOWN in pixel space, so flip y.
        var upImg = SIMD2<Double>(up.x, -up.y)
        let n = simd_length(upImg)
        if n > 1e-9 { upImg /= n } else { upImg = SIMD2<Double>(0, -1) }

        // Moon offset from principal point, in pixels, then its vertical component.
        let off = SIMD2<Double>(moonCentrePx.x - principalPx.x,
                                moonCentrePx.y - principalPx.y)
        let vertPx = simd_dot(off, upImg)             // +up in image = higher altitude
        let vertDeg = vertPx * arcsecPerPixel / 3600.0
        return boreElevDeg + vertDeg
    }

    /// Rough 1σ of the altitude sight (deg): the synthetic-horizon tilt floor
    /// (attitude) plus the sub-pixel centroid term.  The tilt floor is the number
    /// the on-device calibration protocol measures; 0.1° is the model's assumption.
    public static func altitudeSigmaDeg(tiltFloorDeg: Double = 0.1,
                                        centroidPx: Double = 0.3,
                                        arcsecPerPixel: Double) -> Double {
        let centroidDeg = centroidPx * arcsecPerPixel / 3600.0
        return (tiltFloorDeg * tiltFloorDeg + centroidDeg * centroidDeg).squareRoot()
    }

    // MARK: Apparent-altitude corrections (refraction + topocentric parallax)
    //
    // Port of imu_fusion/corrections.py — the SINGLE SOURCE OF TRUTH, validated
    // against IAU ERFA / astropy (sub-arcminute) and checked here against the
    // exported golden vectors.  A phone measures the APPARENT altitude of the
    // disk centre (light refracted up; observer topocentric, so a near body is
    // seen lower by parallax).  `predictedAltitude` and PositionFix work in
    // GEOMETRIC GEOCENTRIC altitudes, so a sight is reduced before the fix.
    //
    //   measured apparent --(- refraction)--> topocentric --(+ parallax)--> geocentric
    //
    // Uses the study's Earth radius so the golden vectors reproduce exactly.

    /// Earth radius (km) — matches starfix.EARTH_RADIUS (circumference / 2π).
    static let earthRadiusKm = 6374.574260537331

    /// Equatorial horizontal parallax (deg): sin(HP) = R_earth / distance.
    public static func horizontalParallaxDeg(distanceKm: Double) -> Double {
        deg(asin(earthRadiusKm / distanceKm))
    }

    /// Parallax in altitude (deg) at a geometric altitude: asin(sin HP · cos alt).
    public static func parallaxInAltitudeDeg(geometricAltDeg: Double,
                                             distanceKm: Double) -> Double {
        let hp = rad(horizontalParallaxDeg(distanceKm: distanceKm))
        let c = max(-1.0, min(1.0, sin(hp) * cos(rad(geometricAltDeg))))
        return deg(asin(c))
    }

    /// Atmospheric refraction (deg), Bennett's formula on the APPARENT altitude.
    /// Mirrors starfix.get_refraction (arc-minutes). Defaults: 10°C, 101 kPa, 50%.
    public static func refractionDeg(apparentAltDeg h: Double,
                                     temperatureC: Double = 10.0,
                                     pressureKpa: Double = 101.0,
                                     humidityPct: Double = 50.0) -> Double {
        if h <= -1.0 { return 0.0 }
        let d = (h + 7.31 / (h + 4.4)) * .pi / 180.0
        let humidityFactor = 1.0 - 0.0000008 * humidityPct * temperatureC
        let arcmin = (1.0 / tan(d)) * (pressureKpa / 101.0)
                     * (283.0 / (273.0 + temperatureC)) * humidityFactor
        return arcmin / 60.0
    }

    /// Reduce a MEASURED apparent altitude to the geometric geocentric altitude
    /// that `predictedAltitude` / PositionFix expects. Applied to every sight.
    public static func geometricFromApparent(apparentAltDeg: Double, distanceKm: Double,
                                             temperatureC: Double = 10.0,
                                             pressureKpa: Double = 101.0,
                                             humidityPct: Double = 50.0) -> Double {
        let topo = apparentAltDeg - refractionDeg(apparentAltDeg: apparentAltDeg,
                                                  temperatureC: temperatureC,
                                                  pressureKpa: pressureKpa,
                                                  humidityPct: humidityPct)
        return topo + parallaxInAltitudeDeg(geometricAltDeg: topo, distanceKm: distanceKm)
    }

    /// Forward model: predict the apparent altitude from the geometric geocentric
    /// altitude. Exact inverse of `geometricFromApparent` (used by the golden test).
    public static func apparentFromGeometric(geocentricAltDeg: Double, distanceKm: Double,
                                             temperatureC: Double = 10.0,
                                             pressureKpa: Double = 101.0,
                                             humidityPct: Double = 50.0) -> Double {
        var topo = geocentricAltDeg
            - parallaxInAltitudeDeg(geometricAltDeg: geocentricAltDeg, distanceKm: distanceKm)
        for _ in 0..<3 {
            topo = geocentricAltDeg
                - parallaxInAltitudeDeg(geometricAltDeg: topo, distanceKm: distanceKm)
        }
        var app = topo
        for _ in 0..<4 {
            app = topo + refractionDeg(apparentAltDeg: app, temperatureC: temperatureC,
                                       pressureKpa: pressureKpa, humidityPct: humidityPct)
        }
        return app
    }
}
