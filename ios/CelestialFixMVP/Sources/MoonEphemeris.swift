// MoonEphemeris.swift — compact, dependency-free low-precision Moon position.
//
// Returns the Moon's apparent geographic position (declination + Greenwich hour
// angle) at a UTC instant.  This is a TRUNCATED lunar theory (the largest periodic
// terms only), good to a few arc-minutes — enough to demonstrate the full
// capture → fix pipeline.  The production app replaces this with the validated
// ephemeris ported from `imu_fusion/astro.py` / `starfix` (arc-second accurate,
// cross-checked in `validate_ephemeris.py`).
//
// References: Meeus, *Astronomical Algorithms* (2nd ed.), ch. 47 (truncated) + ch.
// 13/22 for the obliquity and nutation-free apparent place.  © 2026 MIT.

import Foundation

public struct GeographicPosition {
    public let declinationDeg: Double     // Dec, degrees
    public let ghaDeg: Double             // Greenwich hour angle, degrees [0,360)
    public let distanceKm: Double         // Earth–Moon distance (for parallax/size)
}

public enum MoonEphemeris {

    private static func rad(_ d: Double) -> Double { d * .pi / 180.0 }
    private static func deg(_ r: Double) -> Double { r * 180.0 / .pi }
    private static func norm360(_ d: Double) -> Double { let x = d.truncatingRemainder(dividingBy: 360); return x < 0 ? x + 360 : x }

    /// Julian centuries (TT) from J2000.0.  UTC≈TT here (ΔT ignored at this precision).
    private static func julianCenturies(_ date: Date) -> Double {
        // Unix epoch 1970-01-01 = JD 2440587.5
        let jd = date.timeIntervalSince1970 / 86400.0 + 2440587.5
        return (jd - 2451545.0) / 36525.0
    }

    /// Greenwich Mean Sidereal Time (degrees) — Meeus 12.4.
    public static func gmstDeg(_ date: Date) -> Double {
        let jd = date.timeIntervalSince1970 / 86400.0 + 2440587.5
        let t = (jd - 2451545.0) / 36525.0
        let g = 280.46061837 + 360.98564736629 * (jd - 2451545.0)
              + 0.000387933 * t * t - t * t * t / 38710000.0
        return norm360(g)
    }

    /// Apparent Moon geographic position at `date` (UTC).
    public static func position(_ date: Date) -> GeographicPosition {
        let t = julianCenturies(date)

        // Fundamental arguments (degrees) — Meeus 47.
        let Lp = norm360(218.3164477 + 481267.88123421 * t - 0.0015786 * t * t) // mean longitude
        let D  = norm360(297.8501921 + 445267.1114034 * t - 0.0018819 * t * t)  // mean elongation
        let M  = norm360(357.5291092 + 35999.0502909 * t)                       // Sun anomaly
        let Mp = norm360(134.9633964 + 477198.8675055 * t + 0.0087414 * t * t)  // Moon anomaly
        let F  = norm360(93.2720950 + 483202.0175233 * t - 0.0036539 * t * t)   // argument of latitude

        let Dr = rad(D), Mr = rad(M), Mpr = rad(Mp), Fr = rad(F)

        // Longitude (Σl, 1e-6 deg) — largest terms of Meeus table 47.A.
        var sl = 0.0
        sl += 6288774 * sin(Mpr)
        sl += 1274027 * sin(2*Dr - Mpr)
        sl +=  658314 * sin(2*Dr)
        sl +=  213618 * sin(2*Mpr)
        sl -=  185116 * sin(Mr)
        sl -=  114332 * sin(2*Fr)
        sl +=   58793 * sin(2*Dr - 2*Mpr)
        sl +=   57066 * sin(2*Dr - Mr - Mpr)
        sl +=   53322 * sin(2*Dr + Mpr)
        sl +=   45758 * sin(2*Dr - Mr)
        sl -=   40923 * sin(Mr - Mpr)
        sl -=   34720 * sin(Dr)
        sl -=   30383 * sin(Mr + Mpr)
        sl +=   15327 * sin(2*Dr - 2*Fr)
        sl -=   12528 * sin(Mpr + 2*Fr)
        sl +=   10980 * sin(Mpr - 2*Fr)

        // Latitude (Σb, 1e-6 deg) — largest terms of table 47.B.
        var sb = 0.0
        sb += 5128122 * sin(Fr)
        sb +=  280602 * sin(Mpr + Fr)
        sb +=  277693 * sin(Mpr - Fr)
        sb +=  173237 * sin(2*Dr - Fr)
        sb +=   55413 * sin(2*Dr - Mpr + Fr)
        sb +=   46271 * sin(2*Dr - Mpr - Fr)
        sb +=   32573 * sin(2*Dr + Fr)
        sb +=   17198 * sin(2*Mpr + Fr)
        sb +=    9266 * sin(2*Dr + Mpr - Fr)

        // Distance (Σr, 1e-3 km) — largest terms of table 47.A cosines.
        var sr = 0.0
        sr -= 20905355 * cos(Mpr)
        sr -=  3699111 * cos(2*Dr - Mpr)
        sr -=  2955968 * cos(2*Dr)
        sr -=   569925 * cos(2*Mpr)
        sr +=    48888 * cos(Mr)
        sr -=     3149 * cos(2*Fr)
        sr +=   246158 * cos(2*Dr - 2*Mpr)
        sr -=   152138 * cos(2*Dr - Mr - Mpr)
        sr -=   170733 * cos(2*Dr + Mpr)
        sr -=   204586 * cos(2*Dr - Mr)

        let lambda = Lp + sl / 1_000_000.0                 // ecliptic longitude, deg
        let beta   = sb / 1_000_000.0                       // ecliptic latitude, deg
        let dist   = 385000.56 + sr / 1000.0               // km

        // Mean obliquity (Meeus 22.2), nutation ignored at this precision.
        let eps = 23.439291 - 0.0130042 * t
        let lr = rad(lambda), br = rad(beta), er = rad(eps)

        // Ecliptic → equatorial (Meeus 13.3/13.4).
        let ra = atan2(sin(lr) * cos(er) - tan(br) * sin(er), cos(lr))
        let dec = asin(sin(br) * cos(er) + cos(br) * sin(er) * sin(lr))

        let raDeg = norm360(deg(ra))
        let gha = norm360(gmstDeg(date) - raDeg)           // apparent GHA (geocentric)
        return GeographicPosition(declinationDeg: deg(dec), ghaDeg: gha, distanceKm: dist)
    }

    /// Moon's geocentric angular radius (arc-seconds) at the given distance.
    public static func angularRadiusArcsec(distanceKm: Double) -> Double {
        let moonRadiusKm = 1737.4
        return atan2(moonRadiusKm, distanceKm) * 180.0 / .pi * 3600.0
    }
}
