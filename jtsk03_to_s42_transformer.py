#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Coordinate Transformation: JTSK03 (S-JTSK) to S-42/83/03
Based on official Rezortná transformačná služba parameters

This script transforms coordinates from JTSK03 (Slovak Local System - Krovak projection)
to S-42/83/03 (Soviet Unified System, Zone 4)

Source: Rezortná transformačná služba documentation
"""

import math
import numpy as np


def dms2deg(degrees, minutes, seconds):
    """Convert degrees, minutes, seconds to decimal degrees"""
    return degrees + minutes / 60.0 + seconds / 3600.0


class CoordinateTransformer:
    """Transform coordinates between JTSK03 and S-42/83/03 systems"""
    
    # Bessel ellipsoid parameters (used for JTSK03)
    BESSEL_A = 6377397.155
    BESSEL_B = 6356078.96290
    BESSEL_E = 0.081696831
    BESSEL_I = 1 / 299.153
    BESSEL_E2 = 1 - (BESSEL_B / BESSEL_A) ** 2
    
    # Krasovsky ellipsoid parameters (used for S-42/83/03)
    KRASOVSKY_A = 6378245.0
    KRASOVSKY_B = 6356863.01887
    KRASOVSKY_E2 = 0.006693421623
    KRASOVSKY_E = 0.0818192702
    
    # GRS80 ellipsoid parameters (intermediate)
    GRS80_A = 6378137.0
    GRS80_B = 6356752.314245
    GRS80_E2 = 0.00669438
    
    # Krovak projection parameters
    KROVAK_ALPHA = 1.000597498372
    KROVAK_K = 1.003419164
    KROVAK_U_K = math.radians(dms2deg(59, 42, 42.69689))
    KROVAK_V_K = math.radians(dms2deg(42, 31, 31.417251))
    KROVAK_A_PARAM = math.pi / 2 - KROVAK_U_K
    KROVAK_S_0 = math.radians(dms2deg(78, 30, 0))
    KROVAK_N = math.sin(KROVAK_S_0)
    KROVAK_R = 6380703.61054
    KROVAK_RHO_0 = 0.9999 * KROVAK_R * (1 / math.tan(KROVAK_S_0))
    
    # Datum transformation parameters: Bessel -> GRS80
    BESSEL_GRS80_PARAMS = {
        "dX": 485.021,
        "dY": 169.465,
        "dZ": 483.839,
        "Rx": math.radians(dms2deg(0, 0, -7.786342)),
        "Ry": math.radians(dms2deg(0, 0, -4.397554)),
        "Rz": math.radians(dms2deg(0, 0, -4.102655)),
        "m": 0
    }
    
    # Datum transformation parameters: GRS80 -> Krasovsky
    GRS80_KRASOVSKY_PARAMS = {
        "dX": 66.171809,
        "dY": 23.995845,
        "dZ": 83.769412,
        "Rx": math.radians(dms2deg(0, 0, 2.50300846)),
        "Ry": math.radians(dms2deg(0, 0, 2.19216164)),
        "Rz": math.radians(dms2deg(0, 0, -2.58393360)),
        "m": 0
    }
    
    # S-42/83/03 Zone 4 parameters
    S42_LON0_ZONE4 = math.radians(27.0)  # Central meridian for zone 4
    S42_FALSE_EASTING = 4500000.0
    S42_FALSE_NORTHING = 0.0
    S42_SCALE = 0.9996  # Scale factor
    
    def __init__(self):
        """Initialize transformer"""
        pass
    
    def jtsk03_to_geographic(self, x, y):
        """
        Convert JTSK03 Krovak projection to geographic (lat/lon) coordinates
        Using proper Krovak inverse transformation
        
        Args:
            x: X coordinate in JTSK03 system (meters)
            y: Y coordinate in JTSK03 system (meters)
            
        Returns:
            tuple: (latitude_rad, longitude_rad) in radians on Bessel ellipsoid
        """
        # Krovak inverse projection
        # Calculate rho and theta from cartesian coordinates
        rho = math.sqrt(x**2 + y**2)
        theta = math.atan2(x, y)
        
        # Calculate D (longitude in the oblique cone system)
        D = theta / self.KROVAK_N
        
        # Calculate psi (conformal latitude in the oblique cone system)
        psi = 2 * (math.atan(math.exp(
            (1.0 / self.KROVAK_N) * math.log(self.KROVAK_RHO_0 / rho)
        )) - math.pi / 4.0)
        
        # Conformal to geodetic latitude transformation
        # First approximation using sine rule
        sin_U = math.sin(self.KROVAK_A_PARAM) * math.sin(psi)
        
        # Clamp to valid range for asin
        sin_U = max(-1.0, min(1.0, sin_U))
        U = math.asin(sin_U)
        
        # Convert conformal latitude to geodetic latitude
        # Using iterative method with Bessel parameters
        e2 = self.BESSEL_E2
        e = math.sqrt(e2)
        
        # Iterative calculation
        lat_rad = U
        for _ in range(10):
            sin_lat = math.sin(lat_rad)
            lat_rad = 2 * math.atan(
                ((1 + e * sin_lat) / (1 - e * sin_lat)) ** (e / 2) * 
                math.tan(math.pi / 4 + U / 2)
            ) - math.pi / 2
        
        # Longitude
        lon_rad = D + self.KROVAK_V_K
        
        return lat_rad, lon_rad
    
    def geographic_to_cartesian(self, lat_rad, lon_rad, ellipsoid='bessel'):
        """
        Convert geographic coordinates to cartesian (X, Y, Z)
        
        Args:
            lat_rad: Latitude in radians
            lon_rad: Longitude in radians
            ellipsoid: 'bessel', 'grs80', or 'krasovsky'
            
        Returns:
            tuple: (X, Y, Z) cartesian coordinates
        """
        if ellipsoid == 'bessel':
            a = self.BESSEL_A
            e2 = self.BESSEL_E2
        elif ellipsoid == 'grs80':
            a = self.GRS80_A
            e2 = self.GRS80_E2
        else:  # krasovsky
            a = self.KRASOVSKY_A
            e2 = self.KRASOVSKY_E2
        
        sin_lat = math.sin(lat_rad)
        cos_lat = math.cos(lat_rad)
        sin_lon = math.sin(lon_rad)
        cos_lon = math.cos(lon_rad)
        
        N = a / math.sqrt(1 - e2 * sin_lat**2)
        
        X = N * cos_lat * cos_lon
        Y = N * cos_lat * sin_lon
        Z = N * (1 - e2) * sin_lat
        
        return X, Y, Z
    
    def cartesian_to_geographic(self, X, Y, Z, ellipsoid='bessel'):
        """
        Convert cartesian (X, Y, Z) to geographic coordinates
        
        Args:
            X, Y, Z: Cartesian coordinates
            ellipsoid: 'bessel', 'grs80', or 'krasovsky'
            
        Returns:
            tuple: (latitude_rad, longitude_rad)
        """
        if ellipsoid == 'bessel':
            a = self.BESSEL_A
            b = self.BESSEL_B
            e2 = self.BESSEL_E2
        elif ellipsoid == 'grs80':
            a = self.GRS80_A
            b = self.GRS80_B
            e2 = self.GRS80_E2
        else:  # krasovsky
            a = self.KRASOVSKY_A
            b = self.KRASOVSKY_B
            e2 = self.KRASOVSKY_E2
        
        lon_rad = math.atan2(Y, X)
        
        p = math.sqrt(X**2 + Y**2)
        e_prime2 = e2 / (1 - e2)
        
        # Iterative method using Bowring's formula
        lat_rad = math.atan2(Z, p * (1 - e2))
        
        for _ in range(10):
            sin_lat = math.sin(lat_rad)
            N = a / math.sqrt(1 - e2 * sin_lat**2)
            lat_rad = math.atan2(Z + e2 * N * sin_lat, p)
        
        return lat_rad, lon_rad
    
    def datum_shift(self, X, Y, Z, params, inverse=False):
        """
        Apply 7-parameter datum transformation (Bursa-Wolf model)
        
        Args:
            X, Y, Z: Input cartesian coordinates
            params: Transformation parameters dictionary
            inverse: If True, apply inverse transformation
            
        Returns:
            tuple: (X_out, Y_out, Z_out) transformed coordinates
        """
        dX = params['dX']
        dY = params['dY']
        dZ = params['dZ']
        Rx = params['Rx']
        Ry = params['Ry']
        Rz = params['Rz']
        m = params['m']
        
        if inverse:
            dX, dY, dZ = -dX, -dY, -dZ
            Rx, Ry, Rz = -Rx, -Ry, -Rz
            m = -m
        
        # Small angles approximation (valid for small rotations in radians)
        # X' = (1+m)(X - Rz*Y + Ry*Z) + dX
        # Y' = (1+m)(Rz*X + Y - Rx*Z) + dY
        # Z' = (1+m)(-Ry*X + Rx*Y + Z) + dZ
        
        scale = 1 + m
        
        X_out = scale * (X - Rz * Y + Ry * Z) + dX
        Y_out = scale * (Rz * X + Y - Rx * Z) + dY
        Z_out = scale * (-Ry * X + Rx * Y + Z) + dZ
        
        return X_out, Y_out, Z_out
    
    def geographic_to_s42_zone4(self, lat_rad, lon_rad):
        """
        Convert geographic coordinates to S-42/83/03 Zone 4 (Transverse Mercator)
        
        Args:
            lat_rad: Latitude in radians (on Krasovsky ellipsoid)
            lon_rad: Longitude in radians (on Krasovsky ellipsoid)
            
        Returns:
            tuple: (easting, northing) in S-42/83/03 Zone 4
        """
        a = self.KRASOVSKY_A
        e2 = self.KRASOVSKY_E2
        e_prime2 = e2 / (1 - e2)
        
        sin_lat = math.sin(lat_rad)
        cos_lat = math.cos(lat_rad)
        tan_lat = math.tan(lat_rad)
        
        N = a / math.sqrt(1 - e2 * sin_lat**2)
        
        # Longitude difference from central meridian
        lon_diff = lon_rad - self.S42_LON0_ZONE4
        
        # Ensure lon_diff is within valid range
        while lon_diff > math.pi:
            lon_diff -= 2 * math.pi
        while lon_diff < -math.pi:
            lon_diff += 2 * math.pi
        
        # Meridional arc
        A0 = 1 - e2 / 4 - 3 * e2**2 / 64 - 5 * e2**3 / 256
        A2 = 3 / 8 * (e2 + e2**2 / 4 - 15 * e2**3 / 128)
        A4 = 15 / 256 * (e2**2 - 3 * e2**3 / 4)
        A6 = -35 / 3072 * e2**3
        
        M = a * (A0 * lat_rad - A2 * math.sin(2 * lat_rad) + 
                 A4 * math.sin(4 * lat_rad) - A6 * math.sin(6 * lat_rad))
        
        T = tan_lat**2
        C = e_prime2 * cos_lat**2
        A = cos_lat * lon_diff
        
        # Transverse Mercator formulas
        easting = (self.S42_SCALE * N * (A + 
                                         A**3 / 6 * (1 - T + C) +
                                         A**5 / 120 * (5 - 18*T + T**2 + 72*C - 58*e_prime2)) +
                   self.S42_FALSE_EASTING)
        
        northing = (self.S42_SCALE * (M + N * tan_lat * 
                                      (A**2 / 2 + 
                                       A**4 / 24 * (5 - T + 9*C + 4*C**2) +
                                       A**6 / 720 * (61 - 58*T + T**2 + 600*C - 330*e_prime2))) +
                    self.S42_FALSE_NORTHING)
        
        return easting, northing
    
    def transform(self, x_jtsk03, y_jtsk03):
        """
        Main transformation method: JTSK03 -> Geographic (Bessel) -> 
        Cartesian (Bessel) -> Cartesian (GRS80) -> Cartesian (Krasovsky) ->
        Geographic (Krasovsky) -> S-42/83/03 Zone 4
        
        Args:
            x_jtsk03: X coordinate in JTSK03 system (meters)
            y_jtsk03: Y coordinate in JTSK03 system (meters)
            
        Returns:
            dict: Complete transformation results
        """
        # Step 1: JTSK03 Krovak -> Geographic (Bessel ellipsoid, radians)
        lat_bessel, lon_bessel = self.jtsk03_to_geographic(x_jtsk03, y_jtsk03)
        
        # Step 2: Geographic (Bessel) -> Cartesian (Bessel)
        X_bessel, Y_bessel, Z_bessel = self.geographic_to_cartesian(
            lat_bessel, lon_bessel, ellipsoid='bessel'
        )
        
        # Step 3: Datum shift Bessel -> GRS80
        X_grs80, Y_grs80, Z_grs80 = self.datum_shift(
            X_bessel, Y_bessel, Z_bessel, 
            self.BESSEL_GRS80_PARAMS
        )
        
        # Step 4: Datum shift GRS80 -> Krasovsky
        X_krasovsky, Y_krasovsky, Z_krasovsky = self.datum_shift(
            X_grs80, Y_grs80, Z_grs80,
            self.GRS80_KRASOVSKY_PARAMS
        )
        
        # Step 5: Cartesian (Krasovsky) -> Geographic (Krasovsky)
        lat_krasovsky, lon_krasovsky = self.cartesian_to_geographic(
            X_krasovsky, Y_krasovsky, Z_krasovsky, ellipsoid='krasovsky'
        )
        
        # Step 6: Geographic (Krasovsky) -> S-42/83/03 Zone 4
        x_s42, y_s42 = self.geographic_to_s42_zone4(lat_krasovsky, lon_krasovsky)
        
        return {
            'x_jtsk03': x_jtsk03,
            'y_jtsk03': y_jtsk03,
            'latitude_bessel': math.degrees(lat_bessel),
            'longitude_bessel': math.degrees(lon_bessel),
            'latitude_krasovsky': math.degrees(lat_krasovsky),
            'longitude_krasovsky': math.degrees(lon_krasovsky),
            'x_s42': x_s42,
            'y_s42': y_s42
        }


def main():
    """Main function - example usage"""
    transformer = CoordinateTransformer()
    
    print("=" * 90)
    print("JTSK03 to S-42/83/03 Zone 4 Coordinate Transformation")
    print("Based on Rezortná transformačná služba official parameters")
    print("=" * 90)
    print()
    
    # Test with the example from Rezortná transformačná služba
    print("Test Case from Rezortná transformačná služba:")
    print("-" * 90)
    
    x_test = 1204350.50
    y_test = 496927.41
    
    print(f"Input (JTSK03):     X={x_test:.2f} m, Y={y_test:.2f} m")
    
    result = transformer.transform(x_test, y_test)
    
    print(f"Geographic (Bessel):   Lat={result['latitude_bessel']:.6f}°, Lon={result['longitude_bessel']:.6f}°")
    print(f"Geographic (Krasovsky): Lat={result['latitude_krasovsky']:.6f}°, Lon={result['longitude_krasovsky']:.6f}°")
    print(f"Output (S-42/83/03 Zone 4): X={result['x_s42']:.3f} m, Y={result['y_s42']:.3f} m")
    print()
    print("Expected from Rezortná transformačná služba:")
    print(f"  Geographic: Lat=48.493798°, Lon=35.619555°")
    print(f"  Output S-42/83/03 Zone 4: X=5422203.220 m, Y=4283434.780 m")
    print()
    print()
    
    # Interactive input
    print("=" * 90)
    print("Interactive Transformation")
    print("=" * 90)
    print("Enter JTSK03 coordinates to transform to S-42/83/03 Zone 4")
    print("(or press Ctrl+C to quit)")
    print()
    
    try:
        while True:
            try:
                x_input = input("Enter X coordinate in JTSK03 (meters): ").strip()
                if x_input.lower() == 'q':
                    break
                    
                x_value = float(x_input)
                y_value = float(input("Enter Y coordinate in JTSK03 (meters): "))
                
                result = transformer.transform(x_value, y_value)
                
                print("\nTransformation Result:")
                print("-" * 90)
                print(f"  Input JTSK03:              X={result['x_jtsk03']:.2f} m, Y={result['y_jtsk03']:.2f} m")
                print(f"  Geographic (Bessel):       Lat={result['latitude_bessel']:.6f}°, Lon={result['longitude_bessel']:.6f}°")
                print(f"  Geographic (Krasovsky):    Lat={result['latitude_krasovsky']:.6f}°, Lon={result['longitude_krasovsky']:.6f}°")
                print(f"  Output S-42/83/03 Zone 4:  X={result['x_s42']:.3f} m, Y={result['y_s42']:.3f} m")
                print()
            except ValueError as e:
                print(f"Invalid input. Please enter valid coordinates: {e}")
                print()
    except (KeyboardInterrupt, EOFError):
        print("\n\nTransformation completed. Goodbye!")


if __name__ == "__main__":
    main()
