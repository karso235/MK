#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Coordinate Transformation: JTSK03 (S-JTSK) to S-42/83/03
Based on mathematical cartography transformation methods from Rezortná transformačná služba

This script transforms coordinates from JTSK03 (Slovak Local System)
to S-42/83/03 (Soviet Unified System, Zone 4)
"""

import math


class CoordinateTransformer:
    """Transform coordinates between JTSK03 and S-42/83/03 systems"""
    
    # Krasovsky ellipsoid parameters (used for both systems)
    A = 6378245.0  # Semi-major axis (meters)
    B = 6356863.019  # Semi-minor axis (meters)
    E2 = 0.006693421622965943  # First eccentricity squared
    E = 0.0818192702210484  # First eccentricity
    
    # JTSK03 (Oblique Conformal Conic - Krovak projection) parameters
    JTSK03_LAT0 = math.radians(49.5)  # Reference latitude
    JTSK03_LON0 = math.radians(24.833333333)  # Reference longitude
    JTSK03_LAT_CONE = math.radians(49.5)  # Cone latitude
    JTSK03_FALSE_EASTING = 0.0
    JTSK03_FALSE_NORTHING = 0.0
    
    # S-42/83/03 Transverse Mercator parameters (Zone 4)
    S42_LON0_ZONE4 = math.radians(27.0)  # Central meridian for zone 4
    S42_FALSE_EASTING = 4500000.0
    S42_FALSE_NORTHING = 0.0
    S42_SCALE = 0.9996  # Scale factor
    
    def __init__(self):
        """Initialize transformer with constants"""
        pass
    
    def jtsk03_to_geographic(self, x, y):
        """
        Convert JTSK03 Krovak projection to geographic (lat/lon) coordinates
        
        Args:
            x: X coordinate in JTSK03 system (meters)
            y: Y coordinate in JTSK03 system (meters)
            
        Returns:
            tuple: (latitude, longitude) in radians
        """
        # JTSK03 uses Krovak projection (oblique conformal conic)
        # Simplified inverse transformation
        
        # These are empirical coefficients for Krovak inverse
        # Based on standard Slovak transformation
        x_adj = x / 1000000.0
        y_adj = y / 1000000.0
        
        # Inverse Krovak transformation coefficients
        # Using polynomial approximation
        lon = 24.8333333333333 - (y_adj * 0.00003365)
        lat = 49.5 - (x_adj * 0.00003365)
        
        # More precise inverse using iterative method
        # Convert to radians for precise calculation
        lat_rad, lon_rad = self._krovak_inverse(x, y)
        
        return lat_rad, lon_rad
    
    def _krovak_inverse(self, x, y):
        """
        Precise inverse Krovak projection
        
        Args:
            x: X coordinate (meters)
            y: Y coordinate (meters)
            
        Returns:
            tuple: (latitude_rad, longitude_rad)
        """
        # Krovak projection inverse
        # Using standard Slovak algorithm
        
        # Constants for Krovak
        alpha = 1.00281635
        k = 0.9999
        
        # Remove false coordinates
        x_adj = x - self.JTSK03_FALSE_EASTING
        y_adj = y - self.JTSK03_FALSE_NORTHING
        
        # Convert to conformal coordinates
        p = math.sqrt(x_adj**2 + y_adj**2)
        theta = math.atan2(y_adj, x_adj)
        
        # Inverse conformal latitude
        lat_c = 2 * math.atan(math.exp(math.log(p / (self.A * k * alpha)) / alpha)) - math.pi / 2
        
        # Longitude from theta
        lon_c = self.JTSK03_LON0 + (theta / alpha)
        
        # Convert conformal coordinates to geodetic
        # Using iterative method
        lat_rad = lat_c
        for _ in range(10):
            sin_lat = math.sin(lat_rad)
            w = math.sqrt(1 - self.E2 * sin_lat**2)
            lat_rad = 2 * math.atan((1 + self.E2 * sin_lat / w) * 
                                    math.tan(lat_c / 2 + math.pi / 4) - 
                                    self.E2 * math.tan(lat_c / 2 + math.pi / 4)) - math.pi / 2
        
        return lat_rad, lon_c
    
    def geographic_to_s42_zone4(self, lat_rad, lon_rad):
        """
        Convert geographic coordinates to S-42/83/03 Zone 4 (Transverse Mercator)
        
        Args:
            lat_rad: Latitude in radians
            lon_rad: Longitude in radians
            
        Returns:
            tuple: (easting, northing) in S-42/83/03 Zone 4
        """
        # S-42/83/03 uses Transverse Mercator projection for Zone 4
        # Central meridian at 27° E
        
        # Calculate parameters
        sin_lat = math.sin(lat_rad)
        cos_lat = math.cos(lat_rad)
        tan_lat = math.tan(lat_rad)
        
        n = self.A / math.sqrt(1 - self.E2 * sin_lat**2)
        e_prime2 = self.E2 / (1 - self.E2)
        
        # Calculate longitude difference from central meridian
        lon_diff = lon_rad - self.S42_LON0_ZONE4
        
        # Ensure lon_diff is within -pi to pi
        while lon_diff > math.pi:
            lon_diff -= 2 * math.pi
        while lon_diff < -math.pi:
            lon_diff += 2 * math.pi
        
        # Calculate meridional arc
        m = (self.A * ((1 - self.E2 / 4 - 3 * self.E2**2 / 64 - 5 * self.E2**3 / 256) * lat_rad -
                       (3 * self.E2 / 8 + 3 * self.E2**2 / 32 - 45 * self.E2**3 / 1024) * math.sin(2 * lat_rad) +
                       (15 * self.E2**2 / 256 - 45 * self.E2**3 / 1024) * math.sin(4 * lat_rad) -
                       (35 * self.E2**3 / 3072) * math.sin(6 * lat_rad)))
        
        # Calculate easting and northing
        t = tan_lat**2
        c = e_prime2 * cos_lat**2
        a = cos_lat * lon_diff
        
        # Easting
        easting = (self.S42_SCALE * n * (a + 
                                         (a**3 / 6) * (1 - t + c) +
                                         (a**5 / 120) * (5 - 18*t + t**2 + 72*c - 58*e_prime2)) +
                   self.S42_FALSE_EASTING)
        
        # Northing
        northing = (self.S42_SCALE * (m + n * tan_lat * 
                                      (a**2 / 2 + 
                                       (a**4 / 24) * (5 - t + 9*c + 4*c**2) +
                                       (a**6 / 720) * (61 - 58*t + t**2 + 600*c - 330*e_prime2))) +
                    self.S42_FALSE_NORTHING)
        
        return easting, northing
    
    def transform(self, x_jtsk03, y_jtsk03):
        """
        Main transformation method: JTSK03 -> Geographic -> S-42/83/03 Zone 4
        
        Args:
            x_jtsk03: X coordinate in JTSK03 system (meters)
            y_jtsk03: Y coordinate in JTSK03 system (meters)
            
        Returns:
            dict: {
                'x_jtsk03': original x,
                'y_jtsk03': original y,
                'latitude': geographic latitude (decimal degrees),
                'longitude': geographic longitude (decimal degrees),
                'x_s42': transformed x in S-42/83/03 Zone 4,
                'y_s42': transformed y in S-42/83/03 Zone 4
            }
        """
        # Step 1: JTSK03 -> Geographic (radians)
        lat_rad, lon_rad = self.jtsk03_to_geographic(x_jtsk03, y_jtsk03)
        
        # Step 2: Geographic -> S-42/83/03 Zone 4
        x_s42, y_s42 = self.geographic_to_s42_zone4(lat_rad, lon_rad)
        
        return {
            'x_jtsk03': x_jtsk03,
            'y_jtsk03': y_jtsk03,
            'latitude': math.degrees(lat_rad),
            'longitude': math.degrees(lon_rad),
            'x_s42': x_s42,
            'y_s42': y_s42
        }


def main():
    """Main function - example usage"""
    transformer = CoordinateTransformer()
    
    print("=" * 80)
    print("JTSK03 to S-42/83/03 Zone 4 Coordinate Transformation")
    print("Based on Rezortná transformačná služba")
    print("=" * 80)
    print()
    
    # Test with the provided example from Rezortná transformačná služba
    print("Test Case from Rezortná transformačná služba:")
    print("-" * 80)
    
    x_test = 1204350.50
    y_test = 496927.41
    
    print(f"Input (JTSK03): X={x_test:.2f} m, Y={y_test:.2f} m")
    
    result = transformer.transform(x_test, y_test)
    
    print(f"  → Geographic: Lat={result['latitude']:.6f}°, Lon={result['longitude']:.6f}°")
    print(f"  → Output (S-42/83/03, Zone 4): X={result['x_s42']:.2f} m, Y={result['y_s42']:.2f} m")
    print()
    print("Expected results from Rezortná transformačná služba:")
    print(f"  → Geographic: Lat=48.493798°, Lon=35.619555°")
    print(f"  → Output (S-42/83/03, Zone 4): X=5422203.220 m, Y=4283434.780 m")
    print()
    print()
    
    # Interactive input
    print("=" * 80)
    print("Interactive Transformation")
    print("=" * 80)
    print("Enter JTSK03 coordinates to transform to S-42/83/03")
    print("(or 'q' to quit)")
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
                print("-" * 80)
                print(f"  Input JTSK03:           X={result['x_jtsk03']:.2f} m, Y={result['y_jtsk03']:.2f} m")
                print(f"  Geographic:             Lat={result['latitude']:.6f}°, Lon={result['longitude']:.6f}°")
                print(f"  Output S-42/83/03 Zone4: X={result['x_s42']:.2f} m, Y={result['y_s42']:.2f} m")
                print()
            except ValueError as e:
                print(f"Invalid input. Please enter valid coordinates: {e}")
                print()
    except (KeyboardInterrupt, EOFError):
        print("\n\nTransformation completed. Goodbye!")


if __name__ == "__main__":
    main()
