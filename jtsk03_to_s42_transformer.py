#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Coordinate Transformation: JTSK03 (S-JTSK) to S-42/83/03
Based on mathematical cartography transformation methods

This script transforms coordinates from JTSK03 (Slovak Local System)
to S-42/83/03 (Soviet Unified System)
"""

import math


class CoordinateTransformer:
    """Transform coordinates between JTSK03 and S-42/83/03 systems"""
    
    # JTSK03 System Parameters (Slovak Local Coordinate System)
    JTSK03_A = 6378245.0  # Semi-major axis (Krasovsky ellipsoid)
    JTSK03_E2 = 0.006693421622965943  # First eccentricity squared
    
    # Transformation Parameters (JTSK03 to S-42/83/03)
    # These are typical values - adjust based on your specific reference document
    DX = 485.0  # Translation in X direction (meters)
    DY = -169.0  # Translation in Y direction (meters)
    DZ = -483.0  # Translation in Z direction (meters)
    
    # Rotation angles (in radians)
    RX = 7.4e-6  # Rotation around X axis
    RY = -2.1e-6  # Rotation around Y axis
    RZ = 1.37e-5  # Rotation around Z axis
    
    # Scale factor
    SCALE = 1.0
    
    # Reference Meridian for JTSK03
    JTSK03_LON0 = 24.833333333  # degrees
    
    # Oblique Mercator parameters for JTSK03
    JTSK03_LAT0 = 49.5  # degrees
    JTSK03_FALSE_EASTING = 0.0
    JTSK03_FALSE_NORTHING = 0.0
    
    def __init__(self):
        """Initialize transformer with constants"""
        self.a = self.JTSK03_A
        self.e2 = self.JTSK03_E2
        
    def jtsk03_to_geographic(self, x, y):
        """
        Convert JTSK03 local coordinates to geographic (lat/lon)
        
        Args:
            x: X coordinate in JTSK03 system (meters)
            y: Y coordinate in JTSK03 system (meters)
            
        Returns:
            tuple: (latitude, longitude) in decimal degrees
        """
        # Remove false easting/northing
        x_adj = x - self.JTSK03_FALSE_EASTING
        y_adj = y - self.JTSK03_FALSE_NORTHING
        
        # Calculate footpoint latitude
        e = math.sqrt(self.e2)
        e_prime2 = self.e2 / (1 - self.e2)
        
        # Recurrence relations for footpoint latitude
        m = y_adj / self.a
        mu = m / (1 - self.e2 / 4 - 3 * self.e2**2 / 64 - 5 * self.e2**3 / 256)
        
        phi1_rad = (mu + 
                    (3 * self.e2 / 8 + 3 * self.e2**2 / 256 + 45 * self.e2**3 / 16384) * math.sin(2 * mu) +
                    (15 * self.e2**2 / 256 + 45 * self.e2**3 / 1024) * math.sin(4 * mu) +
                    (35 * self.e2**3 / 3072) * math.sin(6 * mu))
        
        cos_phi1 = math.cos(phi1_rad)
        sin_phi1 = math.sin(phi1_rad)
        tan_phi1 = math.tan(phi1_rad)
        n1 = self.a / math.sqrt(1 - self.e2 * sin_phi1**2)
        t1 = tan_phi1**2
        c1 = e_prime2 * cos_phi1**2
        r1 = self.a * (1 - self.e2) / math.sqrt((1 - self.e2 * sin_phi1**2)**3)
        d = x_adj / n1
        
        # Calculate latitude
        lat_rad = (phi1_rad - 
                   (t1 / r1) * (d**2 / 2 - 
                               (d**4 / 24) * (5 + 3*t1 + 10*c1 - 4*c1**2 - 9*e_prime2) +
                               (d**6 / 720) * (61 + 90*t1 + 28*t1**2 + 45*c1 - 252*e_prime2 - 3*c1**2)))
        
        # Calculate longitude
        lon_rad = (math.radians(self.JTSK03_LON0) + 
                   (d - (d**3 / 6) * (1 + 2*t1 + c1) + 
                    (d**5 / 120) * (1 - 2*t1 + c1 + 2*c1**2 - 3*e_prime2)) / cos_phi1)
        
        lat = math.degrees(lat_rad)
        lon = math.degrees(lon_rad)
        
        return lat, lon
    
    def geographic_to_s42(self, lat, lon):
        """
        Convert geographic coordinates to S-42/83/03 system
        
        Args:
            lat: Latitude in decimal degrees
            lon: Longitude in decimal degrees
            
        Returns:
            tuple: (x_s42, y_s42) in S-42/83/03 system
        """
        # S-42 parameters (Soviet Unified System)
        a_s42 = 6378245.0  # Semi-major axis (Krasovsky ellipsoid)
        e2_s42 = 0.006693421622965943
        
        lat_rad = math.radians(lat)
        lon_rad = math.radians(lon)
        
        # For S-42/83/03, use zone 3 with central meridian at 21°
        lon0_s42 = math.radians(21.0)
        
        # Calculate Transverse Mercator projection
        e_s42 = math.sqrt(e2_s42)
        sin_lat = math.sin(lat_rad)
        cos_lat = math.cos(lat_rad)
        tan_lat = math.tan(lat_rad)
        
        n_s42 = a_s42 / math.sqrt(1 - e2_s42 * sin_lat**2)
        t = tan_lat**2
        c = e_s42**2 * cos_lat**2 / (1 - e2_s42)
        a_coeff = cos_lat * (lon_rad - lon0_s42)
        
        # Meridional arc
        n_order = 5
        n = (a_s42 - a_s42 * (1 - e2_s42)) / (a_s42 + a_s42 * (1 - e2_s42))
        n2 = n * n
        n3 = n2 * n
        n4 = n3 * n
        n5 = n4 * n
        
        alpha1 = (3/2 * n - 27/32 * n3 + 269/512 * n5)
        alpha2 = (21/16 * n2 - 55/32 * n4)
        alpha3 = (151/96 * n3 - 417/128 * n5)
        alpha4 = (1097/512 * n4)
        
        m = (a_s42 * (1 - e2_s42) * 
             (1 - alpha1 * math.cos(2 * lat_rad) + 
              alpha2 * math.cos(4 * lat_rad) - 
              alpha3 * math.cos(6 * lat_rad) + 
              alpha4 * math.cos(8 * lat_rad)))
        
        # Easting
        x_s42 = (n_s42 * (a_coeff + 
                         (a_coeff**3 / 6) * (1 - t + c) +
                         (a_coeff**5 / 120) * (5 - 18*t + t**2 + 72*c - 58*e_s42**2)) + 
                500000.0)
        
        # Northing
        y_s42 = (m + n_s42 * tan_lat * 
                (a_coeff**2 / 2 + 
                 (a_coeff**4 / 24) * (5 - t + 9*c + 4*c**2) +
                 (a_coeff**6 / 720) * (61 - 58*t + t**2 + 600*c - 330*e_s42**2)))
        
        return x_s42, y_s42
    
    def transform(self, x_jtsk03, y_jtsk03):
        """
        Main transformation method: JTSK03 -> Geographic -> S-42/83/03
        
        Args:
            x_jtsk03: X coordinate in JTSK03 system (meters)
            y_jtsk03: Y coordinate in JTSK03 system (meters)
            
        Returns:
            dict: {
                'x_jtsk03': original x,
                'y_jtsk03': original y,
                'latitude': geographic latitude,
                'longitude': geographic longitude,
                'x_s42': transformed x in S-42/83/03,
                'y_s42': transformed y in S-42/83/03
            }
        """
        # Step 1: JTSK03 -> Geographic
        lat, lon = self.jtsk03_to_geographic(x_jtsk03, y_jtsk03)
        
        # Step 2: Geographic -> S-42/83/03
        x_s42, y_s42 = self.geographic_to_s42(lat, lon)
        
        return {
            'x_jtsk03': x_jtsk03,
            'y_jtsk03': y_jtsk03,
            'latitude': lat,
            'longitude': lon,
            'x_s42': x_s42,
            'y_s42': y_s42
        }


def main():
    """Main function - example usage"""
    transformer = CoordinateTransformer()
    
    print("=" * 70)
    print("JTSK03 to S-42/83/03 Coordinate Transformation")
    print("=" * 70)
    print()
    
    # Example coordinates (you can modify these)
    examples = [
        (400000, 250000),
        (500000, 300000),
        (600000, 350000),
    ]
    
    for x_jtsk, y_jtsk in examples:
        print(f"Input (JTSK03): X={x_jtsk:.2f} m, Y={y_jtsk:.2f} m")
        
        result = transformer.transform(x_jtsk, y_jtsk)
        
        print(f"  → Geographic: Lat={result['latitude']:.6f}°, Lon={result['longitude']:.6f}°")
        print(f"  → Output (S-42/83/03): X={result['x_s42']:.2f} m, Y={result['y_s42']:.2f} m")
        print()
    
    # Interactive input
    print("=" * 70)
    print("Interactive Transformation")
    print("=" * 70)
    
    try:
        while True:
            try:
                x_input = float(input("Enter X coordinate in JTSK03 (or 'q' to quit): "))
                y_input = float(input("Enter Y coordinate in JTSK03: "))
                
                result = transformer.transform(x_input, y_input)
                
                print("\nResult:")
                print(f"  Input JTSK03:      X={result['x_jtsk03']:.2f} m, Y={result['y_jtsk03']:.2f} m")
                print(f"  Geographic:        Lat={result['latitude']:.6f}°, Lon={result['longitude']:.6f}°")
                print(f"  Output S-42/83/03: X={result['x_s42']:.2f} m, Y={result['y_s42']:.2f} m")
                print()
            except ValueError as e:
                print(f"Invalid input: {e}")
                print()
    except (KeyboardInterrupt, EOFError):
        print("\nTransformation completed. Goodbye!")


if __name__ == "__main__":
    main()
