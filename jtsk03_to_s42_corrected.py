#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Coordinate Transformation: JTSK03 (S-JTSK) to S-42/83/03 Zone 4
Corrected implementation based on Krovak inverse transformation formulas
Pages 25-30 from 10_prednaska_MK_vzajomne_transformacie

This script transforms coordinates from JTSK03 (Slovak Local System - Krovak projection)
to S-42/83/03 (Soviet Unified System, Zone 4)

Author: Corrected from lecture notes
"""

import math
import numpy as np
from scipy.integrate import quad


def dms2deg(degrees, minutes, seconds):
    """Convert degrees, minutes, seconds to decimal degrees"""
    return degrees + minutes / 60.0 + seconds / 3600.0


def print_dms(radians, precision=6):
    """Convert radians to DMS string format"""
    degrees = math.degrees(radians)
    d = int(degrees)
    m = int((abs(degrees) - abs(d)) * 60)
    s = ((abs(degrees) - abs(d)) * 60 - m) * 60
    
    if precision == 6:
        return f"{d}°{m}'{s:.6f}''"
    else:
        return f"{d}°{m}'{s:.2f}''"


def meridian(ellipsoid, phi):
    """Calculate meridional arc length"""
    e2 = ellipsoid["e"]**2 if "e" in ellipsoid else ellipsoid["e2"]
    
    A0 = 1 - e2/4 - 3*e2**2/64 - 5*e2**3/256
    A2 = (3/8) * (e2 + e2**2/4 - 15*e2**3/128)
    A4 = (15/256) * (e2**2 - 3*e2**3/4)
    A6 = (-35/3072) * e2**3
    
    return ellipsoid["a"] * (A0*phi - A2*math.sin(2*phi) + A4*math.sin(4*phi) - A6*math.sin(6*phi))


def lateral(ellipsoid, phi):
    """Calculate lateral (radius of curvature in prime vertical)"""
    e2 = ellipsoid["e"]**2 if "e" in ellipsoid else ellipsoid["e2"]
    return ellipsoid["a"] / math.sqrt(1 - e2 * math.sin(phi)**2)


class CoordinateTransformer:
    """Transform coordinates between JTSK03 and S-42/83/03 systems"""
    
    def __init__(self):
        """Initialize ellipsoid parameters"""
        # Bessel ellipsoid (for JTSK03)
        self.Bessel = {
            "a": 6377397.155,
            "b": 6356078.96290,
            "e": 0.081696831,
            "i": 1/299.153
        }
        self.Bessel["e2"] = np.sqrt(((self.Bessel["a"]**2) - (self.Bessel["b"]**2)) / (self.Bessel["b"]**2))
        
        # Krasovsky ellipsoid (for S-42/83/03)
        self.Krasovsky = {
            "a": 6378245.000,
            "b": 6356863.01877,
            "e": 0.081813333,
            "i": 1/298.300
        }
        self.Krasovsky["e2"] = np.sqrt(((self.Krasovsky["a"]**2) - (self.Krasovsky["b"]**2)) / (self.Krasovsky["b"]**2))
        
        # GRS80 ellipsoid (intermediate)
        self.GRS80 = {
            "a": 6378137.000,
            "b": 6356752.31425,
            "e": 0.081819191,
            "i": 1/298.257223563
        }
        self.GRS80["e2"] = np.sqrt(((self.GRS80["a"]**2) - (self.GRS80["b"]**2)) / (self.GRS80["b"]**2))
        
        # Krovak projection parameters
        self.krovak_alpha = 1.000597498372
        self.krovak_k = 1.003419164
        self.krovak_U_k = np.radians(dms2deg(59, 42, 42.69689))
        self.krovak_V_k = np.radians(dms2deg(42, 31, 31.417251))
        self.krovak_a = np.pi / 2 - self.krovak_U_k
        self.krovak_S_0 = np.radians(dms2deg(78, 30, 0))
        self.krovak_n = np.sin(self.krovak_S_0)
        self.krovak_R = 6380703.61054
        self.krovak_rho_0 = 0.9999 * self.krovak_R * (1 / np.tan(self.krovak_S_0))
        
        # Bursa-Wolf transformation parameters
        self.bessel_grs80_params = {
            "dX": 485.021,
            "dY": 169.465,
            "dZ": 483.839,
            "Rx": dms2deg(0, 0, -7.786342),
            "Ry": dms2deg(0, 0, -4.397554),
            "Rz": dms2deg(0, 0, -4.102655),
            "m": 0
        }
        
        self.grs80_krasovsky_params = {
            "dX": 66.171809,
            "dY": 23.995845,
            "dZ": 83.769412,
            "Rx": dms2deg(0, 0, 2.50300846),
            "Ry": dms2deg(0, 0, 2.19216164),
            "Rz": dms2deg(0, 0, -2.58393360),
            "m": 0
        }
        
        # S-42/83/03 Zone 4 parameters
        self.s42_lon0_zone4 = np.radians(27.0)  # Central meridian for zone 4
        self.s42_false_easting = 4500000.0
        self.s42_false_northing = 0.0
        self.s42_scale = 0.9996
    
    def krovak_xy2ell(self, bod):
        """
        Inverse Krovak transformation: XY (projected) -> ellipsoidal coordinates
        Based on lecture pages 25-30 formulas
        
        Args:
            bod: Dictionary with 'x' and 'y' coordinates in JTSK03
            
        Returns:
            Dictionary with added 'phi' and 'lambda' (latitude and longitude in radians)
        """
        p = {}
        p["x"], p["y"] = bod["x"], bod["y"]
        
        # Calculate polar coordinates from Cartesian coordinates
        p["rho"] = np.sqrt((p["x"]**2) + (p["y"]**2))
        p["D_eps"] = np.arctan2(p["y"], p["x"])
        
        # Calculate cartographic coordinates from polar coordinates
        p["S"] = 2 * np.arctan(((self.krovak_rho_0 / p["rho"])**(1/self.krovak_n)) * 
                               np.tan((self.krovak_S_0 / 2) + np.pi/4)) - np.pi/2
        p["D"] = p["D_eps"] / np.sin(self.krovak_S_0)
        
        # Calculate spheric coordinates from cartographic coordinates
        p["U"] = np.arcsin(np.cos(self.krovak_a) * np.sin(p["S"]) - 
                          np.sin(self.krovak_a) * np.cos(p["S"]) * np.cos(p["D"]))
        p["V"] = self.krovak_V_k - np.arcsin((np.cos(p["S"]) / np.cos(p["U"])) * np.sin(p["D"]))
        
        # Calculate ellipsoidal coordinates from spheric coordinates
        p["lambda"] = (p["V"] / self.krovak_alpha) - np.radians(dms2deg(17, 40, 0))
        
        Q = np.log(np.tan((p["U"] / 2) + np.pi/4))
        q = (1 / self.krovak_alpha) * (Q - np.log(self.krovak_k))
        
        # Iterative calculation of latitude
        diff = 1
        iteration = 0
        phi_0 = 2 * np.arctan(np.e**q) - (np.pi / 2)
        
        while diff >= 0.0000000000001:
            phi = 2 * np.arctan((np.e**q) / (
                np.sqrt(((1 - self.Bessel["e"] * np.sin(phi_0)) / 
                        (1 + self.Bessel["e"] * np.sin(phi_0)))**self.Bessel["e"]))) - (np.pi / 2)
            iteration = iteration + 1
            diff = abs(phi - phi_0)
            phi_0 = phi
            p["phi"] = phi
        
        bod["lambda"] = p["lambda"]
        bod["phi"] = p["phi"]
        
        return bod
    
    def ell2xyz(self, phi_a, lambda_a, h_a, ell):
        """Convert ellipsoidal to cartesian coordinates"""
        N = lateral(ell, phi_a)
        
        X = (N + h_a) * np.cos(phi_a) * np.cos(lambda_a)
        Y = (N + h_a) * np.cos(phi_a) * np.sin(lambda_a)
        Z = ((1 - ell["e"]**2) * N + h_a) * np.sin(phi_a)
        
        return X, Y, Z
    
    def xyz2ell(self, X, Y, Z, ell):
        """Convert cartesian to ellipsoidal coordinates"""
        p = np.sqrt(X**2 + Y**2)
        
        lambda_a = np.arctan2(Y, X)
        
        diff = 1
        iteration = 0
        t_0 = Z / ((1 - ell["e"]**2) * p)
        phi_0 = np.arctan(t_0)
        
        while diff >= 0.00000000000001:
            t = Z / (p - ((ell["a"] * ell["e"]**2) / (
                np.sqrt(1 + (1 - ell["e"]**2) * (t_0**2)))))
            iteration = iteration + 1
            phi = np.arctan(t)
            diff = abs(phi - phi_0)
            t_0 = t
            phi_0 = phi
            
            phi_a = phi
        
        h_a = np.sqrt(1 + t**2) * (p - ((ell["a"]) / (
            np.sqrt(1 + (1 - ell["e"]**2) * (t**2)))))
        
        return phi_a, lambda_a, h_a
    
    def bursa_wolf_transform(self, Xa, Ya, Za, params, full_matrix=False):
        """Bursa-Wolf 7-parameter transformation"""
        # Convert parameters from degrees to radians
        dX = params['dX']
        dY = params['dY']
        dZ = params['dZ']
        Rx = np.radians(params['Rx'])
        Ry = np.radians(params['Ry'])
        Rz = np.radians(params['Rz'])
        m = params['m'] * 1e-6
        
        if full_matrix:
            # Full rotation matrix
            R = np.array([
                [np.cos(Ry)*np.cos(Rz),
                 np.cos(Rx)*np.sin(Rz) + np.sin(Rx)*np.sin(Ry)*np.cos(Rz),
                 np.sin(Rx)*np.sin(Rz) - np.cos(Rx)*np.sin(Ry)*np.cos(Rz)],
                [-np.cos(Ry)*np.sin(Rz),
                 np.cos(Rx)*np.cos(Rz) - np.sin(Rx)*np.sin(Ry)*np.sin(Rz),
                 np.sin(Rx)*np.cos(Rz) + np.cos(Rx)*np.sin(Ry)*np.sin(Rz)],
                [np.sin(Ry),
                 -np.sin(Rx)*np.cos(Ry),
                 np.cos(Rx)*np.cos(Ry)]
            ])
        else:
            # Small angles approximation
            R = np.array([
                [1, Rz, -Ry],
                [-Rz, 1, Rx],
                [Ry, -Rx, 1]
            ])
        
        A = np.array([Xa, Ya, Za])
        B = np.array([dX, dY, dZ]) + (1 + m) * (R @ A)
        
        return tuple(B)
    
    def GK_ell2xy(self, bod):
        """
        Gauss-Krüger projection: ellipsoidal -> projected coordinates
        For S-42/83/03 Zone 4
        """
        p = {}
        p["phi"] = bod["phi"]
        p["lambda"] = bod["lambda"]
        
        # Zone 4 central meridian: 27°E
        lambda_0 = np.radians(27.0)
        p["lambda"] = p["lambda"] - lambda_0
        
        # Meridional arc
        def fx_M(phi):
            return meridian(self.Krasovsky, phi)
        
        Sp = quad(fx_M, 0, p["phi"])[0]
        
        eta = self.Krasovsky["e2"] * np.cos(p["phi"])
        t = np.tan(p["phi"])
        N = lateral(self.Krasovsky, p["phi"])
        
        # Gauss-Krüger projection formulas
        p["x"] = (Sp + 
                  N * np.sin(p["phi"]) * np.cos(p["phi"]) * (p["lambda"]**2 / 2) +
                  N * np.sin(p["phi"]) * np.cos(p["phi"])**3 * 
                  (5 - t**2 + 9*eta**2 + 4*eta**4) * (p["lambda"]**4 / 24))
        
        p["y"] = (N * p["lambda"] * np.cos(p["phi"]) +
                  N * np.cos(p["phi"])**3 * (1 - t**2 + eta**2) * (p["lambda"]**3 / 6) +
                  N * np.cos(p["phi"])**5 * 
                  (5 - 18*t**2 + t**4 + 14*eta**2 - 58*eta**2*t**2) * (p["lambda"]**5 / 120))
        
        # Add false easting for Zone 4
        p["y"] += self.s42_false_easting
        
        bod["x"] = p["x"]
        bod["y"] = p["y"]
        
        return bod
    
    def transform(self, x_jtsk03, y_jtsk03):
        """
        Main transformation: JTSK03 -> Geographic (Bessel) -> 
        Cartesian (Bessel) -> Cartesian (GRS80) -> Cartesian (Krasovsky) ->
        Geographic (Krasovsky) -> S-42/83/03 Zone 4
        
        Args:
            x_jtsk03: X coordinate in JTSK03 (meters)
            y_jtsk03: Y coordinate in JTSK03 (meters)
            
        Returns:
            dict: Transformation results
        """
        print("\n" + "="*80)
        print("TRANSFORMATION: JTSK03 -> S-42/83/03 Zone 4")
        print("="*80)
        print(f"\nInput JTSK03 coordinates:")
        print(f"  X = {x_jtsk03:.3f} m")
        print(f"  Y = {y_jtsk03:.3f} m")
        print(f"\n{'-'*80}")
        
        # Step 1: Krovak inverse transformation
        print("\nStep 1: Krovak Inverse Transformation (XY -> Geographic on Bessel)")
        print("-"*80)
        P = {"x": x_jtsk03, "y": y_jtsk03}
        P = self.krovak_xy2ell(P)
        print(f"Geographic coordinates (Bessel ellipsoid):")
        print(f"  φ = {print_dms(P['phi'], 6)}")
        print(f"  λ = {print_dms(P['lambda'], 6)}")
        
        # Step 2: Geographic to Cartesian (Bessel)
        print(f"\nStep 2: Geographic -> Cartesian (Bessel ellipsoid)")
        print("-"*80)
        X_bessel, Y_bessel, Z_bessel = self.ell2xyz(P["phi"], P["lambda"], 0, self.Bessel)
        print(f"Cartesian coordinates (Bessel):")
        print(f"  X = {X_bessel:.3f} m")
        print(f"  Y = {Y_bessel:.3f} m")
        print(f"  Z = {Z_bessel:.3f} m")
        
        # Step 3: Datum shift Bessel -> GRS80
        print(f"\nStep 3: Datum Transformation (Bessel -> GRS80)")
        print("-"*80)
        X_grs80, Y_grs80, Z_grs80 = self.bursa_wolf_transform(
            X_bessel, Y_bessel, Z_bessel, 
            self.bessel_grs80_params, full_matrix=False
        )
        print(f"Cartesian coordinates (GRS80):")
        print(f"  X = {X_grs80:.3f} m")
        print(f"  Y = {Y_grs80:.3f} m")
        print(f"  Z = {Z_grs80:.3f} m")
        
        # Step 4: Datum shift GRS80 -> Krasovsky
        print(f"\nStep 4: Datum Transformation (GRS80 -> Krasovsky)")
        print("-"*80)
        X_kras, Y_kras, Z_kras = self.bursa_wolf_transform(
            X_grs80, Y_grs80, Z_grs80,
            self.grs80_krasovsky_params, full_matrix=True
        )
        print(f"Cartesian coordinates (Krasovsky):")
        print(f"  X = {X_kras:.3f} m")
        print(f"  Y = {Y_kras:.3f} m")
        print(f"  Z = {Z_kras:.3f} m")
        
        # Step 5: Cartesian to Geographic (Krasovsky)
        print(f"\nStep 5: Cartesian -> Geographic (Krasovsky ellipsoid)")
        print("-"*80)
        phi_kras, lambda_kras, h_kras = self.xyz2ell(X_kras, Y_kras, Z_kras, self.Krasovsky)
        print(f"Geographic coordinates (Krasovsky ellipsoid):")
        print(f"  φ = {print_dms(phi_kras, 6)}")
        print(f"  λ = {print_dms(lambda_kras, 6)}")
        
        # Step 6: Geographic to S-42/83/03 Zone 4
        print(f"\nStep 6: Gauss-Krüger Projection -> S-42/83/03 Zone 4")
        print("-"*80)
        P_kras = {"phi": phi_kras, "lambda": lambda_kras}
        P_kras = self.GK_ell2xy(P_kras)
        
        print(f"\nFinal S-42/83/03 Zone 4 coordinates:")
        print(f"  X = {P_kras['x']:.3f} m")
        print(f"  Y = {P_kras['y']:.3f} m")
        print(f"\n{'='*80}\n")
        
        return {
            'x_jtsk03': x_jtsk03,
            'y_jtsk03': y_jtsk03,
            'latitude_bessel_deg': math.degrees(P['phi']),
            'longitude_bessel_deg': math.degrees(P['lambda']),
            'latitude_krasovsky_deg': math.degrees(phi_kras),
            'longitude_krasovsky_deg': math.degrees(lambda_kras),
            'x_s42': P_kras['x'],
            'y_s42': P_kras['y']
        }


def main():
    """Main function - example usage"""
    transformer = CoordinateTransformer()
    
    print("\n" + "="*80)
    print("JTSK03 to S-42/83/03 Zone 4 Coordinate Transformation")
    print("Corrected Implementation Based on Lecture Pages 25-30")
    print("="*80)
    
    # Test case from lecture: Trenčiansky hrad (Trenčín Castle)
    print("\n\nTest Case: Point on Trenčiansky hrad")
    print("-"*80)
    
    x_test = 1204350.498
    y_test = 496927.409
    
    result = transformer.transform(x_test, y_test)
    
    print("\nExpected Results (from lecture):")
    print(f"  Geographic (Bessel):    φ ≈ 48°53'40.95\", λ ≈ 18°02'46.38\"")
    print(f"  Geographic (Krasovsky): φ ≈ 48°53'41.48\", λ ≈ 18°02'44.51\"")
    print(f"  S-42/83/03 Zone 4:      X ≈ 5422203.220 m, Y ≈ 4283434.780 m")
    print(f"\nActual Results:")
    print(f"  Geographic (Bessel):    φ = {result['latitude_bessel_deg']:.6f}°, λ = {result['longitude_bessel_deg']:.6f}°")
    print(f"  Geographic (Krasovsky): φ = {result['latitude_krasovsky_deg']:.6f}°, λ = {result['longitude_krasovsky_deg']:.6f}°")
    print(f"  S-42/83/03 Zone 4:      X = {result['x_s42']:.3f} m, Y = {result['y_s42']:.3f} m")
    
    # Interactive mode
    print("\n" + "="*80)
    print("Interactive Mode - Enter Your Own Coordinates")
    print("="*80)
    print("Enter JTSK03 coordinates to transform (or 'q' to quit)")
    print()
    
    try:
        while True:
            try:
                x_input = input("Enter X coordinate (JTSK03, meters) or 'q' to quit: ").strip()
                if x_input.lower() == 'q':
                    break
                
                x_val = float(x_input)
                y_val = float(input("Enter Y coordinate (JTSK03, meters): "))
                
                result = transformer.transform(x_val, y_val)
                
            except ValueError as e:
                print(f"Invalid input: {e}\n")
    except (KeyboardInterrupt, EOFError):
        print("\n\nTransformation completed!")


if __name__ == "__main__":
    main()
