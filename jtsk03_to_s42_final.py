#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
JTSK03 to S-42/83/03 Zone 4 Coordinate Transformation
Based on Vajsablova: Matematicke zaklady kartografie
Correct Krovak inverse formulas using polar coordinates

Transformation chain:
JTSK03 (Krovak/Bessel) -> Geographic (Bessel) -> 
Cartesian (Bessel) -> Cartesian (GRS80) -> 
Cartesian (Krasovsky) -> Geographic (Krasovsky) -> S-42/83/03
"""

import math
import numpy as np


def dms2deg(degrees, minutes, seconds):
    """Convert DMS to decimal degrees"""
    return degrees + minutes / 60.0 + seconds / 3600.0


def deg2dms(decimal_degrees):
    """Convert decimal degrees to DMS format"""
    d = int(decimal_degrees)
    m = int((abs(decimal_degrees) - abs(d)) * 60)
    s = ((abs(decimal_degrees) - abs(d)) * 60 - m) * 60
    return d, m, s


class JTSK03ToS42Transformer:
    """Coordinate transformation JTSK03 -> S-42/83/03"""
    
    def __init__(self):
        # Bessel 1841 ellipsoid (JTSK03)
        self.bessel = {
            'a': 6377397.155,
            'b': 6356078.96290,
            'e': 0.081696831,
            'e2': 0.00667346,
            'n': 0.0823481
        }
        
        # GRS80 ellipsoid
        self.grs80 = {
            'a': 6378137.0,
            'b': 6356752.314245,
            'e': 0.0818191908,
            'e2': 0.00669438
        }
        
        # Krasovsky ellipsoid (S-42)
        self.krasovsky = {
            'a': 6378245.0,
            'b': 6356863.01887,
            'e': 0.0818133,
            'e2': 0.006693421623
        }
        
        # ============ KROVAK PROJECTION PARAMETERS ============
        # From Vajsablova and JTSK standard
        
        # Fundamental constants
        self.alpha = 1.000597498372  # Conformal latitude factor
        self.k = 1.003419164          # Scale factor
        
        # Projection center coordinates
        self.lat0 = np.radians(dms2deg(49, 30, 0))      # Φ₀ = 49°30'
        self.lon0 = np.radians(dms2deg(24, 50, 0))      # Λ₀ = 24°50'
        
        # Auxiliary angles in projection
        self.U0 = np.radians(dms2deg(59, 42, 42.69689)) # U₀
        self.V0 = np.radians(dms2deg(42, 31, 31.417251))# V₀
        self.A = np.pi / 2 - self.U0                    # A = π/2 - U₀
        
        # Cone parameters
        self.S0 = np.radians(dms2deg(78, 30, 0))        # S₀ = 78°30'
        self.n = np.sin(self.S0)                        # n = sin(S₀)
        
        # Earth radius for cone
        self.R = 6380703.61054
        
        # Parallel of standard latitude cone radius
        self.rho0 = 0.9999 * self.R / np.tan(self.S0)
        
        # Bursa-Wolf 7-parameter transformations
        self.bw_bessel_grs80 = {
            'dX': 485.021,
            'dY': 169.465,
            'dZ': 483.839,
            'Rx': np.radians(dms2deg(0, 0, -7.786342)),
            'Ry': np.radians(dms2deg(0, 0, -4.397554)),
            'Rz': np.radians(dms2deg(0, 0, -4.102655)),
            'm': 0
        }
        
        self.bw_grs80_krasovsky = {
            'dX': 66.171809,
            'dY': 23.995845,
            'dZ': 83.769412,
            'Rx': np.radians(dms2deg(0, 0, 2.50300846)),
            'Ry': np.radians(dms2deg(0, 0, 2.19216164)),
            'Rz': np.radians(dms2deg(0, 0, -2.58393360)),
            'm': 0
        }
        
        # S-42 Zone 4 parameters
        self.s42_zone4_lon0 = np.radians(27.0)
        self.s42_false_easting = 4500000.0
        self.s42_scale = 0.9996
    
    def krovak_inverse(self, x, y):
        """
        Krovak inverse transformation using polar coordinates
        Input: x, y (JTSK03 projected coordinates)
        Output: phi, lam (geodetic coordinates on Bessel ellipsoid)
        
        Based on Vajsablova formulas for Krovak projection
        """
        # Step 1: Convert Cartesian to polar coordinates
        # Note: Krovak uses special coordinate system (may have different origin)
        rho = np.sqrt(x**2 + y**2)
        D_eps = np.arctan2(y, x)
        
        # Step 2: Calculate cartographic latitude S
        # From cone standard parallel formula
        S = 2 * np.arctan(
            (self.rho0 / rho) ** (1 / self.n) * 
            np.tan(self.S0 / 2 + np.pi / 4)
        ) - np.pi / 2
        
        # Step 3: Calculate auxiliary angle D (conformal longitude)
        D = D_eps / np.sin(self.S0)
        
        # Step 4: Calculate spheric coordinates U, V
        # (from cartographic to spheric on auxiliary sphere)
        
        # U = asin(cos(A) * sin(S) - sin(A) * cos(S) * cos(D))
        U = np.arcsin(
            np.cos(self.A) * np.sin(S) - 
            np.sin(self.A) * np.cos(S) * np.cos(D)
        )
        
        # V = V0 - asin(cos(S) / cos(U) * sin(D))
        V = self.V0 - np.arcsin(
            (np.cos(S) / np.cos(U)) * np.sin(D)
        )
        
        # Step 5: Convert spheric to ellipsoidal coordinates
        # Lambda calculation
        lam = (V / self.alpha) - self.lon0
        
        # Conformal latitude Q and q
        Q = np.log(np.tan(np.pi / 4 + U / 2))
        q = (1 / self.alpha) * (Q - np.log(self.k))
        
        # Step 6: Iterative calculation of geodetic latitude
        # phi = 2 * arctan((e^q) / sqrt(...)) - π/2
        
        e = self.bessel['e']
        
        # Initial approximation
        phi = 2 * np.arctan(np.e ** q) - np.pi / 2
        
        # Iterate to refine
        for i in range(20):
            sin_phi = np.sin(phi)
            denom = np.sqrt(((1 - e * sin_phi) / (1 + e * sin_phi)) ** e)
            phi_new = 2 * np.arctan((np.e ** q) / denom) - np.pi / 2
            
            if abs(phi_new - phi) < 1e-14:
                break
            phi = phi_new
        
        return phi, lam
    
    def geodetic_to_cartesian(self, phi, lam, ellipsoid):
        """Geodetic (phi, lam) -> Cartesian (X, Y, Z)"""
        a = ellipsoid['a']
        e2 = ellipsoid['e2']
        
        N = a / np.sqrt(1 - e2 * np.sin(phi)**2)
        
        X = N * np.cos(phi) * np.cos(lam)
        Y = N * np.cos(phi) * np.sin(lam)
        Z = N * (1 - e2) * np.sin(phi)
        
        return X, Y, Z
    
    def cartesian_to_geodetic(self, X, Y, Z, ellipsoid):
        """Cartesian (X, Y, Z) -> Geodetic (phi, lam)"""
        a = ellipsoid['a']
        e2 = ellipsoid['e2']
        
        lam = np.arctan2(Y, X)
        p = np.sqrt(X**2 + Y**2)
        
        # Iterative method (Bowring)
        phi = np.arctan2(Z, p * (1 - e2))
        
        for _ in range(20):
            sin_phi = np.sin(phi)
            N = a / np.sqrt(1 - e2 * sin_phi**2)
            phi_new = np.arctan2(Z + e2 * N * sin_phi, p)
            
            if abs(phi_new - phi) < 1e-14:
                break
            phi = phi_new
        
        return phi, lam
    
    def bursa_wolf(self, X, Y, Z, params, full_matrix=False):
        """7-parameter Bursa-Wolf datum transformation"""
        dX = params['dX']
        dY = params['dY']
        dZ = params['dZ']
        Rx = params['Rx']
        Ry = params['Ry']
        Rz = params['Rz']
        m = params['m'] * 1e-6
        
        if full_matrix:
            # Full rotation matrix (all elements computed)
            Rot = np.array([
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
            # Small angle approximation
            Rot = np.array([
                [1, Rz, -Ry],
                [-Rz, 1, Rx],
                [Ry, -Rx, 1]
            ])
        
        point = np.array([X, Y, Z])
        result = np.array([dX, dY, dZ]) + (1 + m) * (Rot @ point)
        
        return result
    
    def gauss_kruger_forward(self, phi, lam):
        """
        Gauss-Krüger forward projection for S-42 Zone 4
        Using Krasovsky ellipsoid
        """
        a = self.krasovsky['a']
        e2 = self.krasovsky['e2']
        e_prime2 = e2 / (1 - e2)
        
        lon_diff = lam - self.s42_zone4_lon0
        
        # Meridional arc (series expansion)
        n = np.sqrt(e2) / np.sqrt(1 - e2)
        n2, n3 = n**2, n**3
        
        A0 = 1 - e2/4 - 3*e2**2/64 - 5*e2**3/256
        A2 = 3/8 * (e2 + e2**2/4 - 15*e2**3/128)
        A4 = 15/256 * (e2**2 - 3*e2**3/4)
        A6 = -35/3072 * e2**3
        
        M = a * (A0*phi - A2*np.sin(2*phi) + A4*np.sin(4*phi) - A6*np.sin(6*phi))
        
        # Latitude-dependent terms
        N = a / np.sqrt(1 - e2 * np.sin(phi)**2)
        T = np.tan(phi)**2
        C = e_prime2 * np.cos(phi)**2
        A = np.cos(phi) * lon_diff
        
        # Gauss-Krüger series
        x = self.s42_scale * (
            M + 
            N * np.tan(phi) * A**2 / 2 +
            N * np.tan(phi) * A**4 / 24 * (5 - T + 9*C + 4*C**2) +
            N * np.tan(phi) * A**6 / 720 * (61 - 58*T + T**2 + 600*C - 330*e_prime2)
        )
        
        y = self.s42_scale * (
            N * A +
            N * A**3 / 6 * (1 - T + C) +
            N * A**5 / 120 * (5 - 18*T + T**2 + 72*C - 58*e_prime2)
        )
        
        y += self.s42_false_easting
        
        return x, y
    
    def transform(self, x_jtsk03, y_jtsk03, verbose=True):
        """Complete JTSK03 -> S-42/83/03 transformation"""
        
        if verbose:
            print("\n" + "="*80)
            print("JTSK03 (Krovak) to S-42/83/03 Zone 4 Transformation")
            print("="*80)
            print(f"\n[INPUT] JTSK03 Coordinates:")
            print(f"  X = {x_jtsk03:14.3f} m")
            print(f"  Y = {y_jtsk03:14.3f} m\n")
        
        # STEP 1: Krovak inverse
        if verbose:
            print("[STEP 1] Krovak Inverse (XY -> Geodetic on Bessel)")
            print("-"*80)
        
        phi_b, lam_b = self.krovak_inverse(x_jtsk03, y_jtsk03)
        
        if verbose:
            d, m, s = deg2dms(np.degrees(phi_b))
            print(f"  φ_Bessel  = {d:3d}°{m:02d}'{s:09.6f}\"")
            d, m, s = deg2dms(np.degrees(lam_b))
            print(f"  λ_Bessel  = {d:3d}°{m:02d}'{s:09.6f}\"")
        
        # STEP 2: Bessel Geodetic -> Cartesian
        if verbose:
            print("\n[STEP 2] Geodetic -> Cartesian (Bessel)")
            print("-"*80)
        
        X_b, Y_b, Z_b = self.geodetic_to_cartesian(phi_b, lam_b, self.bessel)
        
        if verbose:
            print(f"  X_Bessel  = {X_b:14.3f} m")
            print(f"  Y_Bessel  = {Y_b:14.3f} m")
            print(f"  Z_Bessel  = {Z_b:14.3f} m")
        
        # STEP 3: Datum Shift Bessel -> GRS80
        if verbose:
            print("\n[STEP 3] Datum Transformation (Bessel -> GRS80)")
            print("-"*80)
        
        X_g, Y_g, Z_g = self.bursa_wolf(X_b, Y_b, Z_b, self.bw_bessel_grs80, full_matrix=False)
        
        if verbose:
            print(f"  X_GRS80   = {X_g:14.3f} m")
            print(f"  Y_GRS80   = {Y_g:14.3f} m")
            print(f"  Z_GRS80   = {Z_g:14.3f} m")
        
        # STEP 4: Datum Shift GRS80 -> Krasovsky
        if verbose:
            print("\n[STEP 4] Datum Transformation (GRS80 -> Krasovsky)")
            print("-"*80)
        
        X_k, Y_k, Z_k = self.bursa_wolf(X_g, Y_g, Z_g, self.bw_grs80_krasovsky, full_matrix=True)
        
        if verbose:
            print(f"  X_Krasov  = {X_k:14.3f} m")
            print(f"  Y_Krasov  = {Y_k:14.3f} m")
            print(f"  Z_Krasov  = {Z_k:14.3f} m")
        
        # STEP 5: Krasovsky Cartesian -> Geodetic
        if verbose:
            print("\n[STEP 5] Cartesian -> Geodetic (Krasovsky)")
            print("-"*80)
        
        phi_k, lam_k = self.cartesian_to_geodetic(X_k, Y_k, Z_k, self.krasovsky)
        
        if verbose:
            d, m, s = deg2dms(np.degrees(phi_k))
            print(f"  φ_Krasov  = {d:3d}°{m:02d}'{s:09.6f}\"")
            d, m, s = deg2dms(np.degrees(lam_k))
            print(f"  λ_Krasov  = {d:3d}°{m:02d}'{s:09.6f}\"")
        
        # STEP 6: Gauss-Krüger -> S-42 Zone 4
        if verbose:
            print("\n[STEP 6] Gauss-Krüger Projection (S-42 Zone 4)")
            print("-"*80)
        
        x_s42, y_s42 = self.gauss_kruger_forward(phi_k, lam_k)
        
        if verbose:
            print(f"\n[OUTPUT] S-42/83/03 Zone 4:")
            print(f"  X = {x_s42:14.3f} m")
            print(f"  Y = {y_s42:14.3f} m")
            print(f"\n{'='*80}\n")
        
        return {
            'x_s42': x_s42,
            'y_s42': y_s42,
            'phi_bessel_deg': np.degrees(phi_b),
            'lam_bessel_deg': np.degrees(lam_b),
            'phi_krasovsky_deg': np.degrees(phi_k),
            'lam_krasovsky_deg': np.degrees(lam_k),
        }


def main():
    transformer = JTSK03ToS42Transformer()
    
    print("\n" + "="*80)
    print("JTSK03 to S-42/83/03 Zone 4 - Final Corrected Version")
    print("Based on Vajsablova: Matematicke zaklady kartografie")
    print("="*80)
    
    # Test case: Trenčín Castle
    print("\n\n[TEST CASE] Trenčín Castle (Trenčiansky hrad)")
    print("-"*80)
    
    x_test = 1204350.498
    y_test = 496927.409
    
    result = transformer.transform(x_test, y_test, verbose=True)
    
    print("[VERIFICATION]")
    print("-"*80)
    print(f"Expected S-42: X ≈ 5422203.220 m, Y ≈ 4283434.780 m")
    print(f"Actual S-42:   X = {result['x_s42']:12.3f} m, Y = {result['y_s42']:12.3f} m")
    err_x = abs(result['x_s42'] - 5422203.220)
    err_y = abs(result['y_s42'] - 4283434.780)
    print(f"Error:         X = {err_x:12.3f} m, Y = {err_y:12.3f} m")
    
    # Interactive mode
    print("\n" + "="*80)
    print("Interactive Mode")
    print("="*80)
    
    try:
        while True:
            print("\nEnter JTSK03 coordinates (or 'q' to quit):")
            x_in = input("  X [meters]: ").strip()
            
            if x_in.lower() == 'q':
                break
            
            try:
                x_val = float(x_in)
                y_val = float(input("  Y [meters]: "))
                result = transformer.transform(x_val, y_val, verbose=True)
            except ValueError as e:
                print(f"Error: {e}")
    
    except (KeyboardInterrupt, EOFError):
        print("\n\nGoodbye!")


if __name__ == "__main__":
    main()
