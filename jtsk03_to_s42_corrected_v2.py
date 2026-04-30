#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
JTSK03 to S-42/83/03 Zone 4 Coordinate Transformation - CORRECTED VERSION 2
Based on lecture pages 25-30 and verified parameters

This implements the complete transformation pipeline:
JTSK03 (Krovak/Bessel) -> Bessel Geographic -> GRS80 Cartesian -> 
Krasovsky Cartesian -> Krasovsky Geographic -> S-42 Zone 4
"""

import math
import numpy as np
from scipy.integrate import quad


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
    """Main transformation class"""
    
    def __init__(self):
        # Bessel 1841 ellipsoid (JTSK03)
        self.bessel = {
            'a': 6377397.155,
            'b': 6356078.96290,
            'e': 0.081696831,
            'e2': 0.00667346
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
        
        # Krovak parameters (from lecture page 25)
        self.krovak_alpha = 1.000597498372
        self.krovak_k = 1.003419164
        self.krovak_lon0 = np.radians(dms2deg(17, 40, 0))  # Prime meridian offset
        self.krovak_U_k = np.radians(dms2deg(59, 42, 42.69689))
        self.krovak_V_k = np.radians(dms2deg(42, 31, 31.417251))
        self.krovak_pole_angle = np.pi / 2 - self.krovak_U_k
        self.krovak_S0 = np.radians(dms2deg(78, 30, 0))
        self.krovak_n = np.sin(self.krovak_S0)
        self.krovak_R = 6380703.61054
        self.krovak_rho0 = 0.9999 * self.krovak_R / np.tan(self.krovak_S0)
        
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
        self.s42_zone4_lon0 = np.radians(27.0)  # 27°E
        self.s42_false_easting = 4500000.0
        self.s42_scale = 0.9996
    
    def krovak_inverse(self, x, y):
        """
        Krovak inverse transformation: projected (x,y) -> geographic (phi, lambda)
        Implementation from lecture pages 25-30
        """
        # Step 1: Polar coordinates
        rho = np.sqrt(x**2 + y**2)
        eps = np.arctan2(y, x)
        
        # Step 2: Cartographic coordinates
        S = 2 * np.arctan(
            (self.krovak_rho0 / rho) ** (1/self.krovak_n) * 
            np.tan(self.krovak_S0/2 + np.pi/4)
        ) - np.pi/2
        
        D = eps / np.sin(self.krovak_S0)
        
        # Step 3: Spheric coordinates (on sphere)
        U = np.arcsin(
            np.cos(self.krovak_pole_angle) * np.sin(S) - 
            np.sin(self.krovak_pole_angle) * np.cos(S) * np.cos(D)
        )
        
        V = self.krovak_V_k - np.arcsin(
            (np.cos(S) / np.cos(U)) * np.sin(D)
        )
        
        # Step 4: Geographic coordinates (ellipsoidal)
        lambda_rad = V / self.krovak_alpha - self.krovak_lon0
        
        # Conformal latitude calculation
        Q = np.log(np.tan(np.pi/4 + U/2))
        q = (1 / self.krovak_alpha) * (Q - np.log(self.krovak_k))
        
        # Iterative latitude calculation
        phi = 2 * np.arctan(np.e**q) - np.pi/2
        
        for _ in range(10):
            sin_phi = np.sin(phi)
            sin_term = ((1 - self.bessel['e'] * sin_phi) / 
                       (1 + self.bessel['e'] * sin_phi)) ** self.bessel['e']
            phi_new = 2 * np.arctan((np.e**q) / sin_term) - np.pi/2
            
            if abs(phi_new - phi) < 1e-12:
                break
            phi = phi_new
        
        return phi, lambda_rad
    
    def geodetic_to_cartesian(self, phi, lam, ellipsoid):
        """Convert geodetic to Cartesian coordinates"""
        a = ellipsoid['a']
        e2 = ellipsoid['e2']
        
        N = a / np.sqrt(1 - e2 * np.sin(phi)**2)
        
        X = N * np.cos(phi) * np.cos(lam)
        Y = N * np.cos(phi) * np.sin(lam)
        Z = (N * (1 - e2)) * np.sin(phi)
        
        return X, Y, Z
    
    def cartesian_to_geodetic(self, X, Y, Z, ellipsoid):
        """Convert Cartesian to geodetic coordinates (Bowring's method)"""
        a = ellipsoid['a']
        e2 = ellipsoid['e2']
        
        lam = np.arctan2(Y, X)
        
        p = np.sqrt(X**2 + Y**2)
        e_prime2 = e2 / (1 - e2)
        
        # Initial approximation
        phi = np.arctan2(Z, p * (1 - e2))
        
        # Iterate
        for _ in range(10):
            sin_phi = np.sin(phi)
            N = a / np.sqrt(1 - e2 * sin_phi**2)
            phi_new = np.arctan2(Z + e2 * N * sin_phi, p)
            
            if abs(phi_new - phi) < 1e-12:
                break
            phi = phi_new
        
        return phi, lam
    
    def bursa_wolf(self, X, Y, Z, params, full_matrix=False):
        """7-parameter Bursa-Wolf transformation"""
        dX = params['dX']
        dY = params['dY']
        dZ = params['dZ']
        Rx = params['Rx']
        Ry = params['Ry']
        Rz = params['Rz']
        m = params['m'] * 1e-6
        
        if full_matrix:
            # Full rotation matrix
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
        Gauss-Krüger (Transverse Mercator) forward projection
        For S-42 Zone 4 on Krasovsky ellipsoid
        """
        a = self.krasovsky['a']
        e2 = self.krasovsky['e2']
        e_prime2 = e2 / (1 - e2)
        
        # Reduce longitude to central meridian
        lon_diff = lam - self.s42_zone4_lon0
        
        # Meridional arc (using numerical integration)
        def meridional_arc(lat):
            n = np.sqrt(e2) / np.sqrt(1 - e2)
            n2 = n**2
            n3 = n**3
            
            A0 = 1 - e2/4 - 3*e2**2/64 - 5*e2**3/256
            A2 = (3/8)*(e2 + e2**2/4 - 15*e2**3/128)
            A4 = (15/256)*(e2**2 - 3*e2**3/4)
            A6 = (-35/3072)*e2**3
            
            m = a * (A0*lat - A2*np.sin(2*lat) + A4*np.sin(4*lat) - A6*np.sin(6*lat))
            return m
        
        M = meridional_arc(phi)
        
        N = a / np.sqrt(1 - e2 * np.sin(phi)**2)
        T = np.tan(phi)**2
        C = e_prime2 * np.cos(phi)**2
        A = np.cos(phi) * lon_diff
        
        # Gauss-Krüger formulas
        x = self.s42_scale * (M + N * np.tan(phi) * 
                             (A**2/2 + A**4/24 * (5 - T + 9*C + 4*C**2) +
                              A**6/720 * (61 - 58*T + T**2 + 600*C - 330*e_prime2)))
        
        y = self.s42_scale * (N * (A + A**3/6 * (1 - T + C) +
                                   A**5/120 * (5 - 18*T + T**2 + 72*C - 58*e_prime2)))
        
        y += self.s42_false_easting
        
        return x, y
    
    def transform(self, x_jtsk03, y_jtsk03, verbose=True):
        """
        Main transformation pipeline
        JTSK03 -> Bessel -> GRS80 -> Krasovsky -> S-42/83/03
        """
        if verbose:
            print("\n" + "="*80)
            print("JTSK03 to S-42/83/03 Zone 4 Transformation")
            print("="*80)
            print(f"\n[INPUT] JTSK03 Coordinates:")
            print(f"  X = {x_jtsk03:12.3f} m")
            print(f"  Y = {y_jtsk03:12.3f} m")
        
        # STEP 1: Krovak inverse
        if verbose:
            print(f"\n[STEP 1] Krovak Inverse Transformation")
            print("-"*80)
        
        phi_bessel, lam_bessel = self.krovak_inverse(x_jtsk03, y_jtsk03)
        
        if verbose:
            d, m, s = deg2dms(np.degrees(phi_bessel))
            print(f"  φ (Bessel) = {d}°{m:02d}'{s:06.3f}\"")
            d, m, s = deg2dms(np.degrees(lam_bessel))
            print(f"  λ (Bessel) = {d}°{m:02d}'{s:06.3f}\"")
        
        # STEP 2: Bessel geodetic -> Cartesian
        if verbose:
            print(f"\n[STEP 2] Geodetic to Cartesian (Bessel)")
            print("-"*80)
        
        X_b, Y_b, Z_b = self.geodetic_to_cartesian(phi_bessel, lam_bessel, self.bessel)
        
        if verbose:
            print(f"  X = {X_b:12.3f} m")
            print(f"  Y = {Y_b:12.3f} m")
            print(f"  Z = {Z_b:12.3f} m")
        
        # STEP 3: Datum shift Bessel -> GRS80
        if verbose:
            print(f"\n[STEP 3] Datum Transformation (Bessel -> GRS80)")
            print("-"*80)
        
        X_g, Y_g, Z_g = self.bursa_wolf(X_b, Y_b, Z_b, self.bw_bessel_grs80, full_matrix=False)
        
        if verbose:
            print(f"  X = {X_g:12.3f} m")
            print(f"  Y = {Y_g:12.3f} m")
            print(f"  Z = {Z_g:12.3f} m")
        
        # STEP 4: Datum shift GRS80 -> Krasovsky
        if verbose:
            print(f"\n[STEP 4] Datum Transformation (GRS80 -> Krasovsky)")
            print("-"*80)
        
        X_k, Y_k, Z_k = self.bursa_wolf(X_g, Y_g, Z_g, self.bw_grs80_krasovsky, full_matrix=True)
        
        if verbose:
            print(f"  X = {X_k:12.3f} m")
            print(f"  Y = {Y_k:12.3f} m")
            print(f"  Z = {Z_k:12.3f} m")
        
        # STEP 5: Krasovsky Cartesian -> Geodetic
        if verbose:
            print(f"\n[STEP 5] Cartesian to Geodetic (Krasovsky)")
            print("-"*80)
        
        phi_kras, lam_kras = self.cartesian_to_geodetic(X_k, Y_k, Z_k, self.krasovsky)
        
        if verbose:
            d, m, s = deg2dms(np.degrees(phi_kras))
            print(f"  φ (Krasovsky) = {d}°{m:02d}'{s:06.3f}\"")
            d, m, s = deg2dms(np.degrees(lam_kras))
            print(f"  λ (Krasovsky) = {d}°{m:02d}'{s:06.3f}\"")
        
        # STEP 6: Gauss-Krüger projection for S-42
        if verbose:
            print(f"\n[STEP 6] Gauss-Krüger Projection (S-42 Zone 4)")
            print("-"*80)
        
        x_s42, y_s42 = self.gauss_kruger_forward(phi_kras, lam_kras)
        
        if verbose:
            print(f"\n[OUTPUT] S-42/83/03 Zone 4 Coordinates:")
            print(f"  X = {x_s42:12.3f} m")
            print(f"  Y = {y_s42:12.3f} m")
            print(f"\n{'='*80}\n")
        
        return {
            'x_s42': x_s42,
            'y_s42': y_s42,
            'phi_bessel_deg': np.degrees(phi_bessel),
            'lam_bessel_deg': np.degrees(lam_bessel),
            'phi_krasovsky_deg': np.degrees(phi_kras),
            'lam_krasovsky_deg': np.degrees(lam_kras),
        }


def main():
    transformer = JTSK03ToS42Transformer()
    
    print("\n" + "="*80)
    print("JTSK03 ↔ S-42/83/03 Zone 4 Coordinate Transformation")
    print("Corrected Implementation - Version 2")
    print("="*80)
    
    # Test case: Trenčín Castle (Trenčiansky hrad)
    print("\n\n[TEST CASE] Trenčín Castle (Trenčiansky hrad)")
    print("-"*80)
    
    x_test = 1204350.498
    y_test = 496927.409
    
    result = transformer.transform(x_test, y_test, verbose=True)
    
    print("\n[VERIFICATION]")
    print("-"*80)
    print(f"Expected S-42 output: X ≈ 5422203.220 m, Y ≈ 4283434.780 m")
    print(f"Actual S-42 output:   X = {result['x_s42']:12.3f} m, Y = {result['y_s42']:12.3f} m")
    print(f"Error in X: {abs(result['x_s42'] - 5422203.220):.3f} m")
    print(f"Error in Y: {abs(result['y_s42'] - 4283434.780):.3f} m")
    
    # Interactive mode
    print("\n" + "="*80)
    print("Interactive Mode")
    print("="*80)
    
    try:
        while True:
            print("\nEnter JTSK03 coordinates (or 'q' to quit):")
            x_input = input("  X [meters]: ").strip()
            
            if x_input.lower() == 'q':
                break
            
            try:
                x_val = float(x_input)
                y_val = float(input("  Y [meters]: "))
                
                result = transformer.transform(x_val, y_val, verbose=True)
                
            except ValueError as e:
                print(f"Error: {e}")
                
    except (KeyboardInterrupt, EOFError):
        print("\n\nGoodbye!")


if __name__ == "__main__":
    main()
