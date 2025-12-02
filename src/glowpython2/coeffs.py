import numpy as np

COEFF_LEN_NR = 100  # Number of coefficients in GLOW, NR

"""
Classic A coefficients:
!      A1  5200       N(4S-2D)     Wiese et al, 1966
!      A2  6300       O(3P-1D)     Baluja and Zeippen, 1988
!      A3  6364       O(3P-1D)     Baluja and Zeippen, 1988
!      A4  2972       O(3P-1S)     Kernahan & Pang, 1975
!      A5  5577       O(1D-1S)     Kernahan & Pang, 1975
!      A6  3726       O+(4S-2D)    Kernahan & Pang, 1975
!      A7  2470       O+(4S-2P)    Weise et al, 1966
!      A8  7319-30    O+(2D-2P)    Weise et al, 1966
!      A9  (Hertz II) O2(X-c)      Solheim & Llewellyn, 1978
!      A10 (Veg-Kap)  N2(X-A)      Shemansky, 1969
!      A11 3466       N(4S-2P)     Chamberlain, 1961
!      A12 10400      N(2D-2P)     Chamberlain, 1961
"""
classic_acoeff = [
    1.07e-5, 0.00585, 0.00185, 0.0450, 1.0600, 9.7e-5,
    0.0479, 0.1712, 0.0010, 0.7700, 0.00540, 0.07900
]
CLASSIC_ACOEFF = np.array(
    classic_acoeff + [0.0]*(COEFF_LEN_NR-len(classic_acoeff)),
    dtype=np.float32, order='F'
)
assert (len(CLASSIC_ACOEFF) == COEFF_LEN_NR)

"""
Classic B coefficients:
!      B1  O(1S) from O2+ + e      Yee et al, 1988
!      B2  O(1D) from O2+ + e      Abreu et al, 1986
!      B3  N(2D) from NO+ + e      Kley et al, 1976
!      B4  N(2D) from N2+ + e      Queffelec et al, 1985
!      B5  N(2D) from N2+ + O      Frederick & Rusch, 1977
!      B6  O(1D) from N(2D) + O2   Link, 1983; Langford et al, 1985
!      B7  O(1D) from O+(2D) + O   ?
!      B8  O+(2D) from O+(2P) + e  Link, 1982
!      B9  O+(2P) from O + e*      Gerard & Rusch, 1979; Jones, 1975
!      B10 O+(2D) from O + e*      Gerard & Rusch, 1979; Jones, 1975
!      B11 O+(4S) from O + e*      Gerard & Rusch, 1979; Jones, 1975
!      B12 O+(2P) from O2 + e*     Link, 1982; guess of .3/3
!      B13 O+(2D) from O2 + e*     Link, 1982; guess of .3/3
!      B14 O+(4S) from O2 + e*     Link, 1982; guess of .3/3
!      B15 N+ from N2 + e*         Richards & Torr, 1985
!      B16 N(2D) from above        Zipf et al, 1980
!      B17 O(1D) from N+ + O2      Langford et al, 1985
!      B18 O(1S) from N2(A) + O    Sharp & Torr, 1979
!      B19 O(1S) from O2(*) + O    ? (= 0 at present)
!      B20 O(1D) from N(2D) + O    ?
!      B21 NO+ from N+ + O2        Langford et al, 1985
!      B22 O2+ from N+ + O2        Langford et al, 1985
!      B23 N(2P) from N2+ + e      Queffelec et al, 1985
!      B24 N2 + protons -> N + N   ?
!      B25 N(2D) from N2 + e* dis  Zipf et al, 1980
!      B26 N(2P) from N2 + e* dis  Zipf et al, 1980
!      B27 N(2D) from N2 + hv      Richards et al, 1981 (add to B28)
!      B28 N(2P) from N2 + hv      ? (cf Zipf & McGlaughlin, 1978)
!      B29 N(2D) from N(2P) + O    ?
!      B30 O+(2P) from O2 + hv     ?
!      B31 O+(2D)  "               ?
!      B32 O+(4S)  "               ?
!      B33 O(1S) from O2+ + N      Frederick et al, 1976; Kopp ea, 1977
!      B34 O(1S) from N(2D) + NO   Frederick et al, 1976; Kopp ea, 1977
!      B35 O2 + protons -> (O1D)   ?
!      B36 N2+(B) from N2 + e*     Borst & Zipf, 1970; Shemansky & Broadfoot, 1971
!      B37 (0,0) (3914) fr. N2+(B) Shemansky & Broadfoot, 1971
!      B38 (0,1) (4278) fr. N2+(B) Shemansky & Broadfoot, 1971
!      B39 (0,0) (3371) fr. N2(C)  Conway, 1983; Benesch et al, 1966
!      B40 (0,9) (3352) fr. N2(A)  Cartwright, 1978; Shemansky, 1969
!      B41 O+(2Po) fr. O+(2Pe)     Kirby et al, 1979
!      B42 O+(2Do) fr. O+(2Pe)     Kirby et al, 1979
!      B43 N2(C) bound fraction    ?
!      B44 7990 fr. O(3s'3D)       appx. fr. Hecht, p.c.
!      B45                         not currently in use
!      B46 N 1493 fr. N2+hv DI     guess
!      B47 N 1493 fr. N2+e* DI     guess, cf. Mumma and Zipf (1973), Meier (1991)
!      B48 N2(a) from (a,a',w)     estimate from comparison with Ajello & Shemansky, GUVI data, etc.
!      B49 7774, 1356 fr. O-+O+    Melendez, 1999; Qin, 2015 (cf. Tinsley 1973; Julienne, 1974)
"""
classic_bcoeff = [
    0.07, 1.20, 0.76, 1.85, 1.00, 0.10, 0.50, 0.81, 0.20, 0.32,
    0.48, 0.10, 0.10, 0.10, 0.16, 0.50, 0.30, 0.19, 0.00, 0.10,
    0.43, 0.51, 0.10, 0.60, 0.54, 0.44, 0.80, 0.20, 1.00, 0.33,
    0.33, 0.34, 0.21, 0.20, 0.10, 0.11, 0.65, 0.20, 0.24, 2e-2,
    0.18e2, 0.72e2, 0.75e2, 0.1e2, 0.5e-2, 0.5e-2, 0.2e-2, 0.7e2, 
    0.54e2, 0.5e-3
]
CLASSIC_BCOEFF = np.array(
    classic_bcoeff + [0.0]*(COEFF_LEN_NR - len(classic_bcoeff)),
    dtype=np.float32, order='F'
)
assert (len(CLASSIC_BCOEFF) == COEFF_LEN_NR)

modglow_acoeff = [
    6.60e-6, 6.44e-3, 2.15e-3, 0.0750, 1.260,
    8.89e-5, 4.70e-2, 0.1740, .660, 0.3520,
    6.50e-3, 3.45e-3
]
MODGLOW_ACOEFF = np.array(
    modglow_acoeff + [0.0]*(COEFF_LEN_NR - len(modglow_acoeff)),
    dtype=np.float32, order='F'
)
assert (len(MODGLOW_ACOEFF) == COEFF_LEN_NR)

modglow_bcoeff = [
    0.002, 1.389, 0.76, 1.12, 1.0, 0.10, 0.10, 0.737, 0.100, 0.196,
    0.283, 0.033, 0.065, 0.094, 0.10, 0.288, 0.36, 0.36, 0.00, 0.10,
    0.45, 0.50, 0.00, 0.60, 0.54, 0.44, 0.80, 0.20, 0.10, 0.33,
    0.33, 0.34, 0.15, 0.10, 0.10, 0.11, 0.65, 0.2, 0.24, 0.02,
    0.18, 0.72, 0.75, 0.10, .89, 0.31, 0.64, 0.36, 0.9, 0.05,
    0.15, 0.10, 0.609, 0.041
]

MODGLOW_BCOEFF = np.array(
    modglow_bcoeff + [0.0]*(COEFF_LEN_NR - len(modglow_bcoeff)),
    dtype=np.float32, order='F'
)
assert (len(MODGLOW_BCOEFF) == COEFF_LEN_NR)
