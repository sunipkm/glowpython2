import numpy as np

from glowpython2.utils import interpolate_nan
from .base import maxwellian, no_precipitation
from .plots import Plot
from datetime import datetime
from matplotlib.pyplot import show
from dateutil.parser import parse
import sys
import argparse
import json


def Maxwellian():
    parser = argparse.ArgumentParser(
        prog='Maxwellian Precipitation',
        description='Evaluate upper atmosphere VER for a given time and location.',
        epilog='Uses MSIS00 for atmosphere, IRI90 for ionosphere, GLOW for emissions.')

    parser.add_argument('--time', type=str, help='Time in ISO 8601 format.', default='2015-12-13T10:00:00', required=False)
    parser.add_argument('--glat', type=float, help='Geographic latitude in degrees.', default=42.6, required=False)
    parser.add_argument('--glon', type=float, help='Geographic longitude in degrees.', default=-71.2, required=False)
    parser.add_argument('--Q', type=float, help='Flux in erg cm^-2 s^-1.', default=1, required=False)
    parser.add_argument('--Echar', type=float, help='Characteristic energy in eV.', default=100e3, required=False)
    parser.add_argument('--Nbins', type=int, help='Number of energy bins.', default=250, required=False)
    parser.add_argument('--tec', type=float, help='Total electron content in TECU.', default=None, required=False)
    parser.add_argument('--newmodel', action='store_true', help='Use MSIS 2.1 and IRI-2020.', required=False)

    args = parser.parse_args()

    time = parse(args.time)
    glat = args.glat
    glon = args.glon
    # %% flux [erg cm-2 s-1 == mW m-2 s-1]
    Q = args.Q
    # %% characteristic energy [eV]
    Echar = args.Echar
    # %% Number of energy bins
    Nbins = args.Nbins
    tec = args.tec
    magmodel = 'IGRF14' if args.newmodel else 'POGO68'
    version = 'MSIS21_IRI20' if args.newmodel else 'MSIS00_IRI90'
    lstyles = ['-', '--']
    plot = Plot()
    iono = no_precipitation(time, glat, glon, Nbins, tec=tec, version=version, magmodel=magmodel)  # type: ignore
    ne = interpolate_nan(iono["NeIn"].values, inplace=False)
    hmf = iono.attrs.get('hmf2')
    if hmf is not None:
        hmf = json.loads(hmf)['value']
    print(f'[{version}] TEC: {np.trapezoid(ne, iono.alt_km.values*1e5)*1e-12:.2f} TECU, hmf2: {hmf:.1f} km')

    plot.precip(iono["precip"])
    plot.density(iono)

    plot.temperature(iono)

    plot.ver(iono)

    show()


def NoPrecipitation():
    parser = argparse.ArgumentParser(
        prog='No Electron Precipitation',
        description='Evaluate upper atmosphere VER for a given time and location.',
        epilog='Uses MSIS00 for atmosphere, IRI90 for ionosphere, GLOW for emissions.')

    parser.add_argument('--time', type=str, help='Time in ISO 8601 format.', default='2015-12-13T10:00:00', required=False)
    parser.add_argument('--glat', type=float, help='Geographic latitude in degrees.', default=42.6, required=False)
    parser.add_argument('--glon', type=float, help='Geographic longitude in degrees.', default=-71.2, required=False)
    parser.add_argument('--Nbins', type=int, help='Number of energy bins.', default=250, required=False)
    parser.add_argument('--tec', type=float, help='Total electron content in TECU.', default=None, required=False)
    parser.add_argument('--newmodel', action='store_true', help='Use MSIS 2.1 and IRI-2020.', required=False)

    args = parser.parse_args()

    time = parse(args.time)
    glat = args.glat
    glon = args.glon
    Nbins = args.Nbins
    tec = args.tec
    magmodel = 'IGRF14' if args.newmodel else 'POGO68'
    version = 'MSIS21_IRI20' if args.newmodel else 'MSIS00_IRI90'
    plot = Plot()
    iono = no_precipitation(time, glat, glon, Nbins, tec=tec, version=version, magmodel=magmodel)  # type: ignore
    ne = interpolate_nan(iono["NeIn"].values, inplace=False)
    hmf = iono.attrs.get('hmf2')
    if hmf is not None:
        hmf = json.loads(hmf)['value']
    print(f'[{version}] TEC: {np.trapezoid(ne, iono.alt_km.values*1e5)*1e-12:.2f} TECU, hmf2: {hmf:.1f} km')

    plot.density(iono)

    plot.temperature(iono)

    plot.ver(iono)

    show()
