from __future__ import annotations
from typing import SupportsFloat as Numeric
from datetime import datetime
from pathlib import Path
import os

import numpy as np

from .utils import Singleton
from .glowigrf import glowigrf # type: ignore

DATADIR = Path(os.path.dirname(__file__)) / 'data'
DATADIR = DATADIR.resolve()


class Igrf(Singleton):
    def _init(self, logfile: str = ''):
        # Initialize IGRF model parameters
        glowigrf.igrf_init(str(DATADIR), logfile)
    
    def dipangle(self, time: Numeric, lat: Numeric, lon: Numeric, alt: Numeric | np.ndarray) -> float | np.ndarray:
        """Calculate magnetic dip angle in degrees.

        Args:
            time (Numeric): Decimal year.
            lat (Numeric): Geographic latitude in degrees.
            lon (Numeric): Geographic longitude in degrees.
            alt (Numeric | np.ndarray): Altitude in km.
        Returns:
            Numeric | np.ndarray: Magnetic dip angle in degrees.
        """
        if isinstance(alt, np.ndarray):
            dip = np.full(len(alt), 0, dtype=np.float32, order='F')
            glowigrf.diparray_igrf(float(time), float(lat), float(lon), alt.astype(np.float32, order='F'), dip)
            dip = dip.astype(float)
        else:
            dip = float(glowigrf.dipangle_igrf(float(time), float(lat), float(lon), float(alt)))
        return dip

    def fieldstrength(self, time: Numeric, lat: Numeric, lon: Numeric, alt: np.ndarray) -> np.ndarray:
        """Calculate magnetic field strength in T.

        Args:
            time (Numeric): Decimal year.
            lat (Numeric): Geographic latitude in degrees.
            lon (Numeric): Geographic longitude in degrees.
            alt (Numeric | np.ndarray): Altitude in km.
        Returns:
            Numeric | np.ndarray: Magnetic field strength in T.
        """
        bfield = np.full(len(alt), 0, dtype=np.float32, order='F')
        glowigrf.bmagarray_igrf(float(time), float(lat), float(lon), alt.astype(np.float32, order='F'), bfield)
        bfield = bfield.astype(float)
        return bfield