# %%
from __future__ import annotations
from pathlib import Path
from typing import Literal, Sequence, Tuple, Union, SupportsFloat as Numeric
import numpy as np
from .base import FluxSource, DATA_DIR
from .glowfort import cglow as cg # type: ignore

# %%
class Sflux:
    def __init__(self, kind: FluxSource | Tuple[Sequence[float], Sequence[float], Sequence[float]], data_dir: Path = DATA_DIR):
        self._kind = kind
        self._docalc = True
        if kind == 'Hinteregger':
            self._method = 'Hinteregger'
            ifile = data_dir / 'ssflux_hint.dat'
            self._data = np.loadtxt(ifile, skiprows=1)[::-1, :].T # reverse order
            self._fluxmask = self._data[0] < 251 & self._data[0] > 17 # 17-251 nm
            cg.wave1 = self._data[0]
            cg.wave2 = self._data[1]
        elif kind == 'EUVAC':
            self._method = 'EUVAC'
            ifile = data_dir / 'ssflux_euvac.dat'
            self._data = np.loadtxt(ifile, skiprows=1)[::-1, :].T # reverse order
            self._fluxmask = self._data[0] < 51 & self._data[0] > 1 # 1-51 nm
            cg.wave1 = self._data[0]
            cg.wave2 = self._data[1]
        elif isinstance(kind, tuple):
            self._method = 'custom'
            self._data = np.stack(kind, axis=0)
            if self._data.ndim != 2:
                raise ValueError(f"Custom flux data must be 2D, got {self._data.ndim}D.")
            if self._data.shape[1] != 3:
                raise ValueError(f"Custom flux data must have 3 columns, got {self._data.shape[1]}.")
            if self._data.shape[0] != 123:
                raise ValueError(f"Custom flux data must have 123 wavelength bins, got {self._data.shape[0]}.")
            if self._data[1, 0] > self._data[0, 0]:
                self._data = self._data[::-1, :] # reverse order
            cg.wave1 = self._data[:, 0]
            cg.wave2 = self._data[:, 1]
            self._docalc = False
        else:
            raise ValueError(f"Unknown flux source: {kind}. Must be 'Hinteregger', 'EUVAC', or a tuple of values.")
        
    def update(self, f107: Numeric, f107a: Numeric, xuvfac: int = 3) -> None:
        if self._docalc:
            if self._data.shape[1] == 5: # Hinteregger
                r1 = 0.0138*(f107a - 71.5) + 0.005*(f107 - f107a + 3.9) # + 1
                r2 = 0.59425*(f107a - 71.5) + 0.3811*(f107 - f107a + 3.9) # + 1
                sflux = self._data[2] + r1 * self._data[3] + r2 * self._data[4]
                np.clip(sflux, 0, None, out=sflux)
                if xuvfac > 1e-6:
                    sflux[self._fluxmask] *= xuvfac
            elif self._data.shape[1] == 4: # EUVAC
                p107 = (f107 + f107a)*0.5
                sflux = self._data[2] * (1 + self._data[3]*(p107 - 80))
                np.clip(sflux, 0.1*self._data[2], None, out=sflux)
                if xuvfac > 1e-6:
                    sflux[self._fluxmask] *= xuvfac
            cg.sflux = sflux
        else:
            # Set the cg.wave1, cg.wave2, cg.sflux to these values
            cg.sflux = self._data[:, 2]
# %%
