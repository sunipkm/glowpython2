from __future__ import annotations
import os
from pathlib import Path
from time import perf_counter_ns
from datetime import UTC, datetime, timedelta
from typing import Any, Optional, Dict, Tuple, Sequence, SupportsFloat as Numeric
from dataclasses import dataclass
import numpy as np
from xarray import Dataset
from msis21py import Settings as Msis00Settings, ComputedSettings as Msis00ComputedSettings
import geomagdata as gi

from .version import __version__
from .utils import Singleton, glowdate
from .glowatmo import msis00_init, msis00_eval  # type: ignore

DIRNAME = Path(os.path.dirname(__file__)) / 'data'
DATADIR = DIRNAME.resolve()


@dataclass
class Attribute:
    value: Any
    units: Optional[str]
    long_name: str
    description: Optional[str] = None

    def to_json(self) -> str:
        from json import dumps
        from dataclasses import asdict
        return dumps(asdict(self))


class NrlMsis00(Singleton):
    def _init(self):
        self.change_settings(Msis00Settings())
        self._benchmark = False
        self._call = 0
        self._setup = 0
        self._fortran = 0
        self._ds_build = 0
        self._ds_attrib = 0
        self._ds_settings = 0
        self._total = 0

    def change_settings(self, settings: Msis00Settings):
        """Change the current settings.

        Args:
            settings (Settings): New settings to use.
        """
        self.settings = settings
        cs = Msis00ComputedSettings.from_settings(self.settings)
        msis00_init(
            str(DATADIR), sw=cs.switch_legacy
        )

    @property
    def benchmark(self) -> bool:
        return self._benchmark

    @benchmark.setter
    def benchmark(self, value: bool):
        if value != self._benchmark:
            self._call = 0
            self._setup = 0
            self._fortran = 0
            self._ds_build = 0
            self._ds_attrib = 0
            self._ds_settings = 0
            self._total = 0
        self._benchmark = value

    def get_benchmark(self) -> Optional[Dict[str, timedelta]]:
        """Get benchmark data.

        Returns:
            Optional[Dict[str, timedelta]]: Metric and measured time.
        """
        if not self._benchmark or self._call == 0:
            return None
        return {
            'setup': timedelta(milliseconds=self._setup / self._call),
            'fortran': timedelta(milliseconds=self._fortran / self._call),
            'ds_build': timedelta(milliseconds=self._ds_build / self._call),
            'ds_attrib': timedelta(milliseconds=self._ds_attrib / self._call),
            'ds_settings': timedelta(milliseconds=self._ds_settings / self._call),
            'total': timedelta(milliseconds=self._total / self._call),
        }

    def _msiscall(self, lat: Numeric, lon: Numeric, alt: np.ndarray, ydate: int, ut: Numeric, f107: Tuple[Numeric, Numeric], ap: np.ndarray) -> Dataset:
        start = perf_counter_ns()
        fort_densities = np.zeros((10, len(alt)), dtype=np.float32, order='F')
        fort_temps = np.full(len(alt), np.nan, dtype=np.float32, order='F')
        alt = alt.astype(np.float32, order='F')
        setup = perf_counter_ns()
        exot = msis00_eval(ydate, ut, alt, lat, lon,
                           f107[0], f107[1], ap, fort_densities, fort_temps)
        fortran = perf_counter_ns()
        ds = Dataset()
        ds.coords['alt_km'] = (
            ('alt_km',), alt, {'units': 'km', 'long_name': 'Altitude'})
        densities = ['He', 'O', 'N2', 'O2',
                     'Ar', 'H', 'N', 'Anomalous O', 'NO']
        descriptions = ['Helium', 'Atomic Oxygen', 'Molecular Nitrogen',
                        'Molecular Oxygen', 'Argon', 'Hydrogen', 'Nitrogen',
                        'Anomalous Oxygen', 'Nitric Oxide']
        density_idx = list(range(5)) + list(range(6, 10))
        for idx, name, desc in zip(density_idx, densities, descriptions):
            ds[name] = (('alt_km',), fort_densities[idx],
                        {'units': 'cm^-3', 'long_name': f'{desc} Density'})
        ds['mden'] = (('alt_km',), fort_densities[5],
                      {'units': 'g/cm^3', 'long_name': 'Mass Density'})
        ds['Tn'] = (('alt_km',), fort_temps,
                    {'units': 'K', 'long_name': 'Neutral Temperature'})
        ds_build = perf_counter_ns()
        ds.attrs['attributes'] = 'Stored as JSON strings'
        ds.attrs['description'] = 'NRLMSISE-00 model output'
        ds.attrs['version'] = __version__
        ds.attrs['lat'] = Attribute(
            value=float(lat),
            units='degrees',
            long_name='Geographic Latitude',
        ).to_json()
        ds.attrs['lon'] = Attribute(
            value=float(lon),
            units='degrees',
            long_name='Geographic Longitude',
        ).to_json()
        ds.attrs['f107p'] = Attribute(
            value=float(f107[1]),
            units='sfu',
            long_name='Previous Day F10.7 Solar Flux',
        ).to_json()
        ds.attrs['f107a'] = Attribute(
            value=float(f107[0]),
            units='sfu',
            long_name='81-day Average F10.7 Solar Flux',
        ).to_json()
        ds.attrs['Ap'] = Attribute(
            value=float(ap[0]),
            units=None,
            long_name='Daily Ap Geomagnetic Index',
        ).to_json()
        ds.attrs['exot'] = Attribute(
            value=float(exot),
            units='K',
            long_name='Exosphere Temperature',
        ).to_json()
        ds_attrib = perf_counter_ns()
        ds.attrs['settings'] = self.settings.to_json()
        ds_settings = perf_counter_ns()
        if self._benchmark:
            self._call += 1
            self._setup += (setup - start)*1e-6
            self._fortran += (fortran - setup)*1e-6
            self._ds_build += (ds_build - fortran)*1e-6
            self._ds_attrib += (ds_attrib - ds_build)*1e-6
            self._ds_settings += (ds_settings - ds_attrib)*1e-6
            self._total += (ds_settings - start)*1e-6
        return ds

    def evaluate(
        self,
        time: datetime,
        lat: Numeric, lon: Numeric, alt: np.ndarray,
        *,
        geomag_params: Optional[dict[str, Numeric] | Sequence[Numeric]] = None,
    ) -> Dataset:
        """Evaluate the NRLMSISE-00 model.

        Args:
            time (datetime): Datetime object. 
            lat (Numeric): Geographic latitude.
            lon (Numeric): Geographic longitude.
            alt (np.ndarray): Altitude in kilometers.

        Returns:
            Dataset: Computed dataset.
        """
        if time.tzinfo is not None:
            time = time.astimezone(UTC)
        ydate, utsec = glowdate(time)
        if geomag_params is None:
            ip = gi.get_indices([time - timedelta(days=1), time],
                                81, tzaware=True)  # type: ignore
            f107a = float(ip["f107s"].iloc[1])
            f107 = float(ip['f107'].iloc[1])
            f107p = float(ip['f107'].iloc[0])
            ap = float(ip["Ap"].iloc[1])
        elif isinstance(geomag_params, dict):
            f107a = float(geomag_params['f107a'])
            f107 = float(geomag_params['f107'])
            f107p = float(geomag_params['f107p'])
            ap = float(geomag_params['Ap'])
        elif isinstance(geomag_params, Sequence):
            f107a = float(geomag_params[0])
            f107 = float(geomag_params[1])
            f107p = float(geomag_params[2])
            ap = float(geomag_params[3])
        else:
            raise RuntimeError('Invalid type %s for geomag params %s' % (
                str(type(geomag_params), str(geomag_params))))
        lon = lon % 360  # ensure lon is in 0-360 range
        ds = self.lowlevel(
            lat, lon, alt, ydate, utsec,
            f107a, f107p, np.array([ap]*7, dtype=np.float32, order='F')
        )
        ds.attrs['date'] = time.isoformat()
        return ds

    def lowlevel(self, lat: Numeric, lon: Numeric, alt: np.ndarray, ydate: int, ut: Numeric, f107a: Numeric, f107p: Numeric, ap: np.ndarray) -> Dataset:
        """Low level call to evaluate NRLMSISE-00 model.
        Bypasses date and time calculations.

        Args:
            lat (Numeric): Geographic latitude
            lon (Numeric): Geographic longitude
            alt (np.ndarray): Altitude in kilometers
            ydate (int): YYYYDDD date format
            ut (Numeric): Universal time in seconds
            f107a (Numeric): 81-day average F10.7 solar flux
            f107p (Numeric): Previous day F10.7 solar flux
            ap (np.ndarray): Array of 7 daily Ap geomagnetic indices. Length must be 7. The first element is for the current 3-hour period. The remaining 6 elements are only used in 'Storm' mode.

        Returns:
            Dataset: Computed dataset.
        """
        ds = self._msiscall(
            lat, lon, alt, ydate, ut,
            (f107a, f107p),
            ap
        )
        return ds

class Iri90(Singleton):
    def _init(self):
        self._benchmark = False
        self._call = 0
        self._setup = 0
        self._fortran = 0
        self._ds_build = 0
        self._ds_attrib = 0
        self._ds_settings = 0
        self._total = 0

# %%
def test():
    import matplotlib.pyplot as plt
    from pprint import pprint
    from msis21py import NrlMsis21, alt_grid
    msis = NrlMsis21()
    date = datetime(2022, 3, 21, 12, 0, 0, tzinfo=UTC)
    ds1 = msis.evaluate(
        date,
        40.0, 105.0,
        alt_grid(),
    )
    # iri.benchmark = True
    # for idx in range(int(1e4)):
    #     if idx > 0 and idx % 1000 == 0:
    #         print(f"{idx}/{int(1e4)} calculations done")
    #     _ = iri.calculate(
    #         date,
    #         40.0, 105.0,
    #         alt_grid(),
    #         set
    #     )
    # pprint(iri.get_benchmark())
    pprint(ds1)
    ds2 = msis.evaluate(
        datetime(2022, 3, 21, 0, 0, 0, tzinfo=UTC),
        40.0, 105.0,
        alt_grid(),
    )
    fig, ax = plt.subplots(1, 2, sharey=True)
    ds1.Tn.plot(ax=ax[0], y='alt_km')
    ds2.Tn.plot(ax=ax[1], y='alt_km')
    ax[0].set_xscale('log')
    ax[1].set_xscale('log')
    ax[0].set_ylabel('Altitude (km)')
    ax[0].set_title('Neutral temperature at 12:00 UTC')
    ax[1].set_title('Neutral temperature at 00:00 UTC')
    plt.show()
    # ds1.to_netcdf('iri_output.nc')
