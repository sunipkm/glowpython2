from __future__ import annotations
from dataclasses import dataclass, asdict
from datetime import timedelta, datetime, UTC
import os
from pathlib import Path
from time import perf_counter_ns
from typing import Any, Dict, Iterable, Literal, Optional, SupportsFloat as Numeric, Tuple
import numpy as np
from xarray import Dataset
import geomagdata as gi

from .utils import glowdate, Singleton
from .glowatmo import iri90_eval  # type: ignore

DIRNAME = Path(os.path.dirname(__file__)) / 'data'
DATADIR = DIRNAME.resolve()

B0Model = Literal['Table', 'Gulyeava-1987']  # [3], [T, F]
FoF2Model = Literal['CCIR', 'URSI']  # [4] [T, F]
IonComposition = Literal['Standard', 'Danilov-Yaichnikov-1985']  # [5] [T, F]
IriTopside = Literal['IRI-90', 'IRI-79']  # [6] [T, F]
NeModel = Literal['Model', 'Lay-Functions']  # [10] [T, F]


@dataclass
class Settings:
    """Settings for IRI-90 model."""
    b0_model: B0Model = 'Table'
    fof2_model: FoF2Model = 'URSI'
    ion_composition: IonComposition = 'Standard'
    iri_topside: IriTopside = 'IRI-90'
    ne_model: NeModel = 'Model'
    nmf2: Optional[Numeric] = None  # [7], None for model, oarr[0] for custom
    hmF2: Optional[Numeric] = None  # [8], None for model, oarr[1] for custom

    def to_json(self) -> str:
        """Convert settings to JSON string."""
        import json
        return json.dumps(asdict(self))


def literal_to_idx(value: str) -> Iterable[Tuple[int, bool]]:
    if value == 'Table':
        return (3, True),
    elif value == 'Gulyeava-1987':
        return (3, False),
    elif value == 'CCIR':
        return (4, True),
    elif value == 'URSI':
        return (4, False),
    elif value == 'Standard':
        return (5, True),
    elif value == 'Danilov-Yaichnikov-1985':
        return (5, False),
    elif value == 'IRI-90':
        return (6, True),
    elif value == 'IRI-79':
        return (6, False),
    elif value == 'Model':
        return (10, True),
    elif value == 'Lay-Functions':
        return (10, False),
    else:
        raise ValueError(f"Unknown literal value: {value}")


@dataclass
class ComputedSettings:
    """Computed settings for IRI-90 model."""
    jf: np.ndarray
    oarr: np.ndarray

    @staticmethod
    def from_settings(settings: Settings) -> ComputedSettings:
        """Compute the settings from user settings."""
        jf = np.full(12, True, dtype=bool)
        oarr = np.full(30, -1.0, dtype=float)

        for idx, val in literal_to_idx(settings.b0_model):
            jf[idx] = val
        for idx, val in literal_to_idx(settings.fof2_model):
            jf[idx] = val
        for idx, val in literal_to_idx(settings.ion_composition):
            jf[idx] = val
        for idx, val in literal_to_idx(settings.iri_topside):
            jf[idx] = val
        for idx, val in literal_to_idx(settings.ne_model):
            jf[idx] = val
        if settings.nmf2 is not None:
            oarr[0] = float(settings.nmf2)
            jf[7] = False
        if settings.hmF2 is not None:
            oarr[1] = float(settings.hmF2)
            jf[8] = False

        return ComputedSettings(jf=jf, oarr=oarr.astype(np.float32))


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


class Iri90(Singleton):
    def _init(self):
        self.settings: Settings = Settings()
        self._benchmark = False
        self._call = 0
        self._setup = 0
        self._fortran = 0
        self._ds_build = 0
        self._ds_attrib = 0
        self._ds_settings = 0
        self._total = 0

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

    def _iricall(self, lat: Numeric, lon: Numeric, alt: np.ndarray, day: int, ut: Numeric, f107a: Numeric, settings: ComputedSettings) -> Dataset:
        start = perf_counter_ns()
        outf = np.full((11, len(alt)), 0.0, dtype=np.float32, order='F')
        alt = alt.astype(np.float32, order='F')
        setup = perf_counter_ns()
        iri90_eval(
            settings.jf, 0, lat, lon, day, ut, f107a,
            alt, str(DATADIR), outf, settings.oarr
        )
        fortran = perf_counter_ns()
        ds = Dataset()
        ds.coords['alt_km'] = (
            ('alt_km',), alt, {'units': 'km', 'long_name': 'Altitude'})
        densities = ['Ne', 'O+', 'H+', 'He+', 'O2+', 'NO+', 'Cluster', 'N+']
        den_names = [
            'Electron', 'Oxygen Ion', 'Hydrogen Ion', 'Helium Ion', 'Oxygen Molecular Ion',
            'Nitric Oxide Ion', 'Cluster Ion', 'Nitrogen Ion'
        ]
        density_idx = [0, 4, 5, 6, 7, 8, 9, 10]
        temperatures = ['Te', 'Ti']
        temperature_names = ['Electron Temperature', 'Ion Temperature']
        temperature_idx = [3, 2]
        for idx, name, desc in zip(density_idx, densities, den_names):
            ds[name] = (('alt_km',), outf[idx]*1e-6,
                        {'units': 'cm^-3', 'long_name': f'{desc} Density'})
        for idx, name, desc in zip(temperature_idx, temperatures, temperature_names):
            ds[name] = (('alt_km',), outf[idx], {
                        'units': 'K', 'long_name': f'{desc} Temperature'})
        ds_build = perf_counter_ns()
        ds.attrs['attributes'] = 'Stored as JSON strings'
        ds.attrs['description'] = 'IRI-90 model output'
        oarr = settings.oarr.copy().astype(float)
        attrs = {}
        attrs['nmF2'] = Attribute(
            oarr[0]*1e-6, 'cm^-3', 'F2 Peak Density', 'F2 layer peak electron density')
        attrs['hmF2'] = Attribute(
            oarr[1], 'km', 'F2 Peak Height', 'F2 layer peak height')
        attrs['nmF1'] = Attribute(
            oarr[2]*1e-6, 'cm^-3', 'F1 Peak Density', 'F1 layer peak electron density')
        attrs['hmF1'] = Attribute(
            oarr[3], 'km', 'F1 Peak Height', 'F1 layer peak height')
        attrs['nmE'] = Attribute(oarr[4]*1e-6, 'cm^-3',
                                 'E Layer Peak Density', 'E layer peak electron density')
        attrs['hmE'] = Attribute(
            oarr[5], 'km', 'E Layer Peak Height', 'E layer peak height')
        attrs['nmD'] = Attribute(
            oarr[6]*1e-6, 'cm^-3', 'D Layer inflection point density')
        attrs['hmD'] = Attribute(
            oarr[7], 'km', 'D-region inflection point')
        attrs['hhalf'] = Attribute(
            oarr[8], 'km', 'Half Height', 'Height used by Gulyaeva B0 model')
        attrs['B0'] = Attribute(oarr[9], 'km',
                                'B0', 'Bottomside thickness parameter')
        attrs['valley_base'] = Attribute(
            oarr[10]*1e-6, 'cm^-3', 'Density at E-valley base')
        attrs['valley_top'] = Attribute(
            oarr[11], 'km', 'Height of E-valley top')
        attrs['Te-Peak'] = Attribute(oarr[12], 'K', 'Te Peak')
        attrs['hTe-Peak'] = Attribute(oarr[13],
                                      'km', 'hTe Peak', 'Peak Te altitude')
        attrs['Te-MOD(300km)'] = Attribute(oarr[14], 'K',
                                           'Te MOD(300km)', 'Electron temperature at 300 km altitude')
        attrs['Te-MOD(400km)'] = Attribute(oarr[15], 'K',
                                           'Te MOD(400km)', 'Electron temperature at 400 km altitude')
        attrs['Te-MOD(600km)'] = Attribute(oarr[16], 'K',
                                           'Te MOD(600km)', 'Electron temperature at 600 km altitude')
        attrs['Te-MOD(1400km)'] = Attribute(oarr[17], 'K',
                                            'Te MOD(1400km)', 'Electron temperature at 1400 km altitude')
        attrs['Te-MOD(3000km)'] = Attribute(oarr[18], 'K',
                                            'Te MOD(3000km)', 'Electron temperature at 3000 km altitude')
        attrs['Te-MOD(120km)'] = Attribute(oarr[19], 'K', 'Te MOD(120km)',
                                           'Electron temperature at 120 km altitude, Te = Ti = Tn')
        attrs['Ti-MOD(430km)'] = Attribute(oarr[20], 'K',
                                           'Ti MOD(430km)', 'Ion temperature at 430 km altitude')
        attrs['Ti-Te-Eq'] = Attribute(oarr[21], 'km', 'Ti-Te-Eq',
                                      'Height at which ion and electron temperatures are at equilibrium')
        attrs['sza'] = Attribute(oarr[22], 'degrees', 'Solar Zenith Angle',
                                 'Solar zenith angle at the specified location and time')
        attrs['sun_dec'] = Attribute(
            oarr[23], 'degrees', 'Solar Declination', 'Solar declination angle at the specified time')
        attrs['dip'] = Attribute(oarr[24], 'degrees', 'Magnetic Dip Angle',
                                 'Magnetic dip angle at the specified location')
        attrs['dip-lat'] = Attribute(oarr[25], 'degrees',
                                     'Magnetic Dip Latitude', 'Magnetic dip latitude')
        attrs['dip-lat-mod'] = Attribute(oarr[26], 'degrees',
                                         'Magnetic Dip Latitude (Modified)', 'Modified magnetic dip latitude')

        for key, attr in attrs.items():
            if isinstance(attr, Attribute):
                ds.attrs[key] = attr.to_json()
            elif attr is not None:
                ds.attrs[key] = str(attr)
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
        f107a: Optional[Numeric] = None,
        settings: Optional[Settings | ComputedSettings] = None,
        *,
        tzaware: bool = False
    ) -> Tuple[ComputedSettings, Dataset]:
        """Evaluate the IRI-90 model.

        Args:
            time (datetime): Datetime object. 
            lat (Numeric): Geographic latitude.
            lon (Numeric): Geographic longitude.
            alt (np.ndarray): Altitude in kilometers.
            f107a (Optional[Numeric]): 81-day average F10.7 solar flux. If None, it will be retrieved from geomagdata. Defaults to None.
            settings (Optional[Settings  |  ComputedSettings], optional): Settings to use. Defaults to None.
            tzaware (bool, optional): If time is time zone aware. If true, `time` is recast to 'UTC' using `time.astimezone(pytz.utc)`. Defaults to False.

        Returns:
            Tuple[ComputedSettings, Dataset]: Computed settings and dataset. Passing in Settings will return ComputedSettings. For subsequent calls, pass in the returned ComputedSettings to avoid recomputation.
        """
        if tzaware:
            time = time.astimezone(UTC)
        year, utsec = glowdate(time)
        lon = lon % 360  # ensure lon is in 0-360 range
        if f107a is None:
            ret = gi.get_indices([time], smoothdays=81)  # type: ignore
            f107a_ = ret['f107s'].iloc[0]
        else:
            f107a_ = f107a
        (settings, ds) = self.lowlevel(
            lat, lon, alt, year, utsec, f107a_, settings)
        ds.attrs['date'] = time.isoformat()
        return (settings, ds)

    def lowlevel(self, lat: Numeric, lon: Numeric, alt: np.ndarray, day: int, ut: Numeric, f107a: Numeric, settings: Optional[Settings | ComputedSettings] = None) -> Tuple[ComputedSettings, Dataset]:
        """Low level call to evaluate IRI-2020 model.
        Bypasses date and time calculations.

        Args:
            lat (Numeric): Geographic latitude
            lon (Numeric): Geographic longitude
            alt (np.ndarray): Altitude in kilometers
            day (int): Day of the year (1-365/366)
            ut (Numeric): Universal time in seconds
            f107a (Numeric): 81-day average F10.7 solar flux
            settings (Optional[Settings  |  ComputedSettings], optional): Settings to use. Defaults to None.

        Raises:
            TypeError: If settings is not of type Settings or ComputedSettings.

        Returns:
            Tuple[ComputedSettings, Dataset]: Computed settings and dataset. Passing in Settings will return ComputedSettings. For subsequent calls, pass in the returned ComputedSettings to avoid recomputation.
        """
        if settings is None:
            settings = self.settings
        if isinstance(settings, Settings):
            self.settings = settings
            settings = ComputedSettings.from_settings(settings)
        if not isinstance(settings, ComputedSettings):
            raise TypeError(
                "settings must be of type Settings or ComputedSettings")
        ds = self._iricall(lat, lon, alt, day, ut, f107a, settings)
        return settings, ds
