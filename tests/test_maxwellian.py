import numpy as np
import pytest
import glowpython2 as glow


Q = 1          # erg cm-2 s-1
Echar = 100e3  # eV


@pytest.fixture(scope="module")
def iono_precip(sample_time, sample_lat, sample_lon):
    return glow.maxwellian(sample_time, sample_lat, sample_lon, Nbins=100, Q=Q, Echar=Echar)


@pytest.fixture(scope="module")
def iono_no_precip(sample_time, sample_lat, sample_lon):
    return glow.no_precipitation(sample_time, sample_lat, sample_lon, Nbins=100)


# --- output structure ---

def test_has_alt_coord(iono_precip):
    assert "alt_km" in iono_precip.coords


def test_has_energy_coord(iono_precip):
    assert "energy" in iono_precip.coords


def test_has_ver(iono_precip):
    assert "ver" in iono_precip
    assert "wavelength" in iono_precip["ver"].dims


def test_has_precip(iono_precip):
    assert "precip" in iono_precip


# --- physical sanity ---

def test_temperatures_positive(iono_precip):
    for var in ("Te", "Ti", "Tn"):
        assert np.all(iono_precip[var].values > 0)


def test_ver_valid_range(iono_precip):
    ver = iono_precip["ver"].sel(alt_km=slice(75, None)).values
    assert np.all(ver >= 0)


def test_ver_fill_below_boundary(iono_precip):
    # Below ~72 km VER should be NaN or numerically negligible (< 1e-10 cm-3 s-1)
    ver_low = iono_precip["ver"].sel(alt_km=slice(None, 72)).values
    assert np.all(np.isnan(ver_low) | (np.abs(ver_low) < 1e-10))


def test_precip_flux_nonzero(iono_precip):
    assert np.any(iono_precip["precip"].values > 0)


def test_precip_increases_ver(iono_precip, iono_no_precip):
    total_precip = float(iono_precip["ver"].sum())
    total_no_precip = float(iono_no_precip["ver"].sum())
    assert total_precip > total_no_precip


def test_auroral_4278_has_signal(iono_precip):
    # N2+ first-negative (4278 A) is a direct auroral tracer — must be clearly non-zero
    ver = iono_precip["ver"].sel(wavelength="4278", alt_km=slice(80, None))
    assert float(ver.max()) > 1.0


def test_5577_brighter_with_precip(iono_precip, iono_no_precip):
    # OI 5577 A should be at least 10x brighter with precipitation than without
    v_precip = float(iono_precip["ver"].sel(wavelength="5577").max())
    v_no_precip = float(iono_no_precip["ver"].sel(wavelength="5577").max())
    assert v_precip > v_no_precip * 10
