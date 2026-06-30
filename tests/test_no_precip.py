import numpy as np
import pytest
import glowpython2 as glow


@pytest.fixture(scope="module")
def iono(sample_time, sample_lat, sample_lon):
    return glow.no_precipitation(sample_time, sample_lat, sample_lon, Nbins=100)


# --- output structure ---

def test_has_alt_coord(iono):
    assert "alt_km" in iono.coords


def test_alt_range(iono):
    alt = iono.coords["alt_km"].values
    assert alt.min() >= 50
    assert alt.max() <= 1100


def test_has_neutral_species(iono):
    for var in ("O", "O2", "N2", "NO"):
        assert var in iono, f"Missing neutral species: {var}"


def test_has_ion_species(iono):
    for var in ("O+", "O2+", "NO+", "NeIn"):
        assert var in iono, f"Missing ion species: {var}"


def test_has_temperatures(iono):
    for var in ("Te", "Ti", "Tn"):
        assert var in iono, f"Missing temperature: {var}"


def test_has_ver(iono):
    assert "ver" in iono
    assert "wavelength" in iono["ver"].dims


def test_has_precip(iono):
    assert "precip" in iono


# --- physical sanity ---

def test_neutral_densities_positive(iono):
    for var in ("O", "O2", "N2", "NO"):
        assert np.all(iono[var].values >= 0), f"{var} has negative values"


def test_temperatures_positive(iono):
    for var in ("Te", "Ti", "Tn"):
        assert np.all(iono[var].values > 0), f"{var} has non-positive values"


def test_electron_density_positive(iono):
    assert np.all(iono["NeIn"].values >= 0)


def test_ver_valid_range(iono):
    # VER must be non-negative above 75 km
    ver = iono["ver"].sel(alt_km=slice(75, None)).values
    assert np.all(ver >= 0)


def test_ver_fill_below_boundary(iono):
    # Below ~72 km GLOW outputs NaN or 0 for VER.
    # Must be one of those, not arbitrary negatives
    ver_low = iono["ver"].sel(alt_km=slice(None, 72)).values
    assert np.all(np.isnan(ver_low) | (ver_low == 0))


def test_no_precip_flux_is_zero(iono):
    assert np.all(iono["precip"].values == 0)


def test_n2_dominates_at_low_alt(iono):
    low = iono.sel(alt_km=slice(None, 150))
    assert float(low["N2"].mean()) > float(low["O"].mean())
