import numpy as np
import pytest
from glowpython2 import NrlMsis00
from glowpython2.utils import alt_grid
from msis21py import NrlMsis21

SPECIES = ["O", "O2", "N2", "NO"]


@pytest.fixture(scope="module")
def ds00(sample_time, sample_lat, sample_lon):
    return NrlMsis00().evaluate(time=sample_time, lat=sample_lat, lon=sample_lon, alt=alt_grid())


@pytest.fixture(scope="module")
def ds21(sample_time, sample_lat, sample_lon):
    return NrlMsis21().evaluate(time=sample_time, lat=sample_lat, lon=sample_lon, alt=alt_grid())


# --- output structure ---

@pytest.mark.parametrize("ds_fixture", ["ds00", "ds21"])
def test_has_alt_coord(ds_fixture, request):
    ds = request.getfixturevalue(ds_fixture)
    assert "alt_km" in ds.coords


@pytest.mark.parametrize("ds_fixture,species", [
    (f, s) for f in ["ds00", "ds21"] for s in SPECIES
])
def test_has_species(ds_fixture, species, request):
    ds = request.getfixturevalue(ds_fixture)
    assert species in ds, f"{ds_fixture} missing {species}"


@pytest.mark.parametrize("ds_fixture", ["ds00", "ds21"])
def test_has_temperature(ds_fixture, request):
    ds = request.getfixturevalue(ds_fixture)
    assert "Tn" in ds


# --- physical sanity ---

@pytest.mark.parametrize("ds_fixture", ["ds00", "ds21"])
def test_densities_nonnegative(ds_fixture, request):
    ds = request.getfixturevalue(ds_fixture)
    for spec in SPECIES:
        assert np.all(ds[spec].values >= 0), f"{ds_fixture} {spec} has negative values"


@pytest.mark.parametrize("ds_fixture", ["ds00", "ds21"])
def test_temperature_positive(ds_fixture, request):
    ds = request.getfixturevalue(ds_fixture)
    assert np.all(ds["Tn"].values > 0)


@pytest.mark.parametrize("ds_fixture", ["ds00", "ds21"])
def test_n2_dominates_low_alt(ds_fixture, request):
    ds = request.getfixturevalue(ds_fixture)
    low = ds.sel(alt_km=slice(None, 150))
    assert float(low["N2"].mean()) > float(low["O"].mean())


@pytest.mark.parametrize("ds_fixture", ["ds00", "ds21"])
def test_alt_range(ds_fixture, request):
    ds = request.getfixturevalue(ds_fixture)
    alt = ds.coords["alt_km"].values
    assert alt.min() >= 50
    assert alt.max() <= 1100


@pytest.mark.parametrize("ds_fixture", ["ds00", "ds21"])
def test_n2_density_magnitude(ds_fixture, request):
    # N2 at 100 km is well-constrained; should be ~10^12-10^14 cm^-3
    ds = request.getfixturevalue(ds_fixture)
    n2 = float(ds["N2"].sel(alt_km=100, method="nearest"))
    assert 1e11 < n2 < 1e14


@pytest.mark.parametrize("ds_fixture", ["ds00", "ds21"])
def test_thermosphere_warmer_than_mesosphere(ds_fixture, request):
    # Thermospheric heating means Tn at 400 km > Tn at 80 km
    ds = request.getfixturevalue(ds_fixture)
    tn_meso = float(ds["Tn"].sel(alt_km=80, method="nearest"))
    tn_thermo = float(ds["Tn"].sel(alt_km=400, method="nearest"))
    assert tn_thermo > tn_meso


@pytest.mark.parametrize("ds_fixture", ["ds00", "ds21"])
def test_o_dominates_above_200km(ds_fixture, request):
    # Atomic O should be the dominant neutral above ~200 km
    ds = request.getfixturevalue(ds_fixture)
    high = ds.sel(alt_km=slice(200, None))
    assert float(high["O"].mean()) > float(high["N2"].mean())
