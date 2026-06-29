import numpy as np
import pytest
from glowpython2 import Iri90
from glowpython2.utils import alt_grid
from iri20py import Iri2020

SPECIES = ["Ne", "O+", "O2+", "NO+"]
TEMPERATURES = ["Te", "Ti"]


@pytest.fixture(scope="module")
def ds90(sample_time, sample_lat, sample_lon):
    _, ds = Iri90().evaluate(time=sample_time, lat=sample_lat, lon=sample_lon, alt=alt_grid())
    return ds


@pytest.fixture(scope="module")
def ds20(sample_time, sample_lat, sample_lon):
    _, ds = Iri2020().evaluate(time=sample_time, lat=sample_lat, lon=sample_lon, alt=alt_grid())
    return ds


# --- output structure ---

@pytest.mark.parametrize("ds_fixture", ["ds90", "ds20"])
def test_has_alt_coord(ds_fixture, request):
    ds = request.getfixturevalue(ds_fixture)
    assert "alt_km" in ds.coords


@pytest.mark.parametrize("ds_fixture,species", [
    (f, s) for f in ["ds90", "ds20"] for s in SPECIES
])
def test_has_species(ds_fixture, species, request):
    ds = request.getfixturevalue(ds_fixture)
    assert species in ds, f"{ds_fixture} missing {species}"


@pytest.mark.parametrize("ds_fixture,temp", [
    (f, t) for f in ["ds90", "ds20"] for t in TEMPERATURES
])
def test_has_temperature(ds_fixture, temp, request):
    ds = request.getfixturevalue(ds_fixture)
    assert temp in ds, f"{ds_fixture} missing {temp}"


# --- physical sanity ---

@pytest.mark.parametrize("ds_fixture", ["ds90", "ds20"])
def test_electron_density_valid_range(ds_fixture, request):
    # Ne must be positive in the F-region
    ds = request.getfixturevalue(ds_fixture)
    ne = ds["Ne"].sel(alt_km=slice(150, 500)).values
    assert np.all(ne > 0)


@pytest.mark.parametrize("ds_fixture", ["ds90", "ds20"])
def test_electron_density_fill_below_boundary(ds_fixture, request):
    # Below ~80 km IRI is undefined; must be exactly the -1e-6 fill
    ds = request.getfixturevalue(ds_fixture)
    ne_low = ds["Ne"].sel(alt_km=slice(None, 79)).values
    assert np.all(ne_low == pytest.approx(-1e-6, rel=1e-3))


@pytest.mark.parametrize("ds_fixture", ["ds90", "ds20"])
def test_temperatures_valid_range(ds_fixture, request):
    # Te and Ti must be physically reasonable above 150 km
    ds = request.getfixturevalue(ds_fixture)
    for temp in TEMPERATURES:
        vals = ds[temp].sel(alt_km=slice(150, 500)).values
        assert np.all(vals > 100), f"{ds_fixture} {temp} unreasonably cold above 150 km"


def test_temperatures_fill_below_boundary_iri90(ds90):
    # IRI-90 uses -1.0 as fill below 120 km; must be exactly that
    for temp in TEMPERATURES:
        vals = ds90[temp].sel(alt_km=slice(None, 119)).values
        assert np.all(vals == -1.0), f"ds90 {temp} unexpected values below IRI lower boundary"


@pytest.mark.parametrize("ds_fixture", ["ds90", "ds20"])
def test_alt_range(ds_fixture, request):
    ds = request.getfixturevalue(ds_fixture)
    alt = ds.coords["alt_km"].values
    assert alt.min() >= 50
    assert alt.max() <= 2000


@pytest.mark.parametrize("ds_fixture", ["ds90", "ds20"])
def test_te_ti_physically_reasonable(ds_fixture, request):
    # Te and Ti above 200 km should be warm but not absurd
    ds = request.getfixturevalue(ds_fixture)
    for temp in TEMPERATURES:
        vals = ds[temp].sel(alt_km=slice(200, 500)).values
        assert np.all(vals > 500), f"{ds_fixture} {temp} unreasonably cold above 200 km"
        assert np.all(vals < 15000), f"{ds_fixture} {temp} unreasonably hot above 200 km"


@pytest.mark.parametrize("ds_fixture", ["ds90", "ds20"])
def test_ne_peaks_in_f_region(ds_fixture, request):
    # NmF2 must fall in the F-region; anything outside 150-600 km is a model failure
    ds = request.getfixturevalue(ds_fixture)
    ne = ds["Ne"]
    peak_alt = float(ne.alt_km[ne.argmax("alt_km")])
    assert 150 < peak_alt < 600, f"{ds_fixture} Ne peak at {peak_alt:.0f} km — outside F-region"


@pytest.mark.parametrize("ds_fixture", ["ds90", "ds20"])
def test_o_plus_dominates_f_region(ds_fixture, request):
    # O+ must be the dominant ion in the F-region
    ds = request.getfixturevalue(ds_fixture)
    f_region = ds.sel(alt_km=slice(200, 400))
    o_plus = float(f_region["O+"].mean())
    assert o_plus > float(f_region["O2+"].mean())
    assert o_plus > float(f_region["NO+"].mean())
