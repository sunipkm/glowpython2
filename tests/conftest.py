import os
import pytest
from datetime import datetime, UTC

os.environ.setdefault("MPLBACKEND", "Agg")


@pytest.fixture(scope="session")
def sample_time():
    # Define a single time for all tests
    return datetime(2015, 12, 13, 10, 0, 0, tzinfo=UTC)


@pytest.fixture(scope="session")
def sample_lat():
    # Single latitude for tests
    return 65.1


@pytest.fixture(scope="session")
def sample_lon():
    # Longitude for all tests
    return -147.5
