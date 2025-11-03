# %%
from __future__ import annotations
from xarray import DataArray, Variable
from typing import List, Tuple, Union, SupportsFloat as Numeric
import numpy as np

# %%
def rcolumn(
    major_den: DataArray,
    t_n: DataArray,
    chi: Numeric,
    ) -> Tuple[DataArray, DataArray]:
    """## Calculates the column density for each species at zenith angle `chi`.
    First, the vertical column density is calculated, and then a fit to Chapman Grazing
    Incidence Integral [Smith and Smith, JGR 77, 3592, 1972] is used to calculate the
    slant column density.
    - `chi` < 90 deg: Column densities are calculated directly.
    - `chi` > 90 deg: Column density at grazing height for 90 deg is calculated and doubled, and the column density above the grazing height is subtracted.
    - If grazing height is lower than the bottom of the atmosphere supplied, column densities are set to infinity (`np.inf`).

    ### Args:
        - `major_den (DataArray)`: Major species densities.
        - `t_n (DataArray)`: Neutral temperature profile.
        - `chi (Numeric)`: Zenith angle.

    ### Returns:
        - `Tuple[DataArray, DataArray]`: Column densities, vertical column densities.
    """