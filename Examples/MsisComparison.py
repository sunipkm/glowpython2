#!/usr/bin/env python
"""
Compare NRLMSISE-00 and NRLMSIS-2.1 neutral densities and temperatures.
"""
import matplotlib
import matplotlib.pyplot as plt
from datetime import datetime, UTC

from glowpython2 import NrlMsis00
from glowpython2.utils import alt_grid
from msis21py import NrlMsis21

matplotlib.rcParams.update({'mathtext.fontset': 'cm'})
matplotlib.rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})

time = datetime(2022, 1, 13, 10, 0, 0, tzinfo=UTC)
lat = 0.0
lon = 0.0
alt = alt_grid()

ds00 = NrlMsis00().evaluate(time=time, lat=lat, lon=lon, alt=alt)
ds21 = NrlMsis21().evaluate(time=time, lat=lat, lon=lon, alt=alt)

fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
tax = ax.twiny()

species = ['O', 'O2', 'N2', 'NO']
descs  = ['O', 'O$_2$', 'N$_2$', 'NO']
colors = ['r', 'g', 'b', 'm']
lines, labels = [], []

for spec, color, desc in zip(species, colors, descs):
    l21, = ds21[spec].plot(y='alt_km', ax=ax, color=color)
    l00, = ds00[spec].plot(y='alt_km', ax=ax, linestyle='--', color=color)
    lines  += [l21, l00]
    labels += [f'MSIS-21 {desc}', f'MSIS-00 {desc}']

l21, = ds21['Tn'].plot(y='alt_km', ax=tax, color='k')
l00, = ds00['Tn'].plot(y='alt_km', ax=tax, color='k', linestyle='--')
lines  += [l21, l00]
labels += ['MSIS-21 $T_n$', 'MSIS-00 $T_n$']

ax.set_title('NRLMSISE-00 vs NRLMSIS-2.1')
ax.set_xscale('log')
ax.set_xlabel('Number Density [cm$^{-3}$]')
tax.set_xlabel('Neutral Temperature [K]')
ax.legend(lines, labels)

plt.savefig('msis_comparison.png', dpi=300, bbox_inches='tight')
plt.show()
