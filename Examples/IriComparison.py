#!/usr/bin/env python
"""
Compare IRI-1990 and IRI-2020 ion densities and temperatures.
"""
import matplotlib
import matplotlib.pyplot as plt
from datetime import datetime, UTC

from glowpython2 import Iri90
from glowpython2.utils import alt_grid
from iri20py import Iri2020

matplotlib.rcParams.update({'mathtext.fontset': 'cm'})
matplotlib.rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})

time = datetime(2015, 12, 13, 10, 0, 0, tzinfo=UTC)
lat = 42.6
lon = -71.2
alt = alt_grid()

_, ds90 = Iri90().evaluate(time=time, lat=lat, lon=lon, alt=alt)
_, ds20 = Iri2020().evaluate(time=time, lat=lat, lon=lon, alt=alt)

fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
tax = ax.twiny()

species = ['O+', 'O2+', 'N+', 'NO+']
descs  = ['O^+', 'O_2^+', 'N^+', 'NO^+']
colors = ['r', 'g', 'b', 'm']
lines, labels = [], []

for spec, color, desc in zip(species, colors, descs):
    l20, = ds20[spec].plot(y='alt_km', ax=ax, color=color)
    l90, = ds90[spec].plot(y='alt_km', ax=ax, linestyle='--', color=color)
    lines  += [l20, l90]
    labels += [f'IRI-20 ${desc}$', f'IRI-90 ${desc}$']

for ds, ls, label_te, label_ti in [
    (ds20, '-',  'IRI-20 $T_e$', 'IRI-20 $T_i$'),
    (ds90, '--', 'IRI-90 $T_e$', 'IRI-90 $T_i$'),
]:
    l_te, = ds['Te'].plot(y='alt_km', ax=tax, color='k', linestyle=ls)
    l_ti, = ds['Ti'].plot(y='alt_km', ax=tax, color='c', alpha=0.7, linestyle=ls)
    lines  += [l_te, l_ti]
    labels += [label_te, label_ti]

ax.set_title('IRI-90 vs IRI-20')
ax.set_xscale('log')
ax.set_xlabel('Number Density [cm$^{-3}$]')
tax.set_xlabel('Temperature [K]')
ax.set_xlim(1e-3, None)
tax.set_xlim(100, None)
ax.legend(lines, labels, loc='upper left', fontsize='small')

plt.savefig('iri_comparison.png', dpi=300, bbox_inches='tight')
plt.show()
