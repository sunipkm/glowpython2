# %%
from __future__ import annotations
from datetime import datetime, UTC
import matplotlib
from glowpython2 import Iri90
from iri20py import Iri2020
from glowpython2.utils import alt_grid
import matplotlib.pyplot as plt
import numpy as np
# %%
usetex = False
if not usetex:
    # computer modern math text
    matplotlib.rcParams.update({'mathtext.fontset': 'cm'})

matplotlib.rc(
    'font', **{
        'family': 'serif',
        'serif': ['Times' if usetex else 'Times New Roman']
    }
)
# for Palatino and other serif fonts use:
# rc('font',**{'family':'serif','serif':['Palatino']})
matplotlib.rc('text', usetex=usetex)
# %%
MON = 12
iri90 = Iri90()
iri20 = Iri2020()
_, ds90 = iri90.evaluate(
    time=datetime(2015, MON, 13, 10, 0, 0, tzinfo=UTC),
    lat=42.6,
    lon=-71.2,
    alt=alt_grid(),
)
_, ds20 = iri20.evaluate(
    time=datetime(2015, MON, 13, 10, 0, 0, tzinfo=UTC),
    lat=42.6,
    lon=-71.2,
    alt=alt_grid(),
)
# %%
ig, ax = plt.subplots(figsize=(8, 6), dpi=300)
tax = ax.twiny()
species = ['O+', 'O2+', 'N+', 'NO+']
descs = ['O^+', 'O_2^+', 'N^+', 'NO^+']
# species = ['N+']
colors = ['r', 'g', 'b', 'm']
labels = []
lines = []
for spec, color, desc in zip(species, colors, descs):
    l21, = ds20[spec].plot(y='alt_km', ax=ax, color=color) # type: ignore
    l00, = ds90[spec].plot(y='alt_km', ax=ax, linestyle='--', color=color) # type: ignore
    lines.extend([l21, l00])
    labels.extend([f'IRI-20 ${desc}$', f'IRI-90 ${desc}$'])
ax.set_title('IRI-90 vs IRI-20')
l21, = ds20['Te'].plot(y='alt_km', ax=tax, color='k') # type: ignore
l00, = ds90['Te'].plot(y='alt_km', ax=tax, linestyle='--', color='k') # type: ignore
lines.extend([l21, l00])
labels.extend(['IRI-20 $T_e$', 'IRI-90 $T_e$'])
l21, = ds20['Ti'].plot(y='alt_km', ax=tax, color='c', alpha=0.7) # type: ignore
l00, = ds90['Ti'].plot(y='alt_km', ax=tax, linestyle='--', color='c', alpha=0.7) # type: ignore
lines.extend([l21, l00])
labels.extend(['IRI-20 $T_i$', 'IRI-90 $T_i$'])
ax.set_title('IRI-21 vs IRI-90')
ax.set_xscale('log')
ax.set_xlabel('Number Density [cm$^{-3}$]')
tax.set_xlabel('Temperature [K]')
ax.set_xlim(1e-3, None)
tax.set_xlim(100, None)
ax.legend(lines, labels, loc='upper left', fontsize='small')
plt.savefig('iri90_vs_iri20.png', dpi=300)
plt.show()
# %%
