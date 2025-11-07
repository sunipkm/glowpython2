# %%
from datetime import datetime, UTC
from glowpython2 import Iri90
from iri20py import Iri2020
from glowpython2.utils import alt_grid
import matplotlib.pyplot as plt
import numpy as np
# %%
iri90 = Iri90()
iri20 = Iri2020()
_, ds90 = iri90.evaluate(
    time=datetime(2022, 3, 21, 12, 0, 0, tzinfo=UTC),
    lat=0.0,
    lon=0.0,
    alt=alt_grid(),
)
_, ds20 = iri20.evaluate(
    time=datetime(2022, 3, 21, 12, 0, 0, tzinfo=UTC),
    lat=0.0,
    lon=0.0,
    alt=alt_grid(),
)
# %%
ig, ax = plt.subplots(figsize=(8, 6), dpi=300)
tax = ax.twiny()
species = ['O+', 'O2+', 'N+', 'NO+']
# species = ['N+']
colors = ['r', 'g', 'b', 'm']
labels = []
lines = []
for spec, color in zip(species, colors):
    l21, = ds20[spec].plot(y='alt_km', ax=ax, color=color)
    l00, = ds90[spec].plot(y='alt_km', ax=ax, linestyle='--', color=color)
    lines.extend([l21, l00])
    labels.extend([f'IRI-20 {spec}', f'IRI-90 {spec}'])
ax.set_title('IRI-90 vs IRI-20')
l21, = ds20['Te'].plot(y='alt_km', ax=tax, color='k')
l00, = ds90['Te'].plot(y='alt_km', ax=tax, linestyle='--', color='k')
lines.extend([l21, l00])
labels.extend(['IRI-20 T_e', 'IRI-90 T_e'])
l21, = ds20['Ti'].plot(y='alt_km', ax=tax, color='c', alpha=0.7)
l00, = ds90['Ti'].plot(y='alt_km', ax=tax, linestyle='--', color='c', alpha=0.7)
lines.extend([l21, l00])
labels.extend(['IRI-20 T_i', 'IRI-90 T_i'])
ax.set_title('IRI-21 vs IRI-90')
ax.set_xscale('log')
ax.set_xlim(1e-3, None)
tax.set_xlim(100, None)
ax.legend(lines, labels, loc='upper left', fontsize='small')
plt.savefig('iri90_vs_iri20.png', dpi=300)
plt.show()
# %%
