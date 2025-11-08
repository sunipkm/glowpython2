# %%
from datetime import datetime, UTC
from glowpython2 import NrlMsis00, Msis00Settings
from msis21py import NrlMsis21
from glowpython2.utils import alt_grid
import matplotlib
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
MON = 1
msis00 = NrlMsis00()
msis21 = NrlMsis21()
ds00 = msis00.evaluate(
    time=datetime(2022, MON, 13, 10, 0, 0, tzinfo=UTC),
    lat=0.0,
    lon=0.0,
    alt=alt_grid(),
)
ds21 = msis21.evaluate(
    time=datetime(2022, MON, 13, 10, 0, 0, tzinfo=UTC),
    lat=0.0,
    lon=0.0,
    alt=alt_grid(),
)
# %%
fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
tax = ax.twiny()
species = ['O', 'O2', 'N2', 'NO']
descs = ['O', 'O$_2$', 'N$_2$', 'NO']
colors = ['r', 'g', 'b', 'm']
labels = []
lines = []
for spec, color, desc in zip(species, colors, descs):
    l21, = ds21[spec].plot(y='alt_km', ax=ax, color=color)
    l00, = ds00[spec].plot(y='alt_km', ax=ax, linestyle='--', color=color)
    lines.extend([l21, l00])
    labels.extend([f'MSIS21 {desc}', f'MSIS00 {desc}'])
ax.set_title('MSIS00 vs MSIS21')
l21, = ds21['Tn'].plot(y='alt_km', ax=tax, color='k')
l00, = ds00['Tn'].plot(y='alt_km', ax=tax, linestyle='--', color='k')
lines.extend([l21, l00])
labels.extend(['MSIS21 $T_n$', 'MSIS00 $T_n$'])
ax.set_title('MSIS00 vs MSIS21')
ax.set_xscale('log')
ax.legend(lines, labels)
plt.savefig('msis00_vs_msis21.png', dpi=300)
plt.show()
# %%
