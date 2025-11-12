# %%
from datetime import UTC
from matplotlib.gridspec import GridSpec
import pytz
from glowpython import no_precipitation as glownoprecip
from glowpython2 import no_precipitation
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from dateutil.parser import parse
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
time = parse('2022-3-22T18:00:00')
glat = 42.6
glon = -71.2
Nbins = 250
tec = None
versions = ['GLOW', 'MSIS00_IRI90', 'MSIS21_IRI20']
lstyles = ['-', '--', '-.']
alphas = [0.5, 0.7, 0.5]
biglabels = ['GLOW', 'MSIS-2000 + IRI-1990', 'MSIS-2.1 + IRI-2020']
ionos = {}
for version in versions:
    if version == 'GLOW':
        iono = glownoprecip(time, glat, glon, Nbins)
    else:
        magmodel = 'IGRF14' if version == 'MSIS21_IRI20' else 'POGO68'
        iono = no_precipitation(time, glat, glon, Nbins, tec=tec, version=version, magmodel=magmodel)  # type: ignore
    ionos[version] = iono
# %% Density plot
grid = GridSpec(2, 2, height_ratios=[0.05, 1], width_ratios=[1, 1], hspace=0.1, wspace=0.15)
den_figure = plt.figure(figsize=(6.4, 4.8), dpi=300)
den_axs = np.asarray([den_figure.add_subplot(grid[1, 0]), den_figure.add_subplot(grid[1, 1])])
legend_axs = den_figure.add_subplot(grid[0, :])
legend_axs.axis('off')
den_lines = [{}, {}]  # neutrals, ions
descs = {
    "O": "O", "O2": "O$_2$", "N2": "N$_2$", "NO": "NO",
    "O+": "O$^+$", "O2+": "O$_2^+$", "NO+": "NO$^+$", "N(2D)": "N($^2$D)", 'NeIn': 'e$^-$'
}
for version, linestyle, alpha in zip(versions, lstyles, alphas):
    iono = ionos[version]
    ax = den_axs[0]
    ax.set_prop_cycle(None)  # reset color cycle
    for v in ("O", "O2", "N2", "NO"):
        l, = ax.plot(iono[v], iono[v].alt_km, label=f'{version} {v}', linestyle=linestyle, lw=0.75, alpha=alpha)
        if v not in den_lines[0]:
            den_lines[0][v] = l

for version, linestyle, alpha in zip(versions, lstyles, alphas):
    iono = ionos[version]
    ax = den_axs[1]
    ax.set_prop_cycle(None)  # reset color cycle
    for v in ("O+", "O2+", "NO+", "N(2D)", 'NeIn'):
        l, = ax.plot(iono[v], iono[v].alt_km, label=f'{version} {v}', linestyle=linestyle, lw=0.75, alpha=alpha)
        if v not in den_lines[1]:
            den_lines[1][v] = l

den_axs[0].set_xscale("log")
den_axs[0].set_ylim(60, 1000)
den_axs[0].set_ylabel("Altitude [km]")
den_axs[0].set_title("Neutrals", fontsize='medium')
den_axs[0].grid(True)
den_axs[0].set_xlim(1, None)
lines = []
labels = []
for k, v in den_lines[0].items():
    lines.append(v)
    labels.append(descs.get(k))
den_axs[0].legend(lines, labels, loc="best", fontsize='small')
lines = []
labels = []
den_axs[1].set_xscale("log")
den_axs[1].set_title("Ions", fontsize='medium')
den_axs[1].grid(True)
den_axs[1].set_xlim(1, None)
den_axs[1].yaxis.set_ticklabels([])
for k, v in den_lines[1].items():
    lines.append(v)
    labels.append(descs.get(k))
den_axs[1].legend(lines, labels, loc="best", fontsize='small')
for lab, ls, alpha in zip(biglabels, lstyles, alphas):
    legend_axs.plot([], [], label=lab, color='k', linestyle=ls, alpha=alpha)
legend_axs.legend(loc='center', ncol=3, fontsize='small', mode='expand', frameon=False)
den_figure.text(0.5, 0.04, "Number Density [cm$^{-3}$]", ha='center')
den_figure.text(0.5, 0.9, time.astimezone(UTC).astimezone(pytz.timezone('US/Eastern')).isoformat(sep=' '), ha='center')
den_figure.savefig('glow_den.png', dpi=300, bbox_inches='tight')
plt.show()
# %%
# Temperature plot
temp_grid = GridSpec(2, 1, height_ratios=[0.05, 1], hspace=0.1)
temp_figure = plt.figure(figsize=(6.4, 4.8), dpi=300)
temp_ax = temp_figure.add_subplot(temp_grid[1, 0])
temp_legend_ax = temp_figure.add_subplot(temp_grid[0, 0])
temp_legend_ax.axis('off')
temp_lines = {}
for version, linestyle, alpha in zip(versions, lstyles, alphas):
    iono = ionos[version]
    ax = temp_ax
    ax.set_prop_cycle(None)  # reset color cycle
    l, = ax.plot(iono['Te'], iono['Te'].alt_km, label=f'{version} Te', linestyle=linestyle, lw=0.75, alpha=alpha)
    if 'Te' not in temp_lines:
        temp_lines['Te'] = l
    l, = ax.plot(iono['Ti'], iono['Ti'].alt_km, label=f'{version} Ti', linestyle=linestyle, lw=0.75, alpha=alpha)
    if 'Ti' not in temp_lines:
        temp_lines['Ti'] = l
    l, = ax.plot(iono['Tn'], iono['Tn'].alt_km, label=f'{version} Tn', linestyle=linestyle, lw=0.75, alpha=alpha)
    if 'Tn' not in temp_lines:
        temp_lines['Tn'] = l
temp_ax.set_xlabel("Temperature [K]")
temp_ax.set_ylabel("altitude [km]")
temp_ax.set_title("Temperature")
temp_ax.grid(True)
temp_ax.set_ylim(60, 1000)
temp_ax.set_xlim(100, None)
lines = []
labels = []
for k, v in temp_lines.items():
    lines.append(v)
    labels.append(f'${k}$')
temp_ax.legend(lines, labels, loc="best")
for lab, ls, alpha in zip(biglabels, lstyles, alphas):
    temp_legend_ax.plot([], [], label=lab, color='k', linestyle=ls, alpha=alpha)
temp_legend_ax.legend(loc='center', ncol=3, fontsize='small', mode='expand', frameon=False)
temp_figure.text(0.5, 0.9, time.astimezone(UTC).astimezone(pytz.timezone('US/Eastern')).isoformat(sep=' '), ha='center')
temp_figure.savefig('glow_temp.png', dpi=300, bbox_inches='tight')
plt.show()
# %% VER plot
ver_grid = GridSpec(2, 3, height_ratios=[0.05, 1], hspace=0.1, wspace=0.15)
ver_figure = plt.figure(figsize=(6.4, 4.8), dpi=300)
ver_axs = np.asarray([
    ver_figure.add_subplot(ver_grid[1, 0]),
    ver_figure.add_subplot(ver_grid[1, 1]),
    ver_figure.add_subplot(ver_grid[1, 2])
])
ver_legend_ax = ver_figure.add_subplot(ver_grid[0, :])
ver_legend_ax.axis('off')
ver_lines = {
    'Visible': {},
    'Infrared': {},
    'Ultraviolet': {}
}
ver_components = {
    'Visible': ["4278", "5577", "6300", "5200"],
    'Infrared': ["7320", "7774", "8446", "10400"],
    'Ultraviolet': ["1304", "1356", "1493", "3371", "3644", "3726", "LBH"]
}
for version, linestyle, alpha in zip(versions, lstyles, alphas):
    iono = ionos[version]
    for ax, kind in zip(ver_axs, ('Visible', 'Infrared', 'Ultraviolet')):
        ax.set_prop_cycle(None)  # reset color cycle
        for line in ver_components[kind]:
            l, = ax.plot(iono['ver'].sel(wavelength=line), iono['ver'].alt_km, label=f'{version} {line} Å', linestyle=linestyle, lw=0.75, alpha=alpha)
            if line not in ver_lines[kind]:
                ver_lines[kind][line] = l

for ax, kind in zip(ver_axs, ('Visible', 'Infrared', 'Ultraviolet')):
    ax.set_title(f"{kind}", fontsize='medium')
    ax.grid(True)
    ax.set_ylim(60, 1000)
    ax.set_xlim(1e-5, None)
    ax.set_xscale("log")
    lines = []
    labels = []
    for k, v in ver_lines[kind].items():
        lines.append(v)
        if k != 'LBH':
            labels.append(f'{k} Å')
        else:
            labels.append('LBH')
    ax.legend(lines, labels, loc="best", fontsize='small')
for ax in ver_axs[1:]:
    ax.yaxis.set_ticklabels([])
ver_axs[0].set_ylabel("Altitude [km]")
ver_figure.text(0.5, 0.04, "Volume Emission Rate [cm$^{-3}$ s$^{-1}$]", ha='center')
ver_figure.text(0.5, 0.9, time.astimezone(UTC).astimezone(pytz.timezone('US/Eastern')).isoformat(sep=' '), ha='center')
for lab, ls, alpha in zip(biglabels, lstyles, alphas):
    ver_legend_ax.plot([], [], label=lab, color='k', linestyle=ls, alpha=alpha)
ver_legend_ax.legend(loc='center', ncol=3, fontsize='small', mode='expand', frameon=False)
ver_figure.savefig('glow_ver.png', dpi=300, bbox_inches='tight')
plt.show()

# %%
