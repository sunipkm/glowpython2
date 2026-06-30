#!/usr/bin/env python
"""
Compare GLOW model output across atmospheric model configurations:
  - MSIS-2000 + IRI-1990
  - MSIS-2.1 + IRI-2020
  - MSIS-2.1 + IRI-2020 + ModGLOW coefficients
"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from matplotlib.gridspec import GridSpec

from glowpython2 import no_precipitation

matplotlib.rcParams.update({'mathtext.fontset': 'cm'})
matplotlib.rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})

time = datetime(2022, 3, 22, 22, 0, 0)
glat = 42.6
glon = -71.2
Nbins = 250

versions = ['MSIS00_IRI90', 'MSIS21_IRI20', 'MODGLOW']
biglabels = ['MSIS-2000 + IRI-1990', 'MSIS-2.1 + IRI-2020', 'MSIS-2.1 + IRI-2020 + ModGLOW']
lstyles = ['-', '--', '-.']
alphas = [0.5, 0.7, 0.5]

ionos = {}
for version in versions:
    if version == 'MODGLOW':
        iono = no_precipitation(time, glat, glon, Nbins, version='MSIS21_IRI20', magmodel='IGRF14', newcoeffs=True)
    else:
        magmodel = 'IGRF14' if version == 'MSIS21_IRI20' else 'POGO68'
        iono = no_precipitation(time, glat, glon, Nbins, version=version, magmodel=magmodel)
    ionos[version] = iono

# Density plot
grid = GridSpec(2, 2, height_ratios=[0.05, 1], width_ratios=[1, 1], hspace=0.1, wspace=0.15)
den_figure = plt.figure(figsize=(6.4, 4.8), dpi=300)
den_axs = np.asarray([den_figure.add_subplot(grid[1, 0]), den_figure.add_subplot(grid[1, 1])])
legend_ax = den_figure.add_subplot(grid[0, :])
legend_ax.axis('off')

descs = {
    "O": "O", "O2": "O$_2$", "N2": "N$_2$", "NO": "NO",
    "O+": "O$^+$", "O2+": "O$_2^+$", "NO+": "NO$^+$", "N(2D)": "N($^2$D)", 'NeIn': 'e$^-$'
}
den_lines = [{}, {}]
for version, ls, alpha in zip(versions, lstyles, alphas):
    iono = ionos[version]
    den_axs[0].set_prop_cycle(None)
    for v in ("O", "O2", "N2", "NO"):
        l, = den_axs[0].plot(iono[v], iono[v].alt_km, linestyle=ls, lw=0.75, alpha=alpha)
        if v not in den_lines[0]:
            den_lines[0][v] = l
    den_axs[1].set_prop_cycle(None)
    for v in ("O+", "O2+", "NO+", "N(2D)", 'NeIn'):
        l, = den_axs[1].plot(iono[v], iono[v].alt_km, linestyle=ls, lw=0.75, alpha=alpha)
        if v not in den_lines[1]:
            den_lines[1][v] = l

for ax, lines, title in zip(den_axs, den_lines, ['Neutrals', 'Ions']):
    ax.set_xscale('log')
    ax.set_ylim(60, 1000)
    ax.set_xlim(1, None)
    ax.grid(True)
    ax.set_title(title, fontsize='medium')
    ax.legend([v for v in lines.values()], [descs.get(k) for k in lines], loc='best', fontsize='small')
den_axs[0].set_ylabel('Altitude [km]')
den_axs[1].yaxis.set_ticklabels([])
for lab, ls, alpha in zip(biglabels, lstyles, alphas):
    legend_ax.plot([], [], label=lab, color='k', linestyle=ls, alpha=alpha)
legend_ax.legend(loc='center', ncol=3, fontsize='small', mode='expand', frameon=False)
den_figure.text(0.5, 0.04, 'Number Density [cm$^{-3}$]', ha='center')
den_figure.savefig('model_comparison_density.png', dpi=300, bbox_inches='tight')

# Temperature plot
temp_grid = GridSpec(2, 1, height_ratios=[0.05, 1], hspace=0.1)
temp_figure = plt.figure(figsize=(6.4, 4.8), dpi=300)
temp_ax = temp_figure.add_subplot(temp_grid[1, 0])
temp_legend_ax = temp_figure.add_subplot(temp_grid[0, 0])
temp_legend_ax.axis('off')
temp_lines = {}
for version, ls, alpha in zip(versions, lstyles, alphas):
    iono = ionos[version]
    temp_ax.set_prop_cycle(None)
    for v in ('Te', 'Ti', 'Tn'):
        l, = temp_ax.plot(iono[v], iono[v].alt_km, linestyle=ls, lw=0.75, alpha=alpha)
        if v not in temp_lines:
            temp_lines[v] = l
temp_ax.set_xlabel('Temperature [K]')
temp_ax.set_ylabel('Altitude [km]')
temp_ax.set_title('Temperature')
temp_ax.grid(True)
temp_ax.set_ylim(60, 1000)
temp_ax.set_xlim(100, None)
temp_ax.legend(temp_lines.values(), [f'${k}$' for k in temp_lines], loc='best')
for lab, ls, alpha in zip(biglabels, lstyles, alphas):
    temp_legend_ax.plot([], [], label=lab, color='k', linestyle=ls, alpha=alpha)
temp_legend_ax.legend(loc='center', ncol=3, fontsize='small', mode='expand', frameon=False)
temp_figure.savefig('model_comparison_temperature.png', dpi=300, bbox_inches='tight')

# VER plot
ver_grid = GridSpec(2, 3, height_ratios=[0.05, 1], hspace=0.1, wspace=0.15)
ver_figure = plt.figure(figsize=(6.4, 4.8), dpi=300)
ver_axs = [ver_figure.add_subplot(ver_grid[1, i]) for i in range(3)]
ver_legend_ax = ver_figure.add_subplot(ver_grid[0, :])
ver_legend_ax.axis('off')
ver_components = {
    'Visible':     ["4278", "5577", "6300", "5200"],
    'Infrared':    ["7320", "7774", "8446", "10400"],
    'Ultraviolet': ["1304", "1356", "1493", "3371", "3644", "3726", "LBH"],
}
ver_lines = {k: {} for k in ver_components}
for version, ls, alpha in zip(versions, lstyles, alphas):
    iono = ionos[version]
    for ax, kind in zip(ver_axs, ver_components):
        ax.set_prop_cycle(None)
        for line in ver_components[kind]:
            l, = ax.plot(iono['ver'].sel(wavelength=line), iono['ver'].alt_km, linestyle=ls, lw=0.75, alpha=alpha)
            if line not in ver_lines[kind]:
                ver_lines[kind][line] = l
for ax, kind in zip(ver_axs, ver_components):
    ax.set_title(kind, fontsize='medium')
    ax.grid(True)
    ax.set_ylim(60, 1000)
    ax.set_xlim(1e-5, None)
    ax.set_xscale('log')
    labels = [f'{k} Å' if k != 'LBH' else 'LBH' for k in ver_lines[kind]]
    ax.legend(ver_lines[kind].values(), labels, loc='best', fontsize='small')
ver_axs[0].set_ylabel('Altitude [km]')
for ax in ver_axs[1:]:
    ax.yaxis.set_ticklabels([])
for lab, ls, alpha in zip(biglabels, lstyles, alphas):
    ver_legend_ax.plot([], [], label=lab, color='k', linestyle=ls, alpha=alpha)
ver_legend_ax.legend(loc='center', ncol=3, fontsize='small', mode='expand', frameon=False)
ver_figure.text(0.5, 0.04, 'Volume Emission Rate [cm$^{-3}$ s$^{-1}$]', ha='center')
ver_figure.savefig('model_comparison_ver.png', dpi=300, bbox_inches='tight')

plt.show()
