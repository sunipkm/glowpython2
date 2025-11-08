from matplotlib.pyplot import figure
import xarray
import numpy as np

class Plot:
    def __init__(self, fig_kwargs=None):
        self.fig_kwargs = fig_kwargs if fig_kwargs is not None else {}
        self._density_figure = None
        self._density_axs = None
        self._precip_figure = None
        self._precip_ax = None
        self._temperature_figure = None
        self._temperature_ax = None
        self._ver_figure = None
        self._ver_axs = None

    def density(self, iono: xarray.Dataset, linestyle: str = '-'):
        if self._density_figure is None:
            self._density_figure = figure(**self.fig_kwargs)
            self._density_axs = self._density_figure.subplots(1, 2, sharey=True)
        fig = self._density_figure
        axs = self._density_axs
        if axs is not None:
            for ax in axs: ax.set_prop_cycle(None)  # reset color cycle
            fig.suptitle("Number Density")

            ax = axs[0]
            for v in ("O", "N2", "O2", "NO"):
                ax.plot(iono[v], iono[v].alt_km, label=v, linestyle=linestyle)
            ax.set_xscale("log")
            ax.set_ylabel("altitude [km]")
            ax.set_xlabel("Density [cm$^{-3}$]")
            ax.set_title("Neutrals")
            ax.grid(True)
            ax.set_xlim(1, None)
            ax.legend(loc="best")

            ax = axs[1]
            for v in ("O+", "O2+", "NO+", "N(2D)", 'NeIn'):
                ax.plot(iono[v], iono[v].alt_km, label=v, linestyle=linestyle)
            ax.set_xscale("log")
            ax.set_title("Ions")
            ax.grid(True)
            ax.set_xlim(1, None)
            ax.legend(loc="best")


    def precip(self, precip: xarray.DataArray, linestyle: str = '-'):
        if self._precip_figure is None:
            self._precip_figure = figure(**self.fig_kwargs)
            self._precip_ax = self._precip_figure.gca()
        ax = self._precip_ax
        if ax is not None:
            ax.set_prop_cycle(None)  # reset color cycle
            ax.plot(precip["energy"] / 1e3, precip, linestyle=linestyle)
            ax.set_xlabel("Energy bin centers [keV]")
            ax.set_ylabel("hemispherical flux [cm$^{-2}$ s$^{-1}$ eV$^{-1}$]")
            ax.set_title("precipitation: differential number flux")
            ax.grid(True)


    def temperature(self, iono: xarray.Dataset, linestyle: str = '-'):
        if self._temperature_figure is None:
            self._temperature_figure = figure(**self.fig_kwargs)
            self._temperature_ax = self._temperature_figure.gca()
        ax = self._temperature_ax
        time = iono.time
        location = iono.glatlon
        tail = f"\n{time} {location}"
        if ax is not None:
            ax.set_prop_cycle(None)  # reset color cycle
            ax.plot(iono["Ti"], iono["Ti"].alt_km, label="$T_i$", linestyle=linestyle)
            ax.plot(iono["Te"], iono["Te"].alt_km, label="$T_e$", linestyle=linestyle)
            ax.plot(iono["Tn"], iono["Tn"].alt_km, label="$T_n$", linestyle=linestyle)
            ax.set_xlabel("Temperature [K]")
            ax.set_ylabel("altitude [km]")
            ax.set_title("Ion, Electron, Neutral temperatures" + tail)
            ax.grid(True)
            ax.legend()


    def ver(self, iono: xarray.Dataset, linestyle: str = '-'):
        if self._ver_figure is None:
            self.fig_kwargs['constrained_layout'] = True
            self._ver_figure = figure(**self.fig_kwargs)
            self._ver_axs = self._ver_figure.subplots(1, 3, sharey=True)
        fig = self._ver_figure
        axs = self._ver_axs

        time = iono.time
        location = iono.glatlon
        tail = f"\n{time} {location}"

        fig.suptitle(tail)
        if axs is not None:
            for ax in axs: ax.set_prop_cycle(None)  # reset color cycle
            ver_group(iono["ver"].loc[:, ["4278", "5577", "6300", "5200"]], "Visible", axs[0], linestyle=linestyle)
            ver_group(iono["ver"].loc[:, ["7320", "7774", "8446", "10400"]], "Infrared", axs[1], linestyle=linestyle)
            ver_group(
                iono["ver"].loc[:, ["3371", "3644", "3726", "1356", "1493", "1304", "LBH"]],
                "Ultraviolet",
                axs[2],
                linestyle=linestyle
            )
            axs[0].set_ylabel("altitude [km]")
            axs[0].set_xlabel("Volume Emission Rate [Rayleigh]")
            axs[0].set_xlim(1e-3, None)
            axs[1].set_xlim(1e-3, None)
            axs[2].set_xlim(1e-3, None)


def ver_group(iono: xarray.DataArray, ttxt: str, ax, linestyle: str = '-'):
    nm = np.nanmax(iono)
    if nm == 0 or np.isnan(nm):
        return

    colors = {
        "4278": "blue",
        "5577": "xkcd:dark lime green",
        "5200": "xkcd:golden yellow",
        "6300": "red",
    }

    for w in iono.wavelength:
        ax.plot(iono.loc[:, w], iono.alt_km, label=w.item(), color=colors.get(w.item()), linestyle=linestyle)
    ax.set_xscale("log")
    ax.set_ylim(90, 500)
    ax.set_title(ttxt)
    ax.grid(True)
    ax.legend()
