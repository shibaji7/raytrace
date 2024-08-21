import matplotlib as mpl
import matplotlib.dates as mdates
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import utils as utils
from matplotlib.dates import DateFormatter
from loguru import logger

plt.style.use(["science", "ieee"])
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans", "Lucida Grande", "Verdana"]


class RangeTimeIntervalPlot(object):
    """
    Create plots for velocity, width, power, elevation angle, etc.
    """

    def __init__(self, nrang, dates, rad, fig_title="", num_subplots=3, srange_type="slist"):
        self.nrang = nrang
        self.unique_gates = np.linspace(1, nrang, nrang)
        self.rad = rad
        self.dates = dates
        self.num_subplots = num_subplots
        self._num_subplots_created = 0
        self.fig = plt.figure(
            figsize=(8, 3 * num_subplots), dpi=300
        )  # Size for website
        self.fig_title = fig_title
        self.srange_type = srange_type
        utils.setsize(12)
        return

    def addParamPlot(
        self,
        df: pd.DataFrame,
        beam: int,
        title: str,
        p_max: float = 30,
        p_min: float = -30,
        xlabel: str = "Time, UT",
        zparam: str = "v",
        label: str = r"Velocity, ($ms^{-1}$)",
        cmap: str = "Spectral",
        color_bar: bool = True,
        lay_eclipse: bool = False,
    ):
        ax = self._add_axis()
        df = df[df.bmnum == beam]
        X, Y, Z = utils.get_gridded_parameters(
            df, xparam="time", yparam="slist", zparam=zparam, rounding=False
        )
        cmap = mpl.cm.get_cmap(cmap)
        # Configure axes
        ax.xaxis.set_major_formatter(DateFormatter(r"%H^{%M}"))
        hours = mdates.HourLocator(byhour=range(0, 24, 3))
        ax.xaxis.set_major_locator(hours)
        dtime = (
            pd.Timestamp(self.dates[-1]).to_pydatetime()
            - pd.Timestamp(self.dates[0]).to_pydatetime()
        ).total_seconds() / 3600.0
        if dtime > 2.0 and dtime < 4.0:
            hours = mdates.HourLocator(byhour=range(0, 24, 1))
            ax.xaxis.set_minor_locator(hours)
            ax.xaxis.set_minor_formatter(DateFormatter(r"%H^{%M}"))
        elif dtime < 2.0:
            minutes = mdates.MinuteLocator(byminute=range(0, 60, 10))
            ax.xaxis.set_minor_locator(minutes)
            ax.xaxis.set_minor_formatter(DateFormatter(r"%H^{%M}"))
        ax.set_xlabel(xlabel, fontdict={"size": 12, "fontweight": "bold"})
        ax.set_xlim([self.dates[0], self.dates[-1]])
        ax.set_ylim([0, self.nrang])
        ax.set_ylabel("Range gate", fontdict={"size": 12, "fontweight": "bold"})
        im = ax.pcolormesh(
            X, Y, Z.T, lw=0.01, edgecolors="None", cmap=cmap, vmax=p_max, vmin=p_min,
            zorder=3
        )
        if color_bar:
            self._add_colorbar(im, self.fig, ax, label=label)
        ax.text(0.01, 0.95, title, ha="left", va="center", transform=ax.transAxes)
        if lay_eclipse:
            self.overlay_eclipse(ax, beam)
        return

    def overlay_eclipse(self, ax, beam):
        file = f"datasets/eclipse_path/{self.dates[0].strftime('%Y-%m-%d')}/oc.{self.rad}.{beam}.csv"
        logger.info(f"Overlay eclipse path from file: {file}")
        # Add eclipses
        o = pd.read_csv( file, parse_dates=["dates"], )
        o = o[(o.dates >= self.dates[0]) & (o.dates <= self.dates[1])]
        srange = (
            np.arange(101) 
            if self.srange_type == "slist" else
            180 + (45 * np.arange(101))
        )
        tags = [f"gate_{i}" for i in range(101)]
        p = o[tags].values
        p = np.ma.masked_where(p <= 0, p)
        im = ax.contourf(
            o.dates.tolist(),
            srange,
            p.T,
            cmap="gray_r",
            zorder=3,
            alpha=0.3,
            levels=[0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0],
        )
        ax.set_xlim([self.dates[0], self.dates[-1]])
        return

    def addGSIS(self, df, beam, title, xlabel="", ylabel="Range gate", zparam="gflg"):
        # add new axis
        ax = self._add_axis()
        df = df[df.bmnum == beam]
        X, Y, Z = utils.get_gridded_parameters(
            df, xparam="time", yparam="slist", zparam=zparam, rounding=False
        )
        flags = np.array(df[zparam]).astype(int)
        if -1 in flags and 2 in flags:  # contains noise flag
            cmap = mpl.colors.ListedColormap(
                [
                    (0.0, 0.0, 0.0, 1.0),  # black
                    (1.0, 0.0, 0.0, 1.0),  # blue
                    (0.0, 0.0, 1.0, 1.0),  # red
                    (0.0, 1.0, 0.0, 1.0),
                ]
            )  # green
            bounds = [
                -1,
                0,
                1,
                2,
                3,
            ]  # Lower bound inclusive, upper bound non-inclusive
            handles = [
                mpatches.Patch(color="red", label="IS"),
                mpatches.Patch(color="blue", label="GS"),
                mpatches.Patch(color="black", label="US"),
                mpatches.Patch(color="green", label="SAIS"),
            ]
        elif -1 in flags and 2 not in flags:
            cmap = mpl.colors.ListedColormap(
                [
                    (0.0, 0.0, 0.0, 1.0),  # black
                    (1.0, 0.0, 0.0, 1.0),  # blue
                    (0.0, 0.0, 1.0, 1.0),
                ]
            )  # red
            bounds = [-1, 0, 1, 2]  # Lower bound inclusive, upper bound non-inclusive
            handles = [
                mpatches.Patch(color="red", label="IS"),
                mpatches.Patch(color="blue", label="GS"),
                mpatches.Patch(color="black", label="US"),
            ]
        else:
            cmap = mpl.colors.ListedColormap(
                [(1.0, 0.0, 0.0, 1.0), (0.0, 0.0, 1.0, 1.0)]  # blue
            )  # red
            bounds = [0, 1, 2]  # Lower bound inclusive, upper bound non-inclusive
            handles = [
                mpatches.Patch(color="red", label="IS"),
                mpatches.Patch(color="blue", label="GS"),
            ]
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        ax.xaxis.set_major_formatter(DateFormatter(r"%H^{%M}"))
        hours = mdates.HourLocator(byhour=range(0, 24, 4))
        ax.xaxis.set_major_locator(hours)
        dtime = (
            pd.Timestamp(self.dates[-1]).to_pydatetime()
            - pd.Timestamp(self.dates[0]).to_pydatetime()
        ).total_seconds() / 3600.0
        if dtime < 4.0:
            minutes = mdates.MinuteLocator(byminute=range(0, 60, 10))
            ax.xaxis.set_minor_locator(minutes)
            ax.xaxis.set_minor_formatter(DateFormatter(r"%H^{%M}"))
        ax.set_xlabel(xlabel, fontdict={"size": 12, "fontweight": "bold"})
        ax.set_xlim([self.dates[0], self.dates[-1]])
        ax.set_ylim([0, self.nrang])
        ax.set_ylabel(ylabel, fontdict={"size": 12, "fontweight": "bold"})
        ax.pcolormesh(X, Y, Z.T, lw=0.01, edgecolors="None", cmap=cmap, norm=norm)
        ax.set_title(title, loc="left", fontdict={"fontweight": "bold"})
        ax.legend(handles=handles, loc=1)
        return ax

    def add_range_cell_data(
        self,
        df,
        rc,
        title="",
        xlabel="Time, UT",
        ylabel=r"Velocity, $ms^{-1}$",
        bounds=True,
    ):
        x, y = (
            df[(df.bmnum == rc["bmnum"]) & (df.slist == rc["gate"])].time,
            df[(df.bmnum == rc["bmnum"]) & (df.slist == rc["gate"])].v,
        )
        ax = self._add_axis() if self.rc_ax is None else self.rc_ax
        ax.xaxis.set_major_formatter(DateFormatter(r"$%H^{%M}$"))
        hours = mdates.HourLocator(byhour=range(0, 24, 1))
        ax.xaxis.set_major_locator(hours)
        dtime = (
            pd.Timestamp(self.dates[-1]).to_pydatetime()
            - pd.Timestamp(self.dates[0]).to_pydatetime()
        ).total_seconds() / 3600.0
        if dtime < 4.0:
            minutes = mdates.MinuteLocator(byminute=range(0, 60, 10))
            ax.xaxis.set_minor_locator(minutes)
            ax.xaxis.set_minor_formatter(DateFormatter(r"$%H^{%M}$"))
        ax.set_xlabel(xlabel, fontdict={"size": 12, "fontweight": "bold"})
        ax.set_xlim([self.dates[0], self.dates[-1]])
        # ax.set_ylim([0, self.nrang])
        ax.set_ylabel(ylabel, fontdict={"size": 12, "fontweight": "bold"})
        ax.set_title(title, loc="left", fontdict={"fontweight": "bold"})
        ax.plot(x, y, color="k", ls="-", lw=0.6)
        if bounds:
            ci = df[(df.bmnum == rc["bmnum"]) & (df.slist == rc["gate"])]["v.sprd"]
            ax.fill_between(
                x,
                (y - ci),
                (y + ci),
                color=rc["color"],
                alpha=0.3,
                label="Gate=%02d" % rc["gate"],
            )
        self.rc_ax = ax
        return

    def _add_axis(self):
        self._num_subplots_created += 1
        ax = self.fig.add_subplot(self.num_subplots, 1, self._num_subplots_created)
        if self._num_subplots_created == 1:
            ax.text(
                0.01,
                1.05,
                self.fig_title,
                ha="left",
                va="center",
                transform=ax.transAxes,
            )
        return ax

    def _add_colorbar(self, im, fig, ax, label=""):
        """
        Add a colorbar to the right of an axis.
        """

        pos = ax.get_position()
        cpos = [
            pos.x1 + 0.025,
            pos.y0 + 0.0125,
            0.015,
            pos.height * 0.9,
        ]  # this list defines (left, bottom, width, height
        cax = fig.add_axes(cpos)
        cb = fig.colorbar(im, ax=ax, cax=cax)
        cb.set_label(label)
        return

    def save(self, filepath):
        self.fig.savefig(filepath, bbox_inches="tight")
        return

    def close(self):
        self.fig.clf()
        plt.close()
        return
