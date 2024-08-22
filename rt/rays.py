#!/usr/bin/env python3

"""rays.py: Read Rays, 2D/3D"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"


import numpy as np
import pandas as pd
from loguru import logger
from scipy.io import loadmat


class Rays2D(object):

    def __init__(self, date, rad, beam, elvs, folder, sim_fname):
        """
        Read files by dates and elevation angles
        """
        self.date = date
        self.rad = rad
        self.beam = beam
        self.elvs = elvs
        self.folder = (folder,)
        self.sim_fname = sim_fname
        self.read_file()
        self.ray_power = {"gs": None, "is": None, "all": None}
        self.calc_relative_power(type="gs")
        self.calc_relative_power(type="is")
        return

    def read_file(self):
        """
        Read files
        """
        path_data_keys = [
            "ground_range",
            "height",
            "group_range",
            "phase_path",
            "geometric_distance",
            "electron_density",
            "refractive_index",
        ]
        ray_data_keys = [
            "ground_range",
            "group_range",
            "phase_path",
            "geometric_path_length",
            "initial_elev",
            "final_elev",
            "apogee",
            "gnd_rng_to_apogee",
            "plasma_freq_at_apogee",
            "virtual_height",
            "effective_range",
            "deviative_absorption",
            "TEC_path",
            "Doppler_shift",
            "Doppler_spread",
            "frequency",
            "nhops_attempted",
            "ray_label",
        ]
        self.simulation = dict()
        logger.info(f"Log file: {self.sim_fname}")
        sim_dat = loadmat(self.sim_fname)
        self.ray_path_data, self.ray_data = dict(), []
        for i, elv in enumerate(self.elvs):
            ray_data, path_data = dict(), dict()
            for key in path_data_keys:
                path_data[key] = sim_dat["ray_path_data"][0, i][key].ravel()
            for key in ray_data_keys:
                ray_data[key] = sim_dat["ray_data"][0, i][key].ravel()[0]
            self.simulation[elv] = dict(path_data=path_data, ray_data=ray_data)
            self.ray_data.append(ray_data)
            self.ray_path_data[elv] = pd.DataFrame.from_records(path_data)
        self.ray_data = pd.DataFrame.from_records(self.ray_data)
        return

    def calc_relative_power(self, type):
        """
        Calculate relative power,

        Type: GS - [1]
        """
        pwer = pd.DataFrame()
        if type == "gs":
            labels = [1]
        elif type == "is":
            labels = [-1]
        else:
            [1, 0, -1, -2, -3, -4, -5, -6]
        o = self.ray_data.copy()
        o = o[o.ray_label.isin(labels)]
        # By de Larquier, Sebastien [Thesis]
        # o["weights"] = o.plasma_freq_at_apogee**4/(o.group_range**3) \
        #         if type == "is" else 1./(o.group_range**3)
        o["weights"] = 1.0 / (o.group_range**3)
        ranges = 180 + 45 * np.arange(76, dtype=int)
        lag_power, bins = np.histogram(
            o.group_range,
            bins=ranges,
            weights=o.weights,
        )
        pwer["lag_power"], pwer["srange"], pwer["gate"] = (
            lag_power,
            ranges[:-1],
            range(75),
        )
        pwer["date"], pwer["type"] = self.date, type
        self.ray_power[type] = pwer
        return

    @staticmethod
    def read_rays(event, rad, beam, cfg, folder, sim_fname):
        """
        Static method to read all rays
        """
        elvs = np.linspace(
            float(cfg.start_elevation),
            float(cfg.end_elevation),
            int((cfg.end_elevation - cfg.start_elevation) / cfg.elevation_inctiment)
            + 1,
        )
        rays = Rays2D(event, rad, beam, elvs, folder, sim_fname)
        return rays


import matplotlib.pyplot as plt

plt.style.use(["science", "ieee"])
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans", "Lucida Grande", "Verdana"]
import mpl_toolkits.axisartist.floating_axes as floating_axes
import numpy as np
from matplotlib.projections import polar
from matplotlib.transforms import Affine2D
from mpl_toolkits.axisartist.grid_finder import DictFormatter, FixedLocator


class Plots(object):

    def __init__(self, event, cfg, trace_obj, rad, beam, Re=6371.0):
        self.event = event
        self.cfg = cfg
        self.trace_obj = trace_obj
        self.edens = trace_obj.density
        self.pf = np.sqrt(80.6164e-6 * trace_obj.density)
        self.ref_indx = np.sqrt(1 - (self.pf**2 / cfg.frequency**2))
        self.rad = rad
        self.beam = beam
        self.Re = Re if Re is not None else self.cfg.radius_earth
        return

    def to_polar(self, grange, height):
        th = np.array(grange) / self.Re
        r = np.array(height) + self.Re
        return th, r

    def generate_curvedEarthAxes(
        self,
        nyticks=5,
        nxticks=4,
    ):
        self.fig = plt.figure(figsize=(8, 5), dpi=300)
        maxang = self.cfg.max_ground_range_km / self.Re
        minang = 0 / self.Re
        angran = maxang - minang
        angle_ticks = [(0, "{:.0f}".format(0))]
        while angle_ticks[-1][0] < angran:
            tang = angle_ticks[-1][0] + 1.0 / nxticks * angran
            angle_ticks.append((tang, "{:.0f}".format((tang - minang) * self.Re)))
        grid_locator1 = FixedLocator([v for v, s in angle_ticks])
        tick_formatter1 = DictFormatter(dict(angle_ticks))
        altran = float(self.cfg.end_height_km - 0)
        alt_ticks = [
            (
                0 + self.Re,
                "{:.0f}".format(0),
            )
        ]
        while alt_ticks[-1][0] < self.Re + self.cfg.end_height_km:
            alt_ticks.append(
                (
                    altran / float(nyticks) + alt_ticks[-1][0],
                    "{:.0f}".format(
                        altran / float(nyticks) + alt_ticks[-1][0] - self.Re
                    ),
                )
            )
        _ = alt_ticks.pop()
        grid_locator2 = FixedLocator([v for v, s in alt_ticks])
        tick_formatter2 = DictFormatter(dict(alt_ticks))
        tr_rotate = Affine2D().rotate(np.pi / 2 - maxang / 2)
        _ = Affine2D().translate(0, self.Re)
        tr = polar.PolarTransform() + tr_rotate
        grid_helper = floating_axes.GridHelperCurveLinear(
            tr,
            extremes=(
                0,
                angran,
                self.Re,
                self.Re + self.cfg.end_height_km,
            ),
            grid_locator1=grid_locator1,
            grid_locator2=grid_locator2,
            tick_formatter1=tick_formatter1,
            tick_formatter2=tick_formatter2,
        )
        self.ax = floating_axes.FloatingSubplot(self.fig, 111, grid_helper=grid_helper)

        self.ax.axis["left"].label.set_text(r"Height [km]")
        self.ax.axis["bottom"].label.set_text(r"Ground range [km]")
        self.ax.invert_xaxis()
        self.fig.add_subplot(self.ax, transform=tr)
        # create a parasite axes whose transData in RA, cz
        self.aux_ax = self.ax.get_aux_axes(tr)
        # for aux_ax to have a clip path as in ax
        self.aux_ax.patch = self.ax.patch
        # but this has a side effect that the patch is drawn twice, and possibly
        # over some other artists. So, we decrease the zorder a bit to prevent this.
        self.ax.patch.zorder = 0.9
        return self.ax, self.aux_ax

    def lay_rays(self, kind="pf", zoomed_in=[]):
        self.generate_curvedEarthAxes()
        Th, R = self.to_polar(
            self.trace_obj.bearing_object["dist"],
            self.trace_obj.bearing_object["heights"],
        )
        o, cmap, label, norm = self.get_parameter(kind)
        im = self.aux_ax.pcolormesh(
            Th,
            R,
            o,
            norm=norm,
            cmap=cmap,
            alpha=0.8,
        )
        pos = self.ax.get_position()
        cpos = [
            pos.x1 + 0.025,
            pos.y0 + 0.05,
            0.015,
            pos.height * 0.6,
        ]
        cax = self.fig.add_axes(cpos)
        cbax = self.fig.colorbar(
            im, cax, spacing="uniform", orientation="vertical", cmap="plasma"
        )
        _ = cbax.set_label(label)
        rays = self.trace_obj.rays
        for i, elv in enumerate(rays.elvs):
            ray_path_data, ray_data = (
                rays.ray_path_data[elv],
                rays.simulation[elv]["ray_data"],
            )
            th, r = self.to_polar(ray_path_data.ground_range, ray_path_data.height)
            ray_label = ray_data["ray_label"]
            self.aux_ax.plot(th, r, c="k", zorder=3, alpha=0.7, ls="-", lw=0.1)
            col = "k" if ray_label == 1 else "r"
            if ray_label in [-1, 1]:
                self.aux_ax.scatter([th[-1]], [r[-1]], marker="s", s=3, color=col)
        stitle = "%s UT" % self.event.strftime("%Y-%m-%d %H:%M")
        self.aux_ax.text(
            0.95, 1.01, stitle, ha="right", va="center", transform=self.aux_ax.transAxes
        )
        stitle = f"Model: {self.cfg.model.upper()} / {self.rad}-{'%02d'%self.beam}, $f_0$={self.cfg.frequency} MHz"
        self.aux_ax.text(
            0.05, 1.02, stitle, ha="left", va="center", transform=self.aux_ax.transAxes
        )
        # Create Zoomed in panel
        if len(zoomed_in):
            self.__zoomed_in_panel__(Th, R, kind, zoomed_in)
        return

    def __zoomed_in_panel__(self, Th, R, kind, zoomed_in):
        self.zoom_ax = self.ax.inset_axes([0.4, 1.3, 0.3, 0.5])
        o, cmap, _, norm = self.get_parameter(kind)
        self.zoom_ax.pcolormesh(
            Th,
            R,
            o,
            norm=norm,
            cmap=cmap,
            alpha=0.8,
        )
        rays = self.trace_obj.rays
        for i, elv in enumerate(rays.elvs):
            ray_path_data, ray_data = (
                rays.ray_path_data[elv],
                rays.simulation[elv]["ray_data"],
            )
            th, r = self.to_polar(ray_path_data.ground_range, ray_path_data.height)
            self.zoom_ax.plot(th, r, c="k", zorder=3, alpha=0.7, ls="-", lw=0.05)

        th_lim, r_lim = self.to_polar(zoomed_in[0], zoomed_in[1])
        self.zoom_ax.set_xlim(th_lim)
        self.zoom_ax.set_ylim(r_lim)
        th_ticklabels, r_ticklabels = (
            self.zoom_ax.get_xticklabels(),
            self.zoom_ax.get_yticklabels(),
        )
        for t in th_ticklabels:
            new_x = str(int(t._x * self.Re))
            t.set_text("$\mathdefault{" + new_x + "}$")
        self.zoom_ax.set_xticklabels(th_ticklabels)
        for t in r_ticklabels:
            new_y = str(int(t._y - self.Re))
            t.set_text("$\mathdefault{" + new_y + "}$")
        self.zoom_ax.set_yticklabels(r_ticklabels)

        self.zoom_ax.set_xlabel("Ground Range, [km]", fontdict={"size": 8})
        self.zoom_ax.set_ylabel("Height, [km]", fontdict={"size": 8})
        self.aux_ax.indicate_inset_zoom(self.zoom_ax)
        return

    def get_parameter(self, kind):
        import matplotlib.colors as colors

        if kind == "pf":
            o, cmap, label, norm = (
                getattr(self, kind),
                "plasma",
                r"$f_0$ [MHz]",
                colors.Normalize(
                    self.cfg.ray_trace_plot_lim.pf[0],
                    self.cfg.ray_trace_plot_lim.pf[1]
                ),
            )
        if kind == "edens":
            o, cmap, label, norm = (
                getattr(self, kind),
                "plasma",
                r"$N_e$ [$/m^{-3}$]",
                colors.LogNorm(
                    self.cfg.ray_trace_plot_lim.edens[0],
                    self.cfg.ray_trace_plot_lim.edens[1]
                ),
            )
        if kind == "ref_indx":
            o, cmap, label, norm = (
                getattr(self, kind),
                "plasma",
                r"$\eta$",
                colors.Normalize(
                    self.cfg.ray_trace_plot_lim.ref_indx[0],
                    self.cfg.ray_trace_plot_lim.ref_indx[1]
                ),
            )
        return o, cmap, label, norm

    def save(self, filepath):
        self.fig.savefig(filepath, bbox_inches="tight", facecolor=(1, 1, 1, 1))
        return

    def close(self):
        self.fig.clf()
        plt.close()
        return
