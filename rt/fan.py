#!/usr/bin/env python

"""
    fanUtils.py: module to plot Fan plots with various transformation
"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"


import matplotlib.pyplot as plt
import numpy as np
import scienceplots

plt.style.use(["science", "ieee"])
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans", "Lucida Grande", "Verdana"]

import cartopy
import matplotlib.ticker as mticker
import utils
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER

from rt.sdcarto import SDCarto

## from sdcarto import SDCarto


class Fan(object):
    """
    This class holds plots for all radars FoVs
    """

    def __init__(
        self,
        rads,
        date,
        fig_title=None,
        nrows=1,
        ncols=1,
        coord="geo",
    ):
        self.rads = rads
        self.date = date
        self.nrows, self.ncols = nrows, ncols
        self._num_subplots_created = 0
        self.fig = plt.figure(figsize=(3 * ncols, 3 * nrows), dpi=300)
        self.coord = coord
        self.fig_title = fig_title
        utils.setsize(10)
        return

    def setup(
        self,
        plt_lons,
        plt_lats,
        centers=(-90.0, 45.0),
        extent=[-160, -50, 30, 90],
        proj=None,
        terminator=True,
        lay_eclipse=None,
    ):
        self.plt_lons = plt_lons
        self.plt_lats = plt_lats
        self.centers = centers
        self.extent = extent
        self.proj = proj
        self.terminator = terminator
        self.lay_eclipse = lay_eclipse
        return

    def add_axes(
        self,
    ):
        """
        Instatitate figure and axes labels
        """
        self._num_subplots_created += 1
        proj = (
            self.proj
            if self.proj
            else cartopy.crs.Stereographic(
                central_longitude=self.centers[0], central_latitude=self.centers[1]
            )
        )
        ax = self.fig.add_subplot(
            100 * self.nrows + 10 * self.ncols + self._num_subplots_created,
            projection="SDCarto",
            map_projection=proj,
            coords=self.coord,
            plot_date=self.date,
        )
        ax.overaly_coast_lakes(lw=0.4, alpha=0.4)
        ax.set_extent(self.extent, crs=cartopy.crs.PlateCarree())
        plt_lons = self.plt_lons
        plt_lats = self.plt_lats
        gl = ax.gridlines(crs=cartopy.crs.PlateCarree(), linewidth=0.2)
        gl.xlocator = mticker.FixedLocator(plt_lons)
        gl.ylocator = mticker.FixedLocator(plt_lats)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.n_steps = 90
        self.proj = proj
        self.geo = cartopy.crs.PlateCarree()
        if self.lay_eclipse:
            ax.overaly_eclipse_path(self.lay_eclipse, lineWidth=0.2)
        if self.terminator:
            from cartopy.feature.nightshade import Nightshade

            ax.add_feature(Nightshade(self.date, alpha=1))
        return ax

    def annotate_figure(self, axis_num=0):
        # Annotate the generic info in axis
        utils.setsize(10)
        ax = self.fig.get_axes()[axis_num]
        ax.text(
            -0.02,
            0.99,
            "Coord: Geo",
            ha="center",
            va="top",
            transform=ax.transAxes,
            fontsize="small",
            rotation=90,
        )
        ax.text(
            0.01,
            1.05,
            (
                f"{self.date_string()} / {self.fig_title}"
                if self.fig_title
                else f"{self.date_string()}"
            ),
            ha="left",
            va="center",
            transform=ax.transAxes,
        )
        return

    def date_string(self, label_style="web"):
        # Set the date and time formats
        dfmt = "%d %b %Y" if label_style == "web" else "%d %b %Y,"
        tfmt = "%H:%M"
        stime = self.date
        date_str = "{:{dd} {tt}} UT".format(stime, dd=dfmt, tt=tfmt)
        return date_str

    def generate_fov(
        self,
        rad,
        frame,
        beams=[],
        ax=None,
        maxGate=45,
        text_decription=dict(x=0.1, y=0.9, txt="", ha="left", va="center"),
        col="k",
        p_name="vel",
        p_max=30,
        p_min=-30,
        cmap="Spectral",
        label="Velocity [m/s]",
        cbar=True,
        lats=np.linspace(0, 90, num=90 * 2),
        lons=np.linspace(-180, 180, num=91 * 2),
        alts=np.array([100]),
    ):
        """
        Generate plot with dataset overlaid
        """
        utils.setsize(10)
        ax = ax if ax else self.add_axes()
        ax.text(
            text_decription["x"],
            text_decription["y"],
            text_decription["txt"],
            ha=text_decription["ha"],
            va=text_decription["va"],
            transform=ax.transAxes,
        )
        ax.overlay_radar(rad, font_color=col)
        ax.overlay_fov(rad, lineColor=col)
        if len(frame) > 0:
            ax.overlay_data(
                rad,
                frame,
                self.proj,
                maxGate=maxGate,
                p_name=p_name,
                p_max=p_max,
                p_min=p_min,
                cmap=cmap,
                label=label,
                cbar=cbar,
            )
        if beams and len(beams) > 0:
            [
                ax.overlay_fov(
                    rad, beamLimits=[b, b + 1], ls="-", lineColor="r", lineWidth=1.2
                )
                for b in beams
            ]
        if self.lay_eclipse:
            ax.overlay_eclipse(lats, lons, alts)
        return

    def generate_fovs(self, fds, beams=[], laytec=False):
        """
        Generate plot with dataset overlaid
        """
        ax = self.add_axes()
        for rad in self.rads:
            self.generate_fov(
                rad, fds[rad].frame, beams, ax, laytec, col=fds[rad].color
            )
        return ax

    def overlay_fovs(self, rad, beams=[], ax=None, col="k"):
        """
        Generate plot with dataset overlaid
        """
        ax = ax if ax else self.add_axes()
        ax.overlay_radar(rad, font_color=col)
        ax.overlay_fov(rad, lineColor=col)
        if beams and len(beams) > 0:
            [
                ax.overlay_fov(
                    rad,
                    beamLimits=[b, b + 1],
                    ls="-",
                    lineColor="darkred",
                    lineWidth=0.8,
                )
                for b in beams
            ]
        return ax

    def generate_ham_fov(self, ax, source: dict, stations: list = []):
        ax.overlay_point(
            source["lat"],
            source["lon"],
            source["call_sign"],
            marker="D",
            markerColor="b",
            markerSize=4,
        )
        for stn in stations:
            ax.overlay_point(stn["lat"], stn["lon"], stn["call_sign"], markerColor="r")
            xy = self.proj.transform_points(
                self.geo,
                np.array([stn["lon"], source["lon"]]),
                np.array([stn["lat"], source["lat"]]),
            )
            ax.plot(xy[:, 0], xy[:, 1], color="k", lw=0.5, ls="-", alpha=0.8)
        return ax

    def save(self, filepath):
        self.fig.savefig(filepath, bbox_inches="tight", facecolor=(1, 1, 1, 1))
        return

    def close(self):
        self.fig.clf()
        plt.close()
        return
