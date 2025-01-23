#!/usr/bin/env python

"""hamsci.py: module is dedicated to fetch HamSci database."""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = ""
__status__ = "Research"

import datetime as dt
import glob
import os
import sys

import cartopy
import numpy as np
import pytz
from dateutil import parser as dparser
from hamsci_psws import grape1
from loguru import logger

sys.path.extend([".", "rt/", "rt/density/", "Projects/SCurve/"])
import utils
from fan import Fan


class HamSci(object):
    """
    This class is help to extract the dataset from HamSci database and plot.
    """

    def __init__(self, cfg, event, dates, freq=10.0e6):
        """
        Parameters:
        -----------
        base: Base location
        freq: Frequency of operation in MHz (float)
        dates: Start and end dates
        """
        self.freq = freq
        self.event = event
        self.dates = self.parse_dates(dates)
        self.date_range = [
            dates[0].replace(tzinfo=pytz.utc),
            dates[1].replace(tzinfo=pytz.utc),
        ]
        self.cfg = cfg
        self.base_output_folder = os.path.join(
            self.cfg.project_save_location, self.cfg.project_name
        )
        self.node_registry = cfg.psws_node_registry
        self.hamsci_file_loc = cfg.tmp_hamsci_folder
        logger.info("Loging into local files")
        return

    def parse_dates(self, dates):
        """
        Parsing dates
        """
        da = [
            dates[0].replace(minute=0, hour=0, second=0),
            dates[1].replace(minute=0, hour=0, second=0),
        ]
        return da

    def load_files(self, freq=None):
        """
        Load files using grape1 library
        """
        freq = freq if freq else self.freq
        if self.hamsci_file_loc.startswith("~/"):
            self.hamsci_file_loc = self.hamsci_file_loc.replace("~", os.getenv("HOME"))
        logger.info(f"Searching inside: {self.hamsci_file_loc}*.csv")
        files = glob.glob(f"{self.hamsci_file_loc}*.csv")
        logger.info(f"Load files from {self.hamsci_file_loc}; N-files:{len(files)}")
        if len(files) > 0:
            inv = grape1.DataInventory(data_path=self.hamsci_file_loc)
            inv.filter(
                freq=freq,
                sTime=self.date_range[0],
                eTime=self.date_range[1],
            )
            gn = grape1.GrapeNodes(
                fpath=self.node_registry, logged_nodes=inv.logged_nodes
            )
            return inv, gn
        else:
            return None, None

    def setup_pandas_dataset(
        self,
        freq=None,  # Frequency in Hz
    ):
        """
        Plot dataset in multipoint plot
        """
        logger.info("Fecth files and store to pandas")
        freq = freq if freq else self.freq
        self.gds = {}
        inv, gn = self.load_files(freq)

        if inv:
            node_nrs = inv.get_nodes()
            logger.info(f"Nodes {len(node_nrs)}")
            for node in node_nrs:
                try:
                    logger.info(f"Node number: {node}")
                    gd = grape1.Grape1Data(
                        node,
                        freq,
                        self.date_range[0],
                        self.date_range[1],
                        inventory=inv,
                        grape_nodes=gn,
                        data_path=self.hamsci_file_loc,
                    )
                    gd.process_data()
                    self.gds[gd.meta["call_sign"]] = gd
                except:
                    import traceback

                    traceback.print_exc()
                    logger.error("Sorry, issue occured!")
        return self.gds

    def generate_fov(
        self,
        event: dt.datetime,
        cfg_file: str = "cfg/rt2d_iri_2024_iri_hamsci_SCurve.json",
        call_signs: list = ["w2naf"],
        source: dict = dict(call_sign="wwv", lat=40.6776, lon=-105.0461),
        fname: str = None,
        lons: np.array = None,
        lats: np.array = None,
        extent: np.array = None,
        proj=None,
    ):
        event = event if event else self.event
        logger.info(f"Generate FoV of the GC distances, {event}")
        stations = []
        for call_sign in call_signs:
            meta = self.gds[call_sign.upper()].meta
            stations.append(
                dict(
                    lat=meta["lat"],
                    lon=meta["lon"],
                    call_sign=call_sign,
                )
            )

        lons = np.arange(-180, 180, 30) if lons is None else lons
        lats = np.arange(30, 70, 15) if lats is None else lats
        proj = cartopy.crs.Orthographic(-100, 30) if proj is None else proj
        extent = [-130, -60, 20, 70]
        f = Fan(source["call_sign"], event, fig_title="")
        f.setup(
            lons,
            lats,
            extent=extent,
            proj=proj,
            lay_eclipse=None,
        )
        ax = f.add_axes()
        f.generate_ham_fov(ax, source, stations)
        ax.overaly_eclipse_path(cfg_file, year=self.dates[0].year)
        ax.overlay_eclipse(
            np.arange(0, 90, 1), np.arange(-130, -20, 1), [300], xpos=1.01
        )
        fname = (
            fname
            if fname
            else utils.get_hamsci_folder(
                "wwv",
                event,
                self.cfg.model,
                self.base_output_folder,
            )
            + "/fov.png"
        )
        logger.info(f"Save to {fname}")
        f.annotate_figure()
        f.save(fname)
        f.close()
        return


if __name__ == "__main__":
    cfg = utils.read_params_2D("cfg/rt2d_iri_2024_iri_hamsci_SCurve.json")
    h = HamSci(
        cfg,
        dparser.isoparse(cfg.event),
        [dt.datetime(2024, 4, 8), dt.datetime(2024, 4, 9)],
    )
    h.setup_pandas_dataset()
    h.generate_fov(dparser.isoparse(cfg.event))
