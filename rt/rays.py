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
