#!/usr/bin/env python3

"""rt2D.py: simulate python program for RT"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import copy
import datetime as dt
import os

import eclipse
import numpy as np
import utils
from geopy.distance import great_circle as GC
from loguru import logger
from radar import Radar
from rays import Rays2D
from scipy.io import loadmat, savemat


class RadarBeam2dTrace(object):
    """
    Ray trace class to trace all the points
    """

    def __init__(
        self,
        event,
        rad,
        beam,
        cfg,
        model,
        base_output_folder,
    ):
        self.beam = beam
        self.event = event
        self.model = model
        self.rad = rad
        self.base_output_folder = base_output_folder
        self.folder = utils.get_folder(rad, beam, event, model, base_output_folder)
        os.makedirs(self.folder, exist_ok=True)
        self.cfg = cfg
        self.radar = Radar(self.rad, [event, event + dt.timedelta(minutes=5)], cfg)
        self.fig_name = self.folder + "/{time}.png".format(
            time=self.event.strftime("%H%M")
        )
        self.edensity_file = self.folder + "/{dn}.mat".format(
            dn=self.event.strftime("%H.%M")
        )
        self.sim_fname = self.folder + "/{date}_rt.mat".format(
            date=self.event.strftime("%H%M")
        )
        self._estimate_bearing_()
        self._eclipse()
        return

    def _eclipse(self):
        self.p = (
            np.zeros((len(self.bearing_object["lat"]), len(self.bearing_object["ht"])))
            * np.nan
        )
        h = np.mean(self.bearing_object["ht"])
        e = eclipse.Eclipse()
        for i, lat, lon in zip(
            range(len(self.bearing_object["lat"])),
            self.bearing_object["lat"],
            self.bearing_object["lon"],
        ):
            p = e.create_eclipse_shadow(self.event, lat, lon, h)
            if p == 0:
                p = np.nan
            self.p[i, :] = p
        return

    def _estimate_bearing_(self):
        """Estimate laitude and logitude bearings"""
        fname = self.folder + f"/bearing.mat"
        lat, lon = (self.radar.hdw.geographic.lat, self.radar.hdw.geographic.lon)
        p = (lat, lon)
        _, bearing = utils.calculate_bearing(
            lat,
            lon,
            self.radar.fov[0][self.beam, 0],
            self.radar.fov[1][self.beam, 0],
        )
        logger.info(f"Bearing angle of beam {self.beam} is {bearing} deg")
        bearing_object = {}
        lats, lons = (
            self.radar.fov[0][: self.cfg.slant_gate_of_radar, self.beam],
            self.radar.fov[1][: self.cfg.slant_gate_of_radar, self.beam],
        )
        dist = np.array([GC(p, (latx, lonx)).km for latx, lonx in zip(lats, lons)])
        dist, lats, lons = (
            np.insert(dist, 0, 0),
            np.insert(lats, 0, lat),
            np.insert(lons, 0, lon),
        )
        bearing_object["dist"], bearing_object["lat"], bearing_object["lon"] = (
            dist,
            np.array(lats),
            np.array(lons),
        )
        (
            bearing_object["olat"],
            bearing_object["olon"],
            bearing_object["rb"],
            bearing_object["num_range"],
            bearing_object["max_range"],
            bearing_object["range_inc"],
        ) = (
            lat,
            lon,
            bearing,
            float(len(dist)),
            float(self.cfg.max_ground_range_km),
            float(dist[1] - dist[0]),
        )
        bearing_object["ht"] = np.arange(
            self.cfg.start_height_km,
            self.cfg.end_height_km,
            self.cfg.height_incriment_km,
        ).astype(float)
        (
            bearing_object["start_height"],
            bearing_object["height_inc"],
            bearing_object["num_heights"],
            bearing_object["heights"],
        ) = (
            float(self.cfg.start_height_km),
            float(self.cfg.height_incriment_km),
            float(
                len(
                    np.arange(
                        self.cfg.start_height_km,
                        self.cfg.end_height_km,
                        self.cfg.height_incriment_km,
                    )
                )
            ),
            np.arange(
                self.cfg.start_height_km,
                self.cfg.end_height_km,
                self.cfg.height_incriment_km,
            ),
        )

        bearing_object["freq"], bearing_object["tol"], bearing_object["nhops"] = (
            float(self.cfg.frequency),
            float(1e-7),
            float(self.cfg.nhops),
        )
        bearing_object["elev_s"], bearing_object["elev_i"], bearing_object["elev_e"] = (
            float(self.cfg.start_elevation),
            float(self.cfg.elevation_inctiment),
            float(self.cfg.end_elevation),
        )
        bearing_object["radius_earth"] = self.cfg.radius_earth
        savemat(fname, bearing_object)
        self.bearing_object = copy.copy(bearing_object)
        return

    def read_density_rays(self, fname):
        self.density = loadmat(fname)["ne"]
        self.sim_fname = self.folder + "{date}.{bm}_rt.mat".format(
            bm="%02d" % self.beam, date=self.event.strftime("%H%M")
        )
        logger.info("Data-Model comparison: reading rays....")
        self.rays = Rays2D.read_rays(
            self.event, self.rad, self.beam, self.cfg, self.folder, self.sim_fname
        )
        return

    def compile(self, density):
        """Compute RT using Pharlap"""
        self.density = density
        print(self.folder, os.getcwd() + "/pharlap/pharlap_4.5.3/dat")
        pwd = os.getcwd() + "/pharlap/pharlap_4.5.3/dat"
        cmd = "export DIR_MODELS_REF_DAT={pwd};\
                cd pharlap/;\
                matlab -softwareopengl -nodisplay -nodesktop -nosplash -nojvm -r \"UT=[{ut}];rad='{rad}';dic='{dic}';fname='{fname}';bm={bm};\
                rt_2D;exit;\"".format(
            pwd=pwd,
            ut=self.event.strftime("%Y %m %d %H %M"),
            rad=self.rad,
            dic=self.folder,
            bm=self.beam,
            fname=self.sim_fname,
        )
        logger.info(f"Running command: {cmd}")
        os.system(cmd)
        logger.info("Data-Model comparison: reading rays....")
        self.rays = Rays2D.read_rays(
            self.event, self.rad, self.beam, self.cfg, self.folder, self.sim_fname
        )
        return

    def load_rto(self, density):
        self.density = density
        logger.info("Data-Model comparison: reading rays....")
        self.rays = Rays2D.read_rays(
            self.event, self.rad, self.beam, self.cfg, self.folder, self.sim_fname
        )
        return
