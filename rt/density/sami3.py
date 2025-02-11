import datetime as dt
import os

import numpy as np
import pandas as pd
import utils
import xarray as xr
from dateutil import parser as dparser
from loguru import logger
from scipy.io import loadmat, savemat


class SAMI3(object):
    def __init__(
        self,
        cfg,
        event,
        round_to=-1,
    ):
        self.cfg = cfg
        self.event = event
        self.file_name = os.path.join(
            self.cfg.density_file_location, self.cfg.density_file_name
        )
        self.cfg.density_simulated_datetime = dparser.isoparse(
            self.cfg.density_simulated_datetime
        )
        self.round_to = round_to
        self.load_nc_dataset()
        return

    def load_nc_dataset(self):
        self.store = {}
        logger.info(f"Load files -> {self.file_name}")
        ds = xr.open_dataset(self.file_name)
        self.store["time"] = [
            self.cfg.density_simulated_datetime + dt.timedelta(hours=float(d))
            for d in ds.variables["time"][:]
        ]
        if self.round_to > 0:
            logger.info("Round to nearest minutes!!!!!!!!!!!!")
            self.store["time"] = [self.round_time(e) for e in self.store["time"]]
        self.dates = self.store["time"]
        self.store["alt"] = ds.variables["alt0"].values  # in km
        self.store["glat"] = ds.variables["lat0"].values  # in deg
        self.store["glon"] = ds.variables["lon0"].values  # in 0-360 deg
        # self.store["glon"] = (
        #     np.mod(self.store["glon"] - 180.0, 360.0) - 180.0
        # )  # to +/- 180
        self.store["eden"] = ds.variables["dene0"][:]  # in cc
        ds.close()
        del ds
        return

    def round_time(self, t):
        """Round a datetime object to any time lapse in seconds
        t : datetime.datetime object, default now.
        """
        seconds = (t.replace(tzinfo=None) - t.min).seconds
        rounding = (seconds + self.round_to / 2) // self.round_to * self.round_to
        return t + dt.timedelta(0, rounding - seconds, -t.microsecond)

    def find_time_index(self, t):
        """Finds the index of the interval in the array where the number falls.

        Returns:
            The index of the interval where the number falls, or -1 if the number is
            outside the range of the array.
        """
        for i in range(len(self.store["time"]) - 1):
            if self.store["time"][i] <= t < self.store["time"][i + 1]:
                return i, i + 1
        return -1, 0

    def fetch_dataset(
        self,
        time,
        lats,
        lons,
        alts,
        to_file=None,
    ):
        # Selecting based on time index
        if time in self.store["time"]:
            logger.info(f"Into exact timestamp {time}")
            # Select the exact index if timestamp in the simulation
            i = self.store["time"].index(time)
            self.param, self.alts = self.fetch_interpolated_data(lats, lons, alts, i)
        else:
            # Select the two index if timestamp in the simulation
            i, j = self.find_time_index(time)
            logger.info(
                f"{time}/ into between timestamp {self.store['time'][i]} & {self.store['time'][j]}"
            )
            weights = (self.store["time"][1] - self.store["time"][0]).total_seconds()
            px, _ = self.fetch_interpolated_data(lats, lons, alts, i)
            py, self.alts = self.fetch_interpolated_data(lats, lons, alts, j)
            i_wg, j_wg = (
                (time - self.store["time"][i]).total_seconds() / weights,
                (self.store["time"][j] - time).total_seconds() / weights,
            )
            self.param = px * i_wg + py * j_wg

        if to_file:
            savemat(to_file, dict(ne=self.param))
        return self.param, self.alts

    def fetch_interpolated_data(
        self,
        lats,
        lons,
        alts,
        index,
    ):
        n = len(lats)
        D = self.store["eden"][index]
        glat = self.store["glat"]
        glon = self.store["glon"]
        out, ix = np.zeros((len(alts), n)) * np.nan, 0
        for lat, lon in zip(lats, lons):
            lon = np.mod(360 + lon, 360)
            id_lat = np.argmin(np.abs(glat - lat))
            id_lon = np.argmin(np.abs(glon - lon))
            o = D[id_lon, :, id_lat]
            galt = self.store["alt"]
            out[:, ix] = utils.interpolate_by_altitude(
                galt, alts, o, self.cfg.scale, self.cfg.kind, method="extp"
            )  # in cc
            ix += 1
        return out, galt

    def load_from_file(self, to_file: str):
        logger.info(f"Load from file {to_file.split('/')[-1]}")
        self.param = loadmat(to_file)["ne"]
        return self.param
