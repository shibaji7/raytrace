import datetime as dt
import os
import numpy as np
import xarray as xr
from loguru import logger
from scipy import io
from scipy.io import savemat
import utils


class GITM(object):
    def __init__(
        self,
        cfg,
        event,
        file_path,
    ):
        self.cfg = cfg
        self.file_name = file_path
        self.event = event
        self.param = param
        self.load_nc_dataset()
        return

    def load_nc_dataset(self):
        """
        Load netcdf4 dataset available
        """
        self.store = {}
        drop_vars = [
            "wn",
            "vn",
            "un",
            "tn",
            "denn",
            "TEC",
            "wi",
            "vi",
            "ui",
            "SigH",
            "SigP",
            "time",
        ]
        file = self.file_name
        logger.info(f"Load files -> {file}")
        ds = xr.open_dataset(file, drop_variables=drop_vars)
        self.store["time"] = [
            dt.datetime(y, m, d, h, mm)
            for y, m, d, h, mm in zip(
                ds.year.values,
                ds.month.values,
                ds.day.values,
                ds.hour.values,
                ds.minute.values,
            )
        ]
        (
            self.store["glat"],
            self.store["glon"],
            self.store["alt"],
            self.store["eden"],
        ) = (
            ds.glat.values,
            np.mod(360 + ds.glon.values, 360),
            ds.alt.values / 1e3,
            ds.dene.values,
        )
        logger.info(f"Shape of eden: {ds.dene}")
        ds.close()
        del ds
        return

    def fetch_nc_dataset(
        self,
        time,
        lats,
        lons,
        alts,
        to_file=None,
    ):
        i = np.argmin([np.abs((t - time).total_seconds()) for t in self.store["time"]])
        n = len(lats)
        D = self.store["eden"][i]
        glat = self.store["glat"]
        glon = self.store["glon"]
        galt = self.store["alt"]
        out, ix = np.zeros((len(alts), n)) * np.nan, 0
        for lat, lon in zip(lats, lons):
            lon = np.mod(360 + lon, 360)
            idx = np.argmin(np.abs(glon - lon))
            idy = np.argmin(np.abs(glat - lat))
            o = D[:, idx, idy]
            out[:, ix] = (
                utils.interpolate_by_altitude(
                    galt, alts, o, self.cfg.scale, self.cfg.kind, method="extp"
                )
                * 1e-6
            )
            ix += 1
        self.param, self.alts = out, galt
        if to_file:
            savemat(to_file, dict(ne=self.param))
        return out, galt

    def load_dataset(self):
        """
        Load all dataset available
        """
        files = np.array(sorted(glob.glob(f"{self.folder}*.sav")))
        self.file_dates = np.array(
            [
                dt.datetime.strptime(os.path.split(i)[1][7:-4], "%y%m%d_%H%M%S")
                for i in files
            ]
        )
        self.dataset = {}
        for i, f in enumerate(files):
            self.file_dates[i] = self.file_dates[i].replace(year=self.year)
            D = io.readsav(f)
            self.dataset[self.file_dates[i]] = D
            logger.info(f"Loading file: {f}")
        return