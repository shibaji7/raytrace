import datetime as dt
import os
import numpy as np
import xarray as xr
from loguru import logger
from scipy import io
from scipy.io import savemat


class GITM(object):
    def __init__(
        self,
        cfg,
        event,
        file_path,
        param="dene",
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
        file = os.path.join(self.folder, self.file_name)
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

    def fetch_dataset(self, time, latlim, lonlim):
        """
        Fetch data by lat/lon limits
        """
        i = np.argmin([np.abs((t - time).total_seconds()) for t in self.file_dates])
        D = self.dataset[self.file_dates[i]]
        glat = (D["lat1"] * 180 / np.pi)[:, 0]
        glon = (D["lon1"][0, :] * 180) / np.pi
        idx = np.logical_and(glon >= lonlim[0], glon <= lonlim[1])
        idy = np.logical_and(glat >= latlim[0], glat <= latlim[1])
        galt = np.squeeze(D[self.param]) / 1e3
        out = D[self.param][:, idx, :][:, :, idy]
        return out

    def fetch_dataset_by_locations(self, time, lats, lons, alts):
        """
        Fetch data by lat/lon limits
        """
        n = len(lats)
        i = np.argmin([np.abs((t - time).total_seconds()) for t in self.file_dates])
        logger.info(f"Time: {self.file_dates[i]}")
        D = self.dataset[self.file_dates[i]]
        glat = (D["lat1"] * 180 / np.pi)[:, 0]
        glon = (D["lon1"][0, :] * 180) / np.pi
        galt = D["alt1"][:, 0, 0] / 1e3
        out, ix = np.zeros((len(alts), n)) * np.nan, 0
        for lat, lon in zip(lats, lons):
            lon = np.mod(360 + lon, 360)
            idx = np.argmin(np.abs(glon - lon))
            idy = np.argmin(np.abs(glat - lat))
            o = D[self.param][:, idx, idy]
            out[:, ix] = (
                utils.interpolate_by_altitude(
                    galt, alts, o, self.cfg.scale, self.cfg.kind, method="extp"
                )
                * 1e-6
            )
            ix += 1
        self.param, self.alts = out, galt
        return out, galt

    @staticmethod
    def create_object(
        cfg,
        year=2017,
        folder="dataset/GITM/20170821/",
        param="eden",
        time=dt.datetime(2017, 8, 21, 17, 30),
        lats=[],
        lons=[],
        alts=[],
        to_file=None,
        prev=False,
    ):
        mobj = {}
        gitm = GITM(cfg, year, folder, param)
        mobj["ne"], _ = gitm.fetch_dataset_by_locations(time, lats, lons, alts)
        if to_file:
            savemat(to_file, mobj)
        return gitm

    @staticmethod
    def create_density_files(
        cfg,
        year=2017,
        folder="dataset/GITM/20170821/",
        param="eden",
        time=dt.datetime(2017, 8, 21, 17, 30),
        lats=[],
        lons=[],
        alts=[],
        to_file=None,
        prev=None,
    ):
        gitm = GITM.create_object(
            cfg,
            year,
            folder,
            param,
            time,
            lats,
            lons,
            alts,
            to_file,
        )
        if prev:
            GITM.create_object(
                cfg,
                year,
                folder,
                param,
                time - dt.timedelta(minutes=30),
                lats,
                lons,
                alts,
                prev,
            )
        return gitm