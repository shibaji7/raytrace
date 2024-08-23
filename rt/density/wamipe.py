import datetime as dt
import glob
import math
import os

import h5py
import numpy as np
from loguru import logger
from scipy.io import loadmat, savemat


class WAMIPE2d(object):
    def __init__(
        self,
        cfg,
        event,
    ):
        self.cfg = cfg
        self.file_name = self.cfg.density_file_location
        self.event = event
        self.file_name = os.path.join(
            self.cfg.density_file_location, self.cfg.density_file_name
        )
        self.load_grid()
        self._load_IPE_edensity_files_()
        return

    def load_grid(self):
        """
        Load grid file
        """
        logger.info(f"Loading WAM Grid files")
        self.ccord_file = self.cfg.density_file_location + self.cfg.grid_coordinate_file

        # number of grid points in magnetic longitude in the IPE model
        nmp = self.cfg.wam_paramteres.coordinates.nmp
        # number of magnetic flux tubes in the IPE model from the north pole.
        nlp = self.cfg.wam_paramteres.coordinates.nlp
        # number of grid points along a magnetic flux tube in the IPE model from the northern hemisphere at 90km.
        iDIM = self.cfg.wam_paramteres.coordinates.iDIM
        # 1:geographic eastward; 2:northward; 3:upward
        napex = self.cfg.wam_paramteres.coordinates.napex

        m2km = 1.0e-3
        rad2deg = 180.0 / math.pi
        # Retrieve altitude in meters
        self.altitude = (
            self._read_params_from_hdf_file_(
                self.ccord_file,
                self.cfg.wam_paramteres.coordinates.grid_name,
                "altitude",
            )
            * m2km
        )  # Convert meter to km

        # Retrieve geographic longitude in radian
        self.longitude = (
            self._read_params_from_hdf_file_(
                self.ccord_file,
                self.cfg.wam_paramteres.coordinates.grid_name,
                "longitude",
            )
            * rad2deg
        )  # Convert meter to deg
        self.longitude = np.mod((self.longitude + 180), 360) - 180

        # Retrieve geographic colatitude in radian and convert to latitude
        self.latitude = 90.0 - (
            self._read_params_from_hdf_file_(
                self.ccord_file,
                self.cfg.wam_paramteres.coordinates.grid_name,
                "colatitude",
            )
            * rad2deg  # Convert meter to deg
        )
        print(
            self.latitude.min(),
            self.latitude.max(),
            self.altitude.min(),
            self.altitude.max(),
            self.longitude.min(),
            self.longitude.max(),
        )
        print(self.altitude.shape, self.latitude.shape)
        print(self.longitude[:, 0, 0].tolist())
        print(self.longitude[:, 1, 0].tolist())
        raise Exception("Null recieced")
        return

    def find_index(self, lat, lon, alt):
        i_lat, i_lon, i_alt = np.nan, np.nan, np.nan

        return

    def _load_IPE_edensity_files_(self):
        """
        Load ipe-density file
        """
        logger.info(f"Loading ipe density files")
        ccord_file = self.cfg.density_file_location + self.cfg.density_file_name
        self.files = glob.glob(ccord_file)
        self.files.sort()
        self.dates = [
            dt.datetime.strptime((f.split("/")[-1]).split(".")[-2], "%Y%m%d%H%M")
            for f in self.files
        ]
        return

    def _read_params_from_hdf_file_(self, fpath, ds_name, p_name):
        """Load parameters from file"""
        dataset = None
        with h5py.File(fpath, "r") as file:
            # access datasets and attributes in the file
            logger.info(f"Access attribute {p_name} under dataset {ds_name}")
            dataset = file[ds_name][p_name][:]
        return dataset

    def fetch_dataset(
        self,
        time: dt.datetime,
        lats,
        lons,
        alts,
        to_file: str = None,
        dlat: float = 0.2,
        dlon: float = 0.2,
    ):
        self.param = np.zeros((len(alts), len(lats)))
        i = np.argmin([np.abs((t - time).total_seconds()) for t in self.dates])
        file = self.files[i]
        el_den = self.load_data(file) * 1e-6  # Convert to cub-m to cub-c
        print(file, el_den.shape, el_den.min(), el_den.max())
        raise Exception("Null recieced")
        if to_file:
            savemat(to_file, dict(ne=self.param))
        return self.param, alts

    def load_data(self, fname):
        # array dimensions
        el_den = np.zeros(
            (
                self.cfg.wam_paramteres.coordinates.nmp,  # number of grid points in magnetic longitude in the IPE model
                self.cfg.wam_paramteres.coordinates.nlp,  # number of magnetic flux tubes in the IPE model from the north pole.
                self.cfg.wam_paramteres.coordinates.iDIM,  # number of grid points along a magnetic flux tube in the IPE model from the northern hemisphere at 90km.
            )
        )
        for dp in self.cfg.wam_paramteres.density_params:
            el_den[:, :, :] += self._read_params_from_hdf_file_(
                fname, self.cfg.wam_paramteres.dataset_name, dp
            )[:, :, :]
        return el_den

    def load_from_file(self, to_file: str):
        logger.info(f"Load from file {to_file.split('/')[-1]}")
        self.param = loadmat(to_file)["ne"]
        return self.param
