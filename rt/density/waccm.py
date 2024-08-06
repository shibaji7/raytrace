from loguru import logger
import numpy as np
import os
import xarray as xr
import utils
from scipy.io import loadmat, savemat

pconst = {
    "boltz": 1.38066e-23,  # Boltzmann constant  in Jule K^-1
    "h": 6.626e-34,  # Planks constant  in ergs s
    "c": 2.9979e08,  # in m s^-1
    "avo": 6.023e23,  # avogadro's number
    "Re": 6371.0e3,
    "amu": 1.6605e-27,
    "q_e": 1.602e-19,  # Electron charge in C
    "m_e": 9.109e-31,  # Electron mass in kg
    "g": 9.81,  # Gravitational acceleration on the surface of the Earth
    "eps0": 1e-9 / (36 * np.pi),
    "R": 8.31,  # J mol^-1 K^-1
}

class WACCMX2d(object):

    def __init__(
        self,
        cfg,
        event,
    ):
        self.cfg = cfg
        self.event = event
        self.file_name = self.cfg.density_file_location
        self.load_nc_dataset()
        return

    def transform_density(self, P, T, var, unit="cm"):
        """
        Transform electron density from mol/mol to /cc or /cm
        Formula
        -------
        P = nKT, K = Boltzmann constant in SI
        lev: pressure in hPa (1hPa=100Pa)
        T: neutral temperature
        e/O2/...: density in mol/mol

        T in K, P in Pa, n is electron density in /cubic meter or /cc
        """
        logger.info(f"Transform density scale/unit")
        P = P * 1e2
        u = 1 if unit == "cm" else 1e-6
        den = np.zeros_like(var)
        for i, p in enumerate(P):
            na = p * pconst["avo"] / (pconst["R"] * T[:, i, :, :])
            n = u * na * var[:, i, :, :]
            den[:, i, :, :] = n
        return den

    def load_nc_dataset(self):
        """
        Load netcdf4 dataset available
        """
        self.store = {}
        drop_vars = [
            "WI",
            "VI",
            "UI",
            "TIon",
            "TElec",
            "Op_CHMP",
            "Op_CHML",
            "OpDens",
            "Op",
            "NOp",
            "O2p",
            "N2p",
            "SOLAR_MASK",
            "amb_diff",
            "dwind",
            "dfield",
            "op_dt",
        ]
        file = os.path.join(self.folder, self.file_name)
        logger.info(f"Load files -> {file}")
        ds = xr.open_dataset(file, drop_variables=drop_vars)
        self.store["glat"], self.store["glon"] = (
            ds.lat.values,
            np.mod(360 + ds.lon.values, 360),
        )
        self.store["time"] = [
            self.event.replace(hour=0) + dt.timedelta(seconds=float(d))
            for d in ds.variables["time"][:]
        ]
        self.store["alt"] = ds.variables["Z3GM"][:] / 1e3
        self.store["eden"] = ds.variables["e"][:]
        self.store["T"] = ds.variables["T"][:]
        self.store["lev"] = ds.variables["lev"][:]
        self.store["eden"] = self.transform_density(
            self.store["lev"], self.store["T"], self.store["eden"]
        )
        logger.info(f"Shape of eden: {self.store['eden'].shape}")
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
        out, ix = np.zeros((len(alts), n)) * np.nan, 0
        for lat, lon in zip(lats, lons):
            lon = np.mod(360 + lon, 360)
            idx = np.argmin(np.abs(glat - lat))
            idy = np.argmin(np.abs(glon - lon))
            o = D[:, idx, idy]
            galt = self.store["alt"][i, :, idx, idy]
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

    def load_from_file(self, to_file: str):
        logger.info(f"Load from file {to_file.split('/')[-1]}")
        self.param = loadmat(to_file)["ne"]
        return self.param