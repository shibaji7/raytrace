#!/usr/bin/env python3

"""calculate_doppler.py: simulate python program for RT and compute doppler"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import datetime as dt
import glob
import os
from types import SimpleNamespace

import numpy as np
import pandas as pd
from loguru import logger
from scipy.interpolate import RectBivariateSpline
from scipy.io import loadmat, savemat

import rt.utils as utils
from rt.radar import Radar
from rt.rays import Rays2D


class Doppler(object):
    def __init__(
        self,
        cfg: SimpleNamespace,
        start_time: dt.datetime,
        model: str,
        radar: Radar,
        source: dict,
        target: dict,
        beam: int = 0,
        del_time: int = 1,  # Normalized based on this time inteval in minutes
        base: str = "",
    ) -> None:
        self.cfg = cfg
        self.radar = radar
        self.model = model
        self.start_time = start_time
        self.beam = beam
        self.del_time = del_time
        self.base = base
        self.source = source
        self.target = target
        self.folder, self.dop_folder = self._get_folder_()
        return

    def _get_folder_(self):
        folder = (
            utils.get_folder(
                self.radar.rad, self.beam, self.start_time, self.model, self.base
            )
            if self.radar
            else utils.get_hamsci_folder(
                self.source["call_sign"],
                self.start_time,
                self.model,
                self.base,
                call_sign=self.target["call_sign"],
            )
        )
        dop_folder = os.path.join(folder, "Doppler")
        os.system(f"rm -rf {dop_folder}")
        os.makedirs(dop_folder, exist_ok=True)
        return folder, dop_folder

    def _fetch_bearing_rays_(self, date: dt.datetime):
        # Load bearing
        bearing = loadmat(os.path.join(self.folder, "bearing.mat"))
        grange, height = (
            np.array(bearing["dist"]).ravel(),
            np.array(bearing["ht"]).ravel(),
        )
        # Load rays and density
        rfile, efile = (
            os.path.join(self.folder, f"{date.strftime('%H%M')}_rt.mat"),
            os.path.join(self.folder, f"{date.strftime('%H.%M')}.mat"),
        )
        bearing["elvs"] = np.linspace(
            float(bearing["elev_s"]),
            float(bearing["elev_e"]),
            int((bearing["elev_e"] - bearing["elev_s"]) / bearing["elev_i"]) + 1,
        )
        rays = Rays2D(date, bearing["elvs"], self.folder, rfile)
        density = np.ma.masked_invalid(loadmat(efile)["ne"] * 1e6)
        fun = RectBivariateSpline(grange, height, np.log10(density).T)
        return SimpleNamespace(
            **dict(
                rays=rays,
                density=fun,
                time=date,
                bearing=bearing,
            )
        )

    def _compute_doppler_from_prev_time_(
        self, now: dt.datetime, prev: dt.datetime, to_meters: float = 1.0
    ) -> None:
        delt = self.del_time * self.cfg.rise_time_sec
        self.now, self.prev = (now, prev)
        along = (
            f"{self.radar.rad}/{self.beam}"
            if self.radar
            else f"{self.source['call_sign']}/{self.target['call_sign']}"
        )
        logger.info(
            f"Run Doppler compute for {along}, on {self.now}, Model:{self.model}"
        )
        # Fetch bearing and simulation data files
        event, baseline = (
            self._fetch_bearing_rays_(self.now),
            self._fetch_bearing_rays_(self.prev),
        )
        elvs, frequency = (event.bearing["elvs"], event.bearing["freq"])

        doppler = dict(time=now.strftime("%Y-%m-%dT%H:%M"), rays=[])
        base_ne_fn, event_ne_fn = (event.density, baseline.density)
        event_ray, base_ray = event.rays, baseline.rays

        dop_file = os.path.join(self.dop_folder, f"{now.strftime('%H%M')}.mat")
        if not os.path.exists(dop_file):
            for _, elv in enumerate(elvs):  # Loop for each rays
                event_ray_path = event_ray.simulation[elv][
                    "path_data"
                ]  # Extract ray associated to elv
                base_ray_path = base_ray.simulation[elv][
                    "path_data"
                ]  # Extract baseline ray associated to elv

                event_ray_label, base_ray_label = (
                    event_ray.simulation[elv]["ray_data"]["ray_label"],
                    base_ray.simulation[elv]["ray_data"]["ray_label"],
                )
                # Compute change in electron density along the modified ray
                # Need to rethnik about how rays experiance change in Doppler
                # now we implemented as if modified rays observed a difference from
                # baseline. We can also apply difference between old and modified
                # rays, that will need to intepolate height and both rays equal points.

                # Note that we do that only for rays reaching ground i.e., ray_label == 1
                if (event_ray_label == 1) and (base_ray_label == 1):
                    ray_dop = self._solve_doppler_equation_(
                        elv,
                        event_ne_fn,
                        event_ray_path,
                        base_ne_fn,
                        base_ray_path,
                        frequency,
                        to_meters,
                    )
                    # Doppler shift and velocity calculated from .phase_path
                    dp = (
                        event_ray.simulation[elv]["ray_data"]["phase_path"]
                        - base_ray.simulation[elv]["ray_data"]["phase_path"]
                    ) * 1e3  # convert to meters
                    dop_shift = (
                        (-2.0 * frequency * 1e6 / utils.pconst["c"]) * (dp / (delt))
                    ).ravel()[0] * self.cfg.doppler_multiplier
                    dop_vel = (
                        0.5 * dop_shift * utils.pconst["c"] / (frequency * 1e6)
                    ).ravel()[0]
                    setattr(ray_dop, "pharlap_doppler_shift", dop_shift)
                    setattr(ray_dop, "pharlap_doppler_vel", dop_vel)
                    doppler["rays"].append(ray_dop)
            savemat(
                dop_file,
                dict(
                    doppler=doppler,
                ),
            )
        return

    def _solve_doppler_equation_(
        self,
        elv,
        event_ne_fn,
        event_ray_path,
        base_ne_fn,
        base_ray_path,
        frequency,
        to_meters,
    ):
        delt = self.del_time * self.cfg.rise_time_sec
        # Compute change in height for event ray
        d_height = (
            np.diff(event_ray_path["height"], prepend=event_ray_path["height"][0])
            * to_meters
        )  # Convert to meters if not
        # Compute change in electron density along the modified ray
        # Need to rethnik about how rays experiance change in Doppler
        # now we implemented as if modified rays observed a difference from
        # baseline.
        d_ne = (
            (
                10
                ** event_ne_fn(
                    event_ray_path["ground_range"], event_ray_path["height"], grid=False
                )
            )
            - (
                10
                ** base_ne_fn(
                    event_ray_path["ground_range"], event_ray_path["height"], grid=False
                )
            )
        ) * self.cfg.doppler_multiplier
        # Delete all interaction below 50 km (parameterize by config)
        d_ne[event_ray_path["height"] <= 50] = 0.0
        # Compute change in Doppler freqency due to change in refractive index along the ray
        d_frq_dne = (
            (utils.pconst["kconst"] / (utils.pconst["c"] * frequency * 1e6))
            * (d_ne / delt)
            * d_height
            / np.cos(np.deg2rad(90.0 - elv))
        )
        # Compute total change in Doppler frequency for change in refractive index
        frq_dne = np.trapz(
            np.array(d_frq_dne), np.array(event_ray_path["ground_range"])
        )
        # Compute total change in Doppler frequency due to change in reflection height
        dh = (
            (event_ray_path["height"].max() - base_ray_path["height"].max())
            * 1e3
            * self.cfg.doppler_multiplier
        )
        frq_dh = (
            (-2.0 * frequency * 1e6 / utils.pconst["c"])
            * (dh / (delt))
            * np.cos(np.deg2rad(90 - elv))
        )
        vel_dne = 0.5 * frq_dne * utils.pconst["c"] / (frequency * 1e6)
        vel_dh = 0.5 * frq_dh * utils.pconst["c"] / (frequency * 1e6)
        ray_dop = dict(
            elv=elv,  # Elevation
            d_height=d_height.ravel(),  # Differential height
            event_ray_path_height=event_ray_path["height"].ravel(),  # Event ray height
            event_ray_path_ground_range=event_ray_path[
                "ground_range"
            ].ravel(),  # Event ray ground range
            d_frq_dne=d_frq_dne.ravel()[
                0
            ],  # Doppler along the ray path due to refraction
            frq_dne=frq_dne.ravel()[
                0
            ],  # Total doppler along the ray path due to refraction
            vel_dne=vel_dne.ravel()[
                0
            ],  # Total velocity along the ray path due to refraction
            frq_dh=frq_dh.ravel()[0],  # Total doppler at the peak due to reflection
            vel_dh=vel_dh.ravel()[0],  # Total velocity at the peak due to reflection
            geometric_distance=event_ray_path["geometric_distance"].tolist()[
                -1
            ],  # Slant range distance
            vel_tot=(vel_dh + vel_dne).ravel()[0],  # Total cumulative velocity
        )
        return SimpleNamespace(**ray_dop)

    def fetch_ray_dop(self, time: dt.datetime, calculate_dop_stats: bool = True):
        event = self._fetch_bearing_rays_(time)
        folder = os.path.join(self.folder, "Doppler")
        files = glob.glob(folder + f"/{time.strftime('%H%M')}.mat")
        if len(files) > 0:
            files.sort()
            doppler = utils.loadmatlabfiles(files[0])
            setattr(event, "doppler", SimpleNamespace(**doppler))
            if calculate_dop_stats:
                pass
        else:
            setattr(event, "doppler", None)
        return event


class SuperDARNDoppler(Doppler):
    def __init__(
        self,
        cfg: SimpleNamespace,
        start_time: dt.datetime,
        model: str,
        radar: Radar,
        beam: int = 0,
        del_time: int = 1,  # Normalized based on this time inteval in minutes
        base: str = "",
    ) -> None:
        super(SuperDARNDoppler, self).__init__(
            cfg,
            start_time,
            model,
            radar=radar,
            beam=beam,
            source=None,
            target=None,
            del_time=del_time,
            base=base,
        )
        return

    @staticmethod
    def fetch_by_scan_time(
        event,
        rad,
        model,
        beams,
        base,
        frange=180,
        rsep=45,
    ):
        records = []
        for beam in beams:
            folder = utils.get_folder(rad, beam, event, model, base)
            folder = os.path.join(folder, "Doppler")
            file = os.path.join(folder, event.strftime("%H%M.mat"))
            if os.path.exists(file):
                doppler = utils.loadmatlabfiles(file)
                for ray in doppler["doppler"]["rays"]:
                    try:
                        ray = utils._todict_(ray)
                        srange = ray["geometric_distance"]
                        gate = int((srange - frange) / rsep)
                        d = dict(
                            time=doppler["doppler"]["time"],
                            srange=srange,
                            bmnum=beam,
                            slist=gate,
                            vel_tot=ray["vel_tot"],
                            frq_dne=ray["frq_dne"],
                            vel_dne=ray["vel_dne"],
                            frq_dh=ray["frq_dh"],
                            vel_dh=ray["vel_dh"],
                            pharlap_doppler_vel=ray["pharlap_doppler_vel"],
                            pharlap_doppler_shift=ray["pharlap_doppler_shift"],
                            elv=ray["elv"],
                        )
                    except:
                        logger.error(f"Loading {file} error")
                        d = dict(
                            time=doppler["doppler"]["time"],
                            srange=np.nan,
                            bmnum=np.nan,
                            slist=np.nan,
                            vel_tot=np.nan,
                            frq_dne=np.nan,
                            vel_dne=np.nan,
                            frq_dh=np.nan,
                            vel_dh=np.nan,
                            pharlap_doppler_vel=np.nan,
                            pharlap_doppler_shift=np.nan,
                            elv=np.nan,
                        )
                    records.append(d)
        records = pd.DataFrame.from_records(records)
        return records

    @staticmethod
    def fetch_by_beam(
        event,
        rad,
        model,
        beam,
        base,
        frange=180,
        rsep=45,
    ):
        folder = utils.get_folder(rad, beam, event, model, base)
        folder = os.path.join(folder, "Doppler")
        files = glob.glob(folder + "/*.mat")
        records = []
        for file in files:
            logger.info(f"Load doppler file: {file}")
            doppler = utils.loadmatlabfiles(file)
            for ray in doppler["doppler"]["rays"]:
                try:
                    ray = utils._todict_(ray)
                    srange = ray["geometric_distance"]
                    gate = int((srange - frange) / rsep)
                    d = dict(
                        time=doppler["doppler"]["time"],
                        srange=srange,
                        bmnum=beam,
                        slist=gate,
                        vel_tot=ray["vel_tot"],
                        frq_dne=ray["frq_dne"],
                        vel_dne=ray["vel_dne"],
                        frq_dh=ray["frq_dh"],
                        vel_dh=ray["vel_dh"],
                        pharlap_doppler_vel=ray["pharlap_doppler_vel"],
                        pharlap_doppler_shift=ray["pharlap_doppler_shift"],
                        elv=ray["elv"],
                    )
                except:
                    logger.error(f"Loading {file} error")
                    d = dict(
                        time=doppler["doppler"]["time"],
                        srange=np.nan,
                        bmnum=np.nan,
                        slist=np.nan,
                        vel_tot=np.nan,
                        frq_dne=np.nan,
                        vel_dne=np.nan,
                        frq_dh=np.nan,
                        vel_dh=np.nan,
                        pharlap_doppler_vel=np.nan,
                        pharlap_doppler_shift=np.nan,
                        elv=np.nan,
                    )
                records.append(d)
        records = pd.DataFrame.from_records(records)
        return records


class HamSCIDoppler(Doppler):
    def __init__(
        self,
        source: dict,
        target: dict,
        cfg: SimpleNamespace,
        start_time: dt.datetime,
        model: str,
        del_time: int = 1,  # Normalized based on this time inteval in minutes
        base: str = "",
    ) -> None:
        super(HamSCIDoppler, self).__init__(
            cfg,
            start_time,
            model,
            radar=None,
            beam=None,
            source=source,
            target=target,
            del_time=del_time,
            base=base,
        )
        return

    @staticmethod
    def fetch_records(
        source,
        target,
        event,
        model,
        base,
    ):
        folder = utils.get_hamsci_folder(
            source["call_sign"], event, model, base, target["call_sign"]
        )
        folder = os.path.join(folder, "Doppler")
        files = glob.glob(folder + "/*.mat")
        files.sort()
        records = []
        for file in files:
            logger.info(f"Load doppler file: {file}")
            doppler = utils.loadmatlabfiles(file)
            for ray in doppler["doppler"]["rays"]:
                try:
                    ray = utils._todict_(ray)
                    srange = ray["geometric_distance"]
                    d = dict(
                        time=doppler["doppler"]["time"],
                        srange=srange,
                        vel_tot=ray["vel_tot"],
                        frq_dne=ray["frq_dne"],
                        vel_dne=ray["vel_dne"],
                        frq_dh=ray["frq_dh"],
                        vel_dh=ray["vel_dh"],
                        pharlap_doppler_vel=ray["pharlap_doppler_vel"],
                        pharlap_doppler_shift=ray["pharlap_doppler_shift"],
                        elv=ray["elv"],
                    )
                except:
                    logger.error(f"Loading {file} error")
                    d = dict(
                        time=doppler["doppler"]["time"],
                        srange=np.nan,
                        bmnum=np.nan,
                        slist=np.nan,
                        vel_tot=np.nan,
                        frq_dne=np.nan,
                        vel_dne=np.nan,
                        frq_dh=np.nan,
                        vel_dh=np.nan,
                        pharlap_doppler_vel=np.nan,
                        pharlap_doppler_shift=np.nan,
                        elv=np.nan,
                    )
                records.append(d)
        records = pd.DataFrame.from_records(records)
        return records
