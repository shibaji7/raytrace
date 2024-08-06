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

import concurrent.futures
import glob
import os
from types import SimpleNamespace

import numpy as np
import pandas as pd
import utils
from loguru import logger
from rays import Rays2D
from scipy.interpolate import RectBivariateSpline
from scipy.io import loadmat, savemat


class Doppler(object):

    def __init__(
        self,
        event,
        rad,
        beams=[],
        model="gitm",
        rtime=1,
        base="figures/rt/",
        parallel=True,
    ) -> None:
        self.event = event
        self.rad = rad
        self.beams = beams if len(beams) > 0 else np.arange(24)
        self.model = model
        self.rtime = rtime
        self.base = base
        self.parallel = parallel
        if parallel:
            logger.info(f"Run all the doppler parallel calculations...")
            with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
                _ = list(executor.map(self._compute_doppler_from_prev_time, beams))
        else:
            for beam in beams:
                self._compute_doppler_from_prev_time(beam)
                break
        return

    def _get_folder_(self, beam):
        folder = utils.get_folder(
            self.rad, beam, self.event, self.model, True, self.base
        )
        folder = os.path.join("/".join(folder.split("/")[:-1]), "Doppler")
        os.system(f"rm -rf {folder}")
        os.makedirs(folder, exist_ok=True)
        return folder

    def _compute_doppler_from_prev_time(self, beam) -> None:
        folder = self._get_folder_(beam)
        logger.info(
            f"Run Doppler compute for {self.rad}/{beam}, on {self.event}, Model:{self.model}"
        )
        # Fetch bearing and simulation data files
        event = self._fetch_files_bearing_(beam, False)
        elvs, frequency = (event.bearing["elvs"], event.bearing["freq"])
        for i, t in zip(
            range(1, len(event.times)), event.times[1:]
        ):  # Loop through each time
            doppler = dict(time=t.strftime("%Y-%m-%dT%H:%M"), rays=[])
            base_ne_fn, event_ne_fn = (event.densities[i - 1], event.densities[i])
            event_ray, base_ray = event.rays[i], event.rays[i - 1]
            logger.info(f"Solving for time {t}")
            dop_file = os.path.join(folder, f"{t.strftime('%H%M')}.mat")
            if not os.path.exists(dop_file):
                for _, elv in enumerate(elvs):  # Loop for each rays
                    event_ray_path = event_ray.ray_path_data[
                        elv
                    ]  # Extract ray associated to elv
                    base_ray_path = base_ray.ray_path_data[
                        elv
                    ]  # Extract baseline ray associated to elv
                    # Compute change in electron density along the modified ray
                    # Need to rethnik about how rays experiance change in Doppler
                    # now we implemented as if modified rays observed a difference from
                    # baseline. We can also apply difference between old and modified
                    # rays, that will need to intepolate height and both rays equal points.
                    ray_dop = self._solve_doppler_equation_(
                        elv,
                        event_ne_fn,
                        event_ray_path,
                        base_ne_fn,
                        base_ray_path,
                        frequency,
                    )
                    doppler["rays"].append(ray_dop)
                savemat(
                    dop_file,
                    dict(
                        doppler=doppler,
                    ),
                )
            if not self.parallel:
                break
        return

    def _compute_doppler_from_baseline_(self, beam) -> None:
        folder = self._get_folder_(beam)
        logger.info(
            f"Run Doppler compute for {self.rad}/{beam}, on {self.event}, Model:{self.model}"
        )
        # Fetch bearing and simulation data files
        base = self._fetch_files_bearing_(beam, True)
        event = self._fetch_files_bearing_(beam, False)
        elvs, frequency = (base.bearing["elvs"], base.bearing["freq"])
        for i, t in enumerate(base.times):  # Loop through each time
            doppler = dict(time=t.strftime("%Y-%m-%dT%H:%M"), rays=[])
            base_ne_fn, event_ne_fn = (base.densities[i], event.densities[i])
            event_ray, base_ray = event.rays[i], base.rays[i]
            logger.info(f"Solving for time {t}")
            dop_file = os.path.join(folder, f"{t.strftime('%H%M')}.mat")
            if not os.path.exists(dop_file):
                for _, elv in enumerate(elvs):  # Loop for each rays
                    event_ray_path = event_ray.ray_path_data[
                        elv
                    ]  # Extract ray associated to elv
                    base_ray_path = base_ray.ray_path_data[
                        elv
                    ]  # Extract baseline ray associated to elv
                    # Compute change in electron density along the modified ray
                    # Need to rethnik about how rays experiance change in Doppler
                    # now we implemented as if modified rays observed a difference from
                    # baseline. We can also apply difference between old and modified
                    # rays, that will need to intepolate height and both rays equal points.
                    ray_dop = self._solve_doppler_equation_(
                        elv,
                        event_ne_fn,
                        event_ray_path,
                        base_ne_fn,
                        base_ray_path,
                        frequency,
                    )
                    doppler["rays"].append(ray_dop)
                savemat(
                    dop_file,
                    dict(
                        doppler=doppler,
                    ),
                )
            if not self.parallel:
                break
        return

    def _solve_doppler_equation_(
        self,
        elv,
        event_ne_fn,
        event_ray_path,
        base_ne_fn,
        base_ray_path,
        frequency,
    ):
        kconst, cconst, delt = 80.6, 3e8, self.rtime * 60
        # Compute change in height for event ray
        d_height = np.diff(event_ray_path.height, prepend=event_ray_path.height[0])
        # Compute change in electron density along the modified ray
        # Need to rethnik about how rays experiance change in Doppler
        # now we implemented as if modified rays observed a difference from
        # baseline.
        d_ne = 10 ** base_ne_fn(
            event_ray_path.ground_range, event_ray_path.height, grid=False
        ) - 10 ** event_ne_fn(
            event_ray_path.ground_range, event_ray_path.height, grid=False
        )
        # Delete all interaction below 50 km (parameterize by config)
        d_ne[event_ray_path.height <= 50] = 0.0
        # Compute change in Doppler freqency due to change in refractive index along the ray
        d_frq_dne = (
            (kconst / (cconst * frequency * 1e6))
            * (d_ne / delt)
            * d_height
            / np.cos(np.deg2rad(90.0 - elv))
        )
        # Compute total change in Doppler frequency for change in refractive index
        frq_dne = np.trapz(np.array(d_frq_dne), np.array(event_ray_path.ground_range))
        # Compute total change in Doppler frequency due to change in reflection height
        dh = (base_ray_path.height.max() - event_ray_path.height.max()) * 1e3
        frq_dh = (
            (-2.0 * frequency * 1e6 / cconst) * (dh / (delt)) * np.cos(np.deg2rad(elv))
        )
        ray_dop = dict(
            elv=elv,  # Elevation
            d_height=d_height,  # Differential height
            event_ray_path_height=event_ray_path.height,  # Event ray height
            event_ray_path_ground_range=event_ray_path.ground_range,  # Event ray ground range
            d_frq_dne=d_frq_dne,  # Doppler along the ray path
            frq_dne=frq_dne,  # Total doppler along the ray path
            vel_dne=0.5
            * frq_dne
            * cconst
            / (frequency * 1e6),  # Total velocity along the ray path
            frq_dh=frq_dh,  # Total doppler at the peak
            vel_dh=0.5
            * frq_dh
            * cconst
            / (frequency * 1e6),  # Total velocity at the peak
            geometric_distance=event_ray_path.geometric_distance.tolist()[
                -1
            ],  # Slant range distance
        )
        return ray_dop

    def _fetch_files_bearing_(self, beam, control):
        # Get folder
        folder = utils.get_folder(
            self.rad, beam, self.event, self.model, control, self.base
        )
        # Load bearing
        bearing = loadmat(os.path.join(folder, "bearing.mat"))
        grange, height = (
            np.array(bearing["dist"]).ravel(),
            np.array(bearing["ht"]).ravel(),
        )
        # Load rays and density
        rt_files, eden_files = (
            glob.glob(os.path.join(folder, "*_rt.mat")),
            glob.glob(os.path.join(folder, "*.*.mat")),
        )
        bearing["elvs"] = np.linspace(
            float(bearing["elev_s"]),
            float(bearing["elev_e"]),
            int((bearing["elev_e"] - bearing["elev_s"]) / bearing["elev_i"]) + 1,
        )
        rt_files.sort()
        eden_files.sort()
        rays_list, density_list, time_list = [], [], []
        for rfile, dfile in zip(rt_files, eden_files):
            rays = Rays2D(self.event, self.rad, beam, bearing["elvs"], folder, rfile)
            rays_list.append(rays)
            density = loadmat(dfile)["ne"] * 1e6
            fun = RectBivariateSpline(grange, height, np.log10(density).T)
            density_list.append(fun)
            # Load time
            dx = dfile.split("/")[-1].split(".")[:-1]
            time_list.append(self.event.replace(hour=int(dx[0]), minute=int(dx[1])))
        return SimpleNamespace(
            **dict(
                rays=rays_list,
                densities=density_list,
                times=time_list,
                bearing=bearing,
            )
        )

    @staticmethod
    def fetch_by_scan_time(event, rad, model, beams, frange=180, rsep=45):
        records = []
        for beam in beams:
            folder = utils.get_folder(rad, beam, event, model, True)
            folder = os.path.join("/".join(folder.split("/")[:-1]), "Doppler")
            file = os.path.join(folder, event.strftime("%H%M.mat"))
            doppler = utils.loadmatlabfiles(file)
            for ray in doppler["doppler"]["rays"]:
                ray = utils._todict_(ray)
                srange = ray["geometric_distance"]
                vel_tot = ray["vel_dne"] + ray["vel_dh"]
                gate = int((srange - frange) / rsep)
                records.append(
                    dict(
                        time=doppler["doppler"]["time"],
                        srange=srange,
                        bmnum=beam,
                        slist=gate,
                        vel_tot=vel_tot,
                        frq_dne=ray["frq_dne"],
                        v=ray["vel_dne"],
                        frq_dh=ray["frq_dh"],
                        vel_dh=ray["vel_dh"],
                    )
                )
        records = pd.DataFrame.from_records(records)
        return records

    @staticmethod
    def fetch_by_beam(event, rad, model, beam, frange=180, rsep=45):
        folder = utils.get_folder(rad, beam, event, model, True)
        folder = os.path.join("/".join(folder.split("/")[:-1]), "Doppler")
        files = glob.glob(folder + "/*.mat")
        records = []
        for file in files:
            doppler = utils.loadmatlabfiles(file)
            for ray in doppler["doppler"]["rays"]:
                ray = utils._todict_(ray)
                srange = ray["geometric_distance"]
                vel_tot = ray["vel_dne"] + ray["vel_dh"]
                gate = int((srange - frange) / rsep)
                records.append(
                    dict(
                        time=doppler["doppler"]["time"],
                        srange=srange,
                        bmnum=beam,
                        slist=gate,
                        vel_tot=vel_tot,
                        frq_dne=ray["frq_dne"],
                        vel_dne=ray["vel_dne"],
                        frq_dh=ray["frq_dh"],
                        vel_dh=ray["vel_dh"],
                    )
                )
        records = pd.DataFrame.from_records(records)
        return records
