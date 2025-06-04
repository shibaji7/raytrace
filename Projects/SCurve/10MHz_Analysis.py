#### Modify Doppler based on reach ground distance / height and along the path.


__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "chakras4@erau.edu"
__status__ = "Research"

import argparse
import os
import sys

sys.path.append("Projects/SCurve/")
import datetime as dt
import glob

import numpy as np
import pandas as pd
from dateutil import parser as dparser
from geopy.distance import great_circle as GC
from loguru import logger

from rt import eclipse, utils
from rt.rt2d import Rays2D


def get_w2naf_dist():
    source = (40.0150, -105.2705)
    target = (41.335116, -75.600692)
    return GC(source, target).km


def get_eclipse_along_path(dates):
    source = (40.0150, -105.2705)
    target = (41.335116, -75.600692)
    lats, lons = [], []
    gc = GC(source, source)
    dist = np.linspace(0, 4000, 101)
    _, bearing = utils.calculate_bearing(source[0], source[1], target[0], target[1])
    for d in dist:
        x = gc.destination(source, bearing, distance=d)
        lats.append(x[0])
        lons.append(x[1])
    shadow = eclipse.get_rti_eclipse(dates, lats, lons)
    return shadow


def load_doppler(base):
    dist = get_w2naf_dist()
    print(dist)
    records = []
    dop_files = glob.glob(base + "Doppler/*.mat")
    dop_files.sort()
    logger.info(f"Load files: {len(dop_files)}")
    for file in dop_files:
        doppler = utils.loadmatlabfiles(file)
        for ray in doppler["doppler"]["rays"]:
            try:
                ray = utils._todict_(ray)
                d_frq_dne = ray["d_frq_dne"]
                # print(ray.keys(), ray["d_frq_dne"])
                srange, grange = (
                    ray["geometric_distance"],
                    ray["event_ray_path_ground_range"][-1],
                )
                # Select all rays that ended within 100 km to the w2naf
                if np.abs(grange - dist) < 100000:
                    print(grange, dist, ray["elv"])
                    # Select rays based on ground range
                    d = dict(
                        time=doppler["doppler"]["time"],
                        srange=srange,
                        grange=grange,
                        vel_tot=ray["vel_tot"],
                        frq_dne=ray["frq_dne"],
                        vel_dne=ray["vel_dne"],
                        frq_dh=ray["frq_dh"],
                        vel_dh=ray["vel_dh"],
                        pharlap_doppler_vel=ray["pharlap_doppler_vel"],
                        pharlap_doppler_shift=ray["pharlap_doppler_shift"],
                        elv=ray["elv"],
                    )
                    records.append(d)
            except:
                logger.error(f"Loading {file} error")

    records = pd.DataFrame.from_records(records)
    records.dropna(inplace=True)
    return records


def load_rt_files(folder):
    cfg = utils.read_params_2D("cfg/2024GAE/rt2d_2024_eclipse_sami3.json")
    event = dparser.isoparse(cfg.event)
    base = folder + "2024-04-08/wwv/sami3/w2naf/"
    logger.info(f"Load files from : {base}")
    total_files = glob.glob(base + "*.mat")
    logger.info(f"Total files under here: {len(total_files)}")

    rt_files = glob.glob(base + "*_rt.mat")
    rt_files.sort()
    logger.info(f"Load files: {len(rt_files)}")
    rays = []
    for rt in rt_files:
        logger.info(f"File name: {rt}")
        ray = Rays2D.read_rays(event, cfg, base, rt)
        rays.append(ray)
    dop_records = load_doppler(base=base)
    print(dop_records.head())
    dates = pd.to_datetime(dop_records.time.unique())
    shadow = get_eclipse_along_path(dates)
    return dict(rays=rays, dop_records=dop_records, shadow=shadow)


def compute_statistics(folder):
    from analysis_plots import plot_ts

    records = load_rt_files(folder)
    dop = records["dop_records"]
    plot_ts(
        dop,
        [dt.datetime(2024, 4, 8, 17), dt.datetime(2024, 4, 8, 22)],
        shadow=records["shadow"],
        filepath="figures/TS.png",
    )

    # from grapeDRF import GrapeDRF

    # station = "w2naf"
    # sDate = dt.datetime(2024, 4, 8)
    # eDate = dt.datetime(2024, 4, 9)
    # figd = {}
    # # figd["cfreqs"] = [5, 10, 15, 20]
    # figd["cfreqs"] = [5]
    # figd["solar_lat"] = 41.335116
    # figd["solar_lon"] = -75.600692
    # figd["overlaySolarElevation"] = False
    # figd["overlayEclipse"] = False
    # gDRF = GrapeDRF(sDate, eDate, station)
    # fig, axes, png_fpath = gDRF.plot_figure(**figd)

    # fig.savefig(png_fpath, bbox_inches="tight")
    return


if __name__ == "__main__":
    folder = "/home/chakras4/OneDrive/trace/outputs/GAE2024_SAMI3_w2naf_02Hop_10MHz/"
    compute_statistics(folder)
