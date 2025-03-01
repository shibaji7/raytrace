#!/usr/bin/env python3

"""analysis.py: Analyze the eclipse related SCurve runs"""

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
    records = []
    dop_files = glob.glob(base + "Doppler/*.mat")
    dop_files.sort()
    logger.info(f"Load files: {len(dop_files)}")
    for file in dop_files:
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


def load_rt_files(folder):
    cfg = utils.read_params_2D("cfg/rt2d_2024_eclipse_sami3.json")
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
    dates = pd.to_datetime(dop_records.time.unique())
    shadow = get_eclipse_along_path(dates)
    return dict(rays=rays, dop_records=dop_records, shadow=shadow)


def compute_statistics(folders):
    from analysis_plots import plot_ts

    for f in folders:
        records = load_rt_files(f)
        dop = records["dop_records"]
        dop = dop[dop.elv < 40]
        plot_ts(
            dop,
            [dt.datetime(2024, 4, 8, 17), dt.datetime(2024, 4, 8, 22)],
            shadow=records["shadow"],
            filepath="TS.png",
        )
        break
    return


def overlay_data_models():
    from grapeDRF import GrapeDRF

    station = "w2naf"
    sDate = dt.datetime(2024, 4, 8)
    eDate = dt.datetime(2024, 4, 9)
    figd = {}
    figd["cfreqs"] = [20, 15, 10, 5]
    figd["cfreqs"] = [10]
    figd["solar_lat"] = 41.335116
    figd["solar_lon"] = -75.600692
    figd["overlaySolarElevation"] = False
    figd["overlayEclipse"] = False
    gDRF = GrapeDRF(sDate, eDate, station)
    fig, axes, png_fpath = gDRF.plot_figure(**figd)
    folders = [
        # "/home/shibaji/OneDrive/trace/outputs/April2024_SAMI3_eclipse_hamsci_05MHz_SCurve/",
        "/home/shibaji/OneDrive/trace/outputs/April2024_SAMI3_eclipse_hamsci_10MHz_SCurve/",
        # "/home/shibaji/OneDrive/trace/outputs/April2024_SAMI3_eclipse_hamsci_15MHz_SCurve/",
        # "/home/shibaji/OneDrive/trace/outputs/April2024_SAMI3_eclipse_hamsci_20MHz_SCurve/"
    ]

    for i, f in enumerate(folders):
        records = load_rt_files(f)
        dop, shadow = records["dop_records"], records["shadow"]
        dop = dop[dop.elv < 45]
        ax = axes[i]
        ax.set_xlim([sDate.replace(hour=16), eDate])
        # Overlay estimated Doppler
        dop.time = pd.to_datetime(dop.time)
        ax.plot(dop.time, dop.frq_dne + dop.frq_dh, "b.", ms=1.5)
        ax.set_ylim(-3, 3)
        # Plot shadow
        axt = ax.twinx()
        # utime = [dt.datetime.strptime(t, "%Y-%m-%dT%H:%M") for t in dop.time.unique()]
        axt.plot(dop.time.unique(), 1 - np.nanmax(shadow, axis=1), color="k", ls="--")
        axt.plot(dop.time.unique(), 1 - np.nanmin(shadow, axis=1), color="k", ls="--")
        axt.set_ylim(0, 1)
        axt.set_yticks([])
    fig.savefig(png_fpath, bbox_inches="tight")
    return


if __name__ == "__main__":
    overlay_data_models()
    folders = [
        # "/home/shibaji/OneDrive/trace/outputs/April2024_SAMI3_eclipse_hamsci_05MHz_SCurve/",
        "/home/shibaji/OneDrive/trace/outputs/April2024_SAMI3_eclipse_hamsci_10MHz_SCurve/",
        # "/home/shibaji/OneDrive/trace/outputs/April2024_SAMI3_eclipse_hamsci_15MHz_SCurve/",
        # "/home/shibaji/OneDrive/trace/outputs/April2024_SAMI3_eclipse_hamsci_20MHz_SCurve/"
    ]
    # compute_statistics(folders)
    # 3. Plot / work on Simulated Data TS plots,
    #   a. Angle, Slant range relatated stats
    #   b. Run for 5, 10, 15 and 20 MHz
    # 4. Check a how dh vs dn works out. Need E field values.
    pass
