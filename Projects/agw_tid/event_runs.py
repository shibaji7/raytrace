#!/usr/bin/env python3

"""event_runs.py: Running the main body of the simulation"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "chakras4@erau.edu"
__status__ = "Research"

import argparse
import datetime as dt
import os

import numpy as np
from agw_utils import read_all_rays
from dateutil import parser as dparser
from loguru import logger

from rt import radar, utils
from rt.rti import RangeTimeIntervalPlot
from rt.run_sd_simulations import RadarSimulation

if __name__ == "__main__":
    # add_sys_paths()
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--beam", default=11, help="Radar beam number", type=int)
    parser.add_argument(
        "-f",
        "--cfg_file",
        default="cfg/May2017_TIDs/rt2d_gemini_May2017_control_cosmic_ro.json",
        help="Configuration file",
        type=str,
    )
    parser.add_argument("-md", "--method", default="rt", help="Method rt/rti/fan")
    args = parser.parse_args()
    logger.info("\n Parameter list for simulation ")
    for k in vars(args).keys():
        print("     ", k, "->", str(vars(args)[k]))

    cfg = utils.read_params_2D(args.cfg_file)
    cfg.event = dparser.isoparse(cfg.event)

    if args.method == "rt":
        from rt import radar

        rad = utils.read_params_2D(args.cfg_file).rad
        beams = radar.get_beams(rad) if args.beam == -1 else [args.beam]
        for beam in beams:
            rsim = RadarSimulation(args.cfg_file, beam=beam)
            rsim.gerenate_fov_plot()
            rsim.run_2d_simulation()
            rsim.compute_doppler()
            rsim.generate_rti()
    if args.method == "rti":
        radar = radar.Radar(
            cfg.rad,
            [cfg.event, cfg.event + dt.timedelta(minutes=cfg.time_window)],
            cfg,
        )
        print(radar.df.head(), radar.df.columns)
        rays = read_all_rays(cfg, args.beam)
        rays.lag_power = np.log10(rays.lag_power)
        print(rays.lag_power.min(), rays.lag_power.max())
        rays["bmnum"] = args.beam

        rti = RangeTimeIntervalPlot(
            4000,
            [cfg.event, cfg.event + dt.timedelta(hours=5)],
            cfg.rad,
            fig_title="",
            num_subplots=2,
            srange_type="srange",
        )
        rti.addParamPlot(
            rays,
            args.beam,
            title="GEMINI+PHaRLAP",
            xparam="date",
            zparam="lag_power",
            lay_eclipse=cfg.event_type.eclipse,
            kind="scatter",
            label="Lag Power, (dB)",
            p_max=-8,
            p_min=-12,
        )
        rti.addParamPlot(
            radar.df.copy(),
            args.beam,
            zparam="p_l",
            title="Observations",
            xlabel="Time (UT)",
            lay_eclipse=cfg.event_type.eclipse,
            kind="scatter",
            label="Lag Power, (dB)",
            p_max=20,
            p_min=3,
        )
        base = os.path.join(cfg.project_save_location, cfg.project_name)
        filepath = (
            utils.get_folder(
                cfg.rad,
                args.beam,
                cfg.event,
                cfg.model,
                base,
            )
            + "/rti_pl.png"
        )
        logger.info(f"File: {filepath}")
        rti.save(filepath)
        rti.close()
    # if args.method == "fan":
    #     import utils

    #     cfg = utils.read_params_2D(args.cfg_file)
    #     cfg.event = dparser.isoparse(cfg.event)
    #     dates = [cfg.event, cfg.event + dt.timedelta(minutes=cfg.time_window)]
    #     date = dates[0] + dt.timedelta(minutes=cfg.time_gaps)

    #     while date < dates[-1]:
    #         RadarSimulation.genererate_fan(cfg, date, args.cfg_file)
    #         date += dt.timedelta(minutes=cfg.time_gaps)
    utils.clean()
