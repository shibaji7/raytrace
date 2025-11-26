#!/usr/bin/env python3

"""run_eclipse.py: Analyze the eclipse related SCurve runs"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "chakras4@erau.edu"
__status__ = "Research"

import matplotlib.pyplot as plt

plt.rcParams["text.usetex"] = False

import argparse
import datetime as dt
import os

import numpy as np
import pandas as pd
from dateutil import parser as dparser
from loguru import logger

from rt import radar, utils
from rt.rti import RangeTimeIntervalPlot
from rt.run_sd_simulations import RadarSimulation

if __name__ == "__main__":
    # add_sys_paths()
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--beam", default=15, help="Radar beam number", type=int)
    parser.add_argument(
        "-f",
        "--cfg_file",
        default="cfg/Eclipse/rt2d_gitm_2021_eclipse_base.json",
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
