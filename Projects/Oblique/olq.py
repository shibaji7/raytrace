#!/usr/bin/env python3

"""run_w2naf_trace.py: Analyze the eclipse related SCurve runs for w2naf"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "chakras4@erau.edu"
__status__ = "Research"

import argparse
import concurrent.futures
import datetime as dt
import os
import sys

import pandas as pd
from dateutil import parser as dparser
from loguru import logger

sys.path.extend([".", "rt/", "rt/density/", "Projects/SCurve/"])
from hamsci import HamSci

from rt.run_sd_simulations import RadarSimulation

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--cfg_file",
        default="cfg/rt2d_iri_2015_flare.json",
        help="Configuration file",
        type=str,
    )
    parser.add_argument(
        "-b",
        "--beam",
        default=11,
        help="Radar Beam",
        type=int,
    )
    args = parser.parse_args()
    logger.info("\n Parameter list for simulation ")
    for k in vars(args).keys():
        print("     ", k, "->", str(vars(args)[k]))

    sim = RadarSimulation(args.cfg_file, args.beam)
    sim.run_2d_simulation()
    # sim.compute_doppler()
    # sim.generate_ls()

    # TODO
    # 0. Send the movie/.pngs to groups [D]
    # 1. Check & plot where W2NAF is with respect to WWV [D]
    # 2. May need 2-hop GS [D]
    # 3. Plot / work on Simulated Data TS plots,
    #   a. Angle, Slant range relatated stats
    #   b. Run for 5 and 15 MHz
    # 4. Check a how dh vs dn works out. Need E field values.
