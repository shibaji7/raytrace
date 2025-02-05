#!/usr/bin/env python3

"""analysis.py: Running the main body of the simulation"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "chakras4@erau.edu"
__status__ = "Research"

import argparse
import sys

from agw_utils import add_sys_paths
from dateutil import parser as dparser
from loguru import logger

from rt import radar, utils
from rt.doppler import Doppler
from rt.rti import RangeTimeIntervalPlot
from rt.run_simulations import RadarSimulation

if __name__ == "__main__":
    add_sys_paths()
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--beam", default=11, help="Radar beam number", type=int)
    parser.add_argument(
        "-f",
        "--cfg_file",
        default="../../cfg/rt2d_gemini_May2017_tid.json",
        help="Configuration file",
        type=str,
    )
    parser.add_argument("-md", "--method", default="rt", help="Method rt/fan")
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
