#!/usr/bin/env python3

"""run_eclipse_rt.py: Analyze the eclipse related SCurve runs"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "chakras4@erau.edu"
__status__ = "Research"

import os
import sys
import argparse
from loguru import logger
from dateutil import parser as dparser

sys.path.extend([".", "rt/", "rt/density/"])
import utils

def run_eclipse_rt(args):
    from simulate import RadarSimulation
    rs = RadarSimulation(beam=args.beam, cfg_file=args.cfg_file)
    rs.gerenate_fov_plot()
    rs.run_2d_simulation()
    rs.compute_doppler()
    rs.generate_rti()
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--beam", default=11, help="Radar beam number", type=int)
    parser.add_argument(
        "-f",
        "--cfg_file",
        default="cfg/rt2d_iri_2024_SCurve.json",
        help="Configuration file",
        type=str,
    )
    args = parser.parse_args()
    logger.info("\n Parameter list for simulation ")
    for k in vars(args).keys():
        print("     ", k, "->", str(vars(args)[k]))

    cfg = utils.read_params_2D(args.cfg_file)
    cfg.event = dparser.isoparse(cfg.event)
    run_eclipse_rt(args)
