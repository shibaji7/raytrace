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
import datetime as dt
import os
import sys

from dateutil import parser as dparser
from loguru import logger

sys.path.extend([".", "rt/", "rt/density/", "Projects/SCurve/"])
import utils
from hamsci import HamSci


class HamSCISimulation(object):

    def __init__(
        self,
        cfg_file: str,
    ) -> None:

        ## Kill all Matlab Server Hosts for this run
        os.system("killall MathWorksServiceHost.")
        os.system("rm -rf ~/matlab_crash_dump.*")

        self.cfg_file = cfg_file
        self.cfg = utils.read_params_2D(cfg_file)
        self.start_time = dparser.isoparse(self.cfg.event)
        self.end_time = dparser.isoparse(self.cfg.event) + dt.timedelta(
            minutes=self.cfg.time_window
        )
        self.hamsci = HamSci(
            self.cfg, [self.start_time, self.end_time], self.cfg.frequency
        )
        self.hamsci.setup_pandas_dataset()
        return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--cfg_file",
        default="cfg/rt2d_iri_2024_iri_hamsci_SCurve.json",
        help="Configuration file",
        type=str,
    )
    args = parser.parse_args()
    logger.info("\n Parameter list for simulation ")
    for k in vars(args).keys():
        print("     ", k, "->", str(vars(args)[k]))

    sim = HamSCISimulation(args.cfg_file)
