#!/usr/bin/env python3

"""analysis_eclipse.py: Analyze the eclipse related items"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import argparse
import os
import sys

from dateutil import parser as dparser
from loguru import logger

CD_STEPS = ""
_DIR_ = "figures/movies/"


def add_sys_paths():
    """Adding /rt to sys path"""
    global CD_STEPS
    pwd = os.getcwd()
    last_dir = pwd.split("/")[-1]
    index_of_trace_dir = pwd.split("/").index("raytrace")
    logger.info(
        f"Currently in '/{last_dir}', index of trace folder: {index_of_trace_dir}"
    )
    local_libs = [
        os.path.join("/".join((pwd.split("/")[: index_of_trace_dir + 1])), "rt"),
        os.path.join(
            "/".join((pwd.split("/")[: index_of_trace_dir + 1])), "rt", "density"
        ),
    ]
    logger.info(f"Loacl lib-{local_libs}")
    CD_STEPS = "".join(["../"] * int(len(pwd.split("/")) - index_of_trace_dir - 1))
    sys.path.extend(local_libs)
    os.makedirs(os.path.join(CD_STEPS, _DIR_), exist_ok=True)
    return


add_sys_paths()
import utils


def create_movies(cfg, args):
    base_output_folder = os.path.join(
        cfg.project_save_location,
        cfg.project_name,
        cfg.event.strftime("%Y-%m-%d"),
        cfg.rad,
        "%02d" % args.beam,
        cfg.model,
    )
    to_local = os.path.join(CD_STEPS, _DIR_)
    utils.create_movie(
        base_output_folder, f"movie_{'%02d'%args.beam}.avi", "*.png", out_fold=to_local
    )
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--beam", default=11, help="Radar beam number", type=int)
    parser.add_argument(
        "-f",
        "--cfg_file",
        default="cfg/rt2d_waccm_2017_eclipse.json",
        help="Configuration file",
        type=str,
    )
    args = parser.parse_args()
    args.cfg_file = os.path.join(CD_STEPS, args.cfg_file)
    logger.info("\n Parameter list for simulation ")
    for k in vars(args).keys():
        print("     ", k, "->", str(vars(args)[k]))

    cfg = utils.read_params_2D(args.cfg_file)
    cfg.event = dparser.isoparse(cfg.event)
    create_movies(cfg, args)
