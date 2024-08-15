#!/usr/bin/env python3

"""analysis_ped_mode.py: Analyze the pederson mode"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import sys
import os
from loguru import logger
from dateutil import parser as dparser
import argparse

CD_STEPS = ""
ZOOMED_IN = []

def add_sys_paths():
    """Adding /rt to sys path
    """
    global CD_STEPS
    pwd = os.getcwd()
    last_dir = pwd.split("/")[-1]
    index_of_trace_dir = pwd.split("/").index("raytrace")
    logger.info(f"Currently in '/{last_dir}', index of trace folder: {index_of_trace_dir}")
    local_libs = [
        os.path.join(
            "/".join((pwd.split("/")[:index_of_trace_dir+1])),
            "rt"
        ),
        os.path.join(
            "/".join((pwd.split("/")[:index_of_trace_dir+1])),
            "rt", "density"
        )
    ]
    logger.info(f"Loacl lib-{local_libs}")
    CD_STEPS = "".join(["../"]*int(len(pwd.split("/"))-index_of_trace_dir-1))
    sys.path.extend(local_libs)
    return

add_sys_paths()
import utils
from rays import Plots
from rt2d import RadarBeam2dTrace
from gemini import GEMINI2d


def create_regionnal_plots(cfg, beam):
    global CD_STEPS, ZOOMED_IN
    event = dparser.isoparse(cfg.event)
    logger.info(f"Create regional plot for {cfg.rad}/{beam}/{event}")
    base_output_folder = os.path.join(
        CD_STEPS,
        cfg.project_save_location, cfg.project_name, 
    )
    model = GEMINI2d(cfg, event)
    for d in model.dates:
        rto = RadarBeam2dTrace(
            d, cfg.rad, beam,
            cfg, cfg.model,
            base_output_folder,
        )
        eden = model.load_from_file(rto.edensity_file)
        rto.load_rto(eden)
        plot = Plots(event, cfg, rto, cfg.rad, beam)
        plot.lay_rays(kind=cfg.ray_trace_plot_kind, zoomed_in=ZOOMED_IN)
        plot.save("sample.png")
        plot.close()
        break
    
    return

zoomed_in = []
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--beam", default=0, help="Radar beam number", type=int)
    parser.add_argument(
        "-f",
        "--cfg_file",
        default="cfg/rt2d_gemini_May2017_tid.json",
        help="Configuration file",
        type=str,
    )
    # parser.add_argument("-md", "--method", default="rt", help="Method rt/fan")
    args = parser.parse_args()
    args.cfg_file = os.path.join(CD_STEPS, args.cfg_file)
    logger.info("\n Parameter list for simulation ")
    for k in vars(args).keys():
        print("     ", k, "->", str(vars(args)[k]))
    
    cfg = utils.read_params_2D(args.cfg_file)
    create_regionnal_plots(cfg, args.beam)