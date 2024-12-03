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

import argparse
import os
import sys

from dateutil import parser as dparser
from loguru import logger

sys.path.extend([".", "rt/", "rt/density/"])
import utils

CD_STEPS = ""
_DIR_ = "figures/zoomed/"


def add_sys_paths():
    """Adding /rt to sys path"""
    global CD_STEPS, _DIR_
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


def load_files():
    return


def plot_ls(beam, cfg):
    from iri import IRI2d
    from rays import PlotRays
    from rt2d import RadarBeam2dTrace

    base_output_folder = os.path.join(
        CD_STEPS,
        cfg.project_save_location,
        cfg.project_name,
    )
    logger.info(f"Base folder: {base_output_folder}")
    rto = RadarBeam2dTrace(
        cfg.event,
        cfg.rad,
        beam,
        cfg,
        cfg.model,
        base_output_folder,
    )
    model = IRI2d(cfg, cfg.event)
    eden = model.load_from_file(rto.edensity_file)
    rto.load_rto(eden)
    plot = PlotRays(cfg.event, cfg, rto, cfg.rad, beam)
    plot.lay_rays(
        kind=cfg.ray_trace_plot_kind,
    )
    file = os.path.join(CD_STEPS, _DIR_, f"{cfg.event.strftime('%Y%m%d.%H%M')}.png")
    logger.info(f"Saved to file: {file}")
    plot.save(file)
    plot.close()
    return


if __name__ == "__main__":
    add_sys_paths()
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
    plot_ls(args.beam, cfg)
