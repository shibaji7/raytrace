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

import argparse
import datetime as dt
import os
import sys

from dateutil import parser as dparser
from loguru import logger

CD_STEPS = ""
ZOOMED_IN = [[500, 1800], [150, 300]]
ELV_RANGE = []
_DIR_ = "figures/zoomed/"
DATES = [
    dt.datetime(2017, 5, 27, 19, 52),
    dt.datetime(2017, 5, 27, 19, 53),
    dt.datetime(2017, 5, 27, 19, 54),
]


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


add_sys_paths()
import utils
from gemini import GEMINI2d
from rays import Plots
from rt2d import RadarBeam2dTrace

import radar
from doppler import Doppler
from rti import RangeTimeIntervalPlot


def create_regionnal_plots(cfg, beam):
    global CD_STEPS, ZOOMED_IN, ELV_RANGE, _DIR_
    event = dparser.isoparse(cfg.event)
    logger.info(f"Create regional plot for {cfg.rad}/{beam}/{event}")
    base_output_folder = os.path.join(
        CD_STEPS,
        cfg.project_save_location,
        cfg.project_name,
    )
    model = GEMINI2d(cfg, event)
    for d in model.dates:
        if d in DATES:
            rto = RadarBeam2dTrace(
                d,
                cfg.rad,
                beam,
                cfg,
                cfg.model,
                base_output_folder,
            )
            eden = model.load_from_file(rto.edensity_file)
            rto.load_rto(eden)
            plot = Plots(event, cfg, rto, cfg.rad, beam)
            plot.lay_rays(
                kind=cfg.ray_trace_plot_kind, 
                zoomed_in=ZOOMED_IN, 
                elv_range=ELV_RANGE
            )
            plot.save(os.path.join(CD_STEPS, _DIR_, f"{d.strftime('%Y%m%d.%H%M')}.png"))
            plot.close()
    return

def create_zoomed_rti_rays(cfg, beam):
    global CD_STEPS, ELV_RANGE, _DIR_
    start = dparser.isoparse(cfg.event)
    end = start + dt.timedelta(minutes=int(cfg.time_window))
    radr = radar.Radar(
        cfg.rad, [start, end], cfg
    )
    logger.info(f"Create regional plot for {cfg.rad}/{beam}/{start}")
    base_output_folder = os.path.join(
        CD_STEPS,
        cfg.project_save_location,
        cfg.project_name,
    )

    fig_title = f"Model: {cfg.model.upper()} / {cfg.rad.upper()}-{'%02d'%cfg.beam}, {cfg.frequency} MHz \t {start.strftime('%d %b, %Y')}"
    rtint = RangeTimeIntervalPlot(
        100, [start, end], cfg.rad, fig_title=fig_title, num_subplots=2
    )
    rtint.addParamPlot(
        radr.df.copy(),
        beam,
        title="Observations",
        xlabel="",
        lay_eclipse=None,
    )
    records = Doppler.fetch_by_beam(
        start, cfg.rad, cfg.model,
        cfg.beam, base_output_folder,
        frange=cfg.frange, rsep=cfg.rsep,
    )
    if len(ELV_RANGE) == 2:
        records = records[
            (records.elv >= ELV_RANGE[0])
            & (records.elv <= ELV_RANGE[1])
        ]
    if len(records) > 0:
        records.time = records.time.apply(lambda x: dparser.isoparse(x))
        rtint.addParamPlot(
            records,
            beam,
            title="",
            zparam="vel_tot",
            lay_eclipse=cfg.event_type.eclipse,
        )
    rtint.save(os.path.join(CD_STEPS, _DIR_, f"{cfg.rad}.{start.strftime('%Y%m%d')}.png"))
    rtint.close()
    return


zoomed_in = []
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--beam", default=11, help="Radar beam number", type=int)
    parser.add_argument(
        "-f",
        "--cfg_file",
        default="cfg/rt2d_gemini_May2017_tid_high_res.json",
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
    create_zoomed_rti_rays(cfg, args.beam)
