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
import datetime as dt
import os
import sys

from dateutil import parser as dparser
from loguru import logger

sys.path.extend(["../../rt/", "../../rt/density/"])
import radar
import utils
from doppler import Doppler
from rti import RangeTimeIntervalPlot

CD_STEPS = ""
_DIR_ = "figures/zoomed/"
DATES = [
    dt.datetime(2017, 5, 27, 18, 0) + dt.timedelta(minutes=i * 10) for i in range(15)
]
ZOOMED_IN = [[500, 1600], [150, 250]]
ELV_RANGE = [20, 30]


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
    global CD_STEPS, _DIR_, DATES, ZOOMED_IN
    from gemini import GEMINI2d
    from rays import PlotRays
    from rt2d import RadarBeam2dTrace

    base_output_folder = os.path.join(
        CD_STEPS,
        cfg.project_save_location,
        cfg.project_name,
    )
    logger.info(f"Base folder: {base_output_folder}")
    DATES.append(cfg.event)
    for d in DATES:
        rto = RadarBeam2dTrace(
            d,
            cfg.rad,
            beam,
            cfg,
            cfg.model,
            base_output_folder,
        )
        model = GEMINI2d(cfg, d)
        eden = model.load_from_file(rto.edensity_file)
        rto.load_rto(eden)
        plot = PlotRays(d, cfg, rto, cfg.rad, beam, xlim=[0, 1600])
        plot.lay_rays(kind=cfg.ray_trace_plot_kind, zoomed_in=ZOOMED_IN)
        file = os.path.join(CD_STEPS, _DIR_, f"{d.strftime('%Y%m%d.%H%M')}.png")
        logger.info(f"Saved to file: {file}")
        plot.save(file)
        plot.close()
    return


def create_rtis_by_radars(cfg, beam, rads=["fhe", "fhw", "bks"]):
    global CD_STEPS, ELV_RANGE, _DIR_
    start = dt.datetime(2017, 5, 27, 14)
    end = dt.datetime(2017, 5, 28)
    rtint = RangeTimeIntervalPlot(
        60,
        [start, end],
        cfg.rad,
        fig_title="",
        num_subplots=len(rads),
        srange_type="srange",
    )
    for rad in rads:
        radr = radar.Radar(rad, [start, end], cfg)
        fig_title = f"{rad.upper()}-{'%02d'%beam}, {radr.df.tfreq.median()} MHz /\t {start.strftime('%d %b, %Y')}"
        rtint.addParamPlot(
            radr.df.copy(),
            beam,
            title=f"Observation / {fig_title}",
            xlabel="",
            lay_eclipse=None,
        )
    file = os.path.join(CD_STEPS, _DIR_, f"{cfg.rad}.{start.strftime('%Y%m%d')}.png")
    logger.info(f"Save to file: {file}")
    rtint.save(file)
    rtint.close()
    return


def create_zoomed_rti_rays(cfg, beam):
    global CD_STEPS, ELV_RANGE, _DIR_
    start = cfg.event
    end = start + dt.timedelta(minutes=int(cfg.time_window))
    radr = radar.Radar(cfg.rad, [start, end], cfg)
    logger.info(f"Create regional plot for {cfg.rad}/{beam}/{start}")
    base_output_folder = os.path.join(
        CD_STEPS,
        cfg.project_save_location,
        cfg.project_name,
    )

    fig_title = f"Model: {cfg.model.upper()} / {cfg.rad.upper()}-{'%02d'%cfg.beam}, {cfg.frequency} MHz \t {start.strftime('%d %b, %Y')}"
    rtint = RangeTimeIntervalPlot(
        60, [start, end], cfg.rad, fig_title=fig_title, num_subplots=2
    )
    rtint.addParamPlot(
        radr.df.copy(),
        beam,
        title="Observations",
        xlabel="",
        lay_eclipse=None,
    )
    records = Doppler.fetch_by_beam(
        start,
        cfg.rad,
        cfg.model,
        cfg.beam,
        base_output_folder,
        frange=cfg.frange,
        rsep=cfg.rsep,
    )
    if len(ELV_RANGE) == 2:
        records = records[(records.elv >= ELV_RANGE[0]) & (records.elv <= ELV_RANGE[1])]
    if len(records) > 0:
        records.time = records.time.apply(lambda x: dparser.isoparse(x))
        rtint.addParamPlot(
            records,
            beam,
            title="",
            zparam="vel_tot",
            lay_eclipse=cfg.event_type.eclipse,
        )
    file = os.path.join(CD_STEPS, _DIR_, f"{cfg.rad}.{start.strftime('%Y%m%d')}.png")
    logger.info(f"Save to file: {file}")
    rtint.save(file)
    rtint.close()
    return


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
    args = parser.parse_args()
    logger.info("\n Parameter list for simulation ")
    for k in vars(args).keys():
        print("     ", k, "->", str(vars(args)[k]))

    cfg = utils.read_params_2D(args.cfg_file)
    cfg.event = dparser.isoparse(cfg.event)
    # plot_ls(args.beam, cfg)
    # create_zoomed_rti_rays(cfg, args.beam)
    create_rtis_by_radars(cfg, args.beam)
