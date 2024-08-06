#!/usr/bin/env python3

"""simulate.py: simulate python program for RT"""

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
import json
import os
import sys
from types import SimpleNamespace

import numpy as np
from dateutil import parser as dparser
from loguru import logger

sys.path.extend(["rt/", "rt/density/"])


class RadarSimulation(object):

    def __init__(
        self,
        model: str = None,
        start_time: dt.datetime = None,
        rad: str = None,
        beam: int = 0,
        method: str = None,
        time_window: int = 240,
        time_gaps: int = 5,
        eclipse_start_time: dt.datetime = None,
        cfg_file: str = "cfg/rt2D.json",
    ) -> None:
        self.model = model
        self.method = method
        self.rad = rad
        self.beam = beam
        self.time_gaps = time_gaps
        self.time_window = time_window
        self.start_time = start_time
        self.cfg_file = cfg_file
        self.eclipse_start_time = eclipse_start_time
        self.read_params_2D()
        self.base_output_folder = os.path.join(
            self.cfg.project_save_location, self.cfg.project_name
        )
        self.load_radar()
        return

    def load_radar(self) -> None:
        from radar import Radar

        self.radar = Radar(
            self.rad,
            [self.start_time, self.start_time + dt.timedelta(minutes=self.time_window)],
            self.cfg,
        )
        return

    def read_params_2D(self, fname: str = None) -> None:
        fname = fname if fname else self.cfg_file
        logger.info(f"Load config files: {fname}")
        with open(fname, "r") as f:
            self.cfg = json.load(f, object_hook=lambda x: SimpleNamespace(**x))
        return

    def gerenate_fov_plot(
        self,
        lons=None,
        lats=None,
        extent=None,
        proj=None,
        max_gate=75,
    ):
        import cartopy
        import utils
        from fan import Fan

        lons = np.arange(-180, 180, 30) if lons is None else lons
        lats = (
            (
                np.arange(30, 70, 15)
                if self.radar.hdw.geographic.lat > 0
                else np.arange(-70, -30, 15)
            )
            if lats is None
            else lats
        )
        proj = (
            cartopy.crs.Orthographic(
                10 * int(self.radar.hdw.geographic.lon / 10),
                10 * int(self.radar.hdw.geographic.lat / 10 - 1),
            )
            if proj is None
            else proj
        )
        extent = [
            2 * int(self.radar.fov[1][:max_gate, :].min() / 2),
            2 * int(self.radar.fov[1][:max_gate, :].max() / 2),
            3 * int(self.radar.fov[0][:max_gate, :].min() / 3 - 3),
            3 * int(self.radar.fov[0][:max_gate, :].max() / 3 + 3),
        ]
        fan = Fan(
            self.rad,
            self.start_time,
            fig_title=f"{'%02d'%self.beam}",
        )
        fan.setup(
            lons,
            lats,
            extent=extent,
            proj=proj,
            lay_eclipse=None,
        )
        ax = fan.add_axes()
        ax.overaly_eclipse_path(self.cfg_file, year=self.start_time.year)
        fan.generate_fov(self.rad, [], ax=ax, beams=[self.beam])
        fan.save(
            filepath=utils.get_folder(
                self.rad,
                self.beam,
                self.start_time,
                self.model,
                self.base_output_folder,
            )
            + "/fov.png"
        )
        return

    def run_2d_simulation(self):
        logger.info(f"Inside {self.model.upper()} Simulation...")
        import plots
        from iri import IRI2d
        from gitm import GITM2d
        from rt2d import RadarBeam2dTrace

        if self.model == "iri":
            model = IRI2d(self.cfg, self.start_time)
        elif self.model == "gitm":
            model = GITM2d(self.cfg, self.start_time)
        else:
            raise ValueError(
                "Currently supporting following methods: iri, gitm, waccm-x, and wamipe"
            )
        for d in range(int(args.time_window / args.time_gaps)):
            event = self.start_time + dt.timedelta(minutes=d * args.time_gaps)
            logger.info(f"Load e-Density for, {event}")
            rto = RadarBeam2dTrace(
                event,
                self.rad,
                self.beam,
                self.cfg,
                self.model,
                self.base_output_folder,
            )

            # Load all the electron density
            if not os.path.exists(rto.edensity_file):
                eden, _ = model.fetch_dataset(
                    event,
                    rto.bearing_object["lat"],
                    rto.bearing_object["lon"],
                    rto.bearing_object["ht"],
                    to_file=rto.edensity_file,
                )
            else:
                eden = model.load_from_file(rto.edensity_file)

            # Load and compile MATLAB - PHaRLAP
            if not os.path.exists(rto.sim_fname):
                rto.compile(eden)
            else:
                rto.load_rto(eden)

            # Create RT figures
            if not os.path.exists(rto.fig_name):
                plots.plot_rays(
                    rto.folder,
                    rto.fig_name,
                    rto,
                    rf"{self.model.upper()} + {args.rad.upper()}/{str(args.beam)}, $f_0$={str(self.cfg.frequency)} MHz",
                    maxground=self.cfg.max_ground_range_km,
                    eclipse_time=args.eclipse_time,
                )
        return


def genererate_Fan(event, rad, tfreq, frame, model, param="v"):
    sys.path.append("py/")
    import cartopy
    from fan import Fan

    folder = os.path.join("figures", event.strftime("%b%Y"))
    os.makedirs(folder, exist_ok=True)
    p_max, p_min = (10, -10) if param == "v" else (33, 0)
    cmap = "Spectral" if param == "v" else "plasma"
    label = "Velocity [m/s]" if param == "v" else "Power [dB]"
    fig_title = rf"$f_0=${tfreq} MHz / {model.upper()}"
    fan = Fan(
        [rad],
        event,
        fig_title=fig_title,
    )
    fan.setup(
        np.arange(-180, 180, 30),
        np.arange(30, 70, 20),
        extent=[-150, -80, 40, 90],
        proj=cartopy.crs.Orthographic(-120, 45),
    )
    fan.generate_fov(
        rad,
        frame,
        p_name=param,
        p_max=p_max,
        p_min=p_min,
        cmap=cmap,
        label=label,
        lats=np.linspace(0, 90, num=90 * 2),
    )
    fan.save(f"{folder}/fan.{rad}-{date.strftime('%H%M')}.png")
    fan.close()
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m", "--model", default="iri", help="Model name [wam/gitm/waccm/iri]"
    )
    parser.add_argument(
        "-md", "--method", default="rt", help="Method rt/dop/fan/rti/movie"
    )
    parser.add_argument("-r", "--rad", default="cvw", help="Radar code (cvw)")
    parser.add_argument(
        "-bm", "--beam", default=7, type=int, help="Radar beam (default 7)"
    )
    parser.add_argument(
        "-ev",
        "--event",
        default=dt.datetime(2017, 8, 21, 16),
        help="Event date for simulation [YYYY-mm-ddTHH:MM]",
        type=dparser.isoparse,
    )
    parser.add_argument(
        "-tw",
        "--time_window",
        default=5,
        type=int,
        help="Time window to run the models (minutes)",
    )
    parser.add_argument(
        "-tg",
        "--time_gaps",
        default=5,
        type=int,
        help="Time gaps to run the models (1-minutes)",
    )
    parser.add_argument(
        "-et",
        "--eclipse_time",
        default=dt.datetime(2017, 8, 21, 16),
        help="Event date for simulation [YYYY-mm-ddTHH:MM]",
        type=dparser.isoparse,
    )
    parser.add_argument("-cd", "--clear_dop", default=0, type=int, help="clear doppler")
    parser.add_argument(
        "-f", "--cfg_file", default="cfg/rt2d.json", help="Configuration file"
    )
    args = parser.parse_args()
    logger.info("\n Parameter list for simulation ")
    for k in vars(args).keys():
        print("     ", k, "->", str(vars(args)[k]))
    rsim = RadarSimulation(
        model=args.model,
        start_time=args.event,
        rad=args.rad,
        beam=args.beam,
        method=args.method,
        time_window=args.time_window,
        time_gaps=args.time_gaps,
        eclipse_start_time=args.eclipse_time,
        cfg_file=args.cfg_file,
    )
    rsim.gerenate_fov_plot()
    rsim.run_2d_simulation()
    # if args.method == "rt":
    #     sys.path.append("rt/")
    #     from rt2D import (
    #         execute_2DRT_GITM_simulations,
    #         execute_2DRT_IRI_simulations,
    #         execute_2DRT_WACCMX_simulations,
    #         execute_2DRT_WAM_simulations,
    #     )

    #     for b in range(24):
    #         args.beam = b
    #         gerenate_fov_plot(
    #             args.rad,
    #             args.beam,
    #             args.event,
    #             args.model,
    #             not bool(args.event_study),
    #         )
    #         if args.model == "wam":
    #             execute_2DRT_WAM_simulations(args)
    #         if args.model == "waccm":
    #             execute_2DRT_WACCMX_simulations(args)
    #         if args.model == "iri":
    #             execute_2DRT_IRI_simulations(args)
    #         if args.model == "gitm":
    #             execute_2DRT_GITM_simulations(args)
    # if args.method == "dop":
    #     sys.path.append("rt/")
    #     from calculate_doppler import Doppler

    #     beams = np.arange(24)
    #     Doppler(args.event, args.rad, beams, args.model, parallel=True)
    # if args.method == "fan":
    #     sys.path.append("rt/")
    #     from calculate_doppler import Doppler

    #     beams = np.arange(24)
    #     dates = [dt.datetime(2017, 8, 21, 16), dt.datetime(2017, 8, 21, 20)]
    #     date = dates[0] + dt.timedelta(minutes=5)
    #     while date < dates[-1]:
    #         records = Doppler.fetch_by_scan_time(date, args.rad, args.model, beams)
    #         genererate_Fan(date, args.rad, 10.5, records, args.model)
    #         date += dt.timedelta(minutes=5)
    # if args.method == "rti":
    #     records = Doppler.fetch_by_beam(args.event, args.rad, args.model, args.beam)
    # if args.method == "movie":
    #     sys.path.append("py/")
    #     import utils

    #     folder = os.path.join("figures", args.event.strftime("%b%Y"))
    #     utils.create_mp4(folder, "fan*.png", "movie.mp4")
