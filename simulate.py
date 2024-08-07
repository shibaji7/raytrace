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
import concurrent.futures
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
        time_window: int = 240,
        time_gaps: int = 5,
        cfg_file: str = "cfg/rt2d.json",
        worker: int = 6,
        parallel: bool = True,
    ) -> None:
        self.model = model
        self.rad = rad
        self.beam = beam
        self.time_gaps = time_gaps
        self.time_window = time_window
        self.start_time = start_time
        self.cfg_file = cfg_file
        self.worker = worker
        self.parallel = parallel
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
        if self.cfg.iri_param.eclipse:
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
        from gitm import GITM2d
        from iri import IRI2d
        from waccm import WACCMX2d

        if self.model == "iri":
            self.eden_model = IRI2d(self.cfg, self.start_time)
        elif self.model == "gitm":
            self.eden_model = GITM2d(self.cfg, self.start_time)
        elif self.model == "waccm":
            self.eden_model = WACCMX2d(self.cfg, self.start_time)
        else:
            raise ValueError(
                "Currently supporting following methods: iri, gitm, waccm-x, and wamipe"
            )
        events = [
            self.start_time + dt.timedelta(minutes=d * self.time_gaps)
            for d in range(int(self.time_window / self.time_gaps))
        ]
        if self.parallel:
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=self.worker
            ) as executor:
                _ = list(executor.map(self._run_rt_, events))
        else:
            for event in events:
                logger.info(f"Load e-Density for, {event}")
                self._run_rt_(event)
        return

    def _run_rt_(self, event):
        import plots
        from rt2d import RadarBeam2dTrace

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
            eden, _ = self.eden_model.fetch_dataset(
                event,
                rto.bearing_object["lat"],
                rto.bearing_object["lon"],
                rto.bearing_object["ht"],
                to_file=rto.edensity_file,
            )
        else:
            eden = self.eden_model.load_from_file(rto.edensity_file)

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
                rf"{self.model.upper()} + {self.rad.upper()}/{str(self.beam)}, $f_0$={str(self.cfg.frequency)} MHz",
                maxground=self.cfg.max_ground_range_km,
            )
        return

    def compute_doppler(self):
        import utils
        from doppler import Doppler
        from rti import RangeTimeIntervalPlot

        # Initialize Doppler object
        self.dop = Doppler(
            cfg=self.cfg,
            start_time=self.start_time,
            model=self.model,
            radar=self.radar,
            beam=self.beam,
            del_time=self.time_gaps,
            base=self.base_output_folder,
        )
        events = [
            self.start_time + dt.timedelta(minutes=d * self.time_gaps)
            for d in range(int(self.time_window / self.time_gaps))
        ]
        for now, prev in zip(events[1:], events[:-1]):
            logger.info(
                f"Compute doppler for, {now.strftime('%H:%M')}/{prev.strftime('%H:%M')}"
            )
            self.dop._compute_doppler_from_prev_time_(now, prev)

        fig_title = f"Model: {self.model.upper()} / {self.rad.upper()}-{'%02d'%self.beam}, {self.cfg.frequency} MHz"
        rtint = RangeTimeIntervalPlot(
            100, [events[0], events[-1]], self.rad, fig_title=fig_title, num_subplots=1
        )
        records = Doppler.fetch_by_beam(
            self.start_time,
            self.rad,
            self.model,
            self.beam,
            self.base_output_folder,
            frange=self.cfg.frange,
            rsep=self.cfg.rsep,
        )
        if len(records) > 0:
            rtint.addParamPlot(
                records, self.beam, title="", zparam="vel_tot", lay_eclipse=True
            )
        rtint.save(
            filepath=utils.get_folder(
                self.rad,
                self.beam,
                self.start_time,
                self.model,
                self.base_output_folder,
            )
            + "/rti.png"
        )
        rtint.close()
        return

    @staticmethod
    def genererate_fan(event, rad, tfreq, frame, model, param="v"):
        import cartopy
        from fan import Fan

        folder = os.path.join("figures", event.strftime("%b%Y"))
        os.makedirs(folder, exist_ok=True)
        cmap = "Spectral"
        fan = Fan(
            rad,
            event,
            fig_title="",
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
            p_max=10,
            p_min=-10,
            cmap=cmap,
            label="",
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
    parser.add_argument("-r", "--rad", default="fhw", help="Radar code (fhw)")
    parser.add_argument(
        "-bm", "--beam", default=7, type=int, help="Radar beam (default 7)"
    )
    parser.add_argument(
        "-w", "--worker", default=6, type=int, help="Worker numbers [6]"
    )
    parser.add_argument(
        "-ev",
        "--event",
        default=dt.datetime(2017, 8, 21, 17),
        help="Event date for simulation [YYYY-mm-ddTHH:MM]",
        type=dparser.isoparse,
    )
    parser.add_argument(
        "-tw",
        "--time_window",
        default=6,
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
        "-f", "--cfg_file", default="cfg/rt2d.json", help="Configuration file"
    )
    parser.add_argument("-md", "--method", default="rt", help="Method rt/fan")
    args = parser.parse_args()
    logger.info("\n Parameter list for simulation ")
    for k in vars(args).keys():
        print("     ", k, "->", str(vars(args)[k]))
    if args.method == "rt":
        import radar

        beams = radar.get_beams(args.rad) if args.beam == -1 else [args.beam]
        for beam in beams:
            rsim = RadarSimulation(
                model=args.model,
                start_time=args.event,
                rad=args.rad,
                beam=beam,
                time_window=args.time_window,
                time_gaps=args.time_gaps,
                cfg_file=args.cfg_file,
                worker=args.worker,
            )
            rsim.gerenate_fov_plot()
            rsim.run_2d_simulation()
            rsim.compute_doppler()
    if args.method == "fan":
        import radar
        from doppler import Doppler

        beams = radar.get_beams(args.rad)
        dates = [args.event, args.event + dt.timedelta(minutes=args.time_window)]
        date = dates[0] + dt.timedelta(minutes=args.time_gaps)
        with open(args.cfg_file, "r") as f:
            cfg = json.load(f, object_hook=lambda x: SimpleNamespace(**x))
        base = os.path.join(cfg.project_save_location, cfg.project_name)
        while date < dates[-1]:
            records = Doppler.fetch_by_scan_time(
                date, args.rad, args.model, beams, base
            )
            RadarSimulation.genererate_fan(
                date, args.rad, cfg.frequency, records, args.model
            )
            date += dt.timedelta(minutes=args.time_gaps * 3)
            # break
