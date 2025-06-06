#!/usr/bin/env python3

"""run_simulation.py: simulate python program for RT"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import concurrent.futures
import datetime as dt
import os

import numpy as np
from dateutil import parser as dparser
from loguru import logger


class RadarSimulation(object):
    def __init__(
        self,
        cfg_file: str,
        beam: int = 0,
    ) -> None:
        from rt import utils

        ## Kill all Matlab Server Hosts for this run
        os.system("killall MathWorksServiceHost.")

        self.cfg_file = cfg_file
        self.cfg = utils.read_params_2D(cfg_file)

        self.model = self.cfg.model
        self.rad = self.cfg.rad
        self.beam = beam if beam is not None else self.cfg.beam
        self.time_gaps = self.cfg.time_gaps
        self.time_window = self.cfg.time_window
        self.start_time = dparser.isoparse(self.cfg.event)
        self.event_type = self.cfg.event_type

        self.worker = self.cfg.worker
        self.parallel = self.cfg.worker > 0

        self.base_output_folder = os.path.join(
            self.cfg.project_save_location, self.cfg.project_name
        )
        self.load_radar()
        return

    def load_radar(self) -> None:
        from rt.radar import Radar

        self.radar = Radar(
            self.rad,
            [self.start_time, self.start_time + dt.timedelta(minutes=self.time_window)],
            self.cfg,
        )
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

        import rt.utils as utils
        from rt.fan import Fan

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
        if self.cfg.event_type.eclipse:
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
        from rt.density.gemini import GEMINI2d
        from rt.density.gitm import GITM2d
        from rt.density.iri import IRI2d
        from rt.density.sami3 import SAMI3
        from rt.density.waccm import WACCMX2d
        from rt.density.wamipe import WAMIPE2d

        if self.model == "iri":
            self.eden_model = IRI2d(self.cfg, self.start_time)
        elif self.model == "gitm":
            self.eden_model = GITM2d(self.cfg, self.start_time)
        elif self.model == "waccm":
            self.eden_model = WACCMX2d(self.cfg, self.start_time)
        elif self.model == "gemini":
            self.eden_model = GEMINI2d(self.cfg, self.start_time)
        elif self.model == "wamipe":
            self.eden_model = WAMIPE2d(self.cfg, self.start_time)
        elif self.model == "sami3":
            self.eden_model = SAMI3(self.cfg, self.start_time)
        else:
            raise ValueError(
                f"Currently supporting following methods: iri, gitm, waccm-x, and wamipe, and you provided '{self.model}'"
            )
        events = self.get_event_dates()
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
        from rt.rays import PlotRays
        from rt.rt2d import RadarBeam2dTrace

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
                intp_edens_xlim_index=rto.intp_edens_xlim_index,
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
            plot = PlotRays(event, self.cfg, rto, self.rad, self.beam)
            plot.lay_rays(kind=self.cfg.ray_trace_plot_kind, zoomed_in=[])
            plot.save(rto.fig_name)
            plot.close()
        return

    def compute_doppler(self):
        from rt.doppler import SuperDARNDoppler

        # Initialize Doppler object
        self.dop = SuperDARNDoppler(
            cfg=self.cfg,
            start_time=self.start_time,
            model=self.model,
            radar=self.radar,
            beam=self.beam,
            del_time=self.time_gaps,
            base=self.base_output_folder,
        )
        events = self.get_event_dates()
        for now, prev in zip(events[1:], events[:-1]):
            logger.info(
                f"Compute doppler for, {now.strftime('%H:%M')}/{prev.strftime('%H:%M')}"
            )
            self.dop._compute_doppler_from_prev_time_(now, prev)
        return

    def get_event_dates(self):
        events = [
            self.start_time + dt.timedelta(minutes=d * self.time_gaps)
            for d in range(int(self.time_window / self.time_gaps))
        ]
        # events = (
        #     self.eden_model.dates
        #     if self.model == "gemini"
        #     else [
        #         self.start_time + dt.timedelta(minutes=d * self.time_gaps)
        #         for d in range(int(self.time_window / self.time_gaps))
        #     ]
        # )
        # if self.model == "gemini":
        #     events = events[: self.cfg.time_window]
        return events

    def generate_rti(self):
        import rt.utils as utils
        from rt.doppler import SuperDARNDoppler
        from rt.rti import RangeTimeIntervalPlot

        events = self.get_event_dates()
        fig_title = f"Model: {self.model.upper()} / {self.rad.upper()}-{'%02d'%self.beam}, {self.cfg.frequency} MHz \t {self.start_time.strftime('%d %b, %Y')}"
        rtint = RangeTimeIntervalPlot(
            60,
            [events[0], events[-1]],
            self.rad,
            fig_title=fig_title,
            num_subplots=2,
            srange_type="srange",
        )
        rtint.addParamPlot(
            self.radar.df.copy(),
            self.beam,
            title="Observations",
            xlabel="",
            lay_eclipse=self.cfg.event_type.eclipse,
        )
        records = SuperDARNDoppler.fetch_by_beam(
            self.start_time,
            self.rad,
            self.model,
            self.beam,
            self.base_output_folder,
            frange=self.cfg.frange,
            rsep=self.cfg.rsep,
        )
        if len(records) > 0:
            records.time = records.time.apply(lambda x: dparser.isoparse(x))
            rtint.addParamPlot(
                records,
                self.beam,
                title="",
                zparam="vel_tot",
                lay_eclipse=self.cfg.event_type.eclipse,
                kind="scatter",
            )
        filepath = (
            utils.get_folder(
                self.rad,
                self.beam,
                self.start_time,
                self.model,
                self.base_output_folder,
            )
            + "/rti.png"
        )
        logger.info(f"File: {filepath}")
        rtint.save(filepath)
        rtint.close()
        if self.cfg.to_netcdf:
            fname = os.path.join(
                utils.get_folder(
                    self.rad,
                    self.beam,
                    self.start_time,
                    self.model,
                    self.base_output_folder,
                ),
                f"{self.start_time.strftime('%Y%m%d')}-{'%02d'%self.beam}.nc",
            )
            logger.info(f"Save to File: {fname}")
            if not os.path.exists(fname):
                self.radar.beam_to_netCDF(
                    self.beam,
                    fname,
                    model_frame=records,
                    model_params=["vel_tot"],
                )
        return

    @staticmethod
    def genererate_fan(cfg, date, cfg_file):
        import cartopy

        from rt import radar
        from rt.doppler import Doppler
        from rt.fan import Fan

        beams = radar.get_beams(cfg.rad)
        base = os.path.join(cfg.project_save_location, cfg.project_name)
        folder = os.path.join("figures", cfg.event.strftime("%b%Y"))
        os.makedirs(folder, exist_ok=True)
        radar = radar.Radar(
            cfg.rad, [cfg.event, cfg.event + dt.timedelta(minutes=cfg.time_window)], cfg
        )
        records = Doppler.fetch_by_scan_time(
            date, cfg.rad, cfg.model, beams, base, frange=cfg.frange, rsep=cfg.rsep
        )
        obs_records, scan_time, tf = radar.get_scan_by_time(date)

        lons = np.arange(-180, 180, 30)
        lats = (
            np.arange(30, 70, 15)
            if radar.hdw.geographic.lat > 0
            else np.arange(-70, -30, 15)
        )
        proj = cartopy.crs.Orthographic(
            10 * int(radar.hdw.geographic.lon / 10),
            10 * int(radar.hdw.geographic.lat / 10 - 1),
        )
        extent = [
            2 * int(radar.fov[1][: cfg.slant_gate_of_radar, :].min() / 2),
            2 * int(radar.fov[1][: cfg.slant_gate_of_radar, :].max() / 2),
            3 * int(radar.fov[0][: cfg.slant_gate_of_radar, :].min() / 3 - 3),
            3 * int(radar.fov[0][: cfg.slant_gate_of_radar, :].max() / 3 + 3),
        ]
        fan = Fan(
            [cfg.rad],
            date,
            ncols=2,
            fig_title="",
        )
        eclipse_lats, lay_eclipse = np.linspace(0, 90, num=90 * 2), None
        if cfg.event_type.eclipse:
            lay_eclipse = cfg_file
            if radar.hdw.geographic.lat < 0:
                eclipse_lats = np.linspace(-90, 0, num=90 * 2)
        fan.setup(lons, lats, extent=extent, proj=proj, lay_eclipse=lay_eclipse)
        fan.generate_fov(
            cfg.rad,
            obs_records,
            p_name="v",
            p_max=30,
            p_min=-30,
            cmap="Spectral",
            cbar=False,
            text_decription=dict(
                x=0.05,
                y=0.95,
                txt="Observations" + "\n" + rf"$f_0$={tf} MHz / $T_s$={scan_time} m",
                ha="left",
                va="center",
            ),
            lats=eclipse_lats,
        )
        fan.generate_fov(
            cfg.rad,
            records,
            p_name="vel_tot",
            p_max=30,
            p_min=-30,
            cmap="Spectral",
            label=r"Velocity ($ms^-1$)",
            text_decription=dict(
                x=0.05,
                y=0.95,
                txt=rf"Model: {cfg.model} / $f_0$={cfg.frequency} MHz",
                ha="left",
                va="center",
            ),
            lats=eclipse_lats,
        )
        fan.annotate_figure()
        fan.save(f"{folder}/fan.{cfg.rad}-{date.strftime('%H%M')}.png")
        fan.close()

        if cfg.to_netcdf:
            fname = os.path.join(folder, f"{date.strftime('%Y%m%d%H%M')}-{cfg.rad}.nc")
            logger.info(f"Save to File: {fname}")
            if not os.path.exists(fname):
                radar.scan_to_netCDF(
                    obs_records,
                    fname,
                    model_frame=records,
                    model_params=["vel_tot"],
                )
        return
