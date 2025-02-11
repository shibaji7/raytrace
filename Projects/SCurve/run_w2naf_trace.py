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
import concurrent.futures
import datetime as dt
import os
import sys

import pandas as pd
from dateutil import parser as dparser
from loguru import logger

sys.path.extend([".", "rt/", "rt/density/", "Projects/SCurve/"])
from hamsci import HamSci

import rt.utils as utils


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
        self.model = self.cfg.model
        self.target_call_sign = self.cfg.ray_target.station_name
        self.source = dict(
            call_sign=self.cfg.ray_source.station_name,
            lat=self.cfg.ray_source.lat,
            lon=self.cfg.ray_source.lon,
        )
        self.hamsci = HamSci(
            self.cfg,
            self.start_time,
            [self.start_time, self.end_time],
            10 * 1e6,
        )
        self.hamsci.setup_pandas_dataset()
        self.hamsci.generate_fov(
            self.start_time,
            call_signs=[self.cfg.ray_target.station_name],
            source=self.source,
        )
        self.worker = self.cfg.worker
        self.parallel = self.cfg.worker > 0
        self.base_output_folder = os.path.join(
            self.cfg.project_save_location, self.cfg.project_name
        )
        self.target = dict(
            lat=self.hamsci.gds[self.target_call_sign.upper()].meta["lat"],
            lon=self.hamsci.gds[self.target_call_sign.upper()].meta["lon"],
            call_sign=self.target_call_sign,
        )
        return

    def get_event_dates(self):
        if self.cfg.time_gaps <= 0:
            events = self.eden_model.dates
            end_time = self.start_time + dt.timedelta(minutes=self.cfg.time_window)
            events = [e for e in events if e >= self.start_time and e <= end_time]
            self.cfg.time_gaps = (events[1] - events[0]).total_seconds() / 60
        else:
            events = (
                self.eden_model.dates
                if self.model == "gemini"
                else [
                    self.start_time + dt.timedelta(minutes=d * self.cfg.time_gaps)
                    for d in range(int(self.cfg.time_window / self.cfg.time_gaps))
                ]
            )
        if self.model == "gemini":
            events = events[: self.cfg.time_window]
        return events

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
            self.eden_model = SAMI3(self.cfg, self.start_time, 60)
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
            for i, event in enumerate(events):
                logger.info(f"Load e-Density for, {event}")
                self._run_rt_(event)
        return

    def _run_rt_(self, event):
        from rt.rays import PlotRays
        from rt.rt2d import HamSCI2dTrace

        rto = HamSCI2dTrace(
            event,
            self.source,
            self.target,
            self.cfg,
            self.model,
            self.base_output_folder,
        )

        # # Load all the electron density
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

        # # Load and compile MATLAB - PHaRLAP
        if not os.path.exists(rto.sim_fname):
            rto.compile(eden)
        else:
            rto.load_rto(eden)

        # # Create RT figures
        if not os.path.exists(rto.fig_name):
            plot = PlotRays(
                event,
                self.cfg,
                rto,
                f"{self.source['call_sign']}-{self.target['call_sign']}",
                0,
            )
            plot.lay_rays(
                kind=self.cfg.ray_trace_plot_kind, tag_distance=rto.gc_distance
            )
            plot.save(rto.fig_name)
            plot.close()
        return

    def compute_doppler(self):
        from rt.doppler import HamSCIDoppler

        # Initialize Doppler object
        self.dop = HamSCIDoppler(
            self.source,
            self.target,
            cfg=self.cfg,
            start_time=self.start_time,
            model=self.model,
            del_time=self.cfg.time_gaps,
            base=self.base_output_folder,
        )
        events = self.get_event_dates()
        for now, prev in zip(events[1:], events[:-1]):
            logger.info(
                f"Compute doppler for, {now.strftime('%H:%M')}/{prev.strftime('%H:%M')}"
            )
            self.dop._compute_doppler_from_prev_time_(now, prev, to_meters=1e3)
        return

    def bandpass_filter(self, data, lowcut=0.00015, highcut=0.001, fs=1, order=5):
        from scipy.signal import butter, filtfilt

        nyquist = 0.5 * fs  # Nyquist frequency
        low = lowcut / nyquist
        high = highcut / nyquist
        b, a = butter(order, [low, high], btype="band")
        filtered_data = filtfilt(b, a, data)
        return filtered_data

    def generate_ls(self):
        from rt.doppler import HamSCIDoppler
        from rt.rti import TimeSeriesPlot

        records = HamSCIDoppler.fetch_records(
            self.source,
            self.target,
            self.start_time,
            self.model,
            self.base_output_folder,
        )
        events = self.get_event_dates()
        fig_title = f"Model: {self.model.upper()} / {self.source['call_sign']}-{self.target['call_sign']}, {self.cfg.frequency} MHz \t {self.start_time.strftime('%d %b, %Y')}"
        ts = TimeSeriesPlot([events[0], events[-1]], fig_title, num_subplots=2)
        # records = records.groupby(by="time").mean().reset_index()
        records.time = pd.to_datetime(records.time)
        # records = (
        #     records.set_index("time")
        #     .resample("300S")
        #     .asfreq()
        #     .interpolate(method="cubic")
        #     .reset_index()
        # )
        data = self.hamsci.gds[self.target_call_sign.upper()].data["filtered"]["df"]
        from scipy.signal import detrend

        # ts.addParamPlot(
        #     data.UTC,
        #     self.bandpass_filter(detrend(data.Freq)),
        #     lcolor="k",
        #     # ax=ax,
        #     xlabel="",
        #     ylabel="",
        #     kind="scatter",
        # )
        ax = ts.addParamPlot(records.time, records.frq_dne, lcolor="b", kind="scatter")
        ts.addParamPlot(
            records.time,
            records.frq_dh,
            lcolor="r",
            ls="--",
            xlabel="",
            ylabel="",
            kind="scatter",
            ylim=[-5, 5],
        )
        # ts.addParamPlot(
        #     records.time,
        #     records.frq_dh + records.frq_dne,
        #     lcolor="k",
        #     ls="--",
        #     ax=ax,
        #     xlabel="",
        #     ylabel="",
        #     kind="scatter",
        # )
        filepath = (
            utils.get_hamsci_folder(
                self.source["call_sign"],
                self.start_time,
                self.model,
                self.base_output_folder,
                self.target["call_sign"],
            )
            + "/TS.png"
        )
        logger.info(f"File: {filepath}")
        ts.save(filepath)
        ts.close()
        self.save_simulation_TS(records)
        return

    def save_simulation_TS(self, records):
        o = pd.DataFrame()
        o = records[["time"]]
        o["frequency"] = self.cfg.frequency
        o["fd"] = records.frq_dne + records.frq_dh
        o.to_csv(
            utils.get_hamsci_folder(
                self.source["call_sign"],
                self.start_time,
                self.model,
                self.base_output_folder,
                self.target["call_sign"],
            )
            + f"/LoS_doppler_{self.cfg.frequency}.csv",
            index=False,
            header=True,
        )
        return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--cfg_file",
        default="cfg/rt2d_2024_eclipse_sami3.json",
        help="Configuration file",
        type=str,
    )
    args = parser.parse_args()
    logger.info("\n Parameter list for simulation ")
    for k in vars(args).keys():
        print("     ", k, "->", str(vars(args)[k]))

    sim = HamSCISimulation(args.cfg_file)
    sim.run_2d_simulation()
    sim.compute_doppler()
    sim.generate_ls()

    # TODO
    # 0. Send the movie/.pngs to groups [D]
    # 1. Check & plot where W2NAF is with respect to WWV [D]
    # 2. May need 2-hop GS [D]
    # 3. Plot / work on Simulated Data TS plots,
    #   a. Angle, Slant range relatated stats
    #   b. Run for 5 and 15 MHz
    # 4. Check a how dh vs dn works out. Need E field values.
