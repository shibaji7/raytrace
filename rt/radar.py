"""
    This python module is used to read the dataset from fitacf/fitacf3 
    level dataset.
"""

import bz2
import datetime as dt
import glob
import os
from types import SimpleNamespace

import eclipse
import numpy as np
import pandas as pd
import pydarn
from loguru import logger
from tqdm import tqdm


class Radar(object):

    def __init__(
        self,
        rad: str,
        dates: list = None,
        cfg: SimpleNamespace = None,
    ):
        logger.info(f"Initialize radar: {rad}/{dates[0]}")
        self.rad = rad
        self.dates = dates
        self.cfg = cfg
        self.type = cfg.sd_file_type
        self.tmp_sd_folder = cfg.tmp_sd_folder
        tqdm.pandas()
        self.__setup__()
        self.__fetch_data__()
        if cfg.eclipse_shadow:
            self.create_eclipse_shadow()
        return

    def __setup__(self):
        logger.info(f"Setup radar: {self.rad}")
        self.files = glob.glob(
            os.path.join(
                self.cfg.sd_base_file_loc.format(
                    year=self.dates[0].year,
                    rad=self.rad,
                    type=self.type,
                    date=self.dates[0].strftime("%Y%m%d"),
                )
            )
        )
        self.files.sort()
        self.hdw = pydarn.read_hdw_file(self.rad)
        self.fov = pydarn.Coords.GEOGRAPHIC(self.hdw.stid)
        logger.info(f"Files: {len(self.files)}")
        return

    def get_lat_lon_along_beam(self, beam):
        lats, lons = self.fov[0], self.fov[1]
        return lats[:, beam], lons[:, beam]

    def __fetch_data__(self):
        self.fname = os.path.join(
            self.cfg.tmp_sd_folder,
            f"{self.rad}.{self.type}.{self.dates[0].strftime('%Y%m%d')}.csv",
        )
        os.makedirs(self.cfg.tmp_sd_folder, exist_ok=True)
        if os.path.exists(self.fname):
            self.df = pd.read_csv(self.fname, parse_dates=["time"])
        else:
            records = []
            for f in self.files:
                logger.info(f"Reading file: {f}")
                with bz2.open(f) as fp:
                    reader = pydarn.SuperDARNRead(fp.read(), True)
                    records += reader.read_fitacf()
            if len(records) > 0:
                self.__tocsv__(records)
        if "lat" not in self.df.columns:
            self.df.tfreq = np.round(np.array(self.df.tfreq) / 1e3, 1)
            self.update_location_details()
        self.check_the_sounding_mode()
        return

    def create_eclipse_shadow(self):
        logger.info(f"Create shadow!")
        lats, lons = (self.fov[0].T, self.fov[1].T)
        folder = os.path.join(
            self.cfg.base_eclipse_folder, self.dates[0].strftime("%Y-%m-%d")
        )
        os.makedirs(folder, exist_ok=True)
        dates = [
            self.dates[0].replace(hour=0) + dt.timedelta(minutes=i) for i in range(1440)
        ]
        for bm in range(self.hdw.beams):
            file = f"{folder}/oc.{self.rad}.{bm}.csv"
            if not os.path.exists(file):
                p = eclipse.get_rti_eclipse(dates, lats[bm, :], lons[bm, :])
                frame = pd.DataFrame()
                frame["dates"] = dates
                for g in range(p.shape[1]):
                    frame[f"gate_{g}"] = p[:, g]
                frame.to_csv(file, header=True, index=False)
        return

    def __latlon__(self, row):
        lat, lon = self.fov[0].T, self.fov[1].T
        row["lat"], row["lon"] = (lat[row.bmnum, row.slist], lon[row.bmnum, row.slist])
        p = eclipse.get_eclipse(row.time, [300], [row.lat], [row.lon])[0]
        row["occul"] = p[0, 0, 0, 0]
        return row

    def check_the_sounding_mode(self):
        frequency, scan_time, beams = set(), set(), set()
        txt = ""
        for bm in self.df.bmnum.unique():
            beams.add(bm)
            o = self.df[(self.df.bmnum == bm)].groupby(by="time").mean().reset_index()
            x = np.rint((o.time.iloc[1] - o.time.iloc[0]).total_seconds())
            scan_time.add(x)
            tf = set()
            for t in o.tfreq:
                frequency.add(np.rint(t))
                tf.add(str(np.rint(t)))
            txt += f"Beam: {bm}, t={x}, f={','.join(list(tf))}\n"
        return

    def __tocsv__(self, records):
        (
            time,
            v,
            slist,
            p_l,
            frang,
            scan,
            beam,
            w_l,
            gflg,
            elv,
            phi0,
            tfreq,
            rsep,
            skynoise,
            nrang,
        ) = ([], [], [], [], [], [], [], [], [], [], [], [], [], [], [])
        for r in records:
            if "v" in r.keys():
                t = dt.datetime(
                    r["time.yr"],
                    r["time.mo"],
                    r["time.dy"],
                    r["time.hr"],
                    r["time.mt"],
                    r["time.sc"],
                    r["time.us"],
                )
                time.extend([t] * len(r["v"]))
                tfreq.extend([r["tfreq"]] * len(r["v"]))
                rsep.extend([r["rsep"]] * len(r["v"]))
                frang.extend([r["frang"]] * len(r["v"]))
                nrang.extend([r["nrang"]] * len(r["v"]))
                skynoise.extend([r["noise.sky"]] * len(r["v"]))
                scan.extend([r["scan"]] * len(r["v"]))
                beam.extend([r["bmnum"]] * len(r["v"]))
                v.extend(r["v"])
                gflg.extend(r["gflg"])
                slist.extend(r["slist"])
                p_l.extend(r["p_l"])
                w_l.extend(r["w_l"])
                # if "elv" in r.keys(): elv.extend(r["elv"])
                # if "phi0" in r.keys(): phi0.extend(r["phi0"])

        self.df = pd.DataFrame()
        self.df["v"] = v
        self.df["gflg"] = gflg
        self.df["slist"] = slist
        self.df["bmnum"] = beam
        self.df["p_l"] = p_l
        self.df["w_l"] = w_l
        if len(elv) > 0:
            self.df["elv"] = elv
        if len(phi0) > 0:
            self.df["phi0"] = phi0
        self.df["time"] = time
        self.df["tfreq"] = tfreq
        self.df["scan"] = scan
        self.df["rsep"] = rsep
        self.df["frang"] = frang
        self.df["nrang"] = nrang
        self.df["noise.sky"] = skynoise

        if self.dates:
            self.df = self.df[
                (self.df.time >= self.dates[0]) & (self.df.time <= self.dates[1])
            ]
        self.df.to_csv(self.fname, index=False, header=True)
        return

    def to_csv(
        self,
        names=[
            "time",
            "bmnum",
            "gflg",
            "noise.sky",
            "nrang",
            "p_l",
            "scan",
            "slist",
            "tfreq",
            "v",
            "w_l",
        ],
        params={
            "time": "Datetime of the observation (U: datetime)",
            "bmnum": "Beam number (U: 1)",
            "gflg": "Ground scatter flag (U: 1)",
            "noise.sky": "Sky noise (U: 1)",
            "nrang": "Number of range gates (U: 1)",
            "p_l": "Power (U: dB)",
            "scan": "Scan flag (U: 1)",
            "slist": "Range gate (U: 1)",
            "tfreq": "Operating frequency (U: kHz)",
            "v": "LoS Doppler Velocity (U: m/s)",
            "w_l": "Spectral width (U: m/s)",
        },
    ):
        o = self.df[names].copy()
        txt = "=======================================\n Parameter Description:"
        for n in names:
            txt += f" {n} - {params[n]}\n"
        txt += "=======================================\n"
        datetext = self.dates[0].strftime("%Y%m%d")
        fname = f"database/{self.rad}.{self.type}.{datetext}.csv"
        o.to_csv(fname, index=False, header=True)
        with open(fname, "r+") as fp:
            content = fp.read()
            fp.seek(0, 0)
            fp.write(txt.rstrip("\r\n") + "\n" + content)
        return

    def update_location_details(self):
        logger.info(f"Inside Location details!")
        self.df = self.df.progress_apply(self.__latlon__, axis=1)
        self.df.to_csv(self.fname, index=False, header=True)
        return


if __name__ == "__main__":
    import json
    fname = "cfg/rt2D.json"
    with open(fname, "r") as f:
        cfg = json.load(f, object_hook=lambda x: SimpleNamespace(**x))
    dates = [dt.datetime(2024,4,8), dt.datetime(2024,4,9)]
    Radar("bks", dates, cfg)
    Radar("fhe", dates, cfg)
    Radar("fhw", dates, cfg)
    Radar("kap", dates, cfg)
    Radar("gbr", dates, cfg)
    dates = [dt.datetime(2017,8,21), dt.datetime(2017,8,22)]
    Radar("fhe", dates, cfg)
    Radar("fhw", dates, cfg)
    Radar("cve", dates, cfg)
    Radar("cvw", dates, cfg)
    dates = [dt.datetime(2023,10,14), dt.datetime(2023,10,15)]
    Radar("fhe", dates, cfg)
    Radar("fhw", dates, cfg)
    Radar("cve", dates, cfg)
    Radar("cvw", dates, cfg)
    dates = [dt.datetime(2021,12,4), dt.datetime(2021,12,5)]
    Radar("fir", dates, cfg)