#!/usr/bin/env python
import datetime
import logging
import os

logger = logging.getLogger(__name__)

import pickle

import digital_rf as drf
import matplotlib as mpl
import numpy as np
from matplotlib import pyplot as plt
from scipy import signal

# from eclipse_calc import solarContext

mpl.rcParams["font.size"] = 20
mpl.rcParams["font.weight"] = "bold"
mpl.rcParams["axes.grid"] = True
mpl.rcParams["grid.linestyle"] = ":"
mpl.rcParams["figure.figsize"] = np.array([15, 8])
mpl.rcParams["axes.xmargin"] = 0
mpl.rcParams["figure.dpi"] = 240


def load_grape_drf(sDate, eDate, data_dir, channel="ch0"):
    # DATA LOADING #########################
    meta_dir = os.path.join(data_dir, channel, "metadata")
    do = drf.DigitalRFReader(data_dir)
    dmr = drf.DigitalMetadataReader(meta_dir)

    # "center_frequencies": array([ 2.5 ,  3.33,  5.  ,  7.85, 10.  , 14.67, 15.  , 20.  , 25.  ])

    # Get first and last sample index of data.
    # Index is the number of samples since the epoch (time_since_epoch*sample_rate)
    s0, s1 = do.get_bounds(channel)
    first_sample, last_sample = dmr.get_bounds()
    print("metadata bounds are %i to %i" % (first_sample, last_sample))

    start_idx = int(np.uint64(first_sample))
    print("computed start_idx = ", start_idx)

    fields = dmr.get_fields()
    print("Available fields are <%s>" % (str(fields)))

    print("first read - just get one column ")
    data_dict = dmr.read(start_idx, start_idx + 2, "center_frequencies")
    for key in data_dict.keys():
        #  print((key, data_dict[key]))
        freqList = data_dict[key]
        print("freq = ", freqList[0])

    data_dict = dmr.read(start_idx, start_idx + 2, "lat")
    for key in data_dict.keys():
        #   print((key, data_dict[key]))
        theLatitude = data_dict[key]
        print("Latitude: ", theLatitude)

    data_dict = dmr.read(start_idx, start_idx + 2, "long")
    for key in data_dict.keys():
        #  print((key, data_dict[key]))
        theLongitude = data_dict[key]
        print("Longitude: ", theLongitude)

    latest_meta = dmr.read_latest()
    latest_inx = list(latest_meta.keys())[0]
    cntr_freqs = latest_meta[latest_inx]["center_frequencies"]
    properties = do.get_properties(channel)
    fs = properties["samples_per_second"]

    sinx_0 = drf.util.time_to_sample(sDate, fs)
    sinx_1 = drf.util.time_to_sample(eDate, fs)
    blks = do.get_continuous_blocks(sinx_0, sinx_1, channel)
    for sinx, nsamps in blks.items():
        break  # Get the first sample index and number of samples

    bigarray_dct = {}
    for cfreq in cntr_freqs:
        bigarray_dct[cfreq] = np.zeros(nsamps, dtype=complex)

    for cfreq_inx, (cfreq, bigarray) in enumerate(bigarray_dct.items()):
        print()
        print("Working on {!s} MHz....".format(cfreq))
        data = do.read_vector(sinx, nsamps, channel)
        bigarray[:] = data[:, cfreq_inx]

    result = {}
    result["bigarray_dct"] = bigarray_dct
    result["latest_meta"] = latest_meta[latest_inx]
    result["properties"] = properties
    t0 = drf.util.sample_to_datetime(sinx, fs)
    result["timevec_utc"] = [
        t0 + drf.util.samples_to_timedelta(x, fs) for x in range(nsamps)
    ]
    return result


class GrapeDRF(object):
    def __init__(
        self, sDate, eDate, station, output_dir=os.path.join("output", "grapeDRF")
    ):

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        sDate_str = sDate.strftime("%Y%m%d.%H%M")
        eDate_str = eDate.strftime("%Y%m%d.%H%M")

        event_fname = "{!s}-{!s}_{!s}_grapeDRF".format(sDate_str, eDate_str, station)
        png_fname = event_fname + ".png"
        png_fpath = os.path.join(output_dir, png_fname)
        ba_fpath = os.path.join(output_dir, event_fname + ".ba.pkl")
        data_dir = os.path.join("data", "psws_grapeDRF", station)
        self.station = station

        if not os.path.exists(ba_fpath):
            result = load_grape_drf(sDate, eDate, data_dir)
            with open(ba_fpath, "wb") as fl:
                pickle.dump(result, fl)
        else:
            print("Using cached file {!s}...".format(ba_fpath))
            with open(ba_fpath, "rb") as fl:
                result = pickle.load(fl)

        self.result = result
        self.cfreqs = list(result["bigarray_dct"].keys())
        self.fs = result["properties"]["samples_per_second"]
        self.sDate = sDate
        self.eDate = eDate
        self.data_dir = data_dir
        self.event_fname = event_fname
        self.output_dir = output_dir
        self.png_fpath = png_fpath
        self.cmap = mpl.colors.LinearSegmentedColormap.from_list(
            " ", ["black", "darkgreen", "green", "yellow", "red"]
        )
        self.spectrum_timevec = None

    def plot_figure(self, cfreqs=None, png_fpath=None, **kwargs):
        print("Now plotting {!s}...".format(self.event_fname))
        if cfreqs is None:
            cfreqs = self.cfreqs

        ncols = 1
        nrows = len(cfreqs)
        ax_inx = 0
        fig = plt.figure(figsize=(15, 4 * nrows))
        axes = []
        for cfreq in cfreqs:
            print("   {!s} MHz...".format(cfreq))
            ax_inx += 1
            ax = fig.add_subplot(nrows, ncols, ax_inx)
            self.plot_ax(cfreq, ax, fig, **kwargs)
            axes.append(ax)

        fig.tight_layout()
        if png_fpath is None:
            png_fpath = self.png_fpath

        fig.savefig(png_fpath, bbox_inches="tight")
        print(png_fpath)
        return fig, axes, png_fpath

    def plot_ax(
        self,
        cfreq,
        ax,
        fig,
        cmap=None,
        plot_colorbar=False,
        xlim=None,
        solar_lat=None,
        solar_lon=None,
        overlaySolarElevation=False,
        overlayEclipse=False,
    ):

        sDate = self.sDate
        eDate = self.eDate

        ylabel = []
        #        ylabel.append("{!s} MHz".format(cfreq))
        ylabel.append("Doppler Shift (Hz)")
        ax.set_ylabel("\n".join(ylabel))
        ax.set_xlabel("UTC")

        result = self.result
        props = result["properties"]
        bigarray = self.result["bigarray_dct"].get(cfreq)
        if bigarray is None:
            msg = "ERROR: No data for {!s} MHz".format(cfreq)
            ax.text(0.5, 0.5, msg, ha="center", va="center", transform=ax.transAxes)
            print(msg)
            return
        else:
            txt = self.station + f" / $f_0$={cfreq} MHz" + "\n"
            ax.text(0.05, 0.95, txt, ha="left", va="top", transform=ax.transAxes)

        f, t_spec, Sxx = signal.spectrogram(
            bigarray, fs=self.fs, nfft=1024, window="hann", return_onesided=False
        )
        if self.spectrum_timevec is None:
            # TODO: Make this more clean in the future.
            ts0 = min(result["timevec_utc"]).timestamp()
            ts1 = max(result["timevec_utc"]).timestamp()
            ts_vec = np.linspace(ts0, ts1, len(t_spec))
            self.spectrum_timevec = [
                datetime.datetime.utcfromtimestamp(x) for x in ts_vec
            ]

        f = (np.fft.fftshift(f)).astype(
            "float64"
        )  # Frequency needs to be in float64 for some reason...
        Sxx = np.fft.fftshift(Sxx, axes=0)
        Sxx_db = 10 * np.log10(Sxx)
        if cmap is None:
            cmap = self.cmap
        mpbl = ax.pcolormesh(self.spectrum_timevec, f, Sxx_db, cmap=cmap)

        if plot_colorbar:
            cbar = fig.colorbar(mpbl, label="PSD [dB]")

        odct = {"color": "white", "lw": 4, "alpha": 0.75}
        if overlaySolarElevation:
            sts = solarContext.solarTimeseries(sDate, eDate, solar_lat, solar_lon)
            sts.overlaySolarElevation(ax, **odct)

        if overlayEclipse:
            sts = solarContext.solarTimeseries(sDate, eDate, solar_lat, solar_lon)
            sts.overlayEclipse(ax, **odct)

        if xlim is None:
            xlim = (sDate, eDate)

        xticks = ax.get_xticks()
        xtkls = []
        for xtk in xticks:
            dt = mpl.dates.num2date(xtk)
            xtkl = dt.strftime("%H:%M")
            xtkls.append(xtkl)
        ax.set_xticks(xticks)
        ax.set_xticklabels(xtkls)


if __name__ == "__main__":
    station = "w2naf"
    sDate = datetime.datetime(2024, 4, 8)
    eDate = datetime.datetime(2024, 4, 9)

    figd = {}
    figd["cfreqs"] = [20, 15, 10, 5]
    #    figd["cfreqs"]                  = [10]
    figd["solar_lat"] = 41.335116
    figd["solar_lon"] = -75.600692
    figd["overlaySolarElevation"] = False
    figd["overlayEclipse"] = False

    gDRF = GrapeDRF(sDate, eDate, station)
    gDRF.plot_figure(**figd)
