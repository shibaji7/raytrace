#!/usr/bin/env python

"""utils.py: utility module to support other functions."""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import glob
import json
import math
import os
from pathlib import Path
from types import SimpleNamespace

import cv2
import numpy as np
import scipy.io as spio
from loguru import logger
from scipy.interpolate import interp1d


pconst = {
    "kconst": 80.6, # Converstion constant
    "boltz": 1.38066e-23,  # Boltzmann constant  in Jule K^-1
    "h": 6.626e-34,  # Planks constant  in ergs s
    "c": 2.9979e08,  # in m s^-1
    "avo": 6.023e23,  # avogadro's number
    "Re": 6371.0e3,
    "amu": 1.6605e-27,
    "q_e": 1.602e-19,  # Electron charge in C
    "m_e": 9.109e-31,  # Electron mass in kg
    "g": 9.81,  # Gravitational acceleration on the surface of the Earth
    "eps0": 1e-9 / (36 * np.pi),
    "R": 8.31,  # J mol^-1 K^-1
}

def read_params_2D(fname: str = None):
    logger.info(f"Load config files: {fname}")
    with open(fname, "r") as f:
        cfg = json.load(f, object_hook=lambda x: SimpleNamespace(**x))
    return cfg


def get_gridded_parameters(
    q, xparam="beam", yparam="slist", zparam="v", r=0, rounding=True
):
    """
    Method converts scans to "beam" and "slist" or gate
    """
    plotParamDF = q[[xparam, yparam, zparam]]
    if rounding:
        plotParamDF.loc[:, xparam] = np.round(plotParamDF[xparam].tolist(), r)
        plotParamDF.loc[:, yparam] = np.round(plotParamDF[yparam].tolist(), r)
    plotParamDF = plotParamDF.groupby([xparam, yparam]).mean().reset_index()
    plotParamDF = plotParamDF[[xparam, yparam, zparam]].pivot(
        index=xparam, columns=yparam
    )
    x = plotParamDF.index.values
    y = plotParamDF.columns.levels[1].values
    X, Y = np.meshgrid(x, y)
    # Mask the nan values! pcolormesh can't handle them well!
    Z = np.ma.masked_where(
        np.isnan(plotParamDF[zparam].values), plotParamDF[zparam].values
    )
    return X, Y, Z


def get_folder(rad, beam, date, model, base):
    """
    Get folder by date
    """
    import os

    fold = os.path.join(base, date.strftime("%Y-%m-%d"), rad, "%02d" % beam, model)
    os.makedirs(fold, exist_ok=True)
    return fold


def read_params_2D(fname):
    with open(fname, "r") as f:
        param = json.load(f, object_hook=lambda x: SimpleNamespace(**x))
    return param


def clean():
    files = glob.glob(str(Path.home() / "matlab_crash_dump*"))
    for f in files:
        if os.path.exists(f):
            os.remove(f)
    return


def extrap1d(x, y, kind="linear"):
    """This method is used to extrapolate 1D paramteres"""
    interpolator = interp1d(x, y, kind=kind)
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0] + (x - xs[0]) * (ys[1] - ys[0]) / (xs[1] - xs[0])
        elif x > xs[-1]:
            return ys[-1] + (x - xs[-1]) * (ys[-1] - ys[-2]) / (xs[-1] - xs[-2])
        else:
            return interpolator(x)

    def ufunclike(xs):
        return np.array(list(map(pointwise, np.array(xs))))

    return ufunclike


def interpolate_by_altitude(h, hx, param, scale="log", kind="cubic", method="intp"):
    if scale == "linear":
        pnew = (
            interp1d(h, param, kind=kind)(hx)
            if method == "intp"
            else extrap1d(h, param, kind=kind)(hx)
        )
    if scale == "log":
        pnew = (
            10 ** interp1d(h, np.log10(param), kind=kind)(hx)
            if method == "intp"
            else 10 ** extrap1d(h, np.log10(param), kind=kind)(hx)
        )
    return pnew


def smooth(x, window_len=11, window="hanning"):
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")
    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    if window_len < 3:
        return x
    if not window in ["flat", "hanning", "hamming", "bartlett", "blackman"]:
        raise ValueError(
            "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
        )
    s = np.r_[x[window_len - 1 : 0 : -1], x, x[-2 : -window_len - 1 : -1]]
    if window == "flat":
        w = np.ones(window_len, "d")
    else:
        w = eval("np." + window + "(window_len)")
    y = np.convolve(w / w.sum(), s, mode="valid")
    d = window_len - 1
    y = y[int(d / 2) : -int(d / 2)]
    return y


def create_movie(folder, outfile, pat, fps=1):
    """
    Create movies from pngs
    """
    files = glob.glob(f"{folder}/{pat}")
    files.sort()

    fourcc = cv2.VideoWriter_fourcc(*"XVID")
    img = cv2.imread(files[0])
    height, width, layers = img.shape
    size = (int(width / 2) * 2, int(height / 2) * 2)
    video = cv2.VideoWriter(folder + "/" + outfile, fourcc, fps, size)
    for idx in range(len(files)):
        im = cv2.imread(files[idx])
        im = cv2.resize(im, size)
        video.write(im)
        del im
    video.release()
    return


def create_mp4(folder, pat, outfile, fps=3):
    import os

    cmd = f"""
    ffmpeg -framerate {fps} -pattern_type glob -i '{folder}/{pat}' -c:v libx264 {folder}/{outfile}
    """
    os.system(cmd)
    return


def loadmatlabfiles(filename):
    """
    Load a .mat file and convert mat-objects to nested dictionaries.
    """
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_key_(data)


def _check_key_(data_dict):
    """
    Recursively check and convert mat-objects in a dictionary to nested dictionaries.
    """
    for key, value in data_dict.items():
        if isinstance(value, spio.matlab.mio5_params.mat_struct):
            data_dict[key] = _todict_(value)
    return data_dict


def _todict_(matobj):
    """
    Recursively construct nested dictionaries from mat-objects.
    """
    return {
        strg: (
            _todict_(getattr(matobj, strg))
            if isinstance(getattr(matobj, strg), spio.matlab.mio5_params.mat_struct)
            else getattr(matobj, strg)
        )
        for strg in matobj._fieldnames
    }


def setsize(size=8):
    import matplotlib as mpl

    mpl.rcParams.update(
        {"xtick.labelsize": size, "ytick.labelsize": size, "font.size": size}
    )
    return


def calculate_bearing(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])

    dlon = lon2 - lon1
    y = math.sin(dlon) * math.cos(lat2)
    x = math.cos(lat1) * math.sin(lat2) - math.sin(lat1) * math.cos(lat2) * math.cos(
        dlon
    )
    bearing = math.atan2(y, x)

    return (math.degrees(bearing) + 360) % 360, math.degrees(bearing)
