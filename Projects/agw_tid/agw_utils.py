"""agw_utils.py: AGW related utilities"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "chakras4@erau.edu"
__status__ = "Research"

import datetime as dt
import os
import sys

from dateutil import parser as dparser
from loguru import logger

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


def read_all_rays(cfg, beam, ray_type_to_readpower="gs"):
    import glob
    import os

    import pandas as pd

    import rt.utils as utils
    from rt.rays import Rays2D

    base = os.path.join(cfg.project_save_location, cfg.project_name)
    folder = utils.get_folder(
        cfg.rad,
        beam,
        cfg.event,
        cfg.model,
        base,
    )
    sim_fname_tag = folder + "/*_rt.mat"
    files = glob.glob(sim_fname_tag)
    logger.info(f"fetching files {sim_fname_tag}")
    files.sort()
    rays = []
    base_hr = cfg.event.hour
    for f in files:
        hm = f.split("/")[-1].replace("_rt.mat", "")
        date = (
            cfg.event
            + dt.timedelta(hours=int(hm[:2]) - base_hr)
            + dt.timedelta(minutes=int(hm[2:]))
        )
        r = Rays2D.read_rays(date, cfg, folder, f)
        rays.append(r.ray_power[ray_type_to_readpower])
    rays = pd.concat(rays)
    logger.info(f"Num of files loaded ... {len(files)}")
    return rays
