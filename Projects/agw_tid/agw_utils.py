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
