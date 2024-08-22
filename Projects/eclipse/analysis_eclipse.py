#!/usr/bin/env python3

"""analysis_eclipse.py: Analyze the eclipse related items"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import sys
import os
from loguru import logger
from dateutil import parser as dparser
import argparse
import datetime as dt

CD_STEPS = ""
_DIR_ = "figures/movies/"

def add_sys_paths():
    """Adding /rt to sys path
    """
    global CD_STEPS
    pwd = os.getcwd()
    last_dir = pwd.split("/")[-1]
    index_of_trace_dir = pwd.split("/").index("raytrace")
    logger.info(f"Currently in '/{last_dir}', index of trace folder: {index_of_trace_dir}")
    local_libs = [
        os.path.join(
            "/".join((pwd.split("/")[:index_of_trace_dir+1])),
            "rt"
        ),
        os.path.join(
            "/".join((pwd.split("/")[:index_of_trace_dir+1])),
            "rt", "density"
        )
    ]
    logger.info(f"Loacl lib-{local_libs}")
    CD_STEPS = "".join(["../"]*int(len(pwd.split("/"))-index_of_trace_dir-1))
    sys.path.extend(local_libs)
    os.makedirs(
        os.path.join(CD_STEPS, _DIR_), 
        exist_ok=True
    )
    return

add_sys_paths()
import utils

def create_movies(cfg):
    return