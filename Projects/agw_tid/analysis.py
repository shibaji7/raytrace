#!/usr/bin/env python3

"""analysis.py: Running the analsysis"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "chakras4@erau.edu"
__status__ = "Research"

import argparse
import datetime as dt
import os

import numpy as np
from agw_utils import read_all_rays
from dateutil import parser as dparser
from loguru import logger

from rt import radar, utils

if __name__ == "__main__":
    pass