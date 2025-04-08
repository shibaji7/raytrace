#!/usr/bin/env python3

"""ro_analysis.py: Running the main body of the simulation for COSMIC RO"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "chakras4@erau.edu"
__status__ = "Research"

import xarray as xr
from dateutil import parser as dparser
from loguru import logger

if __name__ == "__main__":
    logger.info("Running RO analysis....")
