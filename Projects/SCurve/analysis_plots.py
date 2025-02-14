#!/usr/bin/env python3

"""analysis_plots.py: Analyze the eclipse related SCurve runs"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "chakras4@erau.edu"
__status__ = "Research"

def plot_ts(dop, events, fig_title="",filepath=""):
    from rt.rti import TimeSeriesPlot
    import pandas as pd
    from loguru import logger

    ts = TimeSeriesPlot([events[0], events[-1]], fig_title, num_subplots=2)
    dop.time = pd.to_datetime(dop.time)
    ax = ts.addParamPlot(dop.time, dop.frq_dne+dop.frq_dh, lcolor="b", kind="scatter")
    ax.set_ylim(-3,3)
    logger.info(f"File: {filepath}")
    ts.save(filepath)
    ts.close()
    return