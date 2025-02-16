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

def plot_ts(dop, events, shadow, fig_title="",filepath=""):
    from rt.rti import TimeSeriesPlot
    import pandas as pd
    import numpy as np
    from loguru import logger

    ylim = [-3,3]
    ts = TimeSeriesPlot([events[0], events[-1]], fig_title, num_subplots=3)
    dop.time = pd.to_datetime(dop.time)
    ax = ts.addParamPlot(
        dop.time, dop.frq_dne+dop.frq_dh, lcolor="b", kind="scatter",
        title=r"(a) Total $\Delta$ f", ylim=ylim
    )
    axt = ax.twinx()
    axt.plot(dop.time.unique(), 1-np.nanmax(shadow, axis=1), color="k", ls="--")
    axt.plot(dop.time.unique(), 1-np.nanmin(shadow, axis=1), color="k", ls="--")
    axt.set_ylim(0,1)
    axt.set_yticks([])
    # ax = ts.addParamPlot(
    #     dop.time, dop.frq_dne, lcolor="r", kind="scatter",
    #     title=r"(b) Total $\Delta$ f($\eta$)", ylim=ylim
    # )
    # ax = ts.addParamPlot(
    #     dop.time, dop.frq_dh, lcolor="k", kind="scatter",
    #     title=r"(c) Total $\Delta$ f($h$)", ylim=ylim
    # )
    logger.info(f"File: {filepath}")
    ts.save(filepath)
    ts.close()
    return