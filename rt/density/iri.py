import datetime as dt

import eclipse
import iricore
import numpy as np
from dateutil import parser as dparser
from loguru import logger
from scipy.io import loadmat, savemat


class IRI2d(object):

    def __init__(
        self,
        cfg,
        event: dt.datetime,
    ):
        self.cfg = cfg
        self.event = event
        self.iri_version = self.cfg.iri_param.iri_version
        if self.cfg.event_type.eclipse:
            self.start_mask_time = dparser.isoparse(self.cfg.iri_param.start_mask_time)
        return

    def load_eclipse(self):
        """
        Create Eclipse in the system
        """
        e = eclipse.Eclipse()
        logger.info(f"Inside eclipse mask: {self.time}")
        oclt = np.zeros([len(self.alts), len(self.lons)])
        for i in range(len(self.alts)):
            for j in range(len(self.lons)):
                oclt[i, j] = e.create_eclipse_shadow(
                    self.time, self.lats[j], self.lons[j], self.alts[i]
                )
        oclt[oclt >= 0.96] = 1.0
        logger.info(f"Min/max value O: {np.nanmax(oclt)}/{np.nanmin(oclt)}")
        oclt = np.nan_to_num(oclt)
        oclt = oclt * self.cfg.iri_param.eclipse_shadow_multiplier
        p = 1.0 - oclt
        # if np.nanmax(oclt) >= 0.8:
        #     p = (1.0 - oclt) * self.cfg.iri_param.eclipse_shadow_multiplier
        logger.info(f"Min/max value P: {np.nanmax(p)}/{np.nanmin(p)}")
        self.param = self.param * p
        return

    def fetch_dataset(
        self,
        time: dt.datetime,
        lats,
        lons,
        alts,
        to_file: str = None,
    ):
        self.lats, self.alts, self.lons = (lats, alts, lons)
        self.time = time
        self.param = np.zeros((len(self.alts), len(self.lats)))
        alt_range = [alts[0], alts[-1], alts[1] - alts[0]]
        for i in range(len(self.lats)):
            iriout = iricore.iri(
                self.time,
                alt_range,
                self.lats[i],
                self.lons[i],
                self.iri_version,
            )
            self.param[:, i] = iriout.edens * 1e-6
        if self.cfg.event_type.eclipse and self.time >= self.start_mask_time:
            self.load_eclipse()
        if to_file:
            savemat(to_file, dict(ne=self.param))
        return self.param, self.alts

    def load_from_file(self, to_file: str):
        logger.info(f"Load from file {to_file.split('/')[-1]}")
        self.param = loadmat(to_file)["ne"]
        return self.param
