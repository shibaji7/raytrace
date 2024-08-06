import eclipse
import iricore
import numpy as np
from loguru import logger
from scipy.io import savemat


class IRI2d(object):

    def __init__(
        self,
        cfg,
        event,
    ):
        self.cfg = cfg
        self.event = event
        self.iri_version = self.cfg.iri_version

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
        p = 1.0 - oclt
        logger.info(f"Min/max value P: {np.nanmax(p)}/{np.nanmin(p)}")
        self.param_val = self.param_val * p
        return

    def fetch_dataset(
        self,
        time,
        lats,
        lons,
        alts,
        to_file=None,
    ):
        self.lats, self.alts, self.lons = (lats, alts, lons)
        self.time = time
        iri_version = iri_version if iri_version else self.iri_version
        self.param = np.zeros((len(self.alts), len(self.lats)))
        for i in range(len(self.lats)):
            iriout = iricore.iri(
                self.time,
                self.alt_range,
                self.lats[i],
                self.lons[i],
                iri_version,
            )
            self.param_val[:, i] = iriout.edens * 1e-6
        if self.cfg.eclipse and self.time>=self.cfg.start_mask_time: 
            self.load_eclipse()
        if to_file:
            savemat(to_file, dict(ne=self.param))
        return self.param_val, self.alts