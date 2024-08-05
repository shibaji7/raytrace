import numpy as np
from loguru import logger
from scipy.io import savemat
import eclipse
import iricore

class IRI2D(object):

    def __init__(
        self,
        time,
        lats,
        lons,
        alts,
        to_file=None,
        cfg=None,
        eclipse_prop=dict(),
        iri_version=20,
    ):
        """
        Setup all parameters
        """
        self.time = time
        self.lats = lats
        self.lons = lons
        self.alts = alts
        self.to_file = to_file
        self.eclipse_prop = eclipse_prop
        self.cfg = cfg
        self.alt_range = [alts[0], alts[-1], alts[1] - alts[0]]
        self.iri_version = iri_version
        self.load_data()
        if (
            "start_mask_time" in eclipse_prop
            and eclipse_prop["start_mask_time"] <= time
        ):
            self.load_eclipse()
        self.save_to_file()
        return

    def load_data(self, iri_version=None):
        """
        Load all IRI dataset
        """
        iri_version = iri_version if iri_version else self.iri_version
        self.param_val = np.zeros((len(self.alts), len(self.lats)))
        for i in range(len(self.lats)):
            iriout = iricore.iri(
                self.time,
                self.alt_range,
                self.lats[i],
                self.lons[i],
                iri_version,
            )
            self.param_val[:, i] = iriout.edens * 1e-6
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

    def save_to_file(self, to_file=None):
        """
        Save to file
        """
        to_file = to_file if to_file else self.to_file
        logger.info(f"IRI save to file: {to_file}")
        if to_file:
            mobj = {}
            mobj["ne"] = self.param_val
            savemat(to_file, mobj)
        return