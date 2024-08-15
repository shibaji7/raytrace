import h5py
import numpy as np
import math
import os
from loguru import logger

class WAMIPE2d(object):
    def __init__(
        self,
        cfg,
        event,
    ):
        self.cfg = cfg
        self.file_name = self.cfg.density_file_location
        self.event = event
        self.file_name = os.path.join(
            self.cfg.density_file_location, self.cfg.density_file_name
        )
        self.load_grid()
        return
    
    def load_grid(self):
        """
        Load grid file
        """
        logger.info(f"Loading WAM Grid files")
        self.ccord_file = (
            self.cfg.density_file_location + self.cfg.grid_coordinate_file
        )
        # array dimensions
        # number of grid points in magnetic longitude in the IPE model
        nmp = 80
        # number of magnetic flux tubes in the IPE model from the north pole.
        nlp = 170
        # number of grid points along a magnetic flux tube in the IPE model from the northern hemisphere at 90km.
        iDIM = 1115
        # 1:geographic eastward; 2:northward; 3:upward
        napex = 3
        apex_d1 = np.zeros((nmp,nlp,iDIM,napex))
        apex_d2 = np.zeros((nmp,nlp,iDIM,napex))
        apex_d3 = np.zeros((nmp,nlp,iDIM,napex))
        m2km = 1.0e-3
        rad2deg = 180./math.pi
        return
    
    def _read_params_from_hdf_file_(self):
        """Load parameters from file
        """
        # grid_name = self.cfg.wam_param.grid_name
        # with h5py.File(fpath, 'r') as file:
        #     # access datasets and attributes in the file

        return
    

