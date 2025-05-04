#!/usr/bin/env python

"""collision.py: collision is dedicated to find all types of collision functions."""

__author__ = "Chakraborty, S."
__copyright__ = "Chakraborty, S."
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "chakras4@erau.edu"
__status__ = "Research"

import datetime as dt
from dataclasses import dataclass

import numpy as np
from loguru import logger

from raidpy.constants import pconst


@dataclass
class Collision_en:
    N2: np.array = None
    O2: np.array = None
    O: np.array = None
    H: np.array = None
    He: np.array = None
    total: np.array = None


@dataclass
class Collision_ei:
    O2p: np.array = None
    Op: np.array = None
    total: np.array = None


@dataclass
class Collision_SN:
    en: Collision_en = None
    ei: Collision_ei = None
    total: np.array = None


@dataclass
class Collision:
    nu_ft: np.array = None
    nu_av_cc: np.array = None
    nu_av_mb: np.array = None
    nu_sn: Collision_SN = None


class ComputeCollision(object):
    """
    This class is a global class to estimate any types atmosphreic collision profile.

    msise = msise module
    iri = IRI module
    """

    def __init__(
        self, msise: dict, iri: dict, date: dt.datetime = None, _run_: bool = False
    ):
        self.msise = msise
        self.iri = iri
        self.collision = Collision()
        self.date = date
        if date:
            logger.info(f"Compute the collison for {date}")
        if _run_:
            self.collision.nu_sn = Collision_SN()
            self.collision.nu_sn.en = Collision_en()
            self.collision.nu_sn.ei = Collision_ei()
            self.collision.nu_sn.total = np.zeros_like(iri["etemp"])

            self.collision.nu_ft = self.calculate_FT_collision_frequency()
            self.collision.nu_av_cc = self.calculate_FT_collision_frequency(2.5)
            self.collision.nu_av_mb = self.calculate_FT_collision_frequency(1.5)
            self.calculate_SN_en_collision_frequency()
            self.calculate_SN_ei_collision_frequency()
        return

    def calculate_FT_collision_frequency(self, frac=1.0):
        """
        This method only provides the Friedrich-Tonker electron neutral collision frequency

        nn <float> = neutral density
        Tn <float> = neutral temperature
        Te <float> = electron temperatute


        nu <float> = collision frequency
        https://azformula.com/physics/dimensional-formulae/what-is-dimensional-formula-of-temperature/
        """
        logger.info(
            f"Compute the Friedrich-Tonker electron neutral collision frequency // with a={frac}"
        )
        p = self.msise["t_nn"] * self.msise["Tn"] * pconst["boltz"]
        nu = (2.637e6 / np.sqrt(self.iri["etemp"]) + 4.945e5) * p
        cfr = frac * nu
        return cfr

    def atmospheric_ion_neutral_collision_frequency(self):
        """
        This method only provides the atmsphreic ion neutral collision frequency from collision theory

        nn <float> = neutral density

        nu <float> = collision frequency
        """
        nu = 3.8e-11 * self.msise["t_nn"]
        return nu

    def calculate_SN_ei_collision_frequency(self, gamma=0.5572, zi=2):
        """
        This method provides electron ion collision frequency profile, nu_ei

        Ne <float> = Electron density in m^-3
        Ni <float> = Ion density in m^-3
        Te <float> = Electron temperature in K
        Ti <float> = Ion temperature in K
        zi <integer 1/2> = Ion Z number

        nu <float> = collision frequency
        """
        logger.warning(f"Compute the Schank-Nagy ion (n/e) collision frequency")
        key_maps = dict(O2p="o2", Op="o")
        e = pconst["q_e"]
        k = pconst["boltz"]
        me = pconst["m_e"]
        eps0 = pconst["eps0"]
        k_e = 1 / (4 * np.pi * eps0)
        for key in key_maps.keys():
            Ni = self.iri[key_maps[key]]
            ki2 = 4 * np.pi * Ni * 1e6 * e**2 * zi**2 * k_e / (k * self.iri["itemp"])
            ke2 = 4 * np.pi * self.iri["edens"] * e**2 * k_e / (k * self.iri["etemp"])
            ki = np.sqrt(ki2)
            ke = np.sqrt(ke2)
            lam = np.log(
                4 * k * self.iri["etemp"] / (gamma**2 * zi * e**2 * k_e * ke)
            ) - (((ke2 + ki2) / ki2) * np.log(np.sqrt((ke2 + ki2) / ke2)))
            setattr(
                self.collision.nu_sn.ei,
                key,
                (
                    4
                    * np.sqrt(2 * np.pi)
                    * Ni
                    * (zi * e**2 * k_e) ** 2
                    * lam
                    / (3 * np.sqrt(me) * (k * self.iri["etemp"]) ** (1.5))
                ),
            )
        self.collision.nu_sn.ei.total = (
            self.collision.nu_sn.ei.O2p + self.collision.nu_sn.ei.Op
        )
        self.collision.nu_sn.total += self.collision.nu_sn.ei.total
        return

    def calculate_SN_en_collision_frequency(self):
        """
        This method provides electron neutral collision frequency profile, nu_en
        """
        logger.warning(f"Compute the Schank-Nagy electron neutral collision frequency")
        self.collision.nu_sn.en.N2 = (
            1e-6
            * 2.33e-11
            * self.msise["N2"]
            * (1 - (1.12e-4 * self.iri["etemp"]))
            * self.iri["etemp"]
        )
        self.collision.nu_sn.en.O2 = (
            1e-6
            * 1.82e-10
            * self.msise["O2"]
            * (1 + (3.6e-2 * np.sqrt(self.iri["etemp"])))
            * np.sqrt(self.iri["etemp"])
        )
        self.collision.nu_sn.en.O = (
            1e-6
            * 8.9e-11
            * self.msise["O"]
            * (1 + (5.7e-4 * self.iri["etemp"]))
            * np.sqrt(self.iri["etemp"])
        )
        self.collision.nu_sn.en.He = (
            1e-6 * 4.6e-10 * self.msise["He"] * np.sqrt(self.iri["etemp"])
        )
        self.collision.nu_sn.en.H = (
            1e-6
            * 4.5e-9
            * self.msise["H"]
            * (1 - (1.35e-4 * self.iri["etemp"]))
            * np.sqrt(self.iri["etemp"])
        )
        self.collision.nu_sn.en.total = (
            self.collision.nu_sn.en.N2
            + self.collision.nu_sn.en.O2
            + self.collision.nu_sn.en.O
            + self.collision.nu_sn.en.He
            + self.collision.nu_sn.en.H
        )
        self.collision.nu_sn.total += self.collision.nu_sn.en.total
        return

    def atmospheric_collision_frequency(ni, nn, T):
        """
        This method provides atmospheric collision profile based on ion and neural densities and electron temerature

        ni <float> = Ion density m^-3
        nn <float> = Neutral density m^-3
        T <float> = Electron temerature in K

        nu <float> = collision frequency
        """
        na_profile = lambda T, nn: (1.8 * 1e-8 * nn * np.sqrt(T / 300))
        ni_profile = lambda T, ni: (6.1 * 1e-3 * ni * (300 / T) * np.sqrt(300 / T))
        nu = ni_profile(T, ni) + na_profile(T, nn)
        return nu


if __name__ == "__main__":
    from raidpy.ionosphere.iono import Ionosphere2d

    ion = Ionosphere2d(
        dt.datetime(2024, 4, 8),
        np.array([0, 1, 3]),
        np.array([0, 1, 3]),
        np.array([100, 200, 300]),
    )
    cc = ComputeCollision(
        ion.msise_block.msise, ion.iri_block.iri, date=ion.date, _run_=True
    )
