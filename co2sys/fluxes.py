#!/usr/bin/env python

"""
This script contains all the functions required to Calculate CO2
fluxes.

CO2_flux is the main script (see documentation for that).
Subfunctions for that script are:
    K0
    schmidt_number
    sat_pH2O

There are also several piston velocity parameterisations under k:
    Li86
    Wa92
    Wa99
    Ni00
    Ho06
    Wa09
"""

import numpy as np

# from numpy import array, max, log, exp, zeros_like, copy
# import numexpr as ne


def co2_flux(S, t, U10, dpCO2, k=None):
    """
    Calculates CO2 flux

      F = k . K0 . dpCO2

    INPUT:  Salinity    [psu]
            Temperature [C]
            Wind Speed  [m/s]
            dpCO2       [uatm]
          **k = gas transfer velocity function
                defaults to Wa92
    OUTPUT: CO2 Flux    [mmol/m2/day]

    DIMENSION ANALYSIS
      k   = cm/hr * 1e-2 * 24 = m/day
      K0  = mmol.m-3.atm-1
      pCO2= uatm  * 1e-6 = atm

      F = mmol . m-3 . m . atm-1 . atm . day-1
      F = mmol . m-2 . day-1
    """

    k = kw.Wa92 if k is None else k

    return k(U10, t) * K0(S, t) * dpCO2 * 1e-6 * 1e-2 * 24.


def K0(S, T):
    """
    Returns the Solubility parameter
    INPUT:  T  = temperature [degC] array
            S  = salinity [psu / absolute S] array
    OUTPUT: K0 array
            corrected for moisture
            mmol.m-3.atm-1
    """

    log = np.log
    exp = np.exp

    T += 273.15

    # Solubility Function of CO2
    # Weiss and Price (1980)
    A1, A2, A3, A4 = -160.7333, 215.4152, 89.8920, -1.47759
    B1, B2, B3 = 0.029941, -0.027455, 0.0053407

    # Weiss and Price (1974)
    A1, A2, A3, A4 = -60.2409, 93.4517, 23.3585, 0
    B1, B2, B3 = 0.023517, -0.023656, 0.0047036

    F = A1
    F += A2 * (100. / T)
    F += A3 * log(T / 100.)
    F += A4 * (T / 100.)**2
    F += S * (B1 + B2 * (T / 100.) + B3 * (T / 100.)**2)

    F = exp(F)

    P, pH2O = sat_pH2O(S, T)

    # Solubility parameter

    return (F / (P - pH2O)) * 1e6


def schmidt_number(t):
    """
    returns the Schmidt number for a given temperature [degC]
    """
    t = np.array(t, ndmin=1)

    if max(t) > 200:
        print "Temperature is in Kelvin\nMust be in degrees Celcius"
    else:
        A = 2073.1
        B = 125.62
        C = 3.6276
        D = 0.043219

        return A - B * t + C * (t**2) - D * (t**3)


def sat_pH2O(S, T, **atm):
    """
    Saturation water vapour pressure - Weiss and Price (1980)
    INPUT:    Salinity  [psu]
              Temperature [C]
              ** 'P'= total pressure (default is 1.)
    """

    log = np.log
    exp = np.exp

    # Atm pressure
    P = atm.pop('P', 1.)

    pH2O = exp(24.4543 -
               67.4509 * (100. / T) -
               4.8489 * log(T / 100.) -
               0.000544 * S)

    return P, pH2O


class _gas_transfer_velocities:
    def Li86(self, u, t):
        """
        Calculates k for Liss and Merlivat 1986
        INPUT:  u: wind speed [m.s^-1] array
                t: temp [degC]
        OUTPUT: k600 cm . hr^-1
        """

        Sc = schmidt_number(t)
        k = np.zeros_like(t)
        u = np.array(u)

        i1 = u <= 3.6
        i2 = (u > 3.6) & (u < 13.)
        i3 = u >= 13.

        k[i1] = (0.17 * u[i1]) * (Sc[i1] / 600)**(-2. / 3.)
        k[i2] = ((u[i2] - 3.4) * 2.8) * (600 / Sc[i2])**0.5
        k[i3] = ((u[i3] - 8.4) * 5.9) * (600 / Sc[i3])**0.5
        return k

    def Wa92(self, u, t):
        """
        Calculates k for Wanninkhof 1992
        INPUT:  u: wind speed [m.s^-1] array
                t: temp [degC]
        OUTPUT: k660 cm . hr^-1
        """
        Sc = schmidt_number(t)
        return (0.31 * (u**2)) * (660 / Sc)**0.5

    def Wa99(self, u, t):
        """
        Calculates k for Wanninkhof et al 1999
        INPUT:  u: wind speed [m.s^-1] array
                t: temp [degC]
        OUTPUT: k660 cm . hr^-1
        """
        Sc = schmidt_number(t)
        return (0.0283 * u**3) * (600 / Sc)**0.5

    def Ni00(self, u, t):
        """
        Calculates k for Nightingale et al 2000
        INPUT:  u: wind speed [m.s^-1] array
                t: temp [degC]
        OUTPUT: k600 cm . hr^-1
        """
        Sc = schmidt_number(t)
        return (0.333 * u + 0.222 * u**2)\
            * (600 / Sc)**0.5

    def Ho06(self, u, t):
        """
        Calculates k for Ho et al 2006
        INPUT:  u: wind speed [m.s^-1] array
                t: temp [degC]
        OUTPUT: k600 cm . hr^-1
        """
        Sc = schmidt_number(t)
        return (0.266 * u**2) * (600 / Sc)**0.5

    def Sw07(self, u, t):
        """
        Calculates k for Wanninkhof 1992
        INPUT:  u: wind speed [m.s^-1] array
                t: temp [degC]
        OUTPUT: k660 cm . hr^-1
        """
        Sc = schmidt_number(t)
        return (0.27 * (u**2)) * (660 / Sc)**0.5

    def Wa09(self, u, t):
        """
        Calculates k for Wanninkhof et al 2009
        INPUT:  u: wind speed [m.s^-1] array
                t: temp [degC]
        OUTPUT: k660 cm . hr^-1
        """
        Sc = schmidt_number(t)
        return (3 + (0.1 * u) + (0.064 * (u**2)) + (0.011 * (u**3)))\
            * (660 / Sc)**0.5


gas_transfer_velocity = kw = _gas_transfer_velocities()


def molCmyr_2_TgC(molCyr, aream2):
        mol = molCyr * aream2
        gC = mol * 12.011
        TgC = gC * 1e-12
        return TgC


if __name__ == "__main__":
    CO2_flux(34.5, 20, 10, -50, k=kw.Ni00)