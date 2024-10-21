import numpy as np

# IRF_temperature parameters IPCC AR5, Table 8.SM.9: 1: fast response, 2: slow response, see more info in Olivié and Peters, 2013: doi:10.5194/esd-4-267-2013
# TODO: check if IPCC AR6 has updated these values 
F_1, F_2, TAU_1, TAU_2  = 0.631, 0.429, 8.4, 409.5   # F[K * W-1 m-2], tau [years]


def IRF_temperature(year: int):
    """Temperature impulse response function at year n after a unitary pulse emission of a GHG leading to a radiative forcing of 1 W/m2.
    unit is K."""

    irf = F_1 / TAU_1 * np.exp(-year / TAU_1) + F_2 / TAU_2 * np.exp(-year / TAU_2)

    return irf


