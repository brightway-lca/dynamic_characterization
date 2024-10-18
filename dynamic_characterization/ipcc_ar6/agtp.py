import numpy as np

# Olivié and Peters, 2013: doi:10.5194/esd-4-267-2013

# parameters
F_1 = 0.631 #[K * W-1 m2] Boucher and Reddy (2008) BR08 model, as reported in Olivié and Peters, 2013, Table 4
F_2 = 0.429 #[K * W-1 m2] Boucher and Reddy (2008) BR08 model, as reported in Olivié and Peters, 2013, Table 4

TAU_1 = 8.4 #[years] Boucher and Reddy (2008) BR08 model, as reported in Olivié and Peters, 2013, Table 4
TAU_2 = 409.5 #[years] Boucher and Reddy (2008) BR08 model, as reported in Olivié and Peters, 2013, Table 4

# TODO: check if IPCC has updated these values in AR6

def IRF_temperature(year: int):
    """Temperature impulse response function at year n after a unitary pulse emission of a GHG leading to a radiative forcing of 1 W/m2.
    unit is K."""

    irf = F_1 / TAU_1 * np.exp(-year / TAU_1) + F_2 / TAU_2 * np.exp(-year / TAU_2)

    return irf


def _calculate_irf_temperature_multipliers(time_horizon):
    
    multipliers = []
    integral = 0
    for i in range(time_horizon): # TODO check if exclusion of last year is correct
        integral += IRF_temperature(i)
        multipliers.append(integral)
   
    multipliers = multipliers[::-1]  # reverse order of multipliers to get IRF_T(time_horizon - year)

    return multipliers