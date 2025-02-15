import numpy as np
from pystrict import strict

from PySDM import Formulae
from PySDM.physics.constants import si
from PySDM.initialisation.spectra import Lognormal
from PySDM.initialisation.sampling import spectral_sampling



@strict
class Settings:
    def __init__(
        self,
        n_sd: int,
        w_updraft: float,
        T0: float,
        N_solution_droplet: float,
        r_solution_droplet: float,
        kappa: float,
        rate: str,
    ):

        self.n_sd = n_sd
        self.w_updraft = w_updraft
        self.r_solution_droplet = r_solution_droplet
        self.N_solution_drople = N_solution_droplet
        self.rate = rate

        self.mass_of_dry_air = 1000 * si.kilogram
        self.initial_pressure = 220 * si.hectopascals
        self.initial_ice_supersaturation = 1.
        self.kappa = kappa
        self.initial_temperature  = T0

        self.formulae = Formulae(
            particle_shape_and_density="MixedPhaseSpheres",
            homogeneous_ice_nucleation_rate=rate,
 #           homogeneous_ice_nucleation_rate="Constant",
            constants={"J_HOM": 1.e15},
        )
        const = self.formulae.constants
        pvs_i = self.formulae.saturation_vapour_pressure.pvs_ice(self.initial_temperature)
        self.initial_water_vapour_mixing_ratio = const.eps / (
            self.initial_pressure / self.initial_ice_supersaturation / pvs_i - 1
        )
        dry_air_density =  (self.formulae.trivia.p_d(self.initial_pressure, self.initial_water_vapour_mixing_ratio )
                            / self.initial_temperature
                            / const.Rd )

        spectrum = Lognormal(norm_factor=N_solution_droplet / dry_air_density,  m_mode=r_solution_droplet, s_geom=1.5)
        self.r_dry, self.specific_concentration = spectral_sampling.Logarithmic(spectrum).sample(n_sd)





        self.t_duration = 5400 # total duration of simulation
        self.dt         = 1. 
        self.n_output = 10 # number of output steps


n_sds = ( 1000, )

w_updrafts = (
    10 * si.centimetre / si.second,
)
        
T_starts = ( 220 * si.kelvin, )

N_solution_droplets = ( 2500 / si.centimetre**3, )

r_solution_droplets = ( 0.0555 * si.micrometre, )

kappas = ( 0.64, )

hom_rates = ( "Constant","Koop2000", "Koop_Correction")

setups = []
for n_sd in n_sds:
    for w_updraft in w_updrafts:
        for T0 in T_starts:
            for N_solution_droplet in N_solution_droplets:
                for r_solution_droplet in r_solution_droplets:
                    for kappa in kappas:
                        for rate in hom_rates:
                            setups.append(
                                Settings(
                                    n_sd=n_sd,
                                    w_updraft=w_updraft,
                                    T0=T0,
                                    N_solution_droplet=N_solution_droplet,
                                    r_solution_droplet=r_solution_droplet,
                                    kappa=kappa,
                                    rate=rate,
                                )
                            )
