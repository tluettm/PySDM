import PySDM.products as PySDM_products
from PySDM.backends import CPU
from PySDM.builder import Builder
from PySDM.dynamics import AmbientThermodynamics, Condensation, Freezing, VapourDepositionOnIce
from PySDM.environments import Parcel
from PySDM.physics import constants as const
from PySDM.initialisation import discretise_multiplicities, equilibrate_wet_radii

from PySDM import Formulae
from PySDM.physics.constants import si
from PySDM.initialisation.spectra import Lognormal
from PySDM.initialisation.sampling import spectral_sampling




class Simulation:
    def __init__(self,backend=CPU):

        self.dt = 1.
        self.n_sd = 10

        self.initial_pressure=220 * si.hPa
        self.initial_temperature=220 * si.kelvin
        self.mass_of_dry_air = 1000 * si.kilogram
        self.initial_ice_supersaturation = 1.6

        self.formulae = Formulae(
            particle_shape_and_density="MixedPhaseSpheres",
        )
        const = self.formulae.constants
        pvs_i = self.formulae.saturation_vapour_pressure.pvs_ice(self.initial_temperature)
        self.initial_water_vapour_mixing_ratio = const.eps / (
                self.initial_pressure / self.initial_ice_supersaturation / pvs_i - 1
        )
        dry_air_density = (self.formulae.trivia.p_d(self.initial_pressure, self.initial_water_vapour_mixing_ratio)
                           / self.initial_temperature
                           / const.Rd)

        env = Parcel(
            mixed_phase=True,
            dt=self.dt,
            mass_of_dry_air=self.mass_of_dry_air,
            p0=self.initial_pressure,
            initial_water_vapour_mixing_ratio=self.initial_water_vapour_mixing_ratio,
            T0=self.initial_temperature,
            w=1e-20 * si.centimetre / si.second,
        )

        builder = Builder(
            backend=backend(
                formulae=self.formulae,
                **(
                    {"override_jit_flags": {"parallel": False}}
                    if backend == CPU
                    else {}
                )
            ),
            n_sd=self.n_sd,
            environment=env,
        )

        builder.add_dynamic(AmbientThermodynamics())
        #builder.add_dynamic(Condensation())
        builder.add_dynamic(VapourDepositionOnIce())

        N_dv_solution_droplet =  2500 / si.centimetre ** 3
        r_mean_solution_droplet = 10 * si.micrometre
        sigma_solution_droplet =  1.6
        spectrum = Lognormal(norm_factor=N_dv_solution_droplet / dry_air_density, m_mode=r_mean_solution_droplet,
                             s_geom=sigma_solution_droplet)
        self.radius, self.specific_concentration = spectral_sampling.Logarithmic(spectrum).sample(self.n_sd)
        self.radius = -self.radius
        self.multiplicities = discretise_multiplicities(self.specific_concentration * env.mass_of_dry_air)


        self.signed_mass = self.formulae.particle_shape_and_density.radius_to_mass( self.radius )
        print( self.multiplicities, self.radius, self.signed_mass)

        attributes = {
            "multiplicity": self.multiplicities,
            "signed water mass": self.signed_mass,
        }

        products = [
            PySDM_products.Time(name="t"),
            PySDM_products.AmbientRelativeHumidity(name="RH_ice", unit="%"),
            PySDM_products.AmbientTemperature(name="T"),
            PySDM_products.AmbientWaterVapourMixingRatio(
                name="vapour", var="water_vapour_mixing_ratio"),
            ]

        self.particulator = builder.build(attributes, products)

    def run(self,steps):

            for i in range(steps):
                self.particulator.run(5)
                print(self.particulator.products["t"].get())
                print(self.particulator.products["T"].get())
                print(self.particulator.products["RH_ice"].get())
model = Simulation()

model.run(1000)

