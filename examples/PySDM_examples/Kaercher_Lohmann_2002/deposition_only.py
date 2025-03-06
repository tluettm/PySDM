import PySDM.products as PySDM_products
from PySDM.backends import CPU
from PySDM.builder import Builder
from PySDM.dynamics import AmbientThermodynamics, Condensation, VapourDepositionOnIce
from PySDM.environments import Parcel, Box
from PySDM.initialisation import discretise_multiplicities

from PySDM import Formulae
from PySDM.physics.constants import si
from PySDM.initialisation.spectra import Lognormal
from PySDM.initialisation.sampling import spectral_sampling


class Simulation:
    def __init__(self,backend=CPU):

        env_type_parcel = True

        self.dt = 1.
        self.n_sd = 1

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
        self.dry_air_density = (self.formulae.trivia.p_d(self.initial_pressure, self.initial_water_vapour_mixing_ratio)
                           / self.initial_temperature
                           / const.Rd)
        if env_type_parcel:
            print( "Env type is parcel" )
            env = Parcel(
                mixed_phase=True,
                dt=self.dt,
                mass_of_dry_air=self.mass_of_dry_air,
                p0=self.initial_pressure,
                initial_water_vapour_mixing_ratio=self.initial_water_vapour_mixing_ratio,
                T0=self.initial_temperature,
                w=1e-20 * si.centimetre / si.second,
            )
        else:
            print("Env type is box")
            env = Box(
                dt=self.dt,
                dv=1 * si.m**3
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
        spectrum = Lognormal(norm_factor=N_dv_solution_droplet / self.dry_air_density, m_mode=r_mean_solution_droplet,
                             s_geom=sigma_solution_droplet)
        self.radius, self.specific_concentration = spectral_sampling.Logarithmic(spectrum).sample(self.n_sd)
        self.radius = -self.radius
        self.multiplicities = discretise_multiplicities(self.specific_concentration * self.mass_of_dry_air)
        self.signed_mass = self.formulae.particle_shape_and_density.radius_to_mass( self.radius )


        attributes = {
            "multiplicity": self.multiplicities,
            "signed water mass": self.signed_mass,
        }

        products = [
            PySDM_products.Time(name="t"),
            # PySDM_products.AmbientRelativeHumidity(name="RH_ice", unit="%"),
            # PySDM_products.AmbientTemperature(name="T"),
            # PySDM_products.AmbientWaterVapourMixingRatio(
            #     name="vapour", var="water_vapour_mixing_ratio", unit='kg/kg'),
            ]

        self.particulator = builder.build(attributes, products)

        if not env_type_parcel:
            self.particulator.environment["T"] = self.initial_temperature
            self.particulator.environment["p"] = self.initial_pressure
            self.particulator.environment["water_vapour_mixing_ratio"] =  self.initial_water_vapour_mixing_ratio
            self.particulator.environment["rhod"] = self.dry_air_density
            self.particulator.environment["RH_ice"] = self.initial_ice_supersaturation

    def run(self,steps):
            multiplicity = self.particulator.attributes["multiplicity"].data[0].copy()
            mass = -self.particulator.attributes["signed water mass"].data[0].copy()
            mixing_ratio = multiplicity * mass / self.mass_of_dry_air

            vapour = self.particulator.environment["water_vapour_mixing_ratio"].data[0].copy()

            int_mass = vapour + mixing_ratio

            for i in range(steps):
                self.particulator.run(1000)

                time = self.particulator.products["t"].get()
                temp = self.particulator.environment["T"].data[0]
                RHi = self.particulator.environment["RH_ice"].data[0]

                print(f"{time=},{temp=},{RHi=}")

                multiplicity = self.particulator.attributes["multiplicity"].data[0]
                mass = -self.particulator.attributes["signed water mass"].data[0]
                mixing_ratio = multiplicity * mass / self.mass_of_dry_air

                vapour = self.particulator.environment["water_vapour_mixing_ratio"][0]

                mass_conv = vapour + mixing_ratio - int_mass

                print( f"{vapour=},{mixing_ratio=},{mass_conv=}" )

model = Simulation()

model.run(100)

