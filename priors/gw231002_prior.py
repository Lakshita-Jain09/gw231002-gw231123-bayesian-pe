"""
Prior distribution for GW231002.

This single prior set (chirp_mass in [40, 100] Msun etc.) is shared across
all three waveform-model runs (IMRPhenomXPHM, IMRPhenomTPHM, IMRPhenomPv3)
so that any differences in the recovered posteriors are attributable to the
waveform model and not to the prior choice.
"""
import numpy as np
import bilby
from bilby.core.prior import Uniform, Sine, Cosine
from bilby.gw.prior import UniformSourceFrame


def get_prior(time_of_event):
    prior = bilby.gw.prior.BBHPriorDict()
    prior["chirp_mass"] = Uniform(name="chirp_mass", minimum=40.0, maximum=100.0)
    prior["mass_ratio"] = Uniform(name="mass_ratio", minimum=0.3, maximum=0.99)
    prior["phase"] = Uniform(name="phase", minimum=0, maximum=2 * np.pi)
    prior["geocent_time"] = Uniform(
        name="geocent_time",
        minimum=time_of_event - 0.1,
        maximum=time_of_event + 0.1,
    )
    prior["a_1"] = Uniform(name="a_1", minimum=0.1, maximum=0.99)
    prior["a_2"] = Uniform(name="a_2", minimum=0.1, maximum=0.99)
    prior["tilt_1"] = Sine(name="tilt_1", minimum=0.5, maximum=np.pi)
    prior["tilt_2"] = Sine(name="tilt_2", minimum=0.5, maximum=np.pi)
    prior["phi_12"] = Uniform(name="phi_12", minimum=0.0, maximum=2 * np.pi)
    prior["phi_jl"] = Uniform(name="phi_jl", minimum=0.0, maximum=2 * np.pi)
    prior["dec"] = Cosine(name="dec", minimum=-1.2, maximum=0.5 * np.pi)
    prior["ra"] = Uniform(name="ra", minimum=1.0, maximum=2 * np.pi, boundary="periodic")
    prior["theta_jn"] = Sine(name="theta_jn", minimum=0.0, maximum=np.pi)
    prior["psi"] = Uniform(name="psi", minimum=0.0, maximum=np.pi, boundary="periodic")
    prior["luminosity_distance"] = UniformSourceFrame(
        name="luminosity_distance", minimum=2700.0, maximum=5000.0
    )
    return prior
