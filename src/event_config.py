"""
Central place for everything that differs *between events* (not between
waveform models). Waveform choice is passed separately as a CLI argument to
src/run_pe.py, so a single config entry here covers all three waveform runs
for that event.
"""
from priors import gw231002_prior, gw231123_prior

EVENT_CONFIG = {
    "GW231002": {
        "time_of_event": 1380292774.7,
        "duration": 8,
        "sampling_frequency": 4096,
        "H1_file": "data/H-H1_GWOSC_O4a_16KHZ_R1-1380290560-4096.hdf5",
        "L1_file": "data/L-L1_GWOSC_O4a_16KHZ_R1-1380290560-4096.hdf5",
        "prior_module": gw231002_prior,
        "phase_marginalization": True,
        "time_marginalization": True,
        "distance_marginalization": False,
        "nlive": 3000,
        "dlogz": 0.1,
        "walks": 100,
        "check_point_delta_t": 60.0,
        "maximum_frequency": 1024,
        "minimum_frequency": 20.0,
        "reference_frequency": 20.0,
    },
    "GW231123": {
        "time_of_event": 1384782888.7,
        "duration": 8,
        "sampling_frequency": 4096,
        "H1_file": "data/H-H1_GWOSC_O4a_16KHZ_R1-1384779776-4096.hdf5",
        "L1_file": "data/L-L1_GWOSC_O4a_16KHZ_R1-1384779776-4096.hdf5",
        "prior_module": gw231123_prior,
        "phase_marginalization": False,
        "time_marginalization": True,
        "distance_marginalization": False,
        "nlive": 3000,
        "dlogz": 0.10,
        "walks": 100,
        "check_point_delta_t": 600.0,
        "maximum_frequency": 1024,
        "minimum_frequency": 20.0,
        "reference_frequency": 20.0,
    },
}

WAVEFORMS = ["IMRPhenomXPHM", "IMRPhenomTPHM", "IMRPhenomPv3"]
