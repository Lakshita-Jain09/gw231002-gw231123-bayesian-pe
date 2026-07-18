#!/usr/bin/env python
"""
Bayesian parameter estimation for GW231002 / GW231123 with a choice of
waveform approximant.

This replaces the original six scripts (gw231002_X/T/p.py,
gw231123X/T/pv3.py), which differed only in `waveform_approximant`,
`outdir`, and (for GW231002 vs GW231123) the prior and a couple of
likelihood/sampler settings. Those differences are now captured in
priors/*.py and src/event_config.py, and selected here via CLI flags.

Example
-------
python src/run_pe.py --event GW231002 --waveform IMRPhenomXPHM \
    --h1-file data/H-H1_..._4096.hdf5 --l1-file data/L-L1_..._4096.hdf5

Strain files can be downloaded from GWOSC (https://gwosc.org) for the
relevant GPS segment; they are not included in this repository.
"""
import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import bilby
import lal

from bilby.gw.conversion import (
    convert_to_lal_binary_black_hole_parameters,
    generate_all_bbh_parameters,
)

from data_utils import load_gwosc_strain
from event_config import EVENT_CONFIG, WAVEFORMS
from remnant_fits import compute_remnant_properties, summarize, STANDARD_SUMMARY_KEYS

lal.swig_redirect_standard_output_error(False)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--event", required=True, choices=list(EVENT_CONFIG.keys()))
    p.add_argument("--waveform", required=True, choices=WAVEFORMS)
    p.add_argument("--h1-file", default=None, help="Override default H1 strain file path")
    p.add_argument("--l1-file", default=None, help="Override default L1 strain file path")
    p.add_argument("--outdir", default=None, help="Bilby output directory")
    p.add_argument("--label", default=None, help="Run label (defaults to event name)")
    p.add_argument("--plots-dir", default="results", help="Where to save summary plots")
    p.add_argument("--sigma", type=int, default=2, choices=[1, 2],
                    help="Credible-interval width for printed summaries (1 or 2 sigma)")
    return p.parse_args()


def main():
    args = parse_args()
    cfg = EVENT_CONFIG[args.event]

    time_of_event = cfg["time_of_event"]
    duration = cfg["duration"]
    analysis_start = time_of_event - 6
    psd_duration = duration * 32
    psd_start_time = analysis_start - psd_duration

    h1_file = args.h1_file or cfg["H1_file"]
    l1_file = args.l1_file or cfg["L1_file"]
    outdir = args.outdir or f"outdir_{args.event}_{args.waveform}"
    label = args.label or args.event
    event_plot_dir = os.path.join(args.plots_dir, args.event)
    os.makedirs(event_plot_dir, exist_ok=True)

    # ---------------------------------------------------------------- data
    H1_full = load_gwosc_strain(h1_file)
    L1_full = load_gwosc_strain(l1_file)

    H1_analysis_data = H1_full.crop(analysis_start, analysis_start + duration)
    L1_analysis_data = L1_full.crop(analysis_start, analysis_start + duration)

    H1_analysis_data.plot()
    plt.savefig(os.path.join(event_plot_dir, "H1_strain.png"))
    plt.close()
    L1_analysis_data.plot()
    plt.savefig(os.path.join(event_plot_dir, "L1_strain.png"))
    plt.close()

    H1 = bilby.gw.detector.get_empty_interferometer("H1")
    L1 = bilby.gw.detector.get_empty_interferometer("L1")
    H1.set_strain_data_from_gwpy_timeseries(H1_analysis_data)
    L1.set_strain_data_from_gwpy_timeseries(L1_analysis_data)

    H1_psd_data = H1_full.crop(psd_start_time, psd_start_time + psd_duration)
    L1_psd_data = L1_full.crop(psd_start_time, psd_start_time + psd_duration)

    psd_alpha_h1 = 2 * H1.strain_data.roll_off / duration
    psd_alpha_l1 = 2 * L1.strain_data.roll_off / duration

    H1_psd = H1_psd_data.psd(fftlength=duration, overlap=0,
                              window=("tukey", psd_alpha_h1), method="median")
    L1_psd = L1_psd_data.psd(fftlength=duration, overlap=0,
                              window=("tukey", psd_alpha_l1), method="median")

    H1.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(
        frequency_array=H1_psd.frequencies.value, psd_array=H1_psd.value)
    L1.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(
        frequency_array=L1_psd.frequencies.value, psd_array=L1_psd.value)
    H1.maximum_frequency = cfg["maximum_frequency"]
    L1.maximum_frequency = cfg["maximum_frequency"]

    for ifo, name in [(H1, "H1"), (L1, "L1")]:
        fig, ax = plt.subplots()
        ax.loglog(ifo.strain_data.frequency_array,
                   np.abs(ifo.strain_data.frequency_domain_strain), label="Strain")
        ax.loglog(ifo.power_spectral_density.frequency_array,
                   ifo.power_spectral_density.asd_array, label="PSD")
        ax.set_xlabel("Frequency [Hz]")
        ax.set_ylabel(r"Strain [strain/$\sqrt{Hz}$]")
        ax.legend()
        fig.savefig(os.path.join(event_plot_dir, f"{name}_psd.png"))
        plt.close(fig)

    # -------------------------------------------------------------- prior
    prior = cfg["prior_module"].get_prior(time_of_event)

    # ---------------------------------------------------------- waveform
    interferometers = [H1, L1]
    waveform_arguments = dict(
        waveform_approximant=args.waveform,
        reference_frequency=cfg["reference_frequency"],
        minimum_frequency=cfg["minimum_frequency"],
        catch_waveform_errors=True,
    )
    waveform_generator = bilby.gw.WaveformGenerator(
        duration=duration,
        sampling_frequency=cfg["sampling_frequency"],
        frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
        waveform_arguments=waveform_arguments,
        parameter_conversion=convert_to_lal_binary_black_hole_parameters,
    )

    likelihood = bilby.gw.likelihood.GravitationalWaveTransient(
        interferometers, waveform_generator, priors=prior,
        time_marginalization=cfg["time_marginalization"],
        phase_marginalization=cfg["phase_marginalization"],
        distance_marginalization=cfg["distance_marginalization"],
    )

    # --------------------------------------------------------- sampling
    result = bilby.run_sampler(
        likelihood, prior, sampler="dynesty", outdir=outdir, label=label,
        conversion_function=bilby.gw.conversion.generate_all_bbh_parameters,
        nlive=cfg["nlive"], dlogz=cfg["dlogz"], queue_size=20,
        bound="multi", sample="rwalk", walks=cfg["walks"], resume=True,
        check_point_delta_t=cfg["check_point_delta_t"],
    )

    # ------------------------------------------------------ postprocess
    posterior = generate_all_bbh_parameters(result.posterior)

    fig = result.plot_corner(
        parameters=["chirp_mass", "mass_ratio", "a_1", "a_2", "tilt_1", "tilt_2",
                    "phi_12", "phi_jl", "psi", "theta_jn", "geocent_time",
                    "phase", "ra", "dec", "luminosity_distance"],
        priors=True,
    )
    fig.savefig(os.path.join(event_plot_dir, f"corner_extrinsic_{args.waveform}.png"))

    fig = result.plot_corner(
        parameters=["chirp_mass", "mass_ratio", "chi_eff", "a_1", "a_2",
                    "mass_1", "mass_2", "total_mass"],
        priors=True,
    )
    fig.savefig(os.path.join(event_plot_dir, f"corner_intrinsic_{args.waveform}.png"))

    posterior = compute_remnant_properties(result.posterior.copy())

    print(f"\n=== {args.event} / {args.waveform} : remnant properties ===")
    for col, name in [
        ("final_mass_approx_source", "final_mass_approx_source [Msun]"),
        ("final_spin_approx", "final_spin_approx"),
        ("f_ringdown_220_source_Hz", "f_ringdown_220_source_Hz"),
        ("f_ringdown_220_detector_Hz", "f_ringdown_220_detector_Hz"),
        ("tau_220_s", "tau_220_s"),
        ("Radiated_energy", "E_rad"),
    ]:
        summarize(posterior[col], name, sigma=args.sigma)

    print(f"\n=== {args.event} / {args.waveform} : evidence ===")
    print("log_evidence      =", result.log_evidence)
    print("log_evidence_err  =", result.log_evidence_err)
    print("log_noise_evidence=", result.log_noise_evidence)
    print("log_bayes_factor  =", result.log_bayes_factor)

    print(f"\n=== {args.event} / {args.waveform} : source parameters ===")
    for key in STANDARD_SUMMARY_KEYS:
        if key in result.posterior.columns:
            summarize(result.posterior[key], key, sigma=args.sigma)

    # bilby's built-in NR-fit remnant columns, when available (distinct from
    # the practical fits in remnant_fits.py used above)
    print(f"\n=== {args.event} / {args.waveform} : bilby built-in remnant columns ===")
    for key in ["final_mass", "peak_luminosity", "recoil_velocity"]:
        present = key in result.posterior.columns
        print(f"{key} in posterior: {present}")
        if present:
            summarize(result.posterior[key], key, sigma=args.sigma)


if __name__ == "__main__":
    main()
