"""
Phenomenological remnant-property fits (final mass, final spin, 220 ringdown
frequency/damping time, radiated energy) applied to a bilby posterior, plus a
small summary-statistics helper used throughout the analysis notebooks.

These are practical, NR-calibrated fits (not publication-grade fits with full
error propagation) used to give a quick physical picture of the remnant BH.
"""
import numpy as np

MTSUN_SI = 4.92549095e-6  # solar mass in seconds


def compute_remnant_properties(posterior):
    """
    Given a bilby posterior DataFrame (already run through
    generate_all_bbh_parameters), compute approximate remnant-BH properties
    and append them as new columns. Returns the same DataFrame.
    """
    m1 = posterior["mass_1_source"].values
    m2 = posterior["mass_2_source"].values
    chi1z = posterior["spin_1z"].values
    chi2z = posterior["spin_2z"].values

    eta = (m1 * m2) / (m1 + m2) ** 2
    Mtot = m1 + m2
    chi_eff = (m1 * chi1z + m2 * chi2z) / (m1 + m2)

    # Radiated energy fraction (bounded to a physically sane range)
    erad_frac = 0.048 * (4 * eta) ** 2 * (1 + 0.3 * chi_eff)
    erad_frac = np.clip(erad_frac, 0.01, 0.12)

    final_mass = Mtot * (1 - erad_frac)
    E_rad = erad_frac * (m1 + m2)

    # Final spin (simple practical fit)
    final_spin = np.sqrt(12) * eta - 3.871 * eta ** 2 + 0.5 * chi_eff
    final_spin = np.clip(final_spin, 0.0, 0.99)

    # 220 ringdown mode (Berti-style practical fit)
    fhat_220 = 1.5251 - 1.1568 * (1 - final_spin) ** 0.1292
    Q_220 = 0.7000 + 1.4187 * (1 - final_spin) ** (-0.4990)
    M_sec = final_mass * MTSUN_SI
    f_ringdown_220_Hz = fhat_220 / (2 * np.pi * M_sec)
    tau_220_s = Q_220 / (np.pi * f_ringdown_220_Hz)

    z = posterior["redshift"].values
    f_ringdown_220_det_Hz = f_ringdown_220_Hz / (1 + z)

    posterior["final_mass_approx_source"] = final_mass
    posterior["final_spin_approx"] = final_spin
    posterior["f_ringdown_220_source_Hz"] = f_ringdown_220_Hz
    posterior["f_ringdown_220_detector_Hz"] = f_ringdown_220_det_Hz
    posterior["tau_220_s"] = tau_220_s
    posterior["Radiated_energy"] = E_rad
    return posterior


def summarize(arr, name, sigma=2):
    """
    Print median and credible interval for a 1D array.

    sigma=1 -> 16/50/84 percentiles (68.3% CI)
    sigma=2 -> 2.5/50/97.5 percentiles (95.5% CI)
    """
    arr = np.asarray(arr)
    arr = arr[np.isfinite(arr)]
    if len(arr) == 0:
        print(f"{name}: all NaN")
        return None
    lo, mid, hi = (16, 50, 84) if sigma == 1 else (2.5, 50, 97.5)
    q_lo, q_mid, q_hi = np.percentile(arr, [lo, mid, hi])
    print(f"{name}: {q_mid:.3f} (+{q_hi - q_mid:.3f}, -{q_mid - q_lo:.3f})")
    return q_lo, q_mid, q_hi


def summarize_param(posterior, key, sigma=2):
    if key not in posterior.columns:
        return None
    return summarize(posterior[key].values, key, sigma=sigma)


# Parameters routinely reported in the thesis tables (5.1 / 6.1 / A.1 / A.2)
STANDARD_SUMMARY_KEYS = [
    "chirp_mass", "mass_ratio", "chi_eff", "chi_p",
    "mass_1_source", "mass_2_source", "luminosity_distance", "redshift",
    "theta_jn", "total_mass", "mass_1", "mass_2", "a_1", "a_2",
    "tilt_1", "tilt_2", "phi_12", "phi_jl", "dec", "ra", "psi",
    "phase", "geocent_time", "symmetric_mass_ratio", "spin_1z", "spin_2z",
    "chi_1_in_plane", "chi_2_in_plane", "cos_tilt_1", "cos_tilt_2",
    "comoving_distance", "chirp_mass_source", "total_mass_source",
]
