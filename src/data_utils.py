"""
Utilities for loading GWOSC strain data into gwpy TimeSeries objects.
"""
import h5py
from gwpy.timeseries import TimeSeries


def load_gwosc_strain(filename):
    """
    Load a GWOSC-format strain HDF5 file into a gwpy TimeSeries.

    Parameters
    ----------
    filename : str
        Path to a GWOSC ``*_GWOSC_*.hdf5`` strain file.

    Returns
    -------
    gwpy.timeseries.TimeSeries
    """
    with h5py.File(filename, "r") as f:
        strain = f["strain/Strain"][:]
        gps_start = f["meta/GPSstart"][()]
        duration_file = f["meta/Duration"][()]
        sample_rate = len(strain) / duration_file
    return TimeSeries(strain, t0=gps_start, sample_rate=sample_rate, unit="strain")
