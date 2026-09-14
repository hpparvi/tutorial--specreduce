"""Shared data access for the specreduce tutorial book.

Every chapter loads the OSIRIS frames through this module rather than reading FITS
files directly. Keeping it in one place is what makes "one dataset everywhere"
enforceable, and it keeps trimming -- a detector-level detail the book does not
teach -- out of the chapters' narrative. The frames themselves arrive already
reduced: ``data/full_data/downsample_and_reduce_data.py`` applies the overscan,
bias, gain and flat corrections and attaches the per-pixel uncertainty, so a
chapter sees calibrated charge in electrons from the first cell.

The dataset is described in ``notebooks/01-introduction.ipynb``.
"""

from pathlib import Path

import numpy as np
from astropy.nddata import CCDData

__all__ = [
    "ANCHOR_PIXELS",
    "ANCHOR_WAVELENGTHS",
    "ARC_LAMPS",
    "DATA_DIR",
    "GTC_HGAR",
    "TILT_SOLUTION_FILE",
    "WAVELENGTH_SOLUTION_FILE",
    "fit_tilt_solution",
    "fit_wavelength_solution",
    "parse_section",
    "read_data",
    "read_file",
]

#: Arc lamps in the dataset, in the order ``read_data`` returns them.
ARC_LAMPS = ("HgAr", "Ne", "Xe")


def _find_data_dir():
    """Locate the ``data/`` directory from either this module or the working directory.

    Notebooks execute with their own directory as the working directory, but a book
    build may copy sources elsewhere, so we search upwards from both.
    """
    for start in (Path(__file__).resolve().parent, Path.cwd().resolve()):
        for candidate in (start, *start.parents):
            data_dir = candidate / "data"
            if (data_dir / "osiris_tres_3b.fits").is_file():
                return data_dir
    raise FileNotFoundError(
        "Could not locate the tutorial's data/ directory. It should sit next to "
        "notebooks/ and contain osiris_tres_3b.fits."
    )


DATA_DIR = _find_data_dir()


def parse_section(value):
    """Turn an IRAF-style section string into a pair of numpy slices.

    FITS section keywords such as ``DATASEC`` and ``BIASSEC`` are written
    ``'[x1:x2,y1:y2]'``, in 1-based inclusive FITS convention with the fast axis
    first. Numpy indexes the other way round and is 0-based and half-open, so both
    the order and the bounds have to be converted.
    """
    x, y = value.strip().strip("[]").split(",")
    x1, x2 = (int(v) for v in x.split(":"))
    y1, y2 = (int(v) for v in y.split(":"))
    return np.s_[y1 - 1 : y2, x1 - 1 : x2]


def read_file(fname, trim=True):
    """Read one reduced frame as a `~astropy.nddata.CCDData`, trimmed to ``TRIMSEC``.

    The files are written by ``data/full_data/downsample_and_reduce_data.py`` and are
    already bias-subtracted, flat-fielded and converted to electrons, with a
    `~astropy.nddata.StdDevUncertainty` that combines photon shot noise with the read
    noise. ``CCDData.read`` restores all of that; the only thing left to do here is
    the trim.

    Parameters
    ----------
    fname : str or `~pathlib.Path`
        Path to the FITS file.
    trim : bool, optional
        Trim the frame to ``TRIMSEC``, the region the observatory recommends. Chapters
        1 and 2 pass False to show the detector frame as recorded; every later chapter
        works on the trimmed frame, and Chapter 2 explains why.

    Returns
    -------
    `~astropy.nddata.CCDData`
        Reduced signal in electrons, carrying a `~astropy.nddata.StdDevUncertainty`.
        ``specreduce``'s line-finding and optimal-extraction routines read the
        uncertainty to set detection thresholds and extraction weights, so it is always
        attached. ``meta`` is the FITS header: the observation keywords, the sections
        remapped to this pixel grid, the reduction record, and ``RDNOISE``, the read
        noise per pixel in electrons. A trimmed frame records the section it was cut to
        in ``meta['TRIMMED']``.
    """
    frame = CCDData.read(fname)

    if trim:
        # TRIMSEC is the observatory's recommended usable region. On these frames it
        # drops the prescan, the unilluminated strip below the slit, and two dozen
        # columns at the red end beyond the grism's coverage (Chapter 2 measures all
        # three), so downstream chapters never have to steer around a frame edge.
        section = frame.meta["TRIMSEC"]
        frame = frame[parse_section(section)]
        frame.meta["TRIMMED"] = section

    return frame


def read_data(trim=True):
    """Read the whole tutorial dataset.

    Parameters
    ----------
    trim : bool, optional
        Trim every frame to ``TRIMSEC``; see `read_file`.

    Returns
    -------
    arcs : list of `~astropy.nddata.CCDData`
        The HgAr, Ne, and Xe arc frames, in `ARC_LAMPS` order.
    lamps : tuple of str
        The arc lamp names, i.e. `ARC_LAMPS`.
    obj : `~astropy.nddata.CCDData`
        The TrES-3 science frame.
    """
    arcs = [read_file(DATA_DIR / f"osiris_arc_{lamp}.fits", trim) for lamp in ARC_LAMPS]
    obj = read_file(DATA_DIR / "osiris_tres_3b.fits", trim)
    return arcs, ARC_LAMPS, obj


#: Where Chapter 6 writes the fitted tilt solution, and later chapters read it from.
TILT_SOLUTION_FILE = DATA_DIR / "tilt_solution.asdf"


def fit_tilt_solution(arcs):
    """Fit the book's tilt solution, using the settings Chapter 6 arrives at.

    Chapter 6 derives these numbers and explains what each one does; this function exists
    so that the chapters after it can regenerate the solution when
    `TILT_SOLUTION_FILE` is missing -- if you open Chapter 8 without having run Chapter 6
    first, for instance. In normal use the file is read, not refitted.

    Parameters
    ----------
    arcs : list of `~astropy.nddata.CCDData`
        The arc frames, as returned by `read_data`.

    Returns
    -------
    `~specreduce.tilt_solution.TiltSolution`
        The fitted solution, ready to `~specreduce.tilt_solution.TiltSolution.resample`
        frames or to serialize with ``to_asdf``.
    """
    from specreduce.tilt_correction import TiltCorrection
    from specreduce.tracing import FlatTrace

    tilt = TiltCorrection(arcs, FlatTrace(arcs[0], 135), n_cdisp_samples=7)
    tilt.find_arc_lines(fwhm=2.5, noise_factor=25)
    tilt.fit(degree=4)
    tilt.refine_fit(degree=4, match_distance_bound=1.0)
    return tilt.solution


#: Where Chapter 10 writes the fitted 1D wavelength solution.
WAVELENGTH_SOLUTION_FILE = DATA_DIR / "wavelength_solution.asdf"

#: Air wavelengths of the HgAr lines GTC lists for OSIRIS and the R1000R grism, from the
#: observatory's arc line maps (https://www.gtc.iac.es/instruments/osiris/osiris.php).
GTC_HGAR = np.array(
    [
        5460.735,
        5769.598,
        5790.663,
        6965.431,
        7272.936,
        7635.106,
        7724.207,
        7948.176,
        8115.311,
        8264.522,
        9122.967,
    ]
)

#: Anchor lines identified by hand in Chapter 9: approximate pixel positions read off the
#: detected-line figure, and the GTC reference wavelengths they were matched to. `fit_lines`
#: snaps the pixels to the measured centroids with `match_obs=True`.
ANCHOR_PIXELS = np.array([87, 533, 550, 622, 650, 803], dtype=float)
ANCHOR_WAVELENGTHS = np.array(
    [5460.735, 7635.106, 7724.207, 8115.311, 8264.522, 9122.967]
)


def fit_wavelength_solution(arcs, tilt):
    """Fit the book's 1D wavelength solution, using the settings Chapters 9 and 10 arrive at.

    Chapter 9 identifies the anchor lines by hand and Chapter 10 chooses the refinement
    schedule, and both explain every choice; this function exists so that later chapters can
    regenerate the solution when `WAVELENGTH_SOLUTION_FILE` is missing. In normal use the file
    is read, not refitted.

    Parameters
    ----------
    arcs : list of `~astropy.nddata.CCDData`
        The arc frames, as returned by `read_data`.
    tilt : `~specreduce.tilt_solution.TiltSolution`
        The tilt solution the arcs are rectified with, so that the wavelength solution is
        defined in the rectified space.

    Returns
    -------
    `~specreduce.wavesol1d.WavelengthSolution1D`
        The three-lamp solution.
    """
    import warnings

    import astropy.units as u
    from specreduce.extract import BoxcarExtract
    from specreduce.tracing import FlatTrace
    from specreduce.wavecal1d import WavelengthCalibration1D

    reference = FlatTrace(arcs[0], 135)
    spectra = [
        BoxcarExtract(tilt.resample(frame), reference, width=15).spectrum
        for frame in arcs
    ]

    calibration = WavelengthCalibration1D(
        arc_spectra=spectra,
        line_lists=[GTC_HGAR, ["NeI"], ["XeI"]],
        line_list_bounds=(5000, 10500),
        ref_pixel=512,
        unit=u.angstrom,
        wave_air=True,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # unconverged refits on non-line detections
        calibration.find_lines(fwhm=2.5, noise_factor=15)
    calibration.fit_lines(
        pixels=ANCHOR_PIXELS,
        wavelengths=ANCHOR_WAVELENGTHS,
        degree=3,
        match_obs=True,
        match_cat=True,
        refine_fit=True,
    )
    calibration.refine_fit(max_match_distance=2.0)
    return calibration.solution
