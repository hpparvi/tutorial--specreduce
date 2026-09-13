"""Shared data access for the specreduce tutorial book.

Every chapter loads the OSIRIS frames through this module rather than reading FITS
files directly. Keeping it in one place is what makes "one dataset everywhere"
enforceable, and it keeps bias subtraction, trimming and noise bookkeeping --
detector-level details the book does not teach -- out of the chapters' narrative.

The dataset is described in ``notebooks/01-introduction.ipynb``.
"""

from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.io import fits as pf
from astropy.nddata import CCDData, VarianceUncertainty
from astropy.stats import mad_std

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
            if (data_dir / "osiris_bias.fits.bz2").is_file():
                return data_dir
    raise FileNotFoundError(
        "Could not locate the tutorial's data/ directory. It should sit next to "
        "notebooks/ and contain osiris_bias.fits.bz2."
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


def _read_noise(debiased, header):
    """Read noise of a bias-subtracted frame, in electrons per pixel.

    ``BIASSEC`` points at the prescan: real pixels, clocked out and digitized the
    same way as the rest of the frame, but never exposed to light. After bias
    subtraction the only thing left there is noise, so its robust scatter measures
    the read noise directly -- and measures the *right* quantity, because the
    subtracted bias is a single exposure rather than a master, so its own read noise
    is part of the budget. Estimating it here rather than hard-coding a catalog
    value also handles the extra 2x2 block-summing applied when the tutorial frames
    were made: each pixel here is the sum of four independent reads, so its read
    noise is twice the detector's per-read figure.

    A robust estimator is essential -- the prescan carries a bias-level gradient
    along the slow axis and the occasional hot pixel, and a plain standard deviation
    picks both up.
    """
    if "BIASSEC" not in header:
        raise KeyError(
            "Frame has no BIASSEC keyword, so the read noise cannot be measured "
            "from its prescan."
        )
    return mad_std(debiased[parse_section(header["BIASSEC"])]) * header["GAIN"]


def read_file(fname, bias, trim=True):
    """Read one frame, calibrate it, and wrap it as a `~astropy.nddata.CCDData`.

    The frame is bias-subtracted, converted from ADU to electrons using the ``GAIN``
    keyword, given a per-pixel variance combining photon shot noise with read noise,
    and -- unless ``trim`` is False -- cut down to the region the observatory
    recommends in ``TRIMSEC``. Working in electrons is what makes the shot-noise term
    meaningful: Poisson statistics apply to detected charge, not to the arbitrary
    digitization units the detector reports.

    Parameters
    ----------
    fname : str or `~pathlib.Path`
        Path to the FITS file.
    bias : ndarray
        Bias frame in raw ADU, on the same pixel grid.
    trim : bool, optional
        Trim the frame to ``TRIMSEC``. The read noise is measured from the prescan
        *before* trimming, since trimming removes it. Chapters 1 and 2 pass False
        to show the detector frame as recorded; every later chapter works on the
        trimmed frame, and Chapter 2 explains why.

    Returns
    -------
    `~astropy.nddata.CCDData`
        Bias-subtracted signal in electrons, carrying a
        `~astropy.nddata.VarianceUncertainty`. ``specreduce``'s line-finding and
        optimal-extraction routines read the uncertainty to set detection thresholds
        and extraction weights, so it is always attached. The measured read noise is
        also recorded in ``meta['RDNOISE']``, and survives tilt-correction resampling,
        which carries metadata through. A trimmed frame records the section it was
        cut to in ``meta['TRIMMED']``.
    """
    data, header = pf.getdata(fname, header=True)
    debiased = data.astype("d") - bias

    signal = debiased * header["GAIN"]
    read_noise = _read_noise(debiased, header)

    # Shot noise applies to detected charge only, so the Poisson term is the signal
    # itself -- not the raw counts, which are dominated by the bias pedestal. It is
    # clipped at zero because sky-limited pixels scatter negative, and a negative
    # variance would poison every downstream weight. Dark current is left out: these
    # are 2-12 s exposures on a cooled CCD.
    variance = np.maximum(signal, 0.0) + read_noise**2

    meta = {"RDNOISE": read_noise}
    if trim:
        # TRIMSEC is the observatory's recommended usable region. On these frames it
        # drops the prescan, the unilluminated strip below the slit, and two dozen
        # columns at the red end beyond the grism's coverage (Chapter 2 measures all
        # three), so downstream chapters never have to steer around a frame edge.
        section = parse_section(header["TRIMSEC"])
        signal, variance = signal[section], variance[section]
        meta["TRIMMED"] = header["TRIMSEC"]

    return CCDData(
        signal,
        unit=u.electron,
        uncertainty=VarianceUncertainty(variance),
        meta=meta,
    )


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
    bias = pf.getdata(DATA_DIR / "osiris_bias.fits.bz2").astype("d")
    arcs = [
        read_file(DATA_DIR / f"osiris_arc_{lamp}.fits.bz2", bias, trim)
        for lamp in ARC_LAMPS
    ]
    obj = read_file(DATA_DIR / "osiris_tres_3b.fits.bz2", bias, trim)
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
