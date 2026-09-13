# Reducing Spectroscopic Data with specreduce

This book walks through reducing long-slit spectroscopic data with
[`specreduce`](https://github.com/astropy/specreduce), the Astropy-coordinated package for
spectral tracing, extraction, and calibration. It is written for graduate students and working
astronomers who have a two-dimensional spectral image in hand and want a calibrated
one-dimensional spectrum out of it, and who want to understand what happened in between.

`specreduce` is designed for **both interactive use and automated pipelines**: the same
building blocks work step by step in a notebook (inspecting a trace, tweaking a background
window, checking wavelength residuals by eye) and stitched together into a scripted,
non-interactive reduction. This book teaches both styles, showing each step interactively
before composing the full end-to-end pipeline.

The chapters are ordered as a reduction pipeline, in the order a real reduction runs: the
source is located and the frame is rectified before anything is measured from it, so each
chapter picks up where the previous one finished and ends at a wavelength-calibrated 1D
spectrum.

The book grew out of a set of `specreduce` tutorial notebooks the author wrote for various
instruments and datasets. Those notebooks are now deprecated, but a set of instrument-specific 
case studies will be included later in the future.

```{tableofcontents}
```

:::{note}
**Scope.** This book covers only the steps `specreduce` itself performs: tracing, tilt
correction, background subtraction, extraction, and wavelength calibration. Detector-level
reduction (dark and flat correction, cosmic-ray rejection) is the domain of
[`ccdproc`](https://ccdproc.readthedocs.io) and
[`astroscrappy`](https://astroscrappy.readthedocs.io) and is covered by the
[CCD Data Reduction Guide](https://www.astropy.org/ccd-reduction-and-photometry-guide/).
Flux calibration has a chapter, but it is a placeholder: `specreduce` gains the tools for it
in v1.11.

**Data.** Every chapter works on the same real dataset: a long-slit observation of the TrES-3
system taken with the OSIRIS spectrograph at the Gran Telescopio Canarias, together with its
HgAr, Ne, and Xe arc lamp frames. Chapter 1 introduces it, and explains why it was chosen.
:::

:::{note}
**Tilt correction is instrument-specific.** These OSIRIS frames have a strong tilt, so the
book rectifies them early and every later chapter works on the corrected frame. If your
spectrograph does not need it, skip that chapter: every call in the chapters that follow is
unchanged, and only the numbers differ.
:::

## Credits

### Authors

[Hannu Parviainen](https://orcid.org/0000-0001-5519-1391) (Universidad de La Laguna and
Instituto de Astrofísica de Canarias) wrote the book and chose the dataset. He first reduced
these frames in 2014 with a pipeline of his own, and not well. The result has annoyed him ever
since, and this book is where he finally does it properly.

Contributors who write a chapter, or review one in detail, are added to the author list.

### Funding

% TODO(author): state the funding source(s) for this work, if any.

### Acknowledgements

The book is built on `specreduce`, and would not exist without the people who wrote and
maintain it. Several of the operations taught here, and the API the book is written against,
were developed alongside the book.

The data are based on observations made with the Gran Telescopio Canarias (GTC), installed at
the Spanish Observatorio del Roque de los Muchachos of the Instituto de Astrofísica de
Canarias, on the island of La Palma, under program GTC9-14A. The frames were taken with
OSIRIS for the study published as
[Parviainen et al. (2016), A&A 585, A114](https://ui.adsabs.harvard.edu/abs/2016A%26A...585A.114P/abstract).

The book is published through [Learn Astropy](https://learn.astropy.org), whose build
infrastructure executes every notebook and renders the result. Its structure follows the
[CCD Data Reduction Guide](https://www.astropy.org/ccd-reduction-and-photometry-guide/), which
covers the detector-level steps this book leaves out.

## Software setup

The notebooks need Python 3.11 or later and `specreduce` 1.10 or later, together with the
packages it builds on. Everything is on
[conda-forge](https://conda-forge.org), so a fresh environment is one command. Install Python
through [Miniforge](https://github.com/conda-forge/miniforge) if you do not have conda already,
then:

```sh
conda create -n specreduce -c conda-forge specreduce specutils ccdproc astropy matplotlib jupyterlab
conda activate specreduce
```

If you prefer pip, or already have an environment you want to reuse, the repository's
`requirements.txt` lists the same packages with lower bounds only, so the book always runs
against current releases:

```sh
pip install -r requirements.txt
```

Clone the book's [GitHub repository](https://github.com/hpparvi/tutorial--specreduce), then
launch Jupyter from the `notebooks` directory:

```sh
git clone https://github.com/hpparvi/tutorial--specreduce.git
cd tutorial--specreduce/notebooks
jupyter lab
```

The notebooks import a helper module, `src/common.py`, which loads the frames, subtracts the
bias, converts counts to electrons, and attaches a variance. Every chapter reaches it with
`sys.path.append("../src")`, so if you run a chapter's code outside the notebooks, copy that
file alongside your script.

### Data files

The book's five frames, a science exposure, three arc lamps, and a bias, total about 2.5 MB
as bzip2-compressed FITS and live in the `data` directory at the top of the repository.
`common.read_data` locates that directory itself, so nothing has to be configured. Three
derived products (the tilt solution and the two wavelength solutions, all ASDF) are written into
the same directory by the chapters that fit them, and read back by the chapters that follow;
each of those chapters refits if the file is missing, so any chapter can be run on its own.

## Contributing

Corrections and suggestions are welcome as issues or pull requests on the
[repository](https://github.com/hpparvi/tutorial--specreduce). Notebooks are committed without
outputs; the continuous-integration build executes them and renders the book.
