# Reducing Spectroscopic Data with specreduce

This tutorial book walks through reducing long-slit spectroscopic data with
[`specreduce`](https://github.com/astropy/specreduce), the Astropy-coordinated package
for spectral tracing, extraction, and calibration.

`specreduce` is designed for **both interactive use and automated pipelines**: the same
building blocks work step-by-step in a notebook — inspecting a trace, tweaking a
background window, checking wavelength residuals by eye — and stitched together into a
scripted, non-interactive reduction. This book teaches both styles, showing each step
interactively before composing the full end-to-end pipeline.

The book is organized as a reduction pipeline, in the order a real reduction runs: the
source is located and the frame is rectified before anything is measured from it, so each
chapter picks up where the previous one finished and ends at a wavelength-calibrated 1D
spectrum.

```{tableofcontents}
```

:::{note}
**Scope.** This book covers only the steps `specreduce` itself performs — tracing,
tilt correction, background subtraction, extraction, and wavelength calibration.
Detector-level reduction (dark and flat correction, cosmic-ray rejection) is the
domain of `ccdproc`/`astroscrappy` and is out of scope. Flux calibration has a chapter, but it is
a placeholder: `specreduce` gains the tools for it in v1.11.

**Data.** Every chapter works on the same real dataset: a long-slit observation of the
TrES-3 system taken with the OSIRIS spectrograph at the Gran Telescopio Canarias,
together with its HgAr, Ne, and Xe arc lamp frames. Chapter 1 introduces it.
:::

:::{note}
**Tilt correction is instrument-specific.** These OSIRIS frames have a strong tilt, so the
book rectifies them early and every later chapter works on the corrected frame. If your
spectrograph does not need it, skip that chapter — every call in the chapters that follow
is unchanged, and only the numbers differ.
:::

:::{warning}
Chapter 13 is a placeholder until `specreduce` v1.11 ships its flux-calibration tools.
:::
