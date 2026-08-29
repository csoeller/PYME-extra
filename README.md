# README #

Plugins and associated code for PYME (AKA [python-microscopy](https://python-microscopy.org/)).

This is a mixed bag of extensions/additions to standard PYME and also used as a testbed/platform for

* recipe based processing of SMLM data
* implementing I/O and processing for new formats (especially MINFLUX data in `.npy` and, now preferred, `.zarr.zip` format)
* new ideas for data processing etc

### Installation ###

#### pip installing PYME-extra

A pip based install is currently the most straightforward way to install.

We highly recommend installing into a fresh virtual environment as can be generated with `conda` and related tools:

1. if you don't yet have it, download and install [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda) or [miniforge](https://github.com/conda-forge/miniforge).

2. create and activate a new conda environment with python 3.10 to 3.13 (3.12 and 3.13 are reasonably current at the time of writing - Aug 26 - and therefore recommended)

```
    conda create -n pyme-pip python=3.13
    conda activate pyme-pip
```

Now you are ready to use `pip` to install `python-microscopy` and `PYME-extra` (`#` is the comment sign on macOS/posix shells, just leave out that text when pasting to windows command prompts):

```
	# possibly install python-microscopy first and check that the install succeeds
	pip install python-microscopy
	# then install PYME-extra
	pip install PYME-extra # installation from PyPi
	pymex_install_plugins # important final step: register the plugins systemwide
```

This should be all, at this stage you can launch the main applications, i.e. `visgui` (AKA `PymeVis` or `PYMEVisualize`) and `dh5view` (AKA `PYMEImage`) to open and process microscopy data. See also [Verify Installation](https://www.python-microscopy.org/doc/Installation/Installation.html#verify-installation) in the [PYME docs](https://www.python-microscopy.org/doc).

#### PYME-test-env

For any development install of `PYME-extra` (i.e. using the latest source from github) we recommend installation as part of a [PYME-test-env](https://github.com/csoeller/PYME-test-env) controlled install. All further details please see there.


#### Installing from source

Although we recommend source installs as part of a [PYME-test-env](https://github.com/csoeller/PYME-test-env) controlled install you can also "manually" install directly from  source (i.e. as obtained from github). This is still achieved with pip but now called from within the source directory in which you unpacked PYME-extra (typically done when you cloned the git repository). A plain install from source is done with

```
	pip install .
```

A development install can be achieved by requesting an install in "editable mode". When a package is installed in editable mode, edits to the project source code become effective without the need of a new installation step.

```
	pip install --no-build-isolation -e . # install in in “development mode”
```

In either case (plain or development install), you may need to register the various plugins to implement the extra functionality provided by `PYME-extra`. This is achieved with the plugin installer that will have been installed with PYME-extra. You register with the command

```
	pymex_install_plugins
```

By default it registers the plugins systemwide but you can supply the `--user` option to register only for the current user:

```
	pymex_install_plugins --user # for further details see also pymex_install_plugins -h
```

#### Requirements

External modules required for full functionality currently include (NOTE: this list may fall out of sync with newer releases over time; the fully up-to-date list can always be found in the file `pyproject.toml`):

    python-microscopy
    statsmodels # for FRC smoothing with the lowess filter
    roifile     # to allow using info from ImageJ/Fiji ROIs
    colorcet    # add some colorcet colour tables in PYME
    circle-fit  # needs pip install to get recent version; for 2D NPC analysis
    alphashape # for cluster area and densities in clusters
    zarr>=2,<3 # for MINFLUX I/O
    seaborn # for some prettier plots
    mrcfile # to output 3D data for FSC from a EM FSC server

These should all be installed by the `pip` based install automatically.

We also often use a couple of extra dependencies in notebooks, but strictly speaking no functionality in `PYME-extra` depends directly on these:

    openpyxl
    tabulate

### Issues ###

Note that the showErrs modules rely on a mac installation which uses bash scripts and the [platypus app](https://sveinbjorn.org/platypus) app
to capture STDERR into a temporary file which these modules access. 

Bottom line is that these two error display modules will likely not work on anything but a mac with my PYMEapps wrappers. On other systems they will just generate a message that this functionality is not supported.

The PYME mac app wrappers are available at the [PYME-apps repository](https://github.com/csoeller/PYME-apps).

### Author ###

Christian Soeller
