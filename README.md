# README #

Plugins and associated code for PYME (AKA [python-microscopy](https://python-microscopy.org/)).

This is a mixed bag of extensions/additions to standard PYME and also used as a testbed/platform for

* recipe based processing of SMLM data
* implementing I/O for new formats
* new ideas for data processing

### Installation ###

#### PYME-test-env

These days we recommend installing `PYME-extra` as part of a [PYME-test-env](https://github.com/csoeller/PYME-test-env) controlled install. All further details please see there.

#### Manual install

Installation is achieved with

```
   # cd to PYME-extra subdirectory
   python setup.py install
   python install_plugins.py dist
```

Alternatively, if you want to edit code in place you may want to use a development install

```
   # cd to PYME-extra subdirectory
   python setup.py develop
   python install_plugins.py dist
```

For user installs, plugin files can be installed in the per user config directory by leaving out the `dist` argument to install_plugins.py:

```
   python install_plugins.py
```
#### Requirements

External modules required for full functionality are listed in `requirements.txt`, these currently include

    python-microscopy
    statsmodels # for FRC smoothing with the lowess filter
    roifile     # to allow using info from ImageJ/Fiji ROIs
    colorcet    # add some colorcet colour tables in PYME
    circle-fit  # needs pip install to get recent version; for 2D NPC analysis
    alphashape # for cluster area and densities in clusters
    zarr>=2,<3 # for MINFLUX I/O
    seaborn # for some prettier plots
    mrcfile # to output 3D data for FSC from a EM FSC server

We also often use a couple more dependencies in notebooks, but strictly speaking no functionality in `PYME-extra` depends directly on these:

    openpyxl
    tabulate

### Issues ###

Note that the showErrs modules rely on a mac installation which uses bash scripts and the [platypus app](https://sveinbjorn.org/platypus) app
to capture STDERR into a temporary file which these modules access. 

Bottom line is that these two error display modules will likely not work on anything but a mac with my PYMEapps wrappers. On other systems they will just generate a message that this functionality is not supported.

The PYME mac app wrappers are available at the [PYME-apps repository](https://github.com/csoeller/PYME-apps).

### Contact ###

Christian Soeller (c.soeller _at_ gmail _dot_ com)
