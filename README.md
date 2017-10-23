# README #

Plugins and associated code for PYME.

This is currently used as a testbed for

* the new plugin system
* the recipe based processing of SMLM data
* new ideas for data processing

Once sufficiently mature and deemed of broad enough interest parts could migrate into more permanent core PYME code but that is entirely up to David's views.

As it is it serves as a nice testbed of the newish Python config system. For details on that either consult the docs in your local PYME installation, using a command like:

```
    pydoc PYME.config
```

or look at the [online version](http://www.python-microscopy.org/doc/api/PYME.config.html) at [python-microscopy.org](http://www.python-microscopy.org/) (which can be slightly out of date at times).

### Installation ###

A quick hack, in so far as I do not have a proper setuptools/distutils installation setup, but rather use the PYME.config system to pick up the costom plugins, as well as Python's PYTHONPATH to allow the required modules to be loaded at runtime:

- add the parent directory of this repo to your PYTHONPATH
- add the plugin entries into suitable config files into a directory searched by PYME.config

On my machine, in the per user PYME config file for visgui (~/.PYME/plugins/visgui/myplugins.txt), I currently have


```
#!python

PYMEcs.experimental.clusterTrack
PYMEcs.experimental.fiducials
PYMEcs.experimental.fiducialsNew
PYMEcs.experimental.qPAINT
PYMEcs.experimental.showErrs
PYMEcs.experimental.showShiftMap
PYMEcs.experimental.binEventProperty

```

There is currently only one plugin for the dsviewer (AKA dh5view), I have in a dsviewer plugin file (~/.PYME/plugins/dsviewer/myplugins.txt)

```
#!python

PYMEcs.experimental.showErrsDh5view

```
### Issues ###

Note that the showErrs modules rely on my mac installation which uses bash scripts and the [platypus app](https://sveinbjorn.org/platypus)
to capture STDERR into a temporary file which these modules access. Bottom line is that these modules will likely not work on anything but my local mac. It would be nice to build capturing STDIO/STDERR capture into the basic VisGUI and dh5View scripts as this would enable similar GUI based error inspection.

The PYME mac app wrappers are available at the [python-microscopy-osxapps bitpucket repository](http://bitbucket.org/christian_soeller/python-microscopy-osxapps) if somebody would like to inspect the approach.

### Contact ###

Christian Soeller (c.soeller at exeter.ac.uk)
