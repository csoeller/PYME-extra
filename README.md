# README #

Plugins and associated code for PYME.

This is a mixed bag of extensions/additions to standard PYME and also used as a testbed for

* the new plugin system
* the recipe based processing of SMLM data
* new ideas for data processing

As mentioned it serves as testbed for the newish PYME config system. For details on that either consult the docs in your local PYME installation, using a command like:

```
    pydoc PYME.config
```

or look at the [online version](http://www.python-microscopy.org/doc/api/PYME.config.html) at [python-microscopy.org](http://www.python-microscopy.org/) (which in theory could be slightly out of date if the latest build has been a while ago).

### Installation ###

With the new installation support the installation should be as easy as

```
   # cd to PYME-extra subdirectory
   python setup.py install
   python install_plugins.py dist
```

Alternatively, if you want to edit code in place you may want to use

```
   # cd to PYME-extra subdirectory
   python setup.py develop
   python install_plugins.py dist
```

For user installs, plugin files can be installed in the per user config directory by leaving out the `dist` argument to install_plugins.py:


```
   python install_plugins.py
```


### Issues ###

Note that the showErrs modules rely on a mac installation which uses bash scripts and the [platypus app](https://sveinbjorn.org/platypus) app
to capture STDERR into a temporary file which these modules access. 

Bottom line is that these modules will likely not work on anything but a mac with my PYMEapps wrappers. On other systems they will just generate a message that this functionality is not supported.

The PYME mac app wrappers are available at the [python-microscopy-osxapps bitpucket repository](http://bitbucket.org/christian_soeller/python-microscopy-osxapps) if somebody would like to inspect the approach.

### Contact ###

Christian Soeller (c.soeller _at_ gmail _dot_ com)
