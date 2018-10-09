#!/bin/bash

$PYTHON setup.py install

# Add more build steps here, if they are necessary.

echo "Installing plugin files into PYME config"

$PYTHON install_plugin.py dist


# See
# http://docs.continuum.io/conda/build.html
# for a list of environment variables that are set during the build process.
