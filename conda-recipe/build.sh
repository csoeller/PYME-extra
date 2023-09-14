#!/bin/bash


if [ -n "$OSX_ARCH" ]
    then
        $PYTHON setup.py build
        $PREFIX/python.app/Contents/MacOS/python setup.py install
else
    $PYTHON setup.py install
fi

$PYTHON install_plugins.py dist

# Add more build steps here, if they are necessary.


# See
# http://docs.continuum.io/conda/build.html
# for a list of environment variables that are set during the build process.
