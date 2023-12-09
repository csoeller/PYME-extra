from __future__ import print_function
import sys

from PYME import config
import os
import sys
from distutils.dir_util import copy_tree

# this seems needed so that we can actually see the copying messages from copy_tree
# see also: https://stackoverflow.com/questions/17712722/dir-util-copy-tree-wont-print-the-files-that-it-copies

# note: on some platforms it may use the logging infrastructure, on others the distutils.log framework
# probably easiest to set both
import logging
logging.basicConfig(level=logging.INFO)

from distutils import log
log.set_verbosity(log.INFO)
log.set_threshold(log.INFO)

def main():
    this_dir = os.path.dirname(__file__)

    try:
        if sys.argv[1] == 'dist':
            installdir = config.dist_config_directory
    except IndexError:  # no argument provided, default to user config directory
        installdir = config.user_config_dir

    logging.info("installing plugin files into %s..." % installdir)
    copy_tree(os.path.join(this_dir, 'etc', 'PYME'), installdir, verbose=1)

if __name__ == '__main__':
    main()
