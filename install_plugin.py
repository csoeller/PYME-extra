from __future__ import print_function
import sys

from PYME import config
import os
import sys
from distutils.dir_util import copy_tree

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def main():
    this_dir = os.path.dirname(__file__)

    try:
        if sys.argv[1] == 'dist':
            installdir = config.dist_config_directory
    except IndexError:  # no argument provided, default to user config directory
        installdir = config.user_config_dir

    eprint("installing plugin files into %s..." % installdir)
    copy_tree(os.path.join(this_dir, 'etc', 'PYME'), installdir, verbose=1)

if __name__ == '__main__':
    main()
