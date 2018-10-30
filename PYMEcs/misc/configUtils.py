import PYME.config as config
import os

def get_legacy_scripts_dir():
    return os.path.join(os.path.dirname(config.__file__), 'Acquire/Scripts')

# this tries to replicate what config.get_init_filename does
#  if the config function were changed we would need to change this one as well
def get_init_directories_to_search():
    directories_to_search = [os.path.join(conf_dir, 'init_scripts') for conf_dir in config.config_dirs]
    
    extra_conf_dir = config.config.get('PYMEAcquire-extra_init_dir')
    if not extra_conf_dir is None:
        directories_to_search.insert(0, extra_conf_dir)

    directories_to_search.insert(0, legacy_scripts_directory=get_legacy_scripts_dir())
    
    return directories_to_search

def list_config_dirs():
    print('List of configuration directories:')
    for dir in config.config_dirs:
        print(dir)

def main():
    import sys
    import argparse
    import pprint

    # options parsing
    op = argparse.ArgumentParser(description='inspect PYME config files and settings.')
    op.add_argument('--initdirs', action='store_true',
                    help='list directories searched for init files')
    op.add_argument('-p','--parameters', action='store_true',
                    help='print configuration parameters', dest='config_params')
    op.add_argument('-d','--directories', action='store_true',
                    help='print configuration directories')
    op.add_argument('--protocols', action='store_true',
                    help='print custom protols found')
    op.add_argument('-i','--initfile', default=None,
                    help='locate init file and if found print full path')

    args = op.parse_args()

    if args.initdirs:
        print('List of directories searched for init scripts, legacy path included:')
        for dir in get_init_directories_to_search():
            print(dir)
        sys.exit(0)

    if args.directories:
        list_config_dirs()
        sys.exit(0)

    if args.config_params:
        print('List of configuration parameters:')
        for par in config.config.keys():
            print('%s : %s' % (par,config.config[par]))
        sys.exit(0)

    if args.protocols:
        prots = config.get_custom_protocols()
        print('Custom Protocols:')
        pprint.pprint(prots)
        sys.exit(0)

    if args.initfile is not None:
        inipath = config.get_init_filename(args.initfile, get_legacy_scripts_dir())
        if inipath is None:
            print("Initialisation file %s was not found" % args.initfile)
        else:
            print("Initialisation file %s was resolved as %s" %
                  (args.initfile,os.path.abspath(inipath)))
        sys.exit(0)


    # if we got here carry out a default action
    list_config_dirs()
    
if __name__ == "__main__":
    main()
