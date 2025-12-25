from PYME import config
import os
import logging

DSVIEWER_PLUGINS = \
'''
PYMEcs.experimental.showErrsDh5view
PYMEcs.experimental.mapTools
PYMEcs.experimental.meas2DplotDh5view
PYMEcs.experimental.testChannelByName
PYMEcs.experimental.FRC
PYMEcs.experimental.regExtraCmaps
PYMEcs.experimental.procPoints
PYMEcs.experimental.combine_maps
PYMEcs.experimental.Sofi
PYMEcs.experimental.CalcZfactor
PYMEcs.experimental.ImageJROItools
'''

VISGUI_PLUGINS = \
'''
PYMEcs.experimental.clusterTrack
PYMEcs.experimental.fiducials
PYMEcs.experimental.fiducialsNew
PYMEcs.experimental.qPAINT
PYMEcs.experimental.showErrs
PYMEcs.experimental.showShiftMap
PYMEcs.experimental.binEventProperty
PYMEcs.experimental.onTimes
PYMEcs.experimental.snrEvents
PYMEcs.experimental.randMap
PYMEcs.experimental.mortensen
PYMEcs.experimental.splitRender
PYMEcs.experimental.timedSpecies
PYMEcs.experimental.chaining
PYMEcs.experimental.specLabeling
PYMEcs.experimental.selectROIfilterTable
PYMEcs.experimental.regExtraCmaps
PYMEcs.experimental.pyme2caml
PYMEcs.experimental.Simpler
PYMEcs.experimental.MINFLUX
PYMEcs.experimental.eventProcessing
PYMEcs.experimental.NPCcalcLM
'''

RECIPES = \
'''
PYMEcs.recipes.processing
PYMEcs.recipes.output
PYMEcs.recipes.base
PYMEcs.recipes.localisations
PYMEcs.recipes.simpler
'''

def get_legacy_scripts_dir():
    return os.path.join(os.path.dirname(config.__file__), 'Acquire/Scripts')

# this tries to replicate what config.get_init_filename does
#  if the config function were changed we would need to change this one as well
def get_init_directories_to_search():
    directories_to_search = [os.path.join(conf_dir, 'init_scripts') for conf_dir in config.config_dirs]
    
    extra_conf_dir = config.config.get('PYMEAcquire-extra_init_dir')
    if not extra_conf_dir is None:
        directories_to_search.insert(0, extra_conf_dir)

    directories_to_search.insert(0,get_legacy_scripts_dir())
    
    return directories_to_search

def list_config_dirs():
    print('List of configuration directories:')
    for dir in config.config_dirs:
        print(dir)

def install_plugins():
    # options parsing
    import argparse
    from pathlib import Path

    logging.basicConfig(level=logging.INFO)
    op = argparse.ArgumentParser(description='install PYME-extra plugins')
    op.add_argument('-u','--user', action='store_true',
                    help='install plugin info to user config directory')
    op.add_argument('--dry-run',action="store_true",
                    help='just process options and merely show what would be done')

    args = op.parse_args()
    if args.user:
        installdir = Path(config.user_config_dir)
    else:
        installdir = Path(config.dist_config_directory)

    recfile = installdir / 'plugins' / 'recipes' / 'PYMEcsRecipePlugins.txt'
    visguifile = installdir / 'plugins' / 'visgui' / 'PYMEcsVisguiPlugins.txt'
    dsviewerfile = installdir / 'plugins' / 'dsviewer' / 'PYMEcsDsviewerPlugins.txt'

    logging.info("will install recipes to\n\t%s" % recfile)
    logging.info("will install visgui plugins to\n\t%s" % visguifile)
    logging.info("will install dsviewer plugins to\n\t%s" % dsviewerfile)
    
    if args.dry_run:
        logging.info("dry run, aborting...")
        import sys
        sys.exit(0)

    recfile.write_text(RECIPES)
    visguifile.write_text(VISGUI_PLUGINS)
    dsviewerfile.write_text(DSVIEWER_PLUGINS)

    logging.info("\nPYME-extra recipes and plugins have been registered...")

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
