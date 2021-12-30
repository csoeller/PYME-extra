import PYMEcs.misc.ExtraCmaps

# this module hijacks the plugin system to run colour map registration
# at startup of visgui, dsviewer

# I am exploring the alternative option to do this via the startup console script
# which has the advantage of doing this at startup for dsviewer which is an issue otherwise

def Plug(arg):
    # registerCmaps()
    pass
