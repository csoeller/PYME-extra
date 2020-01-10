from PYMEcs.experimental.ExtraCmaps import cm_ch, cm_cbskb, cm_cbskb2
from PYME.misc.extraCMaps import regCmap
import pylab

registered = False

def Plug(arg):
    global registered
    
    if not registered:
        regCmap(cm_ch)
        regCmap(cm_cbskb)
        regCmap(cm_cbskb2)

        pylab.cm.cmapnames.sort()
        registered = True
