import logging
logger = logging.getLogger(__file__)

from PYME.DSView import ViewIm3D
from PYME.IO.image import ImageStack

class NPCanalysis:

    def __init__(self, dsviewer):
        self.dsviewer = dsviewer
        dsviewer.AddMenuItem('Experimental>NPCs',
                          'Generate NPC mask from Fiji ROISet',
                          self.OnNPCROISet,
                          helpText='Load at FIJI ROISet (e.g. ROISet.zip) and construct a NPC mask from it')

    def OnNPCROISet(self, event=None):
        NPCmaskImg = ImageStack(npcmask, titleStub = 'NPC region mask')
        NPCmaskImg.mdh.copyEntriesFrom(self.dsviewer.image.mdh)
        # add entries for the mask
        NPCmaskImg.mdh['ROISet'] = NPCROI_filename

def Plug(dsviewer):
    """Plugs this module into the gui"""
    NPCanalysis(dsviewer)
