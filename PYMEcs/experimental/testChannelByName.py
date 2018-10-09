import logging
logger = logging.getLogger(__file__)
import wx

from PYMEcs.recipes.base import ExtractChannelByName

def Msg(parent, message, caption = 'Message'):
    dlg = wx.MessageDialog(parent, message, caption, wx.OK)
    dlg.ShowModal()
    dlg.Destroy()

class TestChannelByName:
    """

    """
    def __init__(self, dsviewer):
        self.dsviewer = dsviewer
        self.ec = ExtractChannelByName()
        dsviewer.AddMenuItem('Experimental>Misc',
                          'Test Channel Name Matching by RegEx',
                          self.OnMatch,
                          helpText='test the regex matching of the ExtracChannelByName module')

        
    def OnMatch(self, event=None):        
        if not self.ec.configure_traits(kind='modal'):
            return
        mdh = self.dsviewer.image.mdh
        channelNames = mdh.getOrDefault('ChannelNames',[])
        matches = self.ec._matchChannels(channelNames)
        Msg(None,"ChannelNames: \n\t%s\n\nMatches :\n \t%s" % (channelNames, [channelNames[i] for i in matches]))

        
def Plug(dsviewer):
    """Plugs this module into the gui"""
    dsviewer.tChanByName = TestChannelByName(dsviewer)
