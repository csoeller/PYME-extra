from PYME.recipes.base import register_module, ModuleBase, Filter
from PYME.recipes.traits import HasTraits, Float, List, Bool, Int, CStr, Enum, on_trait_change, Input, Output

from PYME.IO.image import ImageStack
import numpy as np

import logging
logger = logging.getLogger(__name__)

@register_module('ExtractChannelByName')    
class ExtractChannelByName(ModuleBase):
    """Extract one channel from an image using regular expression matching to image channel names - by default this is case insensitive"""
    inputName = Input('input')
    outputName = Output('filtered_image')

    channelNamePattern = CStr('channel0')
    caseInsensitive = Bool(True)
    
    def _matchChannels(self,channelNames):
        # we put this into its own function so that we can call it externally for testing
        import re
        flags = 0
        if self.caseInsensitive:
            flags |= re.I
        idxs = [i for i, c in enumerate(channelNames) if re.search(self.channelNamePattern,c,flags)]
        return idxs
        
    def _pickChannel(self, image):     
        channelNames = image.mdh['ChannelNames']
        idxs = self._matchChannels(channelNames)
        if len(idxs) < 1:
            raise RuntimeError("Expression '%s' did not match any channel names" % self.channelNamePattern)
        if len(idxs) > 1:
            raise RuntimeError(("Expression '%s' did match more than one channel name: " % self.channelNamePattern) +
                               ', '.join([channelNames[i] for i in idxs]))
        idx = idxs[0]
        
        chan = image.data[:,:,:,idx]
        
        im = ImageStack(chan, titleStub = 'Filtered Image')
        im.mdh.copyEntriesFrom(image.mdh)
        im.mdh['ChannelNames'] = [channelNames[idx]]
        im.mdh['Parent'] = image.filename
        
        return im
    
    def execute(self, namespace):
        namespace[self.outputName] = self._pickChannel(namespace[self.inputName])


