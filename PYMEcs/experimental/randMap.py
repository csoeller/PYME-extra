import numpy as np
from traits.api import HasTraits, Str, Int, CStr, List, Enum, Float
from traitsui.api import View, Item, Group
from traitsui.menu import OKButton, CancelButton, OKCancelButtons

class ChannelSelector(HasTraits):
    clist = List([])
    Channel = Enum(values='clist')

    traits_view = View(Group(Item(name = 'Channel')),
                             title = 'Select Channel',
                             buttons = OKCancelButtons
    )

    def add_keys(self,chans):
        for chan in chans:
            if chan not in self.clist:
                self.clist.append(chan)


class TimeSelector(HasTraits):
    clist = List([])
    Channel = Enum(values='clist')
    FromTime = Float()
    ToTime = Float()

    traits_view = View(Group(Item(name = 'Channel'),
                             Item(name = 'FromTime'),
                             Item(name = 'ToTime')),
                       title = 'Select Channel',
                       buttons = OKCancelButtons
    )

    def add_keys(self,chans):
        for chan in chans:
            if chan not in self.clist:
                self.clist.append(chan)



class RandMap:
    """

    """
    def __init__(self, visFr):
        self.visFr = visFr
        self.pipeline = visFr.pipeline

        visFr.AddMenuItem('Experimental>ExtraColumns',
                          'Add random value column',
                          self.OnRandMap,
                          helpText='the added random value column can be used to select a fraction of all events')
        visFr.AddMenuItem('Experimental>ExtraColumns',
                          'Add channel based time selection column',
                          self.OnTimeSelectChannel,
                          helpText='the added column can be used to select a time window of events of a specific colour, filter tsel_channel in range (0.5,2)')
        visFr.AddMenuItem('Experimental>ExtraColumns',
                          'Add channel based random fraction selection column',
                          self.OnRandomSelectChannel,
                          helpText='the added column can be used to select a fraction of events of a specific colour by random sub-sampling, filter rand_channel in range (-.1,fraction)')

    def OnRandMap(self, event=None):
        self.pipeline.selectedDataSource.setMapping('randVal','0*x+np.random.rand(x.size)')
        self.pipeline.Rebuild()

    def OnTimeSelectChannel(self, event=None):
        pipeline = self.pipeline
        if pipeline.selectedDataSource is None:
            return
        if len(pipeline.colourFilter.getColourChans()) < 1:
            return
        timeSelector = TimeSelector()
        timeSelector.add_keys(pipeline.colourFilter.getColourChans())
        if timeSelector.configure_traits(kind='modal'):
            psd = pipeline.selectedDataSource
            tall = np.ones_like(psd['x'], dtype='float')
            tsel = (psd['t'] >= timeSelector.FromTime) * (psd['t'] <= timeSelector.ToTime)

            dispColor = pipeline.colourFilter.currentColour
            pipeline.colourFilter.setColour(timeSelector.Channel)
            idx = pipeline.filter.Index.copy()
            idx[idx] = pipeline.colourFilter.index
            tall[idx] = tsel[idx]
            pipeline.colourFilter.setColour(dispColor)
            pipeline.selectedDataSource.addColumn('tselect_%s' % timeSelector.Channel,tall)
            pipeline.Rebuild()
            self.visFr.CreateFoldPanel()

    def OnRandomSelectChannel(self, event=None):
        pipeline = self.pipeline
        if pipeline.selectedDataSource is None:
            return
        if len(pipeline.colourFilter.getColourChans()) < 1:
            return
        chanSelector = ChannelSelector()
        chanSelector.add_keys(pipeline.colourFilter.getColourChans())
        if chanSelector.configure_traits(kind='modal'):
            psd = pipeline.selectedDataSource
            eall = np.zeros_like(psd['x'], dtype='float')
            erand = np.random.rand(eall.size)

            dispColor = pipeline.colourFilter.currentColour
            pipeline.colourFilter.setColour(chanSelector.Channel)
            idx = pipeline.filter.Index.copy()
            idx[idx] = pipeline.colourFilter.index
            eall[idx] = erand[idx]
            pipeline.colourFilter.setColour(dispColor)
            pipeline.selectedDataSource.addColumn('rand_%s' % chanSelector.Channel,eall)
            pipeline.Rebuild()
            self.visFr.CreateFoldPanel()

            
def Plug(visFr):
    '''Plugs this module into the gui'''
    visFr.randMap = RandMap(visFr)
