import numpy as np
from traits.api import HasTraits, Str, Int, CStr, List, Enum, Float, Bool
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
    FilterEvents = Bool(False)

    traits_view = View(Group(Item(name = 'Channel'),
                             Item(name = 'FromTime'),
                             Item(name = 'ToTime'),
                             Item(name = 'FilterEvents')),
                       title = 'Select Channel',
                       buttons = OKCancelButtons
    )

    def add_keys(self,chans):
        for chan in chans:
            if chan not in self.clist:
                self.clist.append(chan)


class TimeBlock(HasTraits):
    BlockSize = Int(100)   


class ExtraColumns:
    """

    """
    def __init__(self, visFr):
        self.visFr = visFr
        self.pipeline = visFr.pipeline

        visFr.AddMenuItem('Experimental>Deprecated>ExtraColumns',
                          'Add random value column',
                          self.OnRandMap,
                          helpText='the added random value column can be used to select a fraction of all events')
        visFr.AddMenuItem('Experimental>Deprecated>ExtraColumns',
                          'Add time block choice column',
                          self.OnTimeBlocks,
                          helpText='the added column can be used to select blocks of frames as even or odd numbered blocks')
        visFr.AddMenuItem('Experimental>Deprecated>ExtraColumns',
                          'Select 50% of events from even numbered time blocks - EVEN',
                          self.OnSelectTB1,
                          helpText='select half of all events for FRC - even numbered time blocks')
        visFr.AddMenuItem('Experimental>Deprecated>ExtraColumns',
                          'Select 50% of events from from odd numbered time blocks - ODD',
                          self.OnSelectTB2,
                          helpText='select half of all events for FRC - second half')
        visFr.AddMenuItem('Experimental>Deprecated>ExtraColumns',
                          'Select 50% of events from random value column - half 1',
                          self.OnSelectHalf1,
                          helpText='select half of all events for FRC - first half')
        visFr.AddMenuItem('Experimental>Deprecated>ExtraColumns',
                          'Select 50% of events from random value column - half 2',
                          self.OnSelectHalf2,
                          helpText='select half of all events for FRC - second half')
        visFr.AddMenuItem('Experimental>Deprecated>ExtraColumns',
                          'Add channel based time selection column',
                          self.OnTimeSelectChannel,
                          helpText='the added column can be used to select a time window of events of a specific colour, filter tsel_channel in range (0.5,2)')
        visFr.AddMenuItem('Experimental>Deprecated>ExtraColumns',
                          'Add channel based random fraction selection column',
                          self.OnRandomSelectChannel,
                          helpText='the added column can be used to select a fraction of events of a specific colour by random sub-sampling, filter rand_channel in range (-.1,fraction)')

    def OnRandMap(self, event=None):
        self.pipeline.selectedDataSource.setMapping('randVal','0*x+np.random.rand(x.size)')
        self.pipeline.Rebuild()

    def OnTimeBlocks(self, event=None):
        tb = TimeBlock()
        if tb.configure_traits(kind='modal'):
            psd = self.pipeline.selectedDataSource
            blockSel = np.mod((psd['t']/tb.BlockSize).astype('int'),2)
            # self.pipeline.selectedDataSource.setMapping('timeBlock',"np.mod((t/%d).astype('int'),2)" % tb.BlockSize)
            self.pipeline.selectedDataSource.addColumn('timeBlock',blockSel)
            self.pipeline.Rebuild()

    def OnSelectTB1(self, event=None):
        if 'timeBlock' not in self.pipeline.keys():
            print('Need randVal property, call "Add random value column" first')
            return

        self.pipeline.filterKeys['timeBlock'] = (-0.1,0.5)
        self.pipeline.Rebuild()
        self.visFr.CreateFoldPanel()

    def OnSelectTB2(self, event=None):
        if 'timeBlock' not in self.pipeline.keys():
            print('Need randVal property, call "Add random value column" first')
            return

        self.pipeline.filterKeys['timeBlock'] = (0.5,1.5)
        self.pipeline.Rebuild()
        self.visFr.CreateFoldPanel()

            
    def OnSelectHalf1(self, event=None):
        if 'randVal' not in self.pipeline.keys():
            print('Need randVal property, call "Add random value column" first')
            return

        self.pipeline.filterKeys['randVal'] = (-0.1,0.4999)
        self.pipeline.Rebuild()
        self.visFr.CreateFoldPanel()

    def OnSelectHalf2(self, event=None):
        if 'randVal' not in self.pipeline.keys():
            print('Need randVal property, call "Add random value column" first')
            return

        self.pipeline.filterKeys['randVal'] = (0.5,1.1)
        self.pipeline.Rebuild()
        self.visFr.CreateFoldPanel()

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
            if timeSelector.FilterEvents:
                self.pipeline.filterKeys['tselect_%s' % timeSelector.Channel] = (0.5,2.0)
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
    visFr.extraCols = ExtraColumns(visFr)
