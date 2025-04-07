#!/usr/bin/python

##################
# <filename>.py
#
# Copyright David Baddeley, 2012
# d.baddeley@auckland.ac.nz
#
# This file may NOT be distributed without express permision from David Baddeley
#
##################
import wx

import numpy as np
from PYME.contrib.wxPlotPanel import PlotPanel
from PYME.IO import MetaDataHandler
from PYME.DSView import dsviewer as dsviewer
import PYME.IO.image as im
from PYMEcs.pyme_warnings import warn

import os

import logging
logger = logging.getLogger(__name__)

class TrackerPlotPanel(PlotPanel):
    def __init__(self, parent, driftTracker, *args, **kwargs):
        self.dt = driftTracker
        PlotPanel.__init__(self, parent, *args, **kwargs)    

    # add 5th suplot
    # replace 4th plot with offset and
    # new 5th subplot for z-pos (how calculated, z-nominal + dz?, remove offset)
    def draw(self):
        if self.IsShownOnScreen():
            if not hasattr( self, 'subplotxy' ):
                    self.subplotxy = self.figure.add_subplot( 411 )
                    self.subplotz = self.figure.add_subplot( 412 )
                    self.subploto = self.figure.add_subplot( 413 )
                    self.subplotc = self.figure.add_subplot(414)
                    # hopefully this sets it always for this fig, see
                    #   https://matplotlib.org/stable/tutorials/intermediate/tight_layout_guide.html
                    self.figure.set_tight_layout(True)

            try:
                t, dx_nm, dy_nm, dz_nm, corr, corrmax, poffset_nm, pos_um  = np.array(self.dt.get_history(1000)).T
            except ValueError:
                do_plot = False
            else:
                do_plot = True

            if do_plot:
                 # note: we now assume that all history values that are distances are provided in nm
                # this SHOULD match the conventions in driftTracking.py
                
                # a few reused variables
                tolnm = 1e3*self.dt.get_focus_tolerance()
                tdelta = t - self.dt.historyStartTime
                trange = [tdelta.min(), tdelta.max()]
                
                self.subplotxy.cla()
                self.subplotxy.plot(tdelta, dx_nm, 'r')
                self.subplotxy.plot(tdelta, dy_nm, 'g')
                self.subplotxy.set_ylabel('dx/dy (r/g) [nm]')
                self.subplotxy.set_xlim(*trange)
                
                self.subplotz.cla()
                self.subplotz.plot(tdelta, dz_nm, 'b')
                self.subplotz.plot([tdelta[0],tdelta[-1]],[tolnm,tolnm], 'g--')
                self.subplotz.plot([tdelta[0],tdelta[-1]],[-tolnm,-tolnm], 'g--')
                self.subplotz.set_ylabel('dz [nm]')
                self.subplotz.set_xlim(*trange)
                
                self.subploto.cla()
                self.subploto.plot(tdelta, poffset_nm, 'm')
                self.subploto.set_ylabel('offs [nm]')
                self.subploto.set_xlim(*trange)
                
                self.subplotc.cla()
                self.subplotc.plot(tdelta, corr/corrmax, 'r')
                self.subplotc.set_ylabel('C/C_m')
                self.subplotc.set_xlim(*trange)
                self.subplotc.set_xlabel('Time (s)')
    
            self.canvas.draw()


resx = []
resy = []
resz = []

class CalculateZfactorDialog(wx.Dialog):
    def __init__(self):
        self.Zfactorfilename = ''
        wx.Dialog.__init__(self, None, -1, 'Calculate Z-factor')
        sizer1 = wx.BoxSizer(wx.VERTICAL)

        pan = wx.Panel(self, -1)
        vsizermain = wx.BoxSizer(wx.VERTICAL)
        vsizer = wx.BoxSizer(wx.VERTICAL)

        hsizer2 = wx.BoxSizer(wx.HORIZONTAL)
        hsizer2.Add(wx.StaticText(pan, -1, '  Record a Z-stack with 101 slices & 50 nm step. Use the Z-stack to calculate the Z-factor'), 0, wx.ALL, 2)
        vsizer.Add(hsizer2)

        hsizer2 = wx.BoxSizer(wx.HORIZONTAL)
        self.bSelect = wx.Button(pan, -1, 'Select')
        self.bSelect.Bind(wx.EVT_BUTTON, self.OnSelect)
        hsizer2.Add(self.bSelect, 0, wx.ALL, 2)
        self.bPlot = wx.Button(pan, -1, 'Plot')
        self.bPlot.Bind(wx.EVT_BUTTON, self.OnPlot)
        hsizer2.Add(self.bPlot, 0, wx.ALL, 2)
        vsizer.Add(hsizer2)

        self.textZstackFilename = wx.StaticText(pan, -1, 'Z-stack file:    no file selected')
        vsizer.Add(self.textZstackFilename, 0, wx.ALL, 2)

        vsizermain.Add(vsizer, 0, 0, 0)

        self.plotPan = ZFactorPlotPanel(pan, size=(1200,600))
        vsizermain.Add(self.plotPan, 1, wx.EXPAND, 0)

        pan.SetSizerAndFit(vsizermain)
        sizer1.Add(pan, 1,wx.EXPAND, 0)
        self.SetSizerAndFit(sizer1)

    def OnSelect(self, event):
        dlg = wx.FileDialog(self, message="Open a Z-stack Image...", defaultDir=os.getcwd(), 
                            defaultFile="", style=wx.FD_OPEN)

        if dlg.ShowModal() == wx.ID_OK:
            self.Zfactorfilename = dlg.GetPath()

        dlg.Destroy()
        self.textZstackFilename.SetLabel('Z-stack file:   '+self.Zfactorfilename)

    
    def OnPlot(self, event):
        import PYMEcs.Analysis.offlineTracker as otrack

        ds = im.ImageStack(filename=self.Zfactorfilename)
        dataset = ds.data[:,:,:].squeeze()
        refim0 = dataset[:,:,10:91:4]
        calImages0, calFTs0, dz0, dzn0, mask0, X0, Y0 = otrack.genRef(refim0,normalised=False)

        del resx[:]
        del resy[:]
        del resz[:] # empty all these three lists every time before a new plot

        for i in range(dataset.shape[2]):
            image = dataset[:,:,i]
            driftx, drifty, driftz, cm, d = otrack.compare(calImages0, calFTs0, dz0, dzn0, 10, image, mask0, X0, Y0, deltaZ=0.2)
            resx.append(driftx)
            resy.append(drifty)
            resz.append(driftz)

        self.plotPan.draw()
        self.plotPan.Refresh()


class ZFactorPlotPanel(PlotPanel):

    def draw(self):
        dznm = 1e3*np.array(resz)
        dxnm = 110*np.array(resx)
        dynm = 110*np.array(resy)
        t = np.arange(dznm.shape[0])

        dzexp = dznm[50-4:50+5]
        dztheo = np.arange(-200,201,50)
        x = np.arange(-150,151,50)
        y = dznm[50-3:50+4]
        m, b = np.polyfit(x,y,1)
        Zfactor = 1.0/m

        if not hasattr( self, 'subplot' ):
                self.subplot1 = self.figure.add_subplot( 121 )
                self.subplot2 = self.figure.add_subplot( 122 )

        self.subplot1.cla()

        self.subplot1.scatter(t,-dznm,s=5)
        self.subplot1.plot(-dznm, label='z')
        self.subplot1.plot(-dxnm, label='x')
        self.subplot1.plot(-dynm, label='y')

        self.subplot1.grid()
        self.subplot1.legend()

        self.subplot2.cla()

        self.subplot2.plot(dztheo,dzexp,'-o')
        self.subplot2.plot(dztheo,1.0/Zfactor*dztheo,'--')
        font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
        self.subplot2.text(-50, 100, 'Z-factor = %3.1f' % Zfactor, fontdict=font)
        self.subplot2.grid()

        self.canvas.draw()

from PYME.DSView import overlays
import weakref
class DriftROIOverlay(overlays.Overlay):
    def __init__(self, driftTracker):
        self.dt = driftTracker
    
    def __call__(self, view, dc):
        if self.dt.sub_roi is not None:
            dc.SetPen(wx.Pen(colour=wx.CYAN, width=1))
            dc.SetBrush(wx.TRANSPARENT_BRUSH)
            x0, x1, y0, y1 = self.dt.sub_roi
            x0c, y0c = view.pixel_to_screen_coordinates(x0, y0)
            x1c, y1c = view.pixel_to_screen_coordinates(x1, y1)
            sX, sY = x1c-x0c, y1c-y0c
            dc.DrawRectangle(int(x0c), int(y0c), int(sX), int(sY))
            dc.SetPen(wx.NullPen)
        else:
            dc.SetBackground(wx.TRANSPARENT_BRUSH)
            dc.Clear()

from PYME.recipes.traits import HasTraits, Float, Enum, CStr, Bool, Int, List
class DriftTrackConfig(HasTraits):
    zfocusTolerance_nm = Float(50.0,label='ZFocus tolerance in nm',
                               desc="when moving outside of the focus tolerance the z piezo will counteract")
    deltaZ_nm = Float(200,desc='spacing between the planes during recording of the calibration stack in nm',
                      label='Z plane spacing in nm')
    stackHalfSize = Int(35,label='Stack Half Size',
                        desc='number of planes either side of the focal plane that are recorded for calibration')
    minDelay = Int(10,label='minimum delay in frames',
                   desc='minimum delay time in frames between corrections')
    zFactor = Float(1.0,label='z correction factor',
                    desc='z correction factor to match physical movement to sensed z movement')
    plotInterval = Int(10,label='plot interval in frames',
                       desc='interval at which the plots are updated, in units of frames')

from PYMEcs.Acquire.Hardware.driftTracking_n import State
class DriftTrackingControl(wx.Panel):
    def __init__(self, main_frame, driftTracker, winid=-1, showPlots=True):
        ''' This class provides a GUI for controlling the drift tracking system. 
        
        It should be initialised with a reference to the PYMEAcquire main frame, which will stand in as a parent while other GUI items are
        created. Note that the actual parent will be reassigned once the GUI tool panel is created using a Reparent() call.
        '''
        wx.Panel.__init__(self, main_frame, winid)

        main_frame.AddMenuItem('DriftTracking', "Change config parameters", self.OnDriftTrackConfig)
        main_frame.AddMenuItem('DriftTracking', "Save history", self.OnBSaveHist)
        main_frame.AddMenuItem('DriftTracking', "Save calibration stack", self.OnBSaveCalib)
        main_frame.AddMenuItem('DriftTracking', "Calculate z factor", self.OnBCalculateZfactor)
        
        self.dt = driftTracker
        self.dtconfig = DriftTrackConfig(zfocusTolerance_nm=1e3*self.dt.focusTolerance,
                                         deltaZ_nm=1e3*self.dt.deltaZ,
                                         stackHalfSize=self.dt.stackHalfSize)
        self.plotInterval = 10
        self.showPlots = showPlots

        # keep a reference to the main frame. Do this as a weakref to avoid circular references.
        # we need this to be able to access the view to get the current selection and to add overlays.
        self._main_frame = weakref.proxy(main_frame)
        self._view_overlay = None # dummy reference to the overlay so we only create it once

        sizer_1 = wx.BoxSizer(wx.VERTICAL)
        hsizer = wx.BoxSizer(wx.HORIZONTAL)
        
        self.cbTrack = wx.CheckBox(self, -1, 'Track')
        hsizer.Add(self.cbTrack, 0, wx.ALL, 2) 
        self.cbTrack.Bind(wx.EVT_CHECKBOX, self.OnCBTrack)
        self.cbLock = wx.CheckBox(self, -1, 'Lock')
        self.cbLock.Bind(wx.EVT_CHECKBOX, self.OnCBLock)
        hsizer.Add(self.cbLock, 0, wx.ALL, 2)
        self.cbLockActive = wx.CheckBox(self, -1, 'Lock Active')
        self.cbLockActive.Enable(False)
        hsizer.Add(self.cbLockActive, 0, wx.ALL, 2)        
        sizer_1.Add(hsizer, 0, wx.EXPAND, 0)
        
        hsizer = wx.BoxSizer(wx.HORIZONTAL)
        self.bSetPostion = wx.Button(self, -1, 'Set focus to current')
        hsizer.Add(self.bSetPostion, 0, wx.ALL, 2) 
        self.bSetPostion.Bind(wx.EVT_BUTTON, self.OnBSetPostion)
        sizer_1.Add(hsizer, 0, wx.EXPAND, 0)
        
        hsizer = wx.BoxSizer(wx.HORIZONTAL)
        self.tbSubROI = wx.ToggleButton(self, -1, 'Restrict to sub-ROI')
        hsizer.Add(self.tbSubROI, 0, wx.ALL, 2)
        self.tbSubROI.Bind(wx.EVT_TOGGLEBUTTON, self.OnTBToggleSubROI)
        sizer_1.Add(hsizer, 0, wx.EXPAND, 0)

        hsizer = wx.BoxSizer(wx.HORIZONTAL)
        hsizer.Add(wx.StaticText(self, -1, "Calibration:"), 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 2)
        self.gCalib = wx.Gauge(self, -1, 11)
        hsizer.Add(self.gCalib, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 2)
        self.stCalibState = wx.StaticText(self, -1, "UNCALIBRATED")
        hsizer.Add(self.stCalibState, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 2) # second arg 0 or 1?
        sizer_1.Add(hsizer, 0, wx.EXPAND, 0)
        
        hsizer = wx.BoxSizer(wx.HORIZONTAL)
        self.stConfig = wx.StaticText(self, -1,
                                      'tol: %d nm, delZ: %.0f nm, stHsz: %d' %
                                      (self.dtconfig.zfocusTolerance_nm,
                                       self.dtconfig.deltaZ_nm,
                                       self.dtconfig.stackHalfSize), size=[400,-1])
        cfont = self.stConfig.GetFont()
        font = wx.Font(cfont.GetPointSize(), wx.TELETYPE, wx.NORMAL, wx.NORMAL)
        self.stConfig.SetFont(font)
        hsizer.Add(self.stConfig, 0, wx.ALL, 2)        
        sizer_1.Add(hsizer,0, wx.EXPAND, 0)

        hsizer = wx.BoxSizer(wx.HORIZONTAL)
        self.stError = wx.StaticText(self, -1, 'Error:\n\n', size=[200,-1])
        cfont = self.stError.GetFont()
        font = wx.Font(cfont.GetPointSize(), wx.TELETYPE, wx.NORMAL, wx.NORMAL)
        self.stError.SetFont(font)
        hsizer.Add(self.stError, 0, wx.ALL, 2)        
        sizer_1.Add(hsizer,0, wx.EXPAND, 0)
        
        if self.showPlots:
            self.trackPlot = TrackerPlotPanel(self, self.dt, size=[300, 500])
            
            #hsizer.Add(self.stError, 0, wx.ALL, 2)        
            sizer_1.Add(self.trackPlot,0, wx.EXPAND, 0)
        
        self.SetAutoLayout(1)
        self.SetSizer(sizer_1)
        sizer_1.Fit(self)
        sizer_1.SetSizeHints(self)
        self.Layout()
        # end wxGlade

    def OnCBTrack(self, event):
        #print self.cbTrack.GetValue()
        if self.cbTrack.GetValue():
            self.dt.register()
        else:
            self.dt.deregister()
            
    def OnBSetPostion(self, event):
        self.dt.reCalibrate()
        
    def OnTBToggleSubROI(self, event):
        self.toggle_subroi(self.tbSubROI.GetValue())
    
    def toggle_subroi(self, new_state=True):
        ''' Turn sub-ROI tracking on or off, using the current selection in the live image display'''
        if new_state:
            x0, x1, y0, y1, _, _ = self._main_frame.view.do.sorted_selection
            self.dt.set_subroi((x0, x1, y0, y1))
        else:
            self.dt.set_subroi(None)

        if self._view_overlay is None:
            self._view_overlay = self._main_frame.view.add_overlay(DriftROIOverlay(self.dt), 'Drift tracking Sub-ROI')
            
    def OnBSaveCalib(self, event):
        if not hasattr(self.dt, 'state') or (self.dt.state != State.CALIBRATED):
            warn("not in a calibrated state (state is %s), cannot save" % self.dt.state)
        else:
            self.showCalImages()

    def showCalImages(self):
        import numpy as np
        import time
        ds2 = self.dt.refImages

        #metadata handling
        mdh = MetaDataHandler.NestedClassMDHandler()
        mdh.setEntry('StartTime', time.time())
        mdh.setEntry('AcquisitionType', 'Stack')

        #loop over all providers of metadata
        for mdgen in MetaDataHandler.provideStartMetadata:
            mdgen(mdh)
        mdh.setEntry('CalibrationPositions',self.dt.calPositions)

        im = dsviewer.ImageStack(data = ds2, mdh = mdh, titleStub='Unsaved Image')
        if not im.mode == 'graph':
            im.mode = 'lite'

        #print im.mode
        dvf = dsviewer.DSViewFrame(im, mode= im.mode, size=(500, 500))
        dvf.SetSize((500,500))
        dvf.Show()
        
    def OnBSaveHist(self, event):
        if not hasattr(self.dt, 'history') or (len(self.dt.history) <= 0):
            warn("no history that could be saved")
        else:
            dlg = wx.FileDialog(self, message="Save file as...",  
                                defaultFile='history.txt',
                                wildcard='txt File (*.txt)|*.txt',
                                style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)

            if dlg.ShowModal() == wx.ID_OK:
                historyfn = dlg.GetPath()
                np.savetxt(historyfn, self.dt.history, header=' '.join(self.dt.historyColNames))
                dlg = wx.MessageDialog(self._main_frame, "history saved", "Info", wx.OK | wx.ICON_INFORMATION)
                dlg.ShowModal()
                dlg.Destroy()


    def OnDriftTrackConfig(self, event):
        if not self.dtconfig.configure_traits(kind='modal'):
            return
        dtc = self.dtconfig
        self.dt.set_focus_tolerance(1e-3*dtc.zfocusTolerance_nm)
        self.dt.set_delta_Z(1e-3*dtc.deltaZ_nm)
        self.dt.set_stack_halfsize(dtc.stackHalfSize)
        self.dt.minDelay = dtc.minDelay
        self.dt.Zfactor = dtc.zFactor
        self.plotInterval = dtc.plotInterval

        self.stConfig.SetLabel('tol: %d nm, delZ: %.0f nm, stHsz: %d' %
                                      (dtc.zfocusTolerance_nm,
                                       dtc.deltaZ_nm,
                                       dtc.stackHalfSize))
        
    def OnBSetTolerance(self, event):
        self.dt.set_focus_tolerance(float(self.tTolerance.GetValue())/1e3)

    def OnBSetdeltaZ(self, event):
        self.dt.set_delta_Z(float(self.tdeltaZ.GetValue())/1e3)

    def OnBSetHalfsize(self, event):
        self.dt.set_stack_halfsize(int(self.tHalfsize.GetValue()))

    def OnBSetZfactor(self, event):
        self.dt.Zfactor = float(self.tZfactor.GetValue())

    def OnBCalculateZfactor(self, event):
        dlg = CalculateZfactorDialog()
        ret = dlg.ShowModal()
        dlg.Destroy()

    def OnBSetMinDelay(self, event):
        self.dt.minDelay = int(self.tMinDelay.GetValue())

    def OnBSetPlotInterval(self, event):
        self.plotInterval = int(self.tPlotInterval.GetValue())
    
    def OnCBLock(self, event):
        self.dt.set_focus_lock(self.cbLock.GetValue())

    def refresh(self):
        try:
            calibState, NCalibFrames, calibCurFrame = self.dt.get_calibration_progress()
            self.gCalib.SetRange(int(NCalibFrames)) # needs to be int?
            self.gCalib.SetValue(int(calibCurFrame)) # needs to be int?
            self.stCalibState.SetLabel(calibState.name)

            try:
                t, dx_nm, dy_nm, dz_nm, corr, corrmax, poffset_nm, pos_um = self.dt.get_history(1)[-1]
                self.stError.SetLabel(("Error: x = %s nm y = %s nm z = %s nm\n" +
                                       "offs = %s nm c/cm = %4.2f") %
                                      ("{:>+3.2f}".format(dx_nm), "{:>+3.2f}".format(dy_nm),
                                       "{:>+3.1f}".format(dz_nm), "{:>+6.1f}".format(poffset_nm),
                                       corr/corrmax))

            except IndexError:
                pass

            self.cbLock.SetValue(self.dt.get_focus_lock())
            self.cbTrack.SetValue(self.dt.is_tracking())
            self.cbLockActive.SetValue(self.dt.lockActive)
            
            if (len(self.dt.get_history(0)) > 0) and (len(self.dt.get_history(0)) % self.plotInterval == 0) and self.showPlots:
                self.trackPlot.draw()
        except (AttributeError, IndexError):
            logger.exception('error in refresh')
            pass
