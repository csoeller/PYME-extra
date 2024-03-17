# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 15:02:50 2014

@author: David Baddeley
"""

import numpy as np
from numpy.fft import fftn, ifftn, fftshift, ifftshift

import time
from scipy import ndimage
from PYME.Acquire import eventLog

import threading

import logging
logger = logging.getLogger(__name__)
    
from PYME.contrib import dispatch
class StandardFrameSource(object):
    '''This is a simple source which emits frames once per polling interval of the frameWrangler 
    (i.e. corresponding to the onFrameGroup signal of the frameWrangler).
    
    The intention is to reproduce the historical behaviour of the drift tracking code, whilst
    abstracting some of the detailed knowledge of frame handling out of the actual tracking code. 
    
    '''
    def __init__(self, frameWrangler):
        self._fw = frameWrangler
        self._on_frame = dispatch.Signal(['frameData'])
        self._fw.onFrameGroup.connect(self.tick)


    def tick(self, *args, **kwargs):
        self._on_frame.send(sender=self, frameData=self._fw.currentFrame)

    @property
    def shape(self):
        return self._fw.currentFrame.shape
    
    def connect(self, callback):
        self._on_frame.connect(callback)

    def disconnect(self, callback):
        self._on_frame.disconnect(callback)

class OIDICFrameSource(StandardFrameSource):
    """ Emit frames from the camera to the tracking code only for a single OIDIC orientation.

        Currently a straw man / skeleton pending details of OIDIC code.

        TODO - should this reside here, or with the other OIDIC code (which I believe to be in a separate repo)?
    
    """

    def __init__(self, frameWrangler, oidic_controller, oidic_orientation=0):
        super().__init__(frameWrangler)

        self._oidic = oidic_controller
        self._target_orientation = oidic_orientation

    def tick(self, *args, **kwargs):
        # FIXME - change to match actual naming etc ... in OIDIC code.
        # FIXME - check when onFrameGroup is emitted relative to when the OIDIC orientation is set.
        # Is this predictable, or does it depend on the order in which OIDIC and drift tracking are
        # registered with the frameWrangler?
        if self._oidic.orientation == self._target_orientation:
            super().tick(*args, **kwargs)
        else:
            # clobber all frames coming from camera when not in the correct DIC orientation
            pass

from enum import Enum
State = Enum('State', ['UNCALIBRATED', 'CALIBRATING', 'FINISHING_CALIBRATION', 'CALIBRATED'])

class Correlator(object):
    def __init__(self, scope, piezo=None, frame_source=None, sub_roi=None, focusTolerance=.05, deltaZ=0.2, stackHalfSize=35):
        self.piezo = piezo

        if frame_source is None:
            self.frame_source = StandardFrameSource(scope.frameWrangler)
        
        # configuration parameters, accessible via keyword args
        self.focusTolerance = focusTolerance # (in um) how far focus can drift before we correct
        self.deltaZ = deltaZ #z increment (in um) used for calibration
        self.stackHalfSize = stackHalfSize
        # other configuration parameters - not currently accessible via keyword args
        self.skipframes = 1 # number of frames to skip after changing position to let piezo settle
        self.minDelay = 10
        # NOTE: maxTotalCorrection parameter below - we may want to have a reset offset method
        #       to allow resetting the offset during the day as it tends to accumulate;
        #       this has already led to the lock giving up on occasion
        self._maxTotalCorrection = 20.0 # maximum total correction in um
        self.Zfactor = 1.0
        self.logShifts = True

        # we report our tracking info in nm by default
        pixelsize_um = scope.GetPixelSize()
        self.conversion = {'x': 1e3*pixelsize_um[0], 'y':1e3*pixelsize_um[1], 'z':1e3}
        # we may or may not want to use the one below
        # self.trackunits = {'x':'nm', 'y':'nm', 'z':'nm'}

        # 'state-tracking' variables
        self.NCalibFrames = 2*self.stackHalfSize + 1 # this gets recalculated in _prepare_calibration anyway
        self.calibCurFrame = 0
        self.state = State.UNCALIBRATED
        self.tracking = False
        self.lockActive = False
        self.lockFocus = False
        self._last_target_z = -1
  
    def set_subroi(self, bounds):
        """ Set the position of the roi to crop

        Parameters
        ----------

        position : tuple
            The pixel position (x0, x1, y0, y1) in int
        """

        self.sub_roi = bounds
        self.reCalibrate()

    def _crop_frame(self, frame_data):
        if self.sub_roi is None:
            return frame_data.squeeze()   # we may as well do the squeeze here to avoid lots of squeezes elsewhere
        else:
            x0, x1, y0, y1 = self.sub_roi
            return frame_data.squeeze()[x0:x1, y0:y1]

    def set_focus_tolerance(self, tolerance):
        """ Set the tolerance for locking position

        Parameters
        ----------

        tolerance : float
            The tolerance in um
        """

        self.focusTolerance = tolerance

    def get_focus_tolerance(self):
        return self.focusTolerance

    def set_delta_Z(self, delta):
        """ Set the Z increment for calibration stack

        Parameters
        ----------

        delta : float
            The delta in um. This should be the distance over which changes in PSF intensity with depth 
            can be approximated as being linear, with an upper bound of the Nyquist sampling in Z. 
            At Nyquist sampling, the linearity assumption is already getting a bit tenuous. Default = 0.2 um, 
            which is approximately Nyquist sampled at 1.4NA.
        """

        self.deltaZ = delta

    def get_delta_Z(self):
        return self.deltaZ

    def set_stack_halfsize(self, halfsize):
        """ Set the calibration stack half size

        This dictates the maximum size of z-stack you can record whilst retaining focus lock. The resulting 
        calibration range can be calculated as deltaZ*(2*halfsize), and should extend about 1 micron above 
        and below the size of the the largest z-stack to ensure that lock can be maintained at the edges of 
        the stack. The default of 35 gives about 12 um of axial range.

        Parameters
        ----------

        halfsize : int
        """

        self.stackHalfSize = halfsize

    def get_stack_halfsize(self):
        return self.stackHalfSize   

    def set_focus_lock(self, lock=True):
        """ Set locking on or off

        Parameters
        ----------

        lock : bool
            whether the lock should be on
        """

        self.lockFocus = lock

    def get_focus_lock(self):
        return self.lockFocus

    def get_history(self, length=1000):
        try:
            return self.history[-length:]
        except AttributeError:
            return []

    def get_calibration_progress(self):
        """ Returns the current calibration progress:
            calibration is complete when state == State.CALIBRATED
        """

        if self.state == State.UNCALIBRATED:
            curFrame = 0
        elif self.state in [State.CALIBRATING,State.FINISHING_CALIBRATION]:
            curFrame = self.calibCurFrame
        else:
            curFrame = self.NCalibFrames
        return (self.state, self.NCalibFrames, curFrame)

    def is_tracking(self):
        return self.tracking

    def _setRefN(self, frame_data, N):
        d = 1.0*frame_data.squeeze()
        ref = d/d.mean() - 1
        self.refImages[:,:,N] = ref        
        self.calFTs[:,:,N] = ifftn(ref)
        self.calImages[:,:,N] = ref*self.mask

    def _prepare_calibration(self, frame_data):
        # this should not be necessary as we SHOULD only get called if we are uncalibrated
        # self.state = State.UNCALIBRATED #completely uncalibrated
        if self.state != State.UNCALIBRATED:
            raise RuntimeError("this method should only be called when uncalibrated; actual state is %s" % self.state)

        d = 1.0*frame_data.squeeze()        
        
        self.X, self.Y = np.mgrid[0.0:d.shape[0], 0.0:d.shape[1]]
        self.X -= np.ceil(d.shape[0]*0.5)
        self.Y -= np.ceil(d.shape[1]*0.5)
        
        #we want to discard edges after accounting for x-y drift
        self.mask = np.ones_like(d)
        self.mask[:10, :] = 0
        self.mask[-10:, :] = 0
        self.mask[:, :10] = 0
        self.mask[:,-10:] = 0
       
        self.corrRef = 0

        self.calPositions = self.homePos + self.deltaZ*np.arange(-float(self.stackHalfSize), float(self.stackHalfSize + 1))
        self.NCalibFrames = len(self.calPositions)
            
        self.refImages = np.zeros(self.mask.shape[:2] + (self.NCalibFrames,))
        self.calImages = np.zeros(self.mask.shape[:2] + (self.NCalibFrames,))
        self.calFTs = np.zeros(self.mask.shape[:2] + (self.NCalibFrames,), dtype='complex64')
        
        self.lockFocus = False
        self.lockActive = False
        self.logShifts = True
        self.lastAdjustment = 5        

    def _finish_calibration(self):
        if self.state != State.FINISHING_CALIBRATION:
            raise RuntimeError("this method should only be called when finishing the calibration; actual state is %s" % self.state)
        
        # calculate the gradient info (needed in compare calls) from a valid calImages stack
        self.dz = np.gradient(self.calImages)[2].reshape(-1, self.NCalibFrames)
        self.dzn = np.hstack([1./np.dot(self.dz[:,i], self.dz[:,i]) for i in range(self.NCalibFrames)])
        
    def compare(self, frame_data):
        d = 1.0*frame_data.squeeze()
        dm = d/d.mean() - 1
        
        #where is the piezo suppposed to be
        #nomPos = self.piezo.GetPos(0)
        nomPos = self.piezo.GetTargetPos(0)
        
        #find closest calibration position
        posInd = np.argmin(np.abs(nomPos - self.calPositions))
        
        #retrieve calibration information at this location        
        calPos = self.calPositions[posInd]
        FA = self.calFTs[:,:,posInd]
        refA = self.calImages[:,:,posInd] 

        ddz = self.dz[:,posInd]
        dzn = self.dzn[posInd]
        
        #what is the offset between our target position and the calibration position         
        posDelta = nomPos - calPos
        
        #print('%s' % [nomPos, posInd, calPos, posDelta])
        
        #find x-y drift
        C = ifftshift(np.abs(ifftn(fftn(dm)*FA)))
        
        Cm = C.max()    
        
        Cp = np.maximum(C - 0.5*Cm, 0)
        Cpsum = Cp.sum()
        
        dx = (self.X*Cp).sum()/Cpsum
        dy = (self.Y*Cp).sum()/Cpsum
        
        ds = ndimage.shift(dm, [-dx, -dy])*self.mask
        
        #print A.shape, As.shape
        
        self.ds_A = (ds - refA)
        
        #calculate z offset between actual position and calibration position
        dz = self.Zfactor*self.deltaZ*np.dot(self.ds_A.ravel(), ddz)*dzn
        
        #posInd += np.round(dz / self.deltaZ)
        #posInd = int(np.clip(posInd, 0, self.NCalibStates))
            
        #add the offset back to determine how far we are from the target position
        dz = dz - posDelta
        
        return dx, dy, dz, Cm, nomPos, posInd, calPos, posDelta
        
    def compare_log_and_correct(self,frameData):
        # compare returns pixel coordinates for dx, dy, um for dz
        dx, dy, dz, cCoeff, nomPos, posInd, calPos, posDelta = self.compare(frameData)
        self.corrRef = max(self.corrRef, cCoeff) # keep track of historically maximal correlation amplitude
            
        #print dx, dy, dz
        dx_nm, dy_nm, dz_nm = (self.conversion['x']*dx, self.conversion['y']*dy, self.conversion['z']*dz)

        offset = self.piezo.GetOffset()
        offset_nm = 1e3*offset

        pos_um = self.piezo.GetPos(0)

        #FIXME: logging shouldn't call piezo.GetOffset() etc ... for performance reasons
        #       (is this still true, we keep the values cached in memory??)
        # this is the local logging, not to the actual localisation data acquiring instance of PYMEAcquire
        self.history.append((time.time(), dx_nm, dy_nm, dz_nm, cCoeff, self.corrRef, offset_nm, pos_um))
        eventLog.logEvent('PYME2ShiftMeasure', '%3.1f, %3.1f, %3.1f' % (dx_nm, dy_nm, dz_nm))
            
        self.lockActive = self.lockFocus and (cCoeff > .5*self.corrRef) # we release the lock when the correlation becomes too weak
        if self.lockActive:
            if abs(offset) > self._maxTotalCorrection:
                self.lockFocus = False
                logger.info("focus lock released, maximal Offset value exceeded (%.1f um)" % self._maxTotalCorrection)
            if abs(dz) > self.focusTolerance and self.lastAdjustment >= self.minDelay:
                # this sets the correction on the connected piezo
                self.piezo.SetOffset(offset - dz)                    
                self.lastAdjustment = 0
            else:
                self.lastAdjustment += 1
            
            if self.logShifts:
                # this logs to the connected copy of PYMEAcquire via the RESTServer
                if hasattr(self.piezo, 'LogShiftsCorrelAmp'):
                    self.piezo.LogShiftsCorrelAmp(dx_nm, dy_nm, dz_nm, self.lockActive, coramp=cCoeff/self.corrRef)
                else:
                    self.piezo.LogShifts(dx_nm, dy_nm, dz_nm, self.lockActive)

    def tick(self, frameData = None, **kwargs):
        if frameData is None:
            raise ValueError('frameData must be specified')
        else:
            frameData = self._crop_frame(frameData)
        
        targetZ = self.piezo.GetTargetPos(0)
        
        if not 'mask' in dir(self) or not self.frame_source.shape[:2] == self.mask.shape[:2]:
            # this just sets the UNCALIBRATED state and leaves the rest to the _prepare_calibration call
            self.state = State.UNCALIBRATED
            
        #called on a new frame becoming available
        if self.state == State.UNCALIBRATED:
            #print "cal init"

            #redefine our positions for the calibration
            self.homePos = self.piezo.GetPos(0)
            self._prepare_calibration(frameData)
            self.calibCurFrame = 0
            self.skipcounter = self.skipframes
            # move to our first position in the calib stack
            self.piezo.MoveTo(0, self.calPositions[0])
            
            # preps done, now switch state to calibrating
            self.state = State.CALIBRATING
        elif self.state == State.CALIBRATING:
            # print "cal proceed"
            if self.skipcounter >= 1:
                self.skipcounter -= 1 # we let the last move settle...
            else:
                # piezo step completed - record current image and move on to next position
                self._setRefN(frameData, int(self.calibCurFrame))
                if self.calibCurFrame == self.NCalibFrames-1: # we are mostly done, this was our last plane
                    self.state = State.FINISHING_CALIBRATION
                else: # not done yet
                    self.skipcounter = self.skipframes # reset skip counter
                    self.calibCurFrame += 1            # and go to next frame position
                    self.piezo.MoveTo(0, self.calPositions[int(self.calibCurFrame)])                
                
        elif self.state == State.FINISHING_CALIBRATION:
            # print "cal finishing"
            self._finish_calibration()
            self.piezo.MoveTo(0, self.homePos) # move back to where we started
            
            #reset our history log
            self.history = []
            self.historyColNames = ['time','dx_nm','dy_nm','dz_nm','corrAmplitude','corrAmpMax','piezoOffset_nm','piezoPos_um']
            self.historyStartTime = time.time()
            
            self.state = State.CALIBRATED # now we are fully calibrated
            
        elif self.state == State.CALIBRATED:
            # print "fully calibrated"
            if np.allclose(self._last_target_z, targetZ): # check we are on target in z
                self.compare_log_and_correct(frameData)

        else:
            raise RuntimeError("unknown calibration state %s, giving up" % self.state)
        
        self._last_target_z = targetZ                    

    def reCalibrate(self):
        self.state = State.UNCALIBRATED
        self.corrRef = 0
        self.lockActive = False
        
    def register(self):
        self.frame_source.connect(self.tick)
        self.tracking = True
        
    def deregister(self):
        self.frame_source.disconnect(self.tick)
        self.tracking = False
