from PYME.Acquire.Hardware.Piezos import offsetPiezoREST

from PYME.Acquire import eventLog
import time

from PYME.util import webframework

import logging
logger = logging.getLogger(__name__)

class OffsetPiezoCorrelLog(offsetPiezoREST.OffsetPiezo):

    @webframework.register_endpoint('/LogShiftsCA', output_is_json=False)
    def LogShiftsCA(self, dx, dy, dz, active=True, coramp=-1.0):
        import wx
        #eventLog.logEvent('ShiftMeasure', '%3.4f, %3.4f, %3.4f' % (dx, dy, dz))
        wx.CallAfter(eventLog.logEvent, 'ShiftMeasure', '%3.4f, %3.4f, %3.4f' % (float(dx), float(dy), float(dz)), time.time())
        wx.CallAfter(eventLog.logEvent, 'PiezoOffset', '%3.4f, %d' % (self.GetOffset(), int(active)), time.time())
        wx.CallAfter(eventLog.logEvent, 'CorrelationAmplitude', '%3.4f' % (float(coramp)), time.time())

class OffsetPiezoCorrelLogClient(offsetPiezoREST.OffsetPiezoClient):

    def LogShiftsCorrelAmp(self, dx, dy, dz, active=True, coramp=-1.0):
        res = self._session.get(self.urlbase +
                                '/LogShiftsCA?dx=%3.3f&dy=%3.3f&dz=%3.3f&active=%d&coramp=%3.3f' % (dx, dy, dz, active,coramp))


def getClient():
    #TODO - move away from hard-coded ports!!!
    return OffsetPiezoCorrelLogClient()

def getServer():
    return offsetPiezoREST.server_class(offset_piezo_base_class=OffsetPiezoCorrelLog)
