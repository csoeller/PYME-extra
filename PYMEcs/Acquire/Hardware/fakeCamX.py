from PYME.Acquire.Hardware.Simulator.fakeCam import FakeCamera
from PYME.IO import MetaDataHandler

# the only difference of this class is to supply extra metadata to allow testing the
# camera map processing code
class FakeCameraX(FakeCamera):
    
    def GenStartMetadata(self, mdh):
        super().GenStartMetadata(mdh)
        # these are useful when we want to test camera map making code
        mdh.setEntry('Camera.SerialNumber', 'FAKE-007')
        mdh.setEntry('Camera.SensorWidth',self.GetCCDWidth())
        mdh.setEntry('Camera.SensorHeight',self.GetCCDHeight())
        mdh.setEntry('Camera.Model', 'FakeCameraX')

