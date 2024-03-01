
# implementation using the JoystickBase base class
from PYME.Acquire.Hardware.Piezos.piezo_pipython_gcs import JoystickBase
class DigitalJoystick(JoystickBase):
    def __init__(self):
        super().__init__()

    # this method needs to be tweaked to the specific controller and joystick model
    def init(self,gcspiezo):
        # register gcspiezo - this needs to always happen
        self.gcspiezo = gcspiezo

        # configure joystick axes and velocity
        #
        # associate controller axis 1 (1st arg) using velocity control mode (2nd arg, value 3)
        # with HID device_ID 2 (third arg) using its 1st axis ('Axis_1') (1: 4th arg)
        # NOTE: only the HID axis number (e.g. 1) is accepted as argument,
        #       NOT the name returned by qHIS (e.g. 'Axis_1')
        self.gcspiezo.pi.HIA(1,3,2,1)
        # associate controller axis 2 (1st arg) using velocity control mode (2nd arg, value 3)
        # with HID device_ID 2 (third arg) using its 2nd axis ('Axis_2') (2: 4th arg)
        self.gcspiezo.pi.HIA(2,3,2,2)
        # set closed loop velocity to a reasonable but  smallish value (e.g. here 1.5 apparently means 1.5 mm/s)
        self.gcspiezo.pi.VEL(self.gcspiezo.axes,[1.5 for axis in self.gcspiezo.axes])

        self._initialised = True
 
    # GCS commands to enable the joystick
    def enablecommands(self): # this method should only be used from the parent gcspiezo object
        # testing showed that we need to give a list of booleans, one for each axis
        self.gcspiezo.pi.HIN(self.gcspiezo.axes,[True for axis in self.gcspiezo.axes])

    # GCS commands to enable the joystick
    def disablecommands(self): # this method should only be used from the parent gcspiezo object
        # testing showed that we need to give a list of booleans, one for each axis
        self.gcspiezo.pi.HIN(self.gcspiezo.axes,[False for axis in self.gcspiezo.axes])
