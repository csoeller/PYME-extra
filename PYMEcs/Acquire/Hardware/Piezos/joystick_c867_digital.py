# a piezo_pipython_gcs compatible joystick class needs
# - Enable() and IsEnabled() methods for use by PYME position ui and scope object
# - an init() method to set the joystick parameters via GCS commands and register the gcspiezo "parent" instance
# - an enablecommands() method of GCS commands that will be used by the gcspiezo "parent" object to enable the joystick
# - a disablecommands() method of GCS commands that will be used by the gcspiezo "parent" object to disable the joystick
# example usage:
#    from PYME.Acquire.Hardware.Piezos.piezo_pipython_gcs import GCSPiezoThreaded
#    from PYMEcs.Acquire.Hardware.Piezos.joystick_c867_digital import digitalJoystick
#    scope.stage = GCSPiezoThreaded('PI C-867 Piezomotor Controller SN 0122013807', axes=['1', '2'],
#                                   refmodes='FRF',joystick=digitalJoystick())
class digitalJoystick:
    def __init__(self):
        self.gcspiezo = None
        self._initialised = False

    # this method needs to be tweaked to the specific controller and joystick model
    def init(self,gcspiezo):
        self.gcspiezo = gcspiezo
        # associate controller axis 1 (1st arg) using velocity control mode (2nd arg, value 3)
        # with HID device_ID 2 (third arg) using its 'Axis_1' (4th arg)
        self.gcspiezo.pi.gcscommands.HIA(1,3,2,'Axis_1')
        # equivalent for 2nd axis
        self.gcspiezo.pi.gcscommands.HIA(2,3,2,'Axis_2')
        # set closed loop velocity to a certain value (e.g. 1 meaning 1 mm/s or what?)
        self.gcspiezo.pi.gcscommands.VEL(self.gcspiezo.axes,1.0)

        self._initialised = True
        
    def Enable(self, enabled = True):
        self.check_initialised()
        if not self.IsEnabled() == enabled:
            if enabled:
                self.gcspiezo.enable_joystick()
            else:
                self.gcspiezo.disable_joystick()

    def IsEnabled(self):
        self.check_initialised()
        return self.gcspiezo._joystick_enabled

    # GCS commands to enable the joystick
    def enablecommands(self,pidevice): # this method should only be used from the parent gcspiezo object
        # possible we should use self.gcspiezo.axes here, e.g. pidevice.gcscommands.HIN(self.gcspiezo.axes,True)
        pidevice.gcscommands.HIN(self.gcspiezo.axes,True)

    # GCS commands to enable the joystick
    def disablecommands(self,pidevice): # this method should only be used from the parent gcspiezo object
        # possible we should use self.gcspiezo.axes here
        pidevice.gcscommands.HIN(self.gcspiezo.axes,False)

    def check_initialised(self):
        if not self._initialised:
            raise RuntimeError("joystick must have been initialised before using these methods")
        
