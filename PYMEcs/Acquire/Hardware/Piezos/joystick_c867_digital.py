# a gcspiezo compatible joystick class needs
# - Enable() and IsEnabled() methods for use by position ui and scope object
# - an init() method to set the joystick parameters via GCS commands and register the gcspiezo "parent" object
# - an enablecommands() method that will be used by the gcspiezo "parent" object to enable the joystick
# - a disablecommands() method that will be used by the gcspiezo "parent" object to disable the joystick

class digitalJoystick:
    def __init__(self):
        self.gcspiezo = None
        self._initialised = False
        
    def init(self,gcspiezo):
        self.gcspiezo = gcspiezo
        # associate controller axis 1 (1st arg) using velocity control mode (2nd arg, value 3)
        # with HID device_ID 2 (third arg) using its 'Axis_1' (4th arg)
        self.gcspiezo.pi.gcscommands.HIA(1,3,2,'Axis_1')
        # equivalent for 2nd axis
        self.gcspiezo.pi.gcscommands.HIA(2,3,2,'Axis_2')
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

    def enablecommands(self,pidevice): # these methods should only be used from the parent gcspiezo object
        # possible we should use self.gcspiezo.axes here
        pidevice.gcscommands.HIN(['1','2'],True)
        
    def disablecommands(self,pidevice): # these methods should only be used from the parent gcspiezo object
        # possible we should use self.gcspiezo.axes here
        pidevice.gcscommands.HIN(['1','2'],False)

    def check_initialised(self):
        if not self._initialised:
            raise RuntimeError("joystick must have been initialised before using these methods")
        
