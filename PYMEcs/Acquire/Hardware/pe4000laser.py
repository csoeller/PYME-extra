from PYME.Acquire.Hardware.lasers import Laser
from PYMEcs.Acquire.Hardware.coolled_pe4000 import PE4000

class pe4000laser(Laser):
   def __init__(self, name,turnOn=False, portname='COM3', maxpower=1.0, power_fudge=1.0,**kwargs):
      self.device = PE4000(portname)
      self.device.load_wavelength(560)
      self.name = name
      #self.ser_port = serial.Serial(portname, 500000, timeout=.1, writeTimeout=2)
      self.powerControlable = True
      self.power = 0.02

   def IsOn(self):
      return self.isOn

   def TurnOn(self):
      self.device.set_channel('C', intensity=100.0*self.power, on=True)
      self.isOn = True
      
   def TurnOff(self):
      self.device.all_off()
      self.isOn = False

   def SetPower(self, power):
      if power < 0 or power > 1:
         raise RuntimeError('Error setting laser power: Power must be between 0 and 1')
      self.power = power

   def GetPower(self):
      return self.power
