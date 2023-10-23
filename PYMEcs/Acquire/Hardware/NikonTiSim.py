class LightPath:
    def __init__(self, names = ['EYE', 'L100', 'R100', 'L80']):
        self.simulatedPosition = 1
        self.names = names
        self.wantChangeNotification = []
        
        self.lastPosition = self.GetPosition()
        
    def SetPosition(self, pos):
        self.simulatedPosition = (pos + 1)
        self.lastPosition = pos
        self.OnChange()
        
    def GetPosition(self):
        return int(self.simulatedPosition) - 1
        
    def SetPort(self, port):
        self.SetPosition(self.names.index(port))
        self.OnChange()
        
    def GetPort(self):
        return self.names[self.GetPosition()]
    
    def ProvideMetadata(self,mdh):
        mdh.setEntry('NikonTi.LightPath', self.GetPort())
        
    def OnChange(self):
        for a in self.wantChangeNotification:
            a()
            
    def Poll(self):
        pos = self.GetPosition()
        if not self.lastPosition == pos:
            self.lastPosition = pos
            self.OnChange()
