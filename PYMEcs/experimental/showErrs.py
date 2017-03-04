import wx

# this is a quick hack
# it needs to be properly wrapped in a frame etc
# for now a convenient but ugly way to see error messages on OS X

class WindowPopup(wx.PopupWindow):
   """ Pops up a window to provide description for the selection """
   def __init__(self, parent, style, content):
      wx.PopupWindow.__init__(self, parent, style)

      sizer = wx.BoxSizer(wx.VERTICAL)
      self.st = wx.TextCtrl(self, -1, style = wx.TE_MULTILINE | wx.TE_READONLY, size = (900, 400))
      self.st.SetValue(content)
      self.SetSize((900, 400))
      sizer.Add(self.st, 0, wx.EXPAND)
      self.SetSizer(sizer)
      self.Layout()



class ShowErr:
   """

   """
   def __init__(self, visFr):
      self.visFr = visFr
      self.txtwin = None
      visFr.AddMenuItem('Experimental', 'Toggle Error Log Display\tCtrl+E', self.OnToggleErrVisibility)


   def OnToggleErrVisibility(self, event=None):
      import os

      if self.txtwin is None:
         self.txtwin = WindowPopup(self.visFr,wx.SIMPLE_BORDER, '')
      
      if self.txtwin.IsShown():
         self.txtwin.Hide()
      else:
         pid = int(os.getpid()) - 1
         with open(os.path.join('/tmp','visgui-%d.tmp' % pid),"r") as f:
            self.txt = "\n".join(f.readlines())
         self.txtwin.st.SetValue(self.txt)
         self.txtwin.Layout()
         self.txtwin.Show()

def Plug(visFr):
    """Plugs this module into the gui"""
    visFr.showErr = ShowErr(visFr)
