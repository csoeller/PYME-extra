import wx
from  wx.lib.dialogs import ScrolledMessageDialog

def isDarwin():
   import os
   from sys import platform
   return platform == "darwin"

def Warn(parent, message, caption = 'Warning!'):
    dlg = wx.MessageDialog(parent, message, caption, wx.OK | wx.ICON_WARNING)
    dlg.ShowModal()
    dlg.Destroy()

class ShowErr:
   """

   """
   def __init__(self, dsviewer):
      self.dsviewer = dsviewer
      self.txtwin = None
      dsviewer.AddMenuItem('Experimental', 'Errors in scrolled message dialog\tCtrl+E', self.OnErrScrolledDialog)

   def getLogFileName(self):
      import os
      pid = int(os.getpid()) - 1
      return os.path.join('/tmp','dh5view-%d.tmp' % pid)

   def OnErrScrolledDialog(self, event=None):
      if not isDarwin():
         Warn(None,'aborting: functionality only available on mac','Error')
         return
      with open(self.getLogFileName(),"r") as f:
         txt = "\n".join(f.readlines())
      dlg = ScrolledMessageDialog(self.dsviewer, txt, "VisGUI Error Output", size=(900,400),
                                  style=wx.RESIZE_BORDER | wx.DEFAULT_DIALOG_STYLE )
      dlg.ShowModal()
      dlg.Destroy()

def Plug(dsviewer):
    """Plugs this module into the gui"""
    dsviewer.showErr = ShowErr(dsviewer)
