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

   def getLogFileName(self, pidoffset = -1):
      import os
      pid = int(os.getpid()) + pidoffset
      return os.path.join('/tmp','dh5view-%d.tmp' % pid)

   def OnErrScrolledDialog(self, event=None):
      if not isDarwin():
         Warn(None,'aborting: functionality only available on mac','Error')
         return
      ok = False
      import os.path
      fname = self.getLogFileName(pidoffset = -1)
      if os.path.isfile(fname):
         ok = True
      else:
         fname = self.getLogFileName(pidoffset = 0)
         if os.path.isfile(fname):
            ok = True
         else:
            Warn(None,'aborting: cannot find log file, tried %s and %s' % (self.getLogFileName(pidoffset = -1),
                                                                           self.getLogFileName(pidoffset = 0)))
      if ok:
         with open(fname,"r") as f:
            txt = "\n".join(f.readlines())
            dlg = ScrolledMessageDialog(self.dsviewer, txt, "VisGUI Error Output", size=(900,400),
                                        style=wx.RESIZE_BORDER | wx.DEFAULT_DIALOG_STYLE )
            dlg.ShowModal()
            dlg.Destroy()

def Plug(dsviewer):
    """Plugs this module into the gui"""
    dsviewer.showErr = ShowErr(dsviewer)
