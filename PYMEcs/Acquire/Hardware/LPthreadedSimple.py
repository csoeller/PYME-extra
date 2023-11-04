from threading import Thread, Lock
from queue import Queue
from PYME.Acquire.Hardware.NikonTi import LightPath
from http import HTTPStatus

class commandObject(object):
    def __init__(self,cmd,*args):
        self.cmd = cmd
        self.args = args

    def __repr__(self):
        return 'CMDMSG : %s!' % (self.cmd) +  \
            '[' + ','.join(str(x) for x in self.args) + ']'

class returnObject(object):
    def __init__(self,cmd,status,*args):
        self.cmd = cmd
        self.status = status
        self.return_args = args
        
    def __repr__(self):
        return 'RETMSG : %s!%s!' % (self.cmd,self.status) + \
            '[' + ','.join(str(x) for x in self.return_args) + ']'

# ...
# except com_error:
# ...

# note that we use HTTP status codes as we use a REST server anyway and
# thus having status values we can directly pass back via HTTP responses are useful

from pywintypes import com_error

class LPSafe(LightPath):
    def __init__(self):
        super().__init__()

    def GetNamesS(self):
        try:
            names = self.names
            status = HTTPStatus.OK
        except com_error:
            status = HTTPStatus.SERVICE_UNAVAILABLE
            names = None
        return status, names

    def GetPortS(self):
        try:
            port = self.GetPort()
            status = HTTPStatus.OK
        except com_error:
            status = HTTPStatus.SERVICE_UNAVAILABLE
            port = None

        return status, port

    def SetPortS(self,port):
        try:
            self.SetPort(port)
        except com_error:
            status = HTTPStatus.SERVICE_UNAVAILABLE
        else:
            status = HTTPStatus.OK

        return status

    def GetPositionS(self):
        try:
            pos = self.GetPosition()
            status = HTTPStatus.OK
        except com_error:
            status = HTTPStatus.SERVICE_UNAVAILABLE
            pos = None

        return status, pos

    def SetPositionS(self,pos):
        try:
            self.SetPosition(pos)
        except com_error:
            status = HTTPStatus.SERVICE_UNAVAILABLE
        else:
            status = HTTPStatus.OK

        return status
    

class LPThread(Thread):
    def __init__(self, *args, **kwargs):
        super().__init__()
        self.command_queue = Queue()
        self.results_queue = Queue()
        self.lp = LPSafe()
        
        self.daemon = True

    def run(self):
        while True:
            cmd = self.command_queue.get()
            result = self.parse_command(cmd)
            self.results_queue.put(result)

    def parse_command(self,cmd):
        if (not cmd.cmd) or (cmd.cmd is None):
            return returnObject(None, HTTPStatus.NO_CONTENT, 'no command received')
        elif cmd.cmd == 'GetPort':
            status, port = self.lp.GetPortS()
            return returnObject(cmd.cmd, status, port)
        elif cmd.cmd == 'SetPort':
            status = self.lp.SetPortS(cmd.args[0])
            return returnObject(cmd.cmd, status)       
        elif cmd.cmd == 'GetPosition':
            status, pos = self.lp.GetPositionS()
            return returnObject(cmd.cmd, status,pos)
        elif cmd.cmd == 'SetPosition':
            status = self.lp.SetPositionS(cmd.args[0])
            return returnObject(cmd.cmd, status)
        elif cmd.cmd == 'GetNames':
            status, names = self.lp.GetNamesS()
            return returnObject(cmd.cmd, status, names)
        else:
            return returnObject(cmd.cmd, HTTPStatus.NOT_IMPLEMENTED)

    def run_command(self,cmd,*args):
        self.command_queue.put(commandObject(cmd,*args))
        ret = self.results_queue.get()
        if len(ret.return_args) == 0:
            return ret.status
        elif len(ret.return_args) == 1:
            return ret.status, ret.return_args[0]
        else:
            return ret.status, ret.return_args

# test some aspects of the code here
def main():
    lpt = LPThread()
    lpt.start()

    status,port = lpt.run_command('GetPort')
    if status == HTTPStatus.OK:
        print('Port is %s' % port)

    status = lpt.run_command('SetPort','L100')

    status, port = lpt.run_command('GetPort')
    if status == HTTPStatus.OK:
        print('Port is %s' % port)

    status, pos = lpt.run_command('GetPosition')
    if status == HTTPStatus.OK:
        print('Position is %d' % pos)

    status, names = lpt.run_command('GetNames')
    if status == HTTPStatus.OK:
        print('Names:',names)

    status = lpt.run_command(None)   

    print(commandObject('GetIt',1,'L100'))
    print(returnObject('Received',HTTPStatus.OK,1,'L100'))

if __name__ == '__main__':
    main()
