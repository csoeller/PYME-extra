import requests
from urllib.parse import urljoin
import logging
from http import HTTPStatus
from PYMEcs.Acquire.Hardware.NikonTiSim import LightPath as TiLightPath

logger = logging.getLogger(__name__)

# by basing it on the simulator we inherit the "generic" methods, i.e.
# ProvideMetadata, OnChange and Poll
# which are required for the functional object when used with
# PYMEAcquire
class LPClient(TiLightPath):
    def __init__(self, url="http://127.0.0.1:5000", timeout=15):
        self.url = url
        self.timeout = timeout
        self.wantChangeNotification = []

        self.names = self.GetNames() # execute the method once and then keep the values cached
        self.lastPosition = self.GetPosition() # initialise here as later used in Poll method

    def GetPosition(self):
        response = requests.get(urljoin(self.url,'position'),timeout=self.timeout)
        if response.status_code == HTTPStatus.OK:
            position = int(response.content.decode())
            return position
        else:
            logger.warn("HTTP response error: got status %d" % response.status_code)
            return None

    def SetPosition(self,position):
        response = requests.put(urljoin(self.url,'position'),data=str(position),timeout=self.timeout)
        if response.status_code == HTTPStatus.OK:
            newpos = int(response.content.decode())
            if newpos != position:
                logger.warn("asked for position '%d' but got position '%s'" % (newpos,position))
            self.lastPosition = newpos # update so that OnChange works ok
            self.OnChange()
        else:
            logger.warn("HTTP response error: got status %d" % response.status_code)

    def GetPort(self):
        response = requests.get(urljoin(self.url,'port'),timeout=self.timeout)
        if response.status_code == HTTPStatus.OK:
            port = response.content.decode()
            return port
        else:
            logger.warn("HTTP response error: got status %d" % response.status_code)
            return None

    def SetPort(self,port):
        if port in self.names:
            response = requests.put(urljoin(self.url,'port'),data=port,timeout=self.timeout)
            if response.status_code == HTTPStatus.OK:
                newport = response.content.decode()
                if newport != port:
                    logger.warn("asked for port '%s' but got port '%s'" % (newport,port))
                self.lastPosition = self.GetPosition() # update so that OnChange works ok
                self.OnChange()
            else:
                logger.warn("HTTP response error: got status %d" % response.status_code)
        else:
            logger.warn("Asking for unknown port '%s', ignoring" % port)

    def GetNames(self):
        response = requests.get(urljoin(self.url,"names"),timeout=self.timeout)
        if response.status_code == HTTPStatus.OK:
            names = response.json()
            return names
        else:
            logger.warn("HTTP response error: got status %d" % response.status_code)
            return None

