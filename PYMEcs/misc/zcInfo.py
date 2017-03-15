from zeroconf import *
import PYME.misc.pyme_zeroconf as pzc
import time

import logging
logger = logging.getLogger(__name__)

class ZeroconfServiceTypes(object):
    """
    Return all of the advertised services on any local networks
    """
    def __init__(self):
        self.found_services = set()

    def add_service(self, zc, type_, name):
        self.found_services.add(name)

    def remove_service(self, zc, type_, name):
        pass

    @classmethod
    def find(cls, zc=None, timeout=5):
        """
        Return all of the advertised services on any local networks.

        :param zc: Zeroconf() instance.  Pass in if already have an
                instance running or if non-default interfaces are needed
        :param timeout: seconds to wait for any responses
        :return: tuple of service type strings
        """
        local_zc = zc or Zeroconf()
        listener = cls()
        browser = ServiceBrowser(
            local_zc, '_services._dns-sd._udp.local.', listener=listener)

        # wait for responses
        time.sleep(timeout)

        # close down anything we opened
        if zc is None:
            local_zc.close()
        else:
            browser.cancel()

        return tuple(sorted(listener.found_services))

try:
    zt = zeroconf.ZeroconfServiceTypes()
except:
    zt = ZeroconfServiceTypes()

def servicesPresent(timeOut=5):
    services = zt.find(timeout=timeOut)
    print(services)

    return len(services) > 0

def checkServer():
    if servicesPresent():
        logger.info('zeroconf services detected')
    else:
        logger.error('no zeroconf services detected - this should not happen')
        
    ns = pzc.getNS()
    if len(ns.advertised_services) > 0:
        logger.info(repr(ns.advertised_services))
    else:
        logger.error('no advertised pyro services - apparently there is no server running on this network')
    
