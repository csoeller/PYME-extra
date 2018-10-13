from zeroconf import *
import PYME.misc.pyme_zeroconf as pzc
import time

import logging
logging.basicConfig(level=logging.INFO)

# logger = logging.getLogger(__name__)
# logger.setLevel(logging.DEBUG)

# # create console handler and set level to debug
# ch = logging.StreamHandler()
# ch.setLevel(logging.INFO)

# # add ch to logger
# logger.addHandler(ch)


TIMEOUTDEFAULT = 10
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
    def find(cls, zc=None, timeout=TIMEOUTDEFAULT):
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

def servicesPresent(timeOut=TIMEOUTDEFAULT, showServices=False):
    services = zt.find(timeout=timeOut)
    if showServices:
        print("Available Services: %s" % repr(services))

    return len(services) > 0

def checkServer(timeOut=TIMEOUTDEFAULT, showServices=False):
    if servicesPresent(timeOut=timeOut, showServices=showServices):
        logging.info('zeroconf services detected')
    else:
        logging.error('no zeroconf services detected - this should not happen')
        
    ns = pzc.getNS()
    adserv = get_advertised_services(ns)
    if len(adserv) > 0:
        logging.info(repr(adserv))
    else:
        logging.error('no advertised pyro services - apparently there is no server running on this network')
    
def get_advertised_services(ns):
    try:
        services = ns.get_advertised_services()
    except:
        services = ns.advertised_services

    return services
