import PYME.version as PYMEver
from packaging import version as pv

Features = {
    'PluginClass' : {'doc':'Plugin base class that implements weak refs for plugins',
                     'minVersion' : '20.07.08'}
    }


def PYMEversionCheck(feature=None):
    if feature is None:
        return True
    if feature in Features:
        if pv.parse(PYMEver.version) >= pv.parse(Features[feature]['minVersion']):
            return True
        else:
            raise RuntimeError("PYME upgrade required! We need the PYME '%s' feature which is available since PYME version %s, you have version %s" %
                               (feature,Features[feature]['minVersion'],PYMEver.version))
    else:
        raise ValueError("checking for unknown feature '%s'" % feature)
