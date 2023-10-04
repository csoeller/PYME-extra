from PYMEcs.Analysis.fitDarkTimes import TabularRecArrayWrap
import pandas as pd
import numpy as np

# we should probably add this to a module like PYMEcs.IO.tabular
def tabularFromCSV(csvname):
    df = pd.DataFrame.from_csv(csvname)
    rec = df.to_records()
    tb = TabularRecArrayWrap(rec)
    return tb

from PYME.IO.tabular import TabularBase
from PYMEcs.IO.MINFLUX import minflux_npy2pyme
class MinfluxNpySource(TabularBase):
    _name = "MINFLUX NPY File Source"
    def __init__(self, filename):
        """ Input filter for use with NPY data exported from MINFLUX data (typically residing in MSR files)."""


        self.res = minflux_npy2pyme(filename)

        # check for invalid localisations:
        # possible TODO - is this needed/helpful, or should we propagate missing values further?
        # FIXME - this does not seem to work properly yet when invalid data is present
        if np.any(self.res['vld'] < 1):
            self.res = self.res[self.res['vld'] >= 1]

        self.res = self.res.to_records(index=False) # convert into NUMPY recarray
        self._keys = list(self.res.dtype.names)
       

    def keys(self):
        return self._keys

    def __getitem__(self, keys):
        key, sl = self._getKeySlice(keys)
        
        if not key in self._keys:
            raise KeyError('Key (%s) not found' % key)

       
        return self.res[key][sl]

    
    def getInfo(self):
        return 'MINFLUX NPY Data Source\n\n %d points' % len(self.res['x'])
