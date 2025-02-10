from PYMEcs.Analysis.fitDarkTimes import TabularRecArrayWrap
import pandas as pd
import numpy as np

# we should probably add this to a module like PYMEcs.IO.tabular
def tabularFromCSV(csvname):
    df = pd.DataFrame.from_csv(csvname)
    rec = df.to_records()
    tb = TabularRecArrayWrap(rec)
    return tb

