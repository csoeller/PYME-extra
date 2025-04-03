from PYME.warnings import warn
import pickle
from pathlib import Path

def load_NPC_set(fname, ts=None, foreshortening=None, version=None):
    with open(fname,'rb') as fi:
        npcs=pickle.load(fi)
    fpath = Path(fname)
    if ts is not None and ts not in fpath.stem:
        warn("MINFLUX time stamp not in NPC dataset filename; check correct filename chosen: %s vs %s" %
             (ts,fpath.stem))

    if foreshortening is not None:
        try:
            npc_foreshortening = npcs.foreshortening
        except AttributeError:
            npc_foreshortening = 1.0

        if abs(npc_foreshortening-foreshortening) >= 0.01:
            warn("NPC foreshortening is %.2f while dataset foreshortening is %.2f, check this is compatible" %
                 (npc_foreshortening,foreshortening))
        

    if version is not None:
        pass # here add version checking logic

    return npcs

def save_NPC_set(npcs,fname):
    with open(fname, "wb") as file:
            pickle.dump(npcs,file)

def findNPCset(pipeline=None):
    if pipeline is None:
        return None
    if 'npcs' not in dir(pipeline) or pipeline.npcs is None:
        return None
    return pipeline.npcs
