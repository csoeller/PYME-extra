# a few routines to manipulate recipes in yaml format
# this is supposed to make programmatic changing of parameters a little easier

# note that we work on the python data structure which is a mix of dicts and lists
# as it comes in from the yaml reader routines

# we initially implement a minimal set of routines with the goal to uniquely identify
# a module in the recipe by module name and the name of inputs and outputs
# if a unique module was identified then the keyword parameters are used to set values
# in the identified module entry. Currently, this is simplistic and mostly for scalar values.
# TODO: For further tweaking we may need to make the modParams class aware of more complex data structures
# like lists and dicts for certain modules...



import yaml
import warnings

def modname(entry):
        return list(entry.keys())[0]

def modinputs(entry):
    mname = modname(entry)
    inputnames = [entry[mname][key] for key in entry[mname].keys() if key.startswith('input')]
    return inputnames

def modoutputs(entry):
    mname = modname(entry)
    outputnames = [entry[mname][key] for key in entry[mname].keys() if key.startswith('output')]
    return outputnames

def modmatch(entry,name,inputDS,outputDS):
    if modname(entry) != name:
        return False
    if inputDS is not None:
        if inputDS not in modinputs(entry):
            return False
    if outputDS is not None:
        if outputDS not in modoutputs(entry):
            return False
    return True

class modParams(object):
    """A class to modify parameters of identified module entries in a recipe

    Typical usage:

    modified_rec = load_recipe('base_recipe.yaml',
                               recmy.modParams('PYMEcs.NNfilter',
                                               nnMin=10.0, nnMax=30.0),
                               recmy.modParams('PYMEcs.DBSCANClustering2',
                                               inputDS='coalesced_nz',searchRadius=7.0)))
    """
    
    def __init__(self,name,inputDS=None,outputDS=None,**kwargs):
        self.name = name
        self.inputDS = inputDS
        self.outputDS = outputDS
        self.paramdict = kwargs
        
    def modify_recipe(self,rec,inplace=False):        
        if not inplace:
            rec = rec.copy()
        matches = [entry for entry in rec if modmatch(entry,self.name,self.inputDS,self.outputDS)]
        if len(matches) == 0:
            warnings.warn("no matching module entry, ignoring parameters")
        elif len(matches) > 1:
            raise RuntimeError("more than one matched entry, need exactly one match")
        else:
            entry = matches[0]
            for key in self.paramdict:
                entry[self.name][key] = self.paramdict[key]
        return rec

def load_recipe(filename, *modpars, return_yamldumped=True):
    with open(filename) as fi:
        rec = yaml.safe_load(fi)
    for mp in modpars:
        rec = mp.modify_recipe(rec)
    if return_yamldumped:
        return yaml.dump(rec)
    else:
        return rec
