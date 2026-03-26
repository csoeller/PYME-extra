import numpy as np
from scipy import optimize

# lightly edited version of same function in
# https://github.com/ahansenlab/chromatin_dynamics/blob/main/02_data_processing/01_MINFLUX/loading.py
def find_timestep(data,use_clumpIndex=False):
    tim = data['tim']
    if use_clumpIndex:
        tid = data['clumpIndex']
    else:
        tid = data['tid']

    dt = np.diff(tim)[np.diff(tid) == 0] # exclude steps across different trajectories
    if np.any(dt == 0):
        raise ValueError("Found zero time lag!")

    # Numerics might be better if everything is O(1)
    scale = np.min(dt)
    dt = dt / scale
    mindt = np.min(dt)
    
    # Step 1: rough estimate through MSE
    def mse(step):
        ints = np.round(dt/step).astype(int)
        return np.sum((dt-step*ints)**2)

    res = optimize.minimize(mse, mindt,
                            bounds=[(mindt, np.inf)],
                           )
    if not res.success:
        print(res)
        raise RuntimeError

    step = res.x

    # Step 2: identify real integer steps
    udts = []
    Ns = []
    cur = 0.5*step
    while cur < np.max(dt):
        ind = (dt > cur) & (dt < cur+step)
        Ns.append(np.sum(ind))
        if Ns[-1] > 0:
            udts.append(np.mean(dt[ind]))
            cur = udts[-1] + 0.5*step
        else:
            udts.append(np.nan)
            cur += step
    udts = np.array(udts)
    Ns   = np.array(Ns)

    # Step 3: fit actual best lag time
    ind = ~np.isnan(udts)
    with np.errstate(divide='ignore'):
        sigma = 1/np.sqrt(Ns[ind]-1)
    res = optimize.curve_fit(lambda x, a: a*x,
                             np.arange(len(udts))[ind]+1,
                             udts[ind],
                             sigma=sigma,
                            )

    return res[0][0]*scale

def gen_tconsensus(data,use_clumpIndex=False,scaleDT=1.0):
    dtc = scaleDT * float(find_timestep(data,use_clumpIndex=use_clumpIndex))
  
    tim = data['tim']
    if use_clumpIndex:
        tid = data['clumpIndex']
    else:
        tid = data['tid']    
    utid = np.unique(tid)

    tt_ms = np.zeros_like(tim)
    tdt_ms = np.zeros_like(tim)
    
    for my_tid in utid:
        ind = tid == my_tid
        my_tim = tim[ind]
        
        lag = np.round(np.diff(my_tim)/dtc).astype(int) * dtc
        t_trace_ms = 1e3*np.insert(np.cumsum(lag), 0, 0)
        t_trace_dtms = 1e3*np.insert(lag, 0, 0)
        tt_ms[ind] = t_trace_ms
        tdt_ms[ind] = t_trace_dtms

    return tt_ms, tdt_ms
