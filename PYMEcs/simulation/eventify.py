from PYME.simulation.locify import locify, eventify2
import numpy as np

# we essentially copy the version of eventify from david's code but modify the way that intensities are modelled
# it may be be ok to draw these from an exponential but here we used a Poisson

fresultdtype=[('tIndex', '<i4'),
              ('fitResults', [('A', '<f4'),('x0', '<f4'),('y0', '<f4'), ('z0', '<f4'),('sigma', '<f4')]),
              ('fitError', [('A', '<f4'),('x0', '<f4'),('y0', '<f4'), ('z0', '<f4'),('sigma', '<f4')])]


_s2 = 110.**2
_a2 = 70.**2

_s2_a2_12 = _s2 + _a2/12
_8_pi_s2_2_a2 = (8*np.pi*_s2**2)/_a2
_r_2_pi = 1.0/2*np.pi

def eventify3(x, y, meanIntensity, meanDuration, backGroundIntensity, meanEventNumber, sf=2, tm=10000, z=0,
             z_err_scale=1.0, paint_mode=True):
    """ PAINT version of eventify """
    
    
    #Is =
    Ns = np.random.poisson(meanEventNumber, x.shape)

    if np.isscalar(z):
        z = z * np.ones_like(x)

    evts = []
    #t = 0

    for x_i, y_i, z_i, N_i in zip(x, y, z, Ns):
        duration = np.random.exponential(meanDuration, size=N_i)
        n_frames = np.ceil(duration).astype('i')
        
        if paint_mode:
            # David's code uses np.random.exponential here!
            # check what we would assume here on the basis of realistic description of the process
            Is = np.random.poisson(meanIntensity, size=N_i)
            # possibly we should do something different still:
            #  for a bust draw a single burstintensity from an exponential dist
            #  then WITHIN the burst vary the frame intensity around the burstintensity via Poisson
            ts = 2*tm*np.random.uniform(size=N_i)
        else:
            Is = np.random.exponential(meanIntensity)*np.ones(N_i)
            ts = np.random.exponential(tm, size=N_i)
        
        evts_i = np.empty(n_frames.sum(), dtype=fresultdtype)
        
        evts_i['fitResults']['x0'] = x_i
        evts_i['fitResults']['y0'] = y_i
        evts_i['fitResults']['z0'] = z_i
        
        k = 0
        
        for j in range(N_i):
            k2 = k + n_frames[j]
            evts_i['fitResults']['A'][k:k2] = Is[j]
            evts_i['fitResults']['A'][k2-1] = Is[j]*(duration[j] %1)
            evts_i['tIndex'][k:k2] = ts[j]+np.arange(n_frames[j])
            
            k = k2
            
            #n_frames = int(np.ceil(duration[j]))
            #t = 2*tm*np.random.uniform()
            
            #I_i = np.random.exponential(meanIntensity)
            
            #t_i = ts[j]+np.arange(n_frames[j])
            
            #I_i = Is[j]*np.ones(n_frames[j])
            #I_i[-1] = I_i[-1]*(duration[j] %1)
        
            #evts += [(x_i, y_i, I_i, t+k) for k in range(duration)] + [(x_i, y_i, I_i*(duration%1), t+floor(duration))]
            #evts.extend([FitResultR(x_i, y_i, z_i, I_i, t + k, backGroundIntensity, z_err_mult=z_err_scale) for k in
            #             range(int(np.floor(duration)))])
            #evts.append(FitResultR(x_i, y_i, z_i, I_i, t_i, backGroundIntensity, z_err_mult=z_err_scale))
        evts.append(evts_i)

    evts = np.hstack(evts)

    #xn, yn, In = evts[:,0], evts[:,1], evts[:,2]
    
    

    In = evts['fitResults']['A']*_r_2_pi

    detect = np.exp(-(In) ** 2 / (2 * sf ** 2 * backGroundIntensity)) < np.random.uniform(size=In.shape)

    #xn = xn[detect]
    #yn = yn[detect]
    #In = In[detect]

    evts = evts[detect]
    
    I = evts['fitResults']['A']

    err_x = np.sqrt(_s2_a2_12 / I + _8_pi_s2_2_a2 * backGroundIntensity / (I * I))
    evts['fitError']['x0'] = err_x

    evts['fitResults']['A'] = I*_r_2_pi
    evts['fitError']['A'] = np.sqrt(I) * _r_2_pi
    
    
    evts['fitError']['x0'] =err_x
    evts['fitError']['y0'] = err_x
    evts['fitError']['z0'] = z_err_scale*err_x

    #fill in the things we don't really need.
    evts['fitResults']['sigma'] = 110.
    evts['fitError']['sigma'] = err_x

    s = evts['fitResults']['x0'].shape

    evts['fitResults']['x0'] = evts['fitResults']['x0'] + evts['fitError']['x0'] * np.random.normal(size=s)
    evts['fitResults']['y0'] = evts['fitResults']['y0'] + evts['fitError']['y0'] * np.random.normal(size=s)
    evts['fitResults']['z0'] = evts['fitResults']['z0'] + evts['fitError']['z0'] * np.random.normal(size=s)

    #filter

    return evts

    #return xn, yn, In
