import PYMEcs.misc.rectsFromSVG as rs
import numpy as np

m1vec = np.array([[1.7,-2.7,3.2]])
m2xy = rs.transform_rot(m1vec[0,0],m1vec[0,1],-70.0)
m2vec = np.array([[m2xy[0],m2xy[1],-m1vec[0,2]]])

d1vec = np.array([[1.2,0.5,1.0]])
d2vec = np.array([[-0.4,-1.0,-1.2]])

d3xy = rs.transform_rot(d1vec[0,0],d1vec[0,1],-70.0)
d4xy = rs.transform_rot(d2vec[0,0],d2vec[0,1],-70.0)
d3vec = np.array([[d3xy[0],d3xy[1],d1vec[0,2]]])
d4vec = np.array([[d4xy[0],d4xy[1],d2vec[0,2]]])

defconfig = {
    'm1':m1vec, # displacement vector to marker 1
    'm2':m2vec, # displacement vector to marker 2
    'd1':d1vec, # displacement vector to dye1 on marker 1
    'd2':d2vec, # displacement vector to dye2 on marker 1
    'd3':d3vec, # displacement vector to dye3 on marker 2
    'd4':d4vec # displacement vector to dye4 on marker 2
    # these are all specified in a coordinate system for quadrant one (lower left quadrant?)
}

#callps = np.vstack([rs.points_from_sprect(sprect) for sprect in sprects])
import pandas as pd

def points_from_sprectarray_3d(sprects, mode='default', config=None, return_all=False):
    pointdf_array = []
    sulocdf_array = []
    m1df_array = []
    m2df_array = []
    for count, rect in enumerate(sprects):
        # as currently done the points are constructed each time afresh for each rect (i.e. RyR)
        # could cache this in a suitable way
        dyedf,m1df,m2df,sulocdf = points_from_sprect_3d(rect,mode=mode,config=config,return_locs=True)
        ryrid = count + 1 # give each RyR a unique id (starting at 1)
        dyedf['ryrid'] = ryrid
        pointdf_array.append(dyedf)
        sulocdf['ryrid'] = ryrid
        sulocdf_array.append(sulocdf)
        m1df['ryrid'] = ryrid
        m2df['ryrid'] = ryrid
        m1df_array.append(m1df)
        m2df_array.append(m2df)
    alldyesdf = pd.concat(pointdf_array,axis=0,ignore_index=True)
    allsulocdf = pd.concat(sulocdf_array,axis=0,ignore_index=True)
    allm1df = pd.concat(m1df_array,axis=0,ignore_index=True)
    allm2df = pd.concat(m2df_array,axis=0,ignore_index=True)
    if return_all:
        return alldyesdf, allm1df, allm2df, allsulocdf
    else:
        return alldyesdf, allsulocdf

def points_from_sprect_3d(sprect,mode='default', config=None, return_locs=False):
    dyedf, m1df, m2df, sulocdf = get_all_points(mode,config)
    dyedft = transform_dye_locs(dyedf,sprect)
    sulocdft = transform_dye_locs(sulocdf,sprect)
    m1dft = transform_dye_locs(m1df,sprect)
    m2dft = transform_dye_locs(m2df,sprect)
    # do we also want to return m1, m2, sulocs?
    if return_locs:
        return (dyedft, m1dft, m2dft, sulocdft)
    else:
        return dyedft
    

# this could in future be cached
def get_all_points(mode,config):
    # get the q1 point for this mode
    # typically these have been defined with the square origin at 0,0
    q1 = get_q1_point(mode)
    # transform this point to a coordinate system with the square center at 0,0, i.e. left corner at -13.5,-13.5
    q1ctrd = transform_to_centered_coords(q1)
    # now make points for the other quadrants and also the other derived points
    # make one hash/pandas dataframe to have different entries for all the eventual dye coordinates
    # one array for all m1 positions
    # one array for all m2 positions
    # one array for the subunit centroids
    dyedf, m1df, m2df, sulocdf = derive_all_others_from_q1(q1ctrd,config)
    # TODO: possibly transform all back to 0,0 based coordinates
    dyedfs = dyedf.copy()
    dyecoords = np.stack((dyedf['x'],dyedf['y'],dyedf['z']),axis=-1)
    dyecoordss = transform_to_centered_coords(dyecoords,forward=False)
    dyedfs['x'] = dyecoordss[:,0]
    dyedfs['y'] = dyecoordss[:,1]
    dyedfs['z'] = dyecoordss[:,2]
    sulocdfs = sulocdf.copy()
    suloc_coords = np.stack((sulocdf['x'],sulocdf['y'],sulocdf['z']),axis=-1)
    suloc_coordss = transform_to_centered_coords(suloc_coords,forward=False)
    sulocdfs['x'] = suloc_coordss[:,0]
    sulocdfs['y'] = suloc_coordss[:,1]
    sulocdfs['z'] = suloc_coordss[:,2]
    # assign transformed coordinates
    m1dfs = m1df.copy()
    m1_coords = np.stack((m1df['x'],m1df['y'],m1df['z']),axis=-1)
    m1_coordss = transform_to_centered_coords(m1_coords,forward=False)
    m1dfs['x'] = m1_coordss[:,0]
    m1dfs['y'] = m1_coordss[:,1]
    m1dfs['z'] = m1_coordss[:,2]

    m2dfs = m2df.copy()
    m2_coords = np.stack((m2df['x'],m2df['y'],m2df['z']),axis=-1)
    m2_coordss = transform_to_centered_coords(m2_coords,forward=False)
    m2dfs['x'] = m2_coordss[:,0]
    m2dfs['y'] = m2_coordss[:,1]
    m2dfs['z'] = m2_coordss[:,2]
    
    return (dyedfs, m1dfs, m2dfs, sulocdfs)
    
def transform_dye_locs(dyedf,sprect):
    dyedfpts = np.stack((dyedf['x'],dyedf['y'],dyedf['z']),axis=-1)
    dyedfptsT = transform_locs(dyedfpts,sprect)
    dyedft = dyedf.copy()
    dyedft['x'] = dyedfptsT[:,0]
    dyedft['y'] = dyedfptsT[:,1]
    dyedft['z'] = dyedfptsT[:,2]
    return dyedft

def transform_locs(coords,sprect):
    angle = sprect.get('rot_angle',0)
    origin = rs.transform_rot(sprect['x'],sprect['y'],angle)
    # first rotate corner points
    R = rs.rmat_ang(angle)
    ptsrxy = coords[:,None,0:2].dot(R).squeeze()
    # now shift corner points to the origin of this rect
    sptsrxy = ptsrxy + origin[None,:]
    if coords.shape[1] > 2:
        coords3d = np.zeros((sptsrxy.shape[0],3))
        coords3d[:,0:2] = sptsrxy
        coords3d[:,2] = coords[:,2]
    else:
        RuntimeError("should not happen") # this should not really ever be called with 2D data
    return coords3d    

def derive_all_others_from_q1(q1reference,config):
    # make sulocs, m1, m2, d1, d2, d3, d4
    # for d1, d2, d3, d4 also make corresponding MID, QID, DID
    # to make dyedf:
    #        append all d1-d4 coordinates to make 'x', 'y' and 'z' coordinates
    #        save corresponding ids in 'mid', 'qid', 'did'
    # also return sulocs, m1, m2 (plain coordinate arrays with shape N,3?)
    m1q1 = q1reference + config['m1']
    m2q1 = q1reference + config['m2']
    
    sulocs = join_with_quads3d(q1reference)
    m1 =     join_with_quads3d(m1q1)
    m2 =     join_with_quads3d(m2q1)
    
    d1q1 = m1q1 + config['d1']
    d2q1 = m1q1 + config['d2']
    d3q1 = m2q1 + config['d3']
    d4q1 = m2q1 + config['d4']
    
    d1 = join_with_quads3d(d1q1)
    d2 = join_with_quads3d(d2q1)
    d3 = join_with_quads3d(d3q1)
    d4 = join_with_quads3d(d4q1)

    sulocsdf = pd.DataFrame.from_dict({'x' : sulocs[:,0],'y' : sulocs[:,1],'z' : sulocs[:,2],})
    sulocsdf['mid'] = -1
    sulocsdf['qid'] = np.arange(4,dtype='i')
    sulocsdf['did'] = -1

    m1df = pd.DataFrame.from_dict({'x' : m1[:,0],'y' : m1[:,1],'z' : m1[:,2],})
    m1df['mid'] = np.zeros((4),dtype='i')
    m1df['qid'] = np.arange(4,dtype='i')
    m1df['did'] = -1
    
    m2df = pd.DataFrame.from_dict({'x' : m2[:,0],'y' : m2[:,1],'z' : m2[:,2],})
    m2df['mid'] = np.ones((4),dtype='i')
    m2df['qid'] = np.arange(4,dtype='i')
    m2df['did'] = -1
    
    dyedf1 = pd.DataFrame.from_dict({'x' : d1[:,0],'y' : d1[:,1],'z' : d1[:,2],})
    dyedf1['mid'] = 0
    dyedf1['qid'] = np.arange(4,dtype='i')
    dyedf1['did'] = 0
    dyedf2 = pd.DataFrame.from_dict({'x' : d2[:,0],'y' : d2[:,1],'z' : d2[:,2],})
    dyedf2['mid'] = 0
    dyedf2['qid'] = np.arange(4,dtype='i')
    dyedf2['did'] = 1
    dyedf3 = pd.DataFrame.from_dict({'x' : d3[:,0],'y' : d3[:,1],'z' : d3[:,2],})
    dyedf3['mid'] = 1
    dyedf3['qid'] = np.arange(4,dtype='i')
    dyedf3['did'] = 2
    dyedf4 = pd.DataFrame.from_dict({'x' : d4[:,0],'y' : d4[:,1],'z' : d4[:,2],})
    dyedf4['mid'] = 1
    dyedf4['qid'] = np.arange(4,dtype='i')
    dyedf4['did'] = 3
    
    dyedf = pd.concat([dyedf1,dyedf2,dyedf3,dyedf4],axis=0,ignore_index=True)
    return dyedf, m1df, m2df, sulocsdf
    
def get_q1_point(mode):
    if mode == 'default':
        q1 = np.array([[3.5,3.5]]) # anti-clock wise
    elif mode == 'RyR-T1366':
        q1 = np.array([[6,5]])
    elif mode == 'RyR-T2023':
        q1 = np.array([[4.6,12.0]])
    elif mode == 'RyR-D4365':
        q1 = np.array([[3.5,12.5]])
    elif mode == 'RyR-D4365-single': # just the single point in the first quadrant
        q1 = np.array([[3.5,12.5]])
    return add_z_from(q1)

def transform_to_centered_coords(q1,forward = True):
    if forward:
        return q1 + [[-13.5,-13.5,0]]
    else:
        return q1 - [[-13.5,-13.5,0]]

def other_quads3d(points):
    outp = []
    for p in points:
        q1 = np.append(rs.transform_rot(p[0],p[1],90.0),p[2])
        q2 = np.append(rs.transform_rot(q1[0],q1[1],90.0),p[2])
        q3 = np.append(rs.transform_rot(q2[0],q2[1],90.0),p[2])
        outp.append(q1)
        outp.append(q2)
        outp.append(q3)
    return np.array(outp)

def join_with_quads3d(points):
    otherps = other_quads3d(points)
    return np.append(points,otherps,axis=0)

def add_z_from(points,zfrom=None):
    if points.shape[1] < 3:
        withz = np.zeros((points.shape[0],3))
        withz[:,0:2] = points
        if zfrom is not None:
            withz[:,2] = zfrom[:,2]
        else:
            withz[:,2] = 0
    else:
        withz = points
    return withz

def pickdyes(dyedf,m1prob,m2prob,dyemean_poisson,locerr=2.0):
    dyedf['muid'] = 100 * dyedf['ryrid'] + 10 * dyedf['qid'] + dyedf['mid']
    muids = np.unique(dyedf['muid'])
    m1uids = muids[muids %2 == 0]
    m2uids = muids[muids %2 == 1]
    m1picked = m1uids[rs.random_pick(m1uids.shape[0],m1prob)]
    m2picked = m2uids[rs.random_pick(m2uids.shape[0],m2prob)]
    allpickeddyes = dyedf.loc[dyedf['muid'].isin(np.append(m1picked,m2picked))]
    rowst = []
    for row in allpickeddyes.itertuples(index=False):
        for count in range(np.random.poisson(dyemean_poisson)):
            rowst.append(row)
    blink=pd.DataFrame(rowst)
    blinkerr = rs.pymedf_add_err(blink,locerr=locerr)
    
    return blinkerr

def join_dye_su(blinkdf,sudf):
    blinkdf['loctype'] = 0 # loctype 0: blink location
    sudf['loctype'] = 1 # loctype 1: subunit location
    sudfe = rs.pymedf_add_err(sudf,locerr=0.01) # we use a very small locerr
    sudfe['muid'] = -1 # no unique marker ID
    jointdfs = pd.concat([blinkdf,sudfe],axis=0,ignore_index=True)
    jointdfs['suid'] = 10 * jointdfs['ryrid'] + jointdfs['qid'] # unique subunit ID

    return jointdfs

def join_dye_all(blinkdf,dyedf,sudf,m1df,m2df):
    blinkdf['loctype'] = 0 # loctype 0: blink location
    sudf['loctype'] = 1 # loctype 1: subunit location
    sudfe = rs.pymedf_add_err(sudf,locerr=0.01) # we use a very small locerr
    sudfe['muid'] = -1 # no unique marker ID
    
    m1df['loctype'] = 2 # loctype 2: marker 1
    m1dfe = rs.pymedf_add_err(m1df,locerr=0.01) # we use a very small locerr
    m1dfe['muid'] = -1 # no unique marker ID
    
    m2df['loctype'] = 3 # loctype 3: marker 2
    m2dfe = rs.pymedf_add_err(m2df,locerr=0.01) # we use a very small locerr
    m2dfe['muid'] = -1 # no unique marker ID

    dyedf['loctype'] = 4 # loctype 4: dye positions
    dyedfe = rs.pymedf_add_err(dyedf,locerr=0.01) # we use a very small locerr
    dyedfe['muid'] = -1 # no unique marker ID
    
    jointdfs = pd.concat([blinkdf,sudfe,m1dfe,m2dfe,dyedfe],axis=0,ignore_index=True)
    jointdfs['suid'] = 10 * jointdfs['ryrid'] + jointdfs['qid'] # unique subunit ID

    return jointdfs

def labelfrac1(blinks,sulocdf):
    bsuid = 10 * blinks['ryrid'] + blinks['qid']
    suid = 10 * sulocdf['ryrid'] + sulocdf['qid']
    return float(np.unique(bsuid).size)/np.unique(suid).size

def labelfrac2(alldf):
    blinks = alldf.loc[alldf['loctype'] == 0]
    sulocdf = alldf.loc[alldf['loctype'] == 1]
    return float(np.unique(blinks['suid']).size)/np.unique(sulocdf['suid']).size

#For error rate simulations
def expected_subunits(alldf):
    blinks = alldf.loc[alldf['loctype'] == 0]
    return (np.unique(blinks['suid']).size)

def find_msug(data):    
    msug_counting = []
    for val in (np.unique(data.ryrid)):
        data_to_check = data.su_num[data.ryrid == val].tolist()
        duplicateFrequencies = {}
        for i in set(data_to_check):
            duplicateFrequencies[i] = data_to_check.count(i)
        count_msugs(duplicateFrequencies, 13.00, msug_counting)
        count_msugs(duplicateFrequencies, 31.00, msug_counting)
        count_msugs(duplicateFrequencies, 101.00, msug_counting)
        count_msugs(duplicateFrequencies, 269.00, msug_counting)

    return(msug_counting)


def count_msugs(freq_dict, val, msug_counting):    
    try:
        if freq_dict[val] > 1:
            msug_counting.append(freq_dict[val])
    except:
        pass
    
    return


def subunit_grouping_number(alldf):
    su_num_mapping = [13.00, 31.00, 101.00, 269.00]
    su_num = [su_num_mapping[value] for value in alldf.qid]
    alldf2 = alldf.assign(su_num=su_num)
    return alldf2


def wrong_su(data):
    acceptable_su_values = [13.00, 31.00, 101.00, 269.00]
    incorrect_grouping_val = [round(item,2) for item in data.su_num if round(item, 2) not in acceptable_su_values]
  
    return (incorrect_grouping_val) 


def distancecalc(x2,y2,z2, x1,y1,z1):
    dist = ((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)**0.5
    
    return dist, x2-x1, y2-y1, z2-z1


def nearest_neighbour_with_deltas(sulocdf, coalesced, NN, dx, dy, dz):
    NNs = []; dxs = []; dys = []; dzs = []
    
    coalesced = coalesced.reset_index(drop=True)
    sulocdf = sulocdf.reset_index(drop=True)
    
    actualsu_locations= [[sulocdf['x'][val], sulocdf['y'][val], sulocdf['z'][val]] for val in range (len(sulocdf['x']))]
    clump_end = [[coalesced.x[val], coalesced.y[val], coalesced.z[val]] for val in range(coalesced.x.size)]
    for val in clump_end:
        dist =[]
        xs =[]
        ys =[]
        zs =[]
        for originalpoint in actualsu_locations:              
            NN_distances, delta_x, delta_y, delta_z = distancecalc(val[0], val[1], val[2], originalpoint[0], originalpoint[1], originalpoint[2])
            dist.append(NN_distances)
            xs.append(delta_x)
            ys.append(delta_y)
            zs.append(delta_z)
        
        NNs.append(np.min(dist))
        dxs.append(xs[np.argmin(dist)])
        dys.append(ys[np.argmin(dist)])
        dzs.append(zs[np.argmin(dist)])
    # In order to keep everything the same size. 
    NN.append(NNs)
    dx.append(dxs)
    dy.append(dys)
    dz.append(dzs)


def subunit_locations(sulocdf, coalesced, su_index, su_val):
    act_su = sulocdf[sulocdf.qid == su_index]
    su = coalesced[coalesced.su_num == su_val]
    
    return act_su, su


def create_df(data):
    coalesced = pd.DataFrame({
        "x": data.x,
        "y" : data.y,
        "z": data.z,
        "su_num": data.su_num,
        "ryr_id":data.ryrid,   
    })
    
    return coalesced


def run_simulations_for_errors(num_repeats, save_location, recipe_location, layout, dyedf, sulocdf, dye1=0.20,dye2=0):
    import PYMEcs.simulation.simulateLabeling as sl
    from PYMEcs.recipes import localisations
    from PYME.LMVis import pipeline

    dbscan_dist = 8
    if layout == "RyR-D4365":
        dbscan_dist = 7
     
    labfraction = []

    NN = [[] for _ in range(4)]
    dx = [[] for _ in range(4)]
    dy = [[] for _ in range(4)]
    dz = [[] for _ in range(4)]

    msug = []; wsug = []

    classifications_made = []; expected_classifications = []
    correct_classifications = []; misclassifications = []

    for j in range (num_repeats):    
        blinks = pickdyes(dyedf, dye1, dye2, 1.0,locerr=2.0) 
        alldf = sl.join_dye_su(blinks,sulocdf) # Join the blinks with the actual subunit locations
        alldf = subunit_grouping_number(alldf) # Add a new column with distinct numbers based on the qid

        labfraction.append(100*labelfrac2(alldf))

        alldf.to_csv(save_location + "/example_data_v{0}.csv".format(j),index=False) # Save each simulation to csv    

        data = pipeline.Pipeline(save_location + "/example_data_v{0}.csv".format(j))  # Load the newly created .csv into PyME
        data.recipe.update_from_yaml(recipe_location + "/Separate actual from event and clump dbscan_{}nm_v2.yaml".format(dbscan_dist)) # Run the DBScan recipe
        data.selectDataSource('merged') #Select the coalesced data

        coalesced = create_df(data)

        for i in range(4):
            act_su, su = subunit_locations(sulocdf, coalesced, i, [13, 31, 101, 269][i]) # Filter the data to only have a single subunit type
            nearest_neighbour_with_deltas(act_su, su, NN[i], dx[i], dy[i], dz[i]) # Calculate NN & deltas for individual subunits 0-3

        msug_initial_counter = find_msug(data) 
        msug_val = sum(msug_initial_counter)-len(msug_initial_counter) # This catches missed subunit groupings by inspecting events with the same RyR and SU #s.
        msug.append(msug_val)

        wsug_val = len(wrong_su(data)) # This catches wrong subunit groupings by inspecting the SU# and checking if the value matches those expected 
        wsug.append(wsug_val)

        misclassifications.append(msug_val+wsug_val) # This is the sum of the 2 possible misclassifications 
        correct_classifications.append(data.x.size-msug_val-wsug_val)     
        classifications_made.append(data.x.size)
        expected_classifications.append(expected_subunits(alldf))

    #Once the simulations are complete, save the data. 
    dyes_found_per = [round(100*(correct_classifications[i]/expected_classifications[i]),2) for i in range(len(correct_classifications))]
    correct_classification_per = [round(100*(correct_classifications[i]/classifications_made[i]),2) for i in range(len(correct_classifications))]
    misclassification_per = [round(100*(misclassifications[i]/classifications_made[i]),2) for i in range(len(correct_classifications))]
    wsug_per = [round(100*(wsug[i]/classifications_made[i]),2) for i in range(len(correct_classifications))]
    msug_per = [round(100*(msug[i]/classifications_made[i]),2) for i in range(len(correct_classifications))]

    result = pd.DataFrame({
        "labelfrac": labfraction,
        "dye_locations" : expected_classifications,
        "classifications_made": classifications_made,
        "correct_classifications": correct_classifications,
        "misclassifications":misclassifications,
        "wrong_subunit_grouping" : wsug,
        "missed_subunit_grouping" : msug,    
        "correct_classification_per": correct_classification_per,
        "misclassification_per": misclassification_per,
        "dyes_found_per": dyes_found_per,
        "msug_per" : msug_per,
        "wsug_per" : wsug_per,
    })

    result.to_csv(save_location + "/results.csv".format(j),index=False)
    
    delta_lists = [
    {
        "dx_list": [item for sublist in dx for item in sublist],
        "dy_list": [item for sublist in dy for item in sublist],
        "dz_list": [item for sublist in dz for item in sublist],
    }
    for dx, dy, dz in zip([dx[0], dx[1], dx[2], dx[3]], [dy[0], dy[1], dy[2], dy[3]], [dz[0], dz[1], dz[2], dz[3]])
    ]

    for i, delta_list in enumerate(delta_lists):
        df = pd.DataFrame(delta_list)
        df.to_csv(save_location + "/delta_lists_su{}.csv".format(i), index=False)