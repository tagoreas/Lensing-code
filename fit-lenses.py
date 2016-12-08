import sys, os, string, emcee
import numpy as np
import astropy.cosmology as cosmo
import astropy.units as apu
import multiprocessing as mp
from helper import setup_libs,get_avd_hmpr,get_lens_calcs
from scipy.optimize import minimize

#####################################
# BEGIN ADJUSTABLE PARAMETERS #######
#####################################

usepixsrc = 1 # flag to either model(1) or not model(0) extended source (requires pixsrc library)
startrow  = 1 # row in 'lenses.dat' file to start modelling at
skiprow   = 1 # model every <skiprow> lenses in 'lenses.dat' file, starting at row <startrow>
symmetry  = 2 # use spherical(0), axisymmetric(1), or triaxial(2) lens models



# initial parameter values, bounds, penalty strengths
# range of parameters for creating initial grid of models
# startp0 is really the minimum value (for creating initial MCMC simplex)
#                     0     1       2     3          4               5            6       7     8       9    10     11  12   13
#                    M_Re slope     q     p        xyang            zang          PA      H0   Reff     R1   R2     RA  Dec  EGAL
startp0 = np.array( [  10.0 , 1.7 , 0.5 , 0.5 ,  0.            ,  0.            , 0.    , 0.47 , 0.  , -2. , 1.25, -10, -10, 0   ] )
rangep0 = np.array( [  2.   , 0.6 , 0.5 , 0.5 , 90.*np.pi/180. , 90.*np.pi/180. , 0.    , 0.4  , 0.  ,  1. , 0.75,  20,  20, 500 ] )

# hard bounds
bounds = np.array( [ [7,15],[1.5,2.5],[0,1],[0,1], \
                         [0,0.5*np.pi],[0,0.5*np.pi],[-20*np.pi/180.,np.pi+20*np.pi/180.], \
                         [0.1,1.5],[0.05,11],[-2,0],[1,2],[-500,500],[-500,500],[0,10000] ] )
# let parameters live outside hard bounds (but only within a certain distance and with a penalty function added).
# parameters are penalized and then reset to bound limits
delbnd = np.array( [ [-0.2,0.2],[-0.2,0.2],[-0.2,0.2],[-0.2,0.2], \
                         [-10.*np.pi/180.,10.*np.pi/180.],[-10.*np.pi/180.,10.*np.pi/180.],[-10.*np.pi/180.,10.*np.pi/180.], \
                         [-0.2,0.2],[-1,1],[-1,1],[-1,1],[-50,50],[-50,50],[-50,50] ] )
sigbnd = np.array( [ [.02,.02],[.02,.02],[.02,.02],[.02,.02], \
                         [2.*np.pi/180.,2.*np.pi/180.],[2.*np.pi/180.,2.*np.pi/180.],[2.*np.pi/180.,2.*np.pi/180.], \
                         [.02,.02],[.02,.02],[.02,.02],[.02,.02],[10,10],[10,10],[5,5] ] )


# number of MCMC walkers and steps (per walker)
nburnwalkers,nburnsteps,nwalkers,nsteps = 300,300,0,0

# number of parallel threads to run
numthreads = mp.cpu_count()
# minimum allowable axis ratio (intermediate or minor)
# velocity dispersion calculations only go down to axis ratio of 0.2 (extrapolation below this will fail!)
minqp = 0.2

# uncertainties on data
ellerr0     = 0.2  # ellipticity
ellangerr0  = 10.0 # position angle (degrees)
avderr0     = 10.  # velocity dispersion (km/s)
poserr0     = 0.05 # position (arcseconds)
tdelerr0    = 0.03 # time delay (fractional uncertainty)
hmprerr0    = 0.10 # half light radius (fractional uncertainty)
tdelerrmin0 = 0.5  # minimum time delay error (days)
hmprerrmin0 = 0.05 # minimum half light radius error (arcseconds)

# seed
#np.random.seed(seed=1)

np.set_printoptions(linewidth=10000000, threshold=10000000, precision=2, suppress=True)

#####################################
# END ADJUSTABLE PARAMETERS #########
#####################################


if 0==symmetry:
    rangep0 = np.delete (rangep0,[2,3,4,5,6])
    startp0 = np.delete (startp0,[2,3,4,5,6])
if 1==symmetry:
    rangep0 = np.delete (rangep0,[3])
    startp0 = np.delete (startp0,[3])   
ndim = len(rangep0)

# global constants and parameters
G          = 4.3022682e-6   # (km/s)^2  kpc/solar masses
arc2rad    = np.pi/(180.0*3600.0)
ckm        = 299792.458     # c in km/s
ckpc       = 9.71561189e-12 # c in kpc/s
c2by4pigarc= 39084.4403     # [c^2/(4 pi G) in solar mass/kpc] /(206265)^2   
c2by4pigkpc= 1.66285417e15  # [c^2/(4 pi G) in solar mass/kpc]





######################################
## FUNCTION PASSED TO MCMC SAMPLER ###
######################################
counter  = 0
countnum = 0
def mcmcfunc (xs00,zlens,zsrc,e1,e2,e12cov,avd,avderr,pos,poserr,tdel,tdelerr,hmpr_arcs,hmprarcserr,egarr):

    global counter,countnum
    counter +=1

    # standardize incoming MCMC parameters array
    if 0==symmetry:
        xs0 = np.array([xs00[0],xs00[1],1.,1.,0.25*np.pi,0.25*np.pi,0,xs00[2],xs00[3],xs00[4],xs00[5],xs00[6],xs00[7],xs00[8]])
    elif 1==symmetry:
        xs0 = np.array([xs00[0],xs00[1],xs00[2],1.,xs00[3],xs00[4],xs00[5],xs00[6],xs00[7],xs00[8],xs00[9],xs00[10],xs00[11],xs00[12]])
    elif 2==symmetry:
        xs0 = np.copy(xs00)
    else:
        print ('Error: Invalid symmetry variable')
        sys.exit()

    # hard priors
    penalty = 0
    for b in range(len(bounds)):
        if xs0[b]<(bounds+delbnd)[b][0]: penalty += (1.+((bounds+delbnd)[b][0]-xs0[b])**2)*1e20
        if xs0[b]>(bounds+delbnd)[b][1]: penalty += (1.+((bounds+delbnd)[b][1]-xs0[b])**2)*1e20
    if 0!=penalty:
        return -0.5*penalty
    
    # soft priors
    val = 0.
    xs = np.copy(xs0)  
    for b in range(len(bounds)):
        if xs0[b]<bounds[b][0]: 
            val += (bounds[b][0]-xs0[b])**2 / sigbnd[b][0]**2
            xs[b] = bounds[b][0]
        if xs0[b]>bounds[b][1]: 
            val += (bounds[b][1]-xs0[b])**2 / sigbnd[b][1]**2
            xs[b] = bounds[b][1]
    val += (xs[11]/100.)**2 + (xs[12]/100.)**2

    # convert from relative axis ratios to absolute axis ratios
    xs[2] = minqp + xs[2]*( 1.0 -minqp)
    xs[3] = minqp + xs[3]*(xs[2]-minqp)


    # some scalings
    xs[0] *=1
    xs[1] *=100.
    xs[2] *=100.
    xs[3] *=100.
    xs[4] *=180./np.pi
    xs[5] *=180./np.pi
    xs[6] *=180./np.pi
    xs[7] *=100.
    xs[8] *=1
    xs[9] *=100.
    xs[10]*=100.
    xs[11]*=0.001
    xs[12]*=0.001
    xs[13]*=0.001


    # get cosmology-dependent parameters and other stuff
    cosmo2 = cosmo.FlatLambdaCDM(H0=xs[7] * apu.km / apu.s / apu.Mpc, Om0=0.307)
    Dod = cosmo2.angular_diameter_distance (zlens).value*1e3
    Dos = cosmo2.angular_diameter_distance (zsrc).value*1e3
    Dds = Dos - (1+zlens) /(1+zsrc) *Dod
    sigmacrkpc = c2by4pigkpc *Dos/Dds/Dod
    sigmacrarc = c2by4pigarc *Dos*Dod/Dds
    kpc2arc = 206264.80624709636/Dod
    arc2kpc    = Dod/206264.80624709636
    hmpr_physical = xs[8]/kpc2arc
    t0 = (1+zlens)/ckpc*(Dod*Dos/Dds)/(3600.*24.)
    numimg = len(pos)
    massRe_physical = 10**xs[0]

    # get half light radius chi2 term
    val0 = ((hmpr_arcs-hmpr_arcs)/hmprarcserr)**2
    val += val0

    # get N-body half light radius and aperture vel. disp.
    avd_nbody,hmpr_nbody = get_avd_hmpr (xs[1],xs[2],xs[3],xs[9],xs[10],xs[4],xs[5])
    # hmpr_nbody<=1 is 2d half light radius in units of 3d half light radius, which is unity
    hmpr_3d_physical = hmpr_physical/hmpr_nbody

    # get lens model predictions
    posoffset = np.concatenate( ( np.array([[xs[11]]*numimg]).T, np.array([[xs[12]]*numimg]).T ), axis=1 )
    lenscalcs    = get_lens_calcs (xs,sigmacrkpc,G,massRe_physical,hmpr_3d_physical,kpc2arc,pos+posoffset,tdel,poserr,tdelerr,t0,arc2rad,1./arc2rad, arc2kpc, egarr)
    ellang1s     = lenscalcs[0]
    ellang2s     = lenscalcs[1]
    massRe_nbody = lenscalcs[2]
    val45        = lenscalcs[3]
    val += val45

    # constrain environment parameter if there is no environment
    if 0==len(egarr):
        mideg = 0.5*(bounds[13][1]+bounds[13][0])*0.001
        sigeg = (bounds[13][1]-bounds[13][0])/20.0
        val45 += (xs[13]-mideg)**2/(sigeg)**2

    # convert Nbody units to physical units and get chi2
    # massRe_nbody<=1 is mass inside Re in units of total mass, which is unity
    avd_model = avd_nbody *np.sqrt(G*(massRe_physical/massRe_nbody)/hmpr_3d_physical)
    val2 = (((avd_model-avd)/avderr)**2)[0]
    val += val2

    
    # ellipticity and PA chi2 terms
    val1 = (np.array([[ellang1s-e1,ellang2s-e2]]).dot(e12cov.dot(np.array([[ellang1s-e1,ellang2s-e2]]).T)))[0,0]
    val += val1

    # half light radius chi2
    val6 = (hmpr_arcs-xs[8])**2/hmprarcserr**2
    val += val6

    # occasionally print out current model
    if counter/10>=countnum: 
        #print (counter,xs,np.array([np.sqrt(val1),np.sqrt(val2),np.sqrt(val45),np.sqrt(val6)]),np.array([val]))
        sys.stdout.write(str(counter)+' ')
        printarr = np.concatenate((xs,[0],np.array([np.sqrt(val1),np.sqrt(val2),np.sqrt(val45),np.sqrt(val6),val])))
        for j in range(len(printarr)):
            sys.stdout.write('{:07.2f} '.format(printarr[j]))
        sys.stdout.write('\n')
        sys.stdout.flush()
        countnum +=1

    # reject models that have failed (hopefully this doesn't ever happen)
    if val!=val:
        return -0.5*1e100
    else:
        return -0.5*val
    


#####################
### MAIN FUNCTION ###
#####################

# load data and tabulated parameters
lenses      = np.loadtxt ('data/lenses.dat')
anggrid     = np.loadtxt ('wts/ang-grid.dat')
numangles   = len(anggrid)
avd_wts  = []
hmpr_wts = []
for jj in range(len(anggrid)):
    avd_wts.append (np.load ('wts/avd-wts-' +str(int(anggrid[jj][0]))+'.npy'))
    hmpr_wts.append(np.load ('wts/hmpr-wts-'+str(int(anggrid[jj][0]))+'.npy'))
for jj in range(len(anggrid)):
    avd_wts[jj]  = np.insert(avd_wts[jj], 0,numangles)
    hmpr_wts[jj] = np.insert(hmpr_wts[jj],0,numangles)
setup_libs (anggrid,avd_wts,hmpr_wts,numthreads,usepixsrc)

# get lenses to analyze
lenses = lenses[startrow-1::skiprow]
lenses = lenses if lenses.ndim==2 else np.array([lenses])

print ('Number of lenses:', len(lenses))
print ('Note: Any warning/error messages will be printed once and then only periodically afterwards.')

# loop over each lens and analyze
for obj0 in lenses:

    np.savetxt('ids.dat',[0]*numthreads,fmt="%d")
    np.savetxt('ids-lock.dat',[])

    outname = str(int(obj0[-2]))+'.'+str(int(obj0[-1]))

    obj = obj0[0:len(obj0)-2]
    objind  = 0

    # extract observables from data file
    zlens           = obj[objind]; objind +=1
    zsrc            = obj[objind]; objind +=1
    r_ein           = obj[objind]; objind +=1
    hmpr_arcs       = obj[objind]; objind +=1
    ell             = obj[objind]; objind +=1
    ellang          = obj[objind]; objind +=1
    avd             = obj[objind]; objind +=1
    numimg          = int(obj[objind]); objind +=1
    posx            = obj[objind+0::5]
    posy            = obj[objind+1::5]
    mag             = obj[objind+2::5]
    tdel            = obj[objind+3::5]
    while ellang<0: ellang += 180
    pos = np.concatenate((np.array([posx]).T,np.array([posy]).T),axis=1)

    # set uncertainties on data
    avderr     = avderr0
    poserr     = np.ones (pos.shape)*poserr0
    tdelerr    = np.clip (tdel*tdelerr0,      tdelerrmin0, np.inf)
    hmprerr    = np.clip (hmpr_arcs*hmprerr0, hmprerrmin0, np.inf)
    # get ellipticity,PA uncertainties
    # convert e,PA uncertainties to e*cos(2PA),e*sin(2PA) covariances
    e1 = ell * np.cos(2*ellang*np.pi/180.)
    e2 = ell * np.sin(2*ellang*np.pi/180.)
    tmprand1 = np.random.normal(0,ellerr0,   1000000) + ell
    tmprand2 = np.random.normal(0,ellangerr0,1000000) + ellang
    ids = np.logical_and(tmprand1>0,tmprand1<1)
    tmprand1=tmprand1[ids]
    tmprand2=tmprand2[ids]
    rand1 = tmprand1 * np.cos(2*tmprand2*np.pi/180.)
    rand2 = tmprand1 * np.sin(2*tmprand2*np.pi/180.)
    e12s  = np.concatenate((np.array([rand1]).T,np.array([rand2]).T),axis=1).T
    e12cov = np.cov(e12s)
    e12cov = np.linalg.inv(e12cov)

    # read in galaxies in environment
    egarr = np.array([])
    if os.path.isfile('data/gal-FJ.dat'):
        egarr = np.loadtxt('data/gal-FJ.dat').flatten()

    # arguments to chi2 function
    mcmcargs = [zlens,zsrc,e1,e2,e12cov,avd,avderr,pos,poserr,tdel,tdelerr,hmpr_arcs,hmprerr,egarr]

    # initialize MCMC sampler
    sampler = emcee.EnsembleSampler (nburnwalkers, ndim, mcmcfunc, threads=numthreads, args=mcmcargs)

    # get grid initial startup models for MCMC
    rangep = np.copy(rangep0)
    startp = np.copy(startp0)
    rsinds = [6,8]
    if 0==symmetry:
        rsinds = [3]
        rangep[rsinds[0]] = 3*hmprerr
        startp[rsinds[0]] = hmpr_arcs - 1.5*hmprerr
    else:
        if 1==symmetry:
            rsinds = [5,7]
        rangep[rsinds[0]] = 3*ellangerr0*np.pi/180.
        rangep[rsinds[1]] = 3*hmprerr
        startp[rsinds[0]] = (ellang   - 1.5*ellangerr0)*np.pi/180.
        startp[rsinds[1]] = hmpr_arcs - 1.5*hmprerr
    p00 = 0.5*rangep + startp
    p0  = np.array([p00 + (np.random.rand(ndim)-0.5)*rangep for i in range(nburnwalkers)])


    # run MCMC sampler
    print ('')
    print ('First 10 rows of initial simplex:')
    print (p0[0:10])
    print ('')
    print ('Running sampler')
    print ('Note: You are likely to get more warnings at the beginning of an MCMC run, when the extremes of parameter ranges are tested.')
    print ('')
    fn = 'chains.'+outname+'.dat'
    f = open(fn, "w")
    f.close()
    for result in sampler.sample(p0, iterations=nburnsteps, storechain=False):
        position = result[0]
        chi2 = np.array([-2*result[1]]).T
        outarr = np.concatenate((position,chi2),axis=1)

        f = open(fn, 'ba')
        np.savetxt(f,outarr)
        f.close()
