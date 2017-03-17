import time
import math
#import h5py
import numpy as np
#from scipy.optimize import bisect
#from scipy.optimize import minimize
#import scipy.integrate as integrate
#from scipy.special import gamma as gfunc
from ctypes import *
import os
from multiprocessing import Lock
#from scipy import interpolate
#from scipy.interpolate import interp1d

thisdir    = os.path.abspath(os.path.dirname(__file__))

stple_libs   = None
stple_f1     = None
stple_f2     = None
stple_f3     = None
stple_f4     = None

libstple_lenscalc = cdll.LoadLibrary(os.path.join(thisdir, 'lib/libstple_lenscalc.so'))
stple_optimchi2cpp = libstple_lenscalc.stple_imageplane_chi2
stple_optimchi2cpp.argtypes = [np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='C_CONTIGUOUS'), \
                                   np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='C_CONTIGUOUS')]
stple_optimchi2cpp.restype  = None
stple_lenscalcsetupcpp = libstple_lenscalc.setup
stple_lenscalcsetupcpp.argtypes = [c_char_p,c_int]
stple_lenscalcsetupcpp.restype  = None
stple_findindcpp = libstple_lenscalc.find_thread_index
stple_findindcpp.argtypes = None
stple_findindcpp.restype  = c_int
stple_releaseindcpp = libstple_lenscalc.release_thread_index
stple_releaseindcpp.argtypes = [c_int]
stple_releaseindcpp.restype  = None

libevalavdhmpr = cdll.LoadLibrary(os.path.join(thisdir, 'lib/libeval-tps-avd-hmpr.so'))
evalavdhmprcpp1 = libevalavdhmpr.setup_common
evalavdhmprcpp1.argtypes = [np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='C_CONTIGUOUS'), \
                            np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='C_CONTIGUOUS')]
evalavdhmprcpp1.restype  = None
evalavdhmprcpp4 = libevalavdhmpr.setup_angle
evalavdhmprcpp4.argtypes = [c_int, \
                            np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='C_CONTIGUOUS'), \
                            np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='C_CONTIGUOUS')]
evalavdhmprcpp4.restype  = None
evalavdhmprcpp2 = libevalavdhmpr.interpolate
evalavdhmprcpp2.argtypes = [np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='C_CONTIGUOUS')]
evalavdhmprcpp2.restype  = c_double
evalavdhmprcpp3 = libevalavdhmpr.interpolatemany
evalavdhmprcpp3.argtypes = [np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='C_CONTIGUOUS')]
evalavdhmprcpp3.restype  = None

def setup_libs (anggrid,avdwts,hmprwts,numthreads,usepixsrc):

    # setup grids and parameter shifting/scaling factors
    arr1 = np.ascontiguousarray(avdwts[0],'float64')
    arr2 = np.ascontiguousarray((anggrid[:,[1,2]]).flatten(),'float64')
    evalavdhmprcpp1 (arr1,arr2)
    
    # setup weights of TPS at different viewing angles
    for jj in range(len(anggrid)):
        arr3 = np.ascontiguousarray(avdwts[jj], 'float64')
        arr4 = np.ascontiguousarray(hmprwts[jj],'float64')
        evalavdhmprcpp4 (int(anggrid[jj][0]),arr3,arr4)
    
    # setup lensing library
    global stple_libs,stple_f1,stple_f2,stple_f3,stple_f4
    stple_libs   = [None]*numthreads
    stple_f1     = [None]*numthreads
    stple_f2     = [None]*numthreads
    stple_f3     = [None]*numthreads
    stple_f4     = [None]*numthreads
    basename = "mylens" if 1==usepixsrc else ""
    basename = basename.encode('utf-8')

    stple_lenscalcsetupcpp (basename,numthreads)    
    for ind in range(numthreads):
        stple_libs[ind] = cdll.LoadLibrary(os.path.join(thisdir, 'lib/libstple_lenscalc.'+str(ind)+'.so'))
        stple_f1[ind] = stple_libs[ind].stple_imageplane_chi2
        stple_f1[ind].argtypes = [np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='C_CONTIGUOUS'), \
                                      np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='C_CONTIGUOUS')]
        stple_f1[ind].restype  = None
        stple_f2[ind] = stple_libs[ind].setup
        stple_f2[ind].argtypes = [c_char_p,c_int]
        stple_f2[ind].restype  = None
        stple_f3[ind] = stple_libs[ind].release_thread_index
        stple_f3[ind].argtypes = [c_int]
        stple_f3[ind].restype  = None
        stple_f4[ind] = stple_libs[ind].find_thread_index
        stple_f4[ind].argtypes = None
        stple_f4[ind].restype  = c_int
        
        stple_f2[ind](basename,numthreads)

def get_avd_hmpr(p1,p2,p3,p4,p5,p6,p7):
    p1 = np.array([p1]).flatten()
    p2 = np.array([p2]).flatten()
    p3 = np.array([p3]).flatten()
    p4 = np.array([p4]).flatten()
    p5 = np.array([p5]).flatten()
    p6 = np.array([p6]).flatten()
    p7 = np.array([p7]).flatten()
    numpts = len(p1)
    parms = np.ascontiguousarray(np.concatenate(([numpts],p1,p2,p3,p4,p5,p6,p7,p1*0.0,p1*0.0)).astype('float64'),'float64')
    try:
        evalavdhmprcpp3 (parms)
    except:
        return 0.,0.
    return np.exp(parms[1+numpts*7:1+numpts*8]), parms[1+numpts*8:1+numpts*9]

def get_lens_calcs (vdparms, sigmacrkpc, G, M, R, kpc2arc, pos, tdel, poserr, tdelerr, t0, arc2rad, rad2arc, arc2kpc, egals):
    
    global stple_in_use,stple_libs,stple_f1,stple_f2

    s      = vdparms[1]
    q      = vdparms[2]
    p      = vdparms[3]
    xyang  = vdparms[4]
    zang   = vdparms[5]
    userPA = vdparms[6]
    r1     = vdparms[9]
    r2     = vdparms[10]

    s      *= 0.01
    q      *= 0.01
    p      *= 0.01
    xyang  *= np.pi/180.
    zang   *= np.pi/180.
    userPA *= np.pi/180.
    r1      = 10**(r1*1.e-2)
    r2      = 10**(r2*1.e-2)

    # formula taken from http://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    cxy    = math.cos (xyang)
    sxy    = math.sin (xyang)
    cz     = math.cos (zang)
    sz     = math.sin (zang)
    axis   = [cxy*sz,sxy*sz,cz] # angle axis
    view   = [1.,0.,0.]         # viewing axis
    cross   = np.cross (axis, view)
    dot     = np.dot(axis,view)
    vss     = np.array([[0,-cross[2],cross[1]],[cross[2],0,-cross[0]],[-cross[1],cross[0],0]])
    rotmat  = np.identity(3) + vss + vss.dot(vss)*(1-dot)/np.dot(cross,cross)
    rotmat  = rotmat.T
    aratio  = np.array([1.,q,p])

    numimg   = len(pos)  
    argarr = np.array([rotmat[0,0],rotmat[0,1],rotmat[0,2], \
                           rotmat[1,0],rotmat[1,1],rotmat[1,2], \
                           rotmat[2,0],rotmat[2,1],rotmat[2,2], \
                           aratio[0],aratio[1],aratio[2],s, \
                           sigmacrkpc,M,R,R*kpc2arc,r1,r2,userPA,t0,numimg])
    posarr     = pos.flatten()     /(R*kpc2arc)
    poserrarr  = poserr.flatten()  /(R*kpc2arc)
    tdelarr    = tdel.flatten()    /t0 *(rad2arc*arc2kpc/R)**2
    tdelerrarr = tdelerr.flatten() /t0 *(rad2arc*arc2kpc/R)**2
    tdelarr   -= tdelarr[0]
    emptyarr   = np.zeros(4)
    sendarr    = np.concatenate((argarr,posarr,tdelarr,poserrarr,tdelerrarr,emptyarr)).astype('float64')
    numpad     = len(posarr)+len(poserrarr)+len(tdelarr)+len(tdelerrarr)+len(argarr)

    egarr  = np.array(np.concatenate(([len(egals)/3,vdparms[13]],egals,np.array([1,vdparms[14],0,0,1.086747931502459])))).astype('float64')
    
    myind = stple_f4[0]()
    stple_f1[myind](sendarr,egarr)
    stple_f3[myind](myind)

    return sendarr[numpad:]


