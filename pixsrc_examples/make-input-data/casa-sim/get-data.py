import numpy as np
np.set_printoptions(threshold='nan')

c = 299792458.0
freq = [245000000000]*4
spws = [0]
pols = [0,1]

projname='merger.uv'

tb.open(projname+'/'+projname+'.alma.out20.noisy.ms', nomodify=F)
n      = tb.nrows()
uvw0   = tb.getcol('UVW',0,n,1)[0:2,:]
data0  = tb.getcol('DATA',0,n,1)
flag0  = tb.getcol('FLAG',0,n,1)
sigma0 = tb.getcol('SIGMA',0,n,1)
out = np.array([]).reshape((0,5))

for i in spws:
    for j in pols:
        flag  = flag0[j,i,:]
        uvw   = uvw0[:,0==flag]/c*freq[i]
        data  = data0[j,i,0==flag]
        sigma = sigma0[j,0==flag]
        
        adder = np.hstack ((uvw.transpose(),np.array([np.real(data)]).transpose(),np.array([np.imag(data)]).transpose(),np.array([sigma]).transpose()))
        out = np.vstack((out,adder))

numtotal = out.shape[0]
np.insert(np.ravel (out,order='F'),0,numtotal).astype ('float64').tofile ('merger.uv.bin')
np.savetxt('merger.uv.vis', out, fmt="%.15e")
tb.close()





tb.open(projname+'/'+projname+'.alma.out20.ms', nomodify=F)
n      = tb.nrows()
uvw0   = tb.getcol('UVW',0,n,1)[0:2,:]
data0  = tb.getcol('DATA',0,n,1)
flag0  = tb.getcol('FLAG',0,n,1)
sigma0 = tb.getcol('SIGMA',0,n,1)
out = np.array([]).reshape((0,5))

for i in spws:
    for j in pols:
        flag  = flag0[j,i,:]
        uvw   = uvw0[:,0==flag]/c*freq[i]
        data  = data0[j,i,0==flag]
        sigma = sigma0[j,0==flag]
        
        adder = np.hstack ((uvw.transpose(),np.array([np.real(data)]).transpose(),np.array([np.imag(data)]).transpose(),np.array([sigma]).transpose()))
        out = np.vstack((out,adder))

numtotal = out.shape[0]
np.insert(np.ravel (out,order='F'),0,numtotal).astype ('float64').tofile ('merger.uv.true.bin')
np.savetxt('merger.uv.true.vis', out, fmt="%.15e")
tb.close()
