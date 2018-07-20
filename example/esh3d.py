#!/usr/bin/env python
import numpy as np
import sys, os

fin = sys.argv[1]
fout = fin.rsplit('.')[0] + '.inp'
half=False
vp=6000.;vs=3464.;rho=2670
Em=rho*vs**2*(3*vp**2-4*vs**2)/(vp**2-vs**2)
vm=(vp**2-2*vs**2)/2/(vp**2-vs**2)
stress=np.array([[0.,0.,0.,0.,0.,0.]],dtype=np.float64)
if half:
    frq=400
    import netCDF4
    nc = netCDF4.Dataset(fin)
    line1=['half']
    line3=['hex 51'] #['tet 29']
    mat=np.array([[Em,vm,rho]],dtype=np.float64)
    t=100; dt=0.01; alpha=0.; beta=0.0025; rfrac=0
    print 'Extracting mesh...'
    coord = np.hstack((nc.variables['coordx'][:].\
            reshape(len(nc.variables['coordx']),1),
            nc.variables['coordy'][:].\
            reshape(len(nc.variables['coordy']),1),
            nc.variables['coordz'][:].\
            reshape(len(nc.variables['coordz']),1)))
    nnd = len(coord)
    hx_node = nc.variables["connect1"][:]
    nel=len(hx_node)
    print '%d nodes, %d elements' %(nnd,nel)
    bnd_el=[]
    for i in nc.variables['ss_prop1'][:]:
        els = nc.variables['elem_ss' + str(i)][:]
        sides = nc.variables['side_ss' + str(i)][:]
        bnd_el.append(np.hstack((els.reshape(len(els),1),sides.reshape(len(sides),1))))
    abs_bc1 = np.hstack((bnd_el[0],   np.ones(shape=(len(bnd_el[0]),1)))) 
    abs_bc2 = np.hstack((bnd_el[1],   np.ones(shape=(len(bnd_el[1]),1))))
    abs_bc3 = np.hstack((bnd_el[2], 2*np.ones(shape=(len(bnd_el[2]),1))))
    abs_bc4 = np.hstack((bnd_el[3], 2*np.ones(shape=(len(bnd_el[3]),1))))
    abs_bc5 = np.hstack((bnd_el[4], 3*np.ones(shape=(len(bnd_el[4]),1))))
    abs_bc  = np.vstack((abs_bc1, abs_bc2, abs_bc3, abs_bc4, abs_bc5))
    top_bc  = bnd_el[5]
    surf=np.zeros(shape=(len(top_bc),1))
    for i in range(len(top_bc)):
        enode=hx_node[top_bc[i,0]]
        x=np.mean(coord[enode-1,0]) 
        y=np.mean(coord[enode-1,1])
        if (x**2+y**2)<100.: surf[i,0]=1
    top_bc = np.hstack((top_bc,surf))
else:
    line1=['full']
    mat=np.array([[Em,vm]],dtype=np.float64)
# ellip(nellip,17): 1-3 ellipsoid centroid coordinate, 4-6 semi-axises, 7-9
# rotation angles around x,y and z axises, 10,11 inclusion Young's modulus 
# and Poisson's ratio, 12-17 eigen strain
ellip=np.array([[0.,0.,-10.,3.,2.E-4,1.,0.,0.,0.,Em,vm,0.,0.,0.,0.1,0.,0.],     
                [0.,0.,-20.,3.,2.E-4,1.,0.,0.,0.,Em,vm,0.,0.,0.,-0.1,0.,0.],
                [0.,0.,-10.,3.,2.,1.,0.,0., np.pi/4.,Em,vm,0.,0.,.001,0.,0.,0.],
                [0.,0.,-20.,3.,2.,1.,0.,0.,-np.pi/4.,Em,vm,0.,0.,-.001,0.,0.,0.]],dtype=np.float64)
nellip=len(ellip)
ocoord=np.array([[0.25,0.25,-30.],
                 [0.25,0.25,-25.],
                 [0.25,0.25,-20.],
                 [0.25,0.25,-15.],
                 [0.25,0.25,-10.],
                 [0.25,0.25, -5.],
                 [0.25,0.25,  0.]],dtype=np.float64)

nobs=len(ocoord)
if os.path.isfile(fout): os.remove(fout)
f = open(fout, 'a')
print 'writing to '+fout+' ...'
np.savetxt(f, line1, fmt='%s')
np.savetxt(f,np.array([[nellip,nobs]],dtype=np.float64),fmt='%g '*2)
if half:
    tol=1E1 # Tolerant traction
    ntol=10 # max iterations
    np.savetxt(f,line3,fmt='%s')
    np.savetxt(f,np.array([[nel,nnd,len(abs_bc),len(top_bc),frq]],dtype=np.int32),'%d '*5)
    np.savetxt(f,np.array([[t,dt,alpha,beta,rfrac,tol,ntol]],dtype=np.float64),fmt='%g '*6+'%d')
    np.savetxt(f,mat,fmt='%g '*3)
    np.savetxt(f,hx_node,delimiter=' ',fmt='%d '*8)
    np.savetxt(f,coord,delimiter=' ',fmt='%g '*3)
    np.savetxt(f,abs_bc,delimiter=' ',fmt='%d '*3)
    np.savetxt(f,top_bc,delimiter=' ',fmt='%d '*3)
else:
    np.savetxt(f,mat,fmt='%g '*2)
np.savetxt(f,stress,fmt='%g '*6)
np.savetxt(f,ellip,delimiter=' ',fmt='%g '*17) 
np.savetxt(f,ocoord,delimiter=' ',fmt='%g '*3)
f.close()
