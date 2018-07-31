#!/usr/bin/env python
import numpy as np
import sys, os

fin = sys.argv[1]
fout = fin.rsplit('.')[0] + '.inp'
half=True
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
ellip=np.array([[0.,0.,-15.25,3.,2.,1.,0.,0.,0.,Em,vm,0.,0.,.001,0.,0.,0.],
                [0.,0.,-15.25,3.,2.,1.E-3,0.,-np.pi/2.,0.,Em,vm,0.,0.,0,0.,0.1*0,0.],
                [-2.,0.,-15.25,3.,2.,1.,0.,0.,0.,Em,vm,.001,.001,.001,0.,0.,0.],
                [2.,2.,-15.25,3.,2.,1.E-3,0.,-np.pi/3.,0.,Em,vm,0.,0.,0.,0.,0.,-.1]],
                 dtype=np.float64)
ellip=ellip[[1],:]
#ellip=np.array([[0.,0.,-15.25,3.,2.,1E-3,0.,0.,0.,Em,vm,0.,0.,.1*0,0.,0.,0.]])   
nellip=len(ellip)

# Okada faults
m=60; n=40; dx=6./m; dy=4./n; x0=-3.+dx/2.; y0=-2+dy/2; a1=3.; a2=2; a3=1E-3; dip=np.pi/2.
rects=np.empty((m*n,9),dtype=np.float64)
k=0
for i in range(m):
    for j in range(n):
        x=x0+i*dx; y=y0+j*dy
        if x**2/a1**2+y**2/a2**2<1:
            thick=2.*a3*np.sqrt(1-x**2/a1**2-y**2/a2**2)
        else:
            thick=0.
        x=np.cos(dip)*x; z=-15.25+np.sin(dip)*x 
        rects[k,:]=[x,y,z,dy,dx,dip,-2.*thick*0.1*1E3,0.,0.]
        k+=1
#rects=np.array([[0.,0.,-15.25,4,6,90,0,0,0.1]]) 
nrect=len(rects)

# Observation grid
m=40; n=30
ocoord=np.empty((m*n,3),dtype=np.float64)
x0=-20.5; x1=20.5; z0=-30.25; z1=0.
x=x0; z=z0; dx=1.; dz=1.
for i in range(m):
    x+=dx; z=z0
    for j in range(n):
        z+=dz
        ocoord[i*n+j,:]=[x,0.,z]
nobs=len(ocoord)
if os.path.isfile(fout): os.remove(fout)
f = open(fout, 'a')
print 'writing to '+fout+' ...'
np.savetxt(f, line1, fmt='%s')
np.savetxt(f,np.array([[nellip,nrect,nobs]],dtype=np.float64),fmt='%d '*3)
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
if nrect>0 and half: np.savetxt(f,rects,delimiter=' ',fmt='%g '*9)
np.savetxt(f,ocoord,delimiter=' ',fmt='%g '*3)
f.close()
