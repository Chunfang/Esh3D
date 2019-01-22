#!/usr/bin/env python
import numpy as np
import sys, os

fin = sys.argv[1]
fout = fin.rsplit('.')[0] + '.inp'
incl=False; inho=True
full=True; half=False; fini=False; 
# Host matrix moduli
vp=6000.;vs=3464.;rho=2670
Em=rho*vs**2*(3*vp**2-4*vs**2)/(vp**2-vs**2)
vm=(vp**2-2*vs**2)/2/(vp**2-vs**2)
mat=np.array([[Em,vm]],dtype=np.float64)

if full: 
    if incl:
        line1=['full-incl']
    elif inho:
        line1=['full-inho']
    # Uniaxial remote stress    
    rstress=np.array([[0.,1E7,0.,0.,0.,0.]],dtype=np.float64)
elif half or fini:
    # Zero remote stress (sources rely on eigenstran or boundary loading)
    rstress = np.zeros(shape=(1,6))
    import netCDF4
    nc = netCDF4.Dataset(fin)
    line3=['hex 51'] #['tet 29']
    # t=100; dt=0.01; alpha=0.; beta=0.0025; rfrac=0
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
    if half:     
        if incl:
            line1=['half-incl']
        elif inho:
            line1=['half-inho']
        # Bc type 1 free, 0 fixed
        bc_typ = np.ones((nnd,3), dtype=np.int8)
        bcx_nodes = np.vstack((nc.variables['node_ns1'][:], nc.variables['node_ns2'][:]))
        bcy_nodes = np.vstack((nc.variables['node_ns3'][:], nc.variables['node_ns4'][:]))
        bcz_nodes = nc.variables['node_ns5'][:] 
        # roller sides and bottom
        for node in bcx_nodes:
            bc_typ[node-1, 0] = 0
        for node in bcy_nodes:    
            bc_typ[node-1, 1] = 0
        for node in bcz_nodes:
            bc_typ[node-1, 2] = 0
        # Uniaxial compression
        trac=np.array([0., 0., 0.])
        trac_bc=bnd_el[5]; tracz = np.dot(np.ones(shape=(len(bnd_el[5]),1)),trac.reshape(1,3))
        # xy plane surface traction match
        surf=np.zeros(shape=(len(trac_bc),1))
        for i in range(len(trac_bc)):
            enode=hx_node[trac_bc[i,0]-1,:]
            x=np.mean(coord[enode-1,0]) 
            y=np.mean(coord[enode-1,1])
            if abs(x)<10. and abs(y)<10.: surf[i,0]=1
        trac_bc = np.hstack((trac_bc,tracz,surf))
    elif fini:
        if incl:
            line1=['fini-incl']
        elif inho:
            line1=['fini-inho']
        trac_bc1=bnd_el[0]; trac1 = np.zeros(shape=(len(bnd_el[0]),3))
        trac_bc2=bnd_el[1]; trac2 = np.zeros(shape=(len(bnd_el[1]),3))
        trac_bc3=bnd_el[2]; trac3 = np.zeros(shape=(len(bnd_el[2]),3))
        trac_bc4=bnd_el[3]; trac4 = np.zeros(shape=(len(bnd_el[3]),3))
        trac_bcx=np.vstack((trac_bc1,trac_bc2)); tracx=np.vstack((trac1,trac2))
        trac_bcy=np.vstack((trac_bc3,trac_bc4)); tracy=np.vstack((trac3,trac4))
        # Uniaxial compression at top
        trac=np.array([0., 0., -1E5])
        trac_bcz=bnd_el[5]; tracz = np.dot(np.ones(shape=(len(bnd_el[5]),1)),trac.reshape(1,3))
        # Surface traction match 
        surfx=np.zeros(shape=(len(trac_bcx),1))
        surfy=np.zeros(shape=(len(trac_bcy),1))
        surfz=np.zeros(shape=(len(trac_bcz),1))
        # yz plane 
        for i in range(len(trac_bcx)):
            enode=hx_node[trac_bcx[i,0]-1,:]
            y=np.mean(coord[enode-1,1])
            z=np.mean(coord[enode-1,2]) 
            if abs(y)<15. and abs(z-15.25)<10.: surfx[i,0]=1
        # xz plane 
        for i in range(len(trac_bcy)):
            enode=hx_node[trac_bcy[i,0]-1,:]
            x=np.mean(coord[enode-1,0])
            z=np.mean(coord[enode-1,2]) 
            if abs(x)<15. and abs(z-15.25)<10.: surfy[i,0]=1
        # xy plane 
        for i in range(len(trac_bcz)):
            enode=hx_node[trac_bcz[i,0]-1,:]
            x=np.mean(coord[enode-1,0])
            y=np.mean(coord[enode-1,1]) 
            if abs(x)<15. and abs(y)<15.: surfz[i,0]=1
        trac_bc=np.vstack((np.hstack((trac_bcx,tracx,surfx)),np.hstack((trac_bcy,tracy,surfy)),np.hstack((trac_bcz,tracz,surfz)))) 
        # Bc type space 1 free, 0 fixed
        bc_typ = np.ones((nnd,3), dtype=np.int8)
        bcz_nodes = nc.variables['node_ns5'][:]
        # clamped bottom
        for node in bcz_nodes:
            bc_typ[node - 1, :] = 0
if incl:        
    # ellip(nellip,15): 1-3 ellipsoid centroid coordinate, 4-6 semi-axises, 7-9
    # rotation angles around x,y and z axises, 10-15 eigen strain
    ellip=np.array([[ 0.,0.,-15.25,3.,2.,1.   ,0.,       0.,0.,  0.,  0.,.001,0., 0.,  0.],
                    [ 0.,0.,-15.25,3.,2.,1.E-3,0.,-np.pi/3.,0.,  0.,  0.,  0.,0., 0., -.1],
                    [ 0.,0.,-15.25,3.,2.,1.E-3,0.,-np.pi/4.,0.,  0.,  0.,  0.,0., 0., -.1],
                    [-2.,0.,-15.25,3.,2.,1.   ,0.,       0.,0.,.001,.001,.001,0., 0.,  0.],
                    [ 2.,2.,-15.25,3.,2.,1.E-3,0.,-np.pi/3.,0.,  0.,  0.,  0.,0., 0., -.1]],
                     dtype=np.float64)
elif inho:
    # ellip(nellip,17): 1-3 ellipsoid centroid coordinate, 4-6 semi-axises, 7-9
    # rotation angles around x,y and z axises, 10,11 inclusion Young's modulus 
    # and Poisson's ratio, 12-17 eigen strain
    Eh=2.*Em; vh=vm
    ellip=np.array([[0.,-10., -15.25, 2.,2.,1. ,0.,       0.,0., Eh, vh, 0., 0., 0., 0., 0., 0.],
                    [0., -5., -15.25, 3.,2.,1. ,0.,-np.pi/6.,0., Eh, vh, 0., 0., 0., 0., 0., 0.],
                    [0.,  0., -15.25, 3.,2.,1. ,0.,       0.,0., Eh, vh, 0., 0., 0., 0., 0., 0.],
                    [0.,  5., -15.24, 3.,2.,1. ,0.,-np.pi/3.,0., Eh, vh, 0., 0., 0., 0., 0., 0.],
                    [0., 10., -15.25, 1.,2.,3. ,0.,       0.,0., Eh, vh, 0., 0., 0., 0., 0., 0.]],
                     dtype=np.float64)
nellip=len(ellip)

if half or fini:
    # Okada faults
    m=60; n=40; dx=6./m; dy=4./n; x0=-3.+dx/2.; y0=-2+dy/2; a1=3.; a2=2; a3=1E-3; dip=np.pi/4.
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
            rects[k,:]=[x,y,z,dy,dx,dip,0,-2*0.099*thick*1E3,0.097*thick*1E3*0]
            k+=1
    nrect=len(rects)
else:
    nrect=0

# Observation grid
#m=40; n=30
#ocoord=np.empty((m*n,3),dtype=np.float64)
#x0=-20.5; x1=20.5; z0=-30.25; z1=0.
#x=x0; z=z0; dx=1.; dz=1.
#for i in range(m):
#    x+=dx; z=z0
#    for j in range(n):
#        z+=dz
#        ocoord[i*n+j,:]=[x,0.,z]
#nobs=len(ocoord)

m=41; n=41
ocoord=np.empty((m*n,3),dtype=np.float64)
x0=-5.; x1=5.; z0=-10.; z1=0.
dx=(x1-x0)/(m-1); dz=(z1-z0)/(n-1); z=z0-dz 
for i in range(n):
    z+=dz; x=x0-dx 
    for j in range(m):
        x+=dx
        ocoord[i*m+j,:]=[0.,x,z]
nobs=len(ocoord)

if os.path.isfile(fout): os.remove(fout)
f = open(fout, 'a')
print 'writing to '+fout+' ...'
np.savetxt(f, line1, fmt='%s')
if full:
    np.savetxt(f,np.array([[nellip,nobs]],dtype=np.float64),fmt='%d '*2)
    np.savetxt(f,mat,fmt='%g '*2)
elif half or fini:
    np.savetxt(f,np.array([[nellip,nrect,nobs]],dtype=np.float64),fmt='%d '*3)
    tol=1E1 # Tolerant traction
    ntol=10 # max iterations
    np.savetxt(f,line3,fmt='%s')
    np.savetxt(f,np.array([[nel,nnd,len(trac_bc)]],dtype=np.int32),'%d '*3)
    np.savetxt(f,np.array([[tol,ntol]],dtype=np.float64),fmt='%g %d')
    np.savetxt(f,mat,fmt='%g '*2)
    np.savetxt(f,hx_node,delimiter=' ',fmt='%d '*8)
    np.savetxt(f,np.hstack((coord,bc_typ)),delimiter=' ',fmt='%g '*3+'%d '*3)
    np.savetxt(f,trac_bc,delimiter=' ',fmt='%d '*2+'%g '*3+'%d')
if inho: np.savetxt(f,rstress,fmt='%g '*6)
if incl:
    np.savetxt(f,ellip,delimiter=' ',fmt='%.17f '*15)
elif inho:
    np.savetxt(f,ellip,delimiter=' ',fmt='%.17f '*17)
if nrect>0 and (half or fini): np.savetxt(f,rects,delimiter=' ',fmt='%g '*9)
np.savetxt(f,ocoord,delimiter=' ',fmt='%g '*3)
f.close()
