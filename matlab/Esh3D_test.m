vm = 0.25;
Em = 8.;
vh = 0.25;
Eh = 8.;
%dimensiona of the ellipsoid. Must always be a1>=a2>=a3 and ortation angles
dim = [3. 2. 1.];
%[alpha beta sigam]
ang = [pi/2 pi/3 pi/4];
% Applied stress ordering is: sigma11, sigma12, sigma13, sigma22,sigma23,sigma33
stressvec = [0 0 0 0 0 0];
eigp= [0. 0. 0. 0. 0. 0.125];
cnt = [0. 0. 0.]; 
% Vertices
% crd = [  900.       -157.8947  -1200;
%     300.       2052.632   -1800;
%     300.       1736.842   -1800;
%     444.6019   1891.354   -1655.398;
%     499.564      49.75378 -1600.436;
%     900.      -1736.842   -1200;
%     300.       -789.4737  -1800;
%     499.5647   -897.6046  -1600.435;
%     499.2963    681.4534  -1600.704;
%     699.3523    574.8531  -1400.648;
%     900.       -473.6842  -1200;
%     699.3906  -1004.055   -1400.609;
%     699.3863   -688.3044  -1400.614;
%     300.       2368.421   -1800;
%     300.        789.4737  -1800;
%     498.2161    997.8018  -1601.784;
%     900.       -789.4737  -1200;
%     300.       3000.      -1800];
crd = [1. 0.  2.;
       2. 0.  0.];
[u, u_abs, strain, strain_abs, stress, stress_abs] = Esh_sol_vert(Em, vm, Eh, vh, dim, ang, stressvec, eigp, cnt, crd);