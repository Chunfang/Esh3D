function [u, u_abs, strain, strain_abs, stress, stress_abs] = Esh_sol_vert(Em, vm, Eh, vh, dim, ang, stressvec, eigp, cnt, crd)
%converts isotropic constants into a matrix
Cm = Ctensord(Em,vm);
Ch = Ctensord(Eh,vh);
stressvec = stressvec';
eigp = eigp';
% fake a1>=a2>=a3
exh = zeros(3, 3);
for i=1:2
    for j=2:3
        if dim(i) < dim(j)
            exh(i, j) = 1;
            tmp = dim(i);
            dim(i) = dim(j);
            dim(j) = tmp;
        end
    end
end
% pre-rotation in order of [z,y,x]
ang_init = pi/2*[exh(2,3) exh(1,3) exh(1,2)];
Rx = [1 0 0;0 cos(ang_init(1)) -sin(ang_init(1));0 sin(ang_init(1)) cos(ang_init(1))];
Ry = [cos(ang_init(2)) 0 sin(ang_init(2));0 1 0;-sin(ang_init(2)) 0 cos(ang_init(2))];
Rz = [cos(ang_init(3)) -sin(ang_init(3)) 0;sin(ang_init(3)) cos(ang_init(3)) 0;0 0 1];
R_init = Rx*Ry*Rz;
angb_init = -ang_init;
Rx = [1 0 0;0 cos(angb_init(1)) -sin(angb_init(1));0 sin(angb_init(1)) cos(angb_init(1))];
Ry = [cos(angb_init(2)) 0 sin(angb_init(2));0 1 0;-sin(angb_init(2)) 0 cos(angb_init(2))];
Rz = [cos(angb_init(3)) -sin(angb_init(3)) 0;sin(angb_init(3)) cos(angb_init(3)) 0;0 0 1];
Rb_init = Rz*Ry*Rx;
delC=Cm-Ch;
% rotation matrices w.r.t the ellipsoid
Rx = [1 0 0;0 cos(ang(1)) -sin(ang(1));0 sin(ang(1)) cos(ang(1))];
Ry = [cos(ang(2)) 0 sin(ang(2));0 1 0;-sin(ang(2)) 0 cos(ang(2))];
Rz = [cos(ang(3)) -sin(ang(3)) 0;sin(ang(3)) cos(ang(3)) 0;0 0 1];
R = Rz*Ry*Rx;
Rx = [1 0 0;0 cos(-ang(1)) -sin(-ang(1));0 sin(-ang(1)) cos(-ang(1))];
Ry = [cos(-ang(2)) 0 sin(-ang(2));0 1 0;-sin(-ang(2)) 0 cos(-ang(2))];
Rz = [cos(-ang(3)) -sin(-ang(3)) 0;sin(-ang(3)) cos(-ang(3)) 0;0 0 1];
Rb = Rx*Ry*Rz;
% rotate stress and initial eigenstrain against oblique ellipsoid
stressten_init = [stressvec(1:3)';stressvec([2 4 5])'; stressvec([3 5 6])'];
stressten_rot = R_init*Rb*stressten_init*(R_init*Rb)';
stressvec_init = stressvec;
stressvec = [stressten_rot(1,:) stressten_rot(2,[2,3]) stressten_rot(3,3)]';
eigpten_init = [eigp(1:3)';eigp([2 4 5])'; eigp([3 5 6])'];
eigpten_rot = R_init*Rb*eigpten_init*(R_init*Rb)';
eigpvec = [eigpten_rot(1,:) eigpten_rot(2,[2,3]) eigpten_rot(3,3)]';
%correspondingly, the applied strain
epsvec=Cm\stressvec;
% pricipal strain and direction
epsv = Cm\stressvec_init; 
epst = [epsv(1:3)';epsv([2 4 5])'; epsv([3 5 6])'];
[epsv_prin epst_prin] = eig(epst);
epst_prinv = [epst_prin(1,1) epst_prin(2,2) epst_prin(3,3)];
%call the internal eshelby tensor.
S4=Eshint(vm,dim);
eigen=(delC*S4-Cm)\(-delC*epsvec-Ch*eigpvec);
for i = 1:size(crd,1)
    vert = (R_init*Rb*(crd(i, :) - cnt)')';
    [D4 disp] = Esh_D4_disp(vm, dim, vert, eigen);
    rD4(i, :, :)=Cmatrix(D4);
    u(i, :) = R*Rb_init*disp';
    %add remote stress
    prin_cos=dot(crd(i, :)'*[1 1 1], epsv_prin, 1);
    ur = epsv_prin*(sqrt((crd(i, 1) - cnt(1))^2+(crd(i, 2) - cnt(2))^2+(crd(i, 3) - cnt(3))^2)*prin_cos.*epst_prinv)';
    u_abs(i, :) = u(i, :)+ur';
    if  vert(1)^2/dim(1)^2+vert(2)^2/dim(2)^2+vert(3)^2/dim(3)^2<=1 % for interior points
        strainr_pert = S4*eigen;
        strainr = epsvec+strainr_pert;
        strainten = [strainr(1:3)';strainr([2 4 5])'; strainr([3 5 6])'];
        strainten = R*Rb_init*strainten*(R*Rb_init)';
        strainten_pert = [strainr_pert(1:3)';strainr_pert([2 4 5])'; strainr_pert([3 5 6])'];
        strainten_pert = R*Rb_init*strainten_pert*(R*Rb_init)';
        strain_abs(i, :) = [strainten(1,:) strainten(2,[2,3]) strainten(3,3)];
        strain(i, :) = [strainten_pert(1,:) strainten_pert(2,[2,3]) strainten_pert(3,3)];
        stressr = stressvec+Cm*(S4*(eigen)-eigen);
        stressr_pert = Cm*(S4*(eigen)-eigen);
        stressrten = [stressr(1:3)';stressr([2 4 5])'; stressr([3 5 6])'];
        stressrten = R*Rb_init*stressrten*(R*Rb_init)';
        stressrten_pert = [stressr_pert(1:3)';stressr_pert([2 4 5])'; stressr_pert([3 5 6])'];
        stressrten_pert = R*Rb_init*stressrten_pert*(R*Rb_init)';
        stress_abs(i, :) = [stressrten(1,:) stressrten(2,[2,3]) stressrten(3,3)];
        stress(i, :) = [stressrten_pert(1,:) stressrten_pert(2,[2,3]) stressrten_pert(3,3)];
    else % exterior point
        strainr = epsvec+squeeze(rD4(i, :, :))*eigen;
        strainr_pert = squeeze(rD4(i, :, :))*eigen;
        strainten = [strainr(1:3)'; strainr([2 4 5])'; strainr([3 5 6])'];
        strainten = R*Rb_init*strainten*(R*Rb_init)';
        strainten_pert = [strainr_pert(1:3)'; strainr_pert([2 4 5])'; strainr_pert([3 5 6])'];
        strainten_pert = R*Rb_init*strainten_pert*(R*Rb_init)';
        strain(i, :) = [strainten_pert(1,:) strainten_pert(2,[2,3]) strainten_pert(3,3)];
        strain_abs(i, :) = [strainten(1,:) strainten(2,[2,3]) strainten(3,3)];
        stressr = stressvec+Cm*squeeze(rD4(i, :, :))*eigen;
        stressr_pert = Cm*squeeze(rD4(i, :, :))*eigen;
        stressrten = [stressr(1:3)';stressr([2 4 5])'; stressr([3 5 6])'];
        stressrten = R*Rb_init*stressrten*(R*Rb_init)';
        stressrten_pert = [stressr_pert(1:3)';stressr_pert([2 4 5])'; stressr_pert([3 5 6])'];
        stressrten_pert = R*Rb_init*stressrten_pert*(R*Rb_init)';
        stress_abs(i, :) = [stressrten(1, :) stressrten(2, [2, 3]) stressrten(3, 3)];
        stress(i, :) = [stressrten_pert(1, :) stressrten_pert(2, [2, 3]) stressrten_pert(3, 3)];
    end
end
return