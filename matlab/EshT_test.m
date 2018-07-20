a=[3,2,1];
vm=0.25;
x=a+1;
[S4,PIVec]=Eshint(vm,a);
eigen=[0;0;0;0;0;0.1];
[D4,u]=Esh_D4_disp(vm,a,x,eigen);