! Copyright (C) 2010-2015 Tabrez Ali 2015-2018 Chunfang Meng. All rights reserved.
! This file is part of Esh3D. See ../COPYING for license information.

module global

#include <petscversion.h>

  use local
  use HDF5
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7)
  implicit none
#include "petscdef.h"
#else
#include <petsc/finclude/petscksp.h>
  use petscksp
  implicit none
#endif
  ! Global variables
  integer :: nnds,nels,nmts,nfrcs,ntrcs,nabs,nsurf,frq,dsp,steps,tstep,nobs,   &
     nobs_loc,nellip,nrect,nsurf_loc,dir,nabs_loc,ntol
  real(8) :: alpha,beta,t,dt,rfrac,val,top,rstress(6),tol
  integer,allocatable :: nodes(:,:),work(:),fnode(:),telsd(:,:),onlst(:,:),    &
     surfel_glb(:),surfel(:),surfid_glb(:),surfid(:),idface(:),absel_glb(:),   &
     absel(:),absid_glb(:),absid(:),absdir_glb(:),absdir(:),oel(:)
  real(8),allocatable :: coords(:,:),mat(:),stress(:,:,:),vvec(:),ocoord(:,:), &
     ocoord_loc(:,:),ellip(:,:),oshape(:,:),surfloc_glb(:,:),surfloc(:,:),     &
     surfdat(:,:),surfmat_glb(:,:),surfmat(:,:),surf_glb(:),surfdat_glb(:,:),  &
     surf(:),resid(:),odat_glb(:,:),odat(:,:),rect(:,:)
  real(8),allocatable,target :: uu(:),tot_uu(:)
  character(4) :: stype
  character(256) :: output_file
  logical :: half
  Vec :: Vec_F,Vec_U,Vec_Um,Vec_Up
  Vec,pointer :: Vec_W(:)
  Mat :: Mat_K,Mat_M,Mat_Minv
  KSP :: Krylov
  PC :: PreCon
  ! Local element/side/node variables
  integer :: el,side,node
  real(8) :: E,nu,dns,H,B
  integer,allocatable :: indx(:),indxp(:),enodes(:)
  real(8),allocatable :: k(:,:),m(:,:),f(:),ecoords(:,:)
  ! Variables for parallel code
  integer :: nprcs,rank,ierr
  integer,allocatable :: epart(:),npart(:) ! Partitioning
  ! L-G Mapping
  integer,allocatable :: nmap(:),emap(:),nl2g(:,:),indxmap(:,:),ol2g(:)
  Vec :: Seq_U
  IS :: From,To
  VecScatter :: Scatter
  real(8),pointer :: pntr(:)

contains

  ! Form local [K]
  subroutine FormLocalK(el,k,indx)
    implicit none
    integer :: el,indx(:)
    real(8) :: k(:,:)
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    E=mat(1); nu=mat(2)
    call FormElK(ecoords,E,nu,k)
    call FormLocalIndx(enodes,indx)
  end subroutine FormLocalK

  ! Form local [M]
  subroutine FormLocalM(el,m,indx)
    implicit none
    integer :: el,indx(:)
    real(8) :: m(:,:)
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    dns=mat(3)
    call FormElM(ecoords,dns,m)
    call FormElIndx(enodes,indx)
  end subroutine FormLocalM

  ! Apply nodal force
  subroutine ApplyNodalForce(node,vvec)
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7)
#include "petsc.h"
#endif
    integer :: node,i,j
    real(8) :: vvec(:)
    do i=1,dmn
       j=dmn*node-dmn+i-1
       val=vvec(i)
       call VecSetValue(Vec_F,j,val,Add_Values,ierr)
    end do
  end subroutine ApplyNodalForce

  ! Apply traction (EbEAve)
  subroutine ApplyTraction(el,side,vvec)
    implicit none
    integer :: el,side,i,snodes(nps)
    real(8) :: vvec(:),area
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    call EdgeAreaNodes(enodes,ecoords,side,area,snodes)
    vvec=vvec*area/dble(nps)
    snodes=nl2g(snodes,2)
    do i=1,nps
       call ApplyNodalForce(snodes(i),vvec)
    end do
  end subroutine ApplyTraction

  subroutine Cross(a,b,r)
    implicit none
    real(8) :: a(3),b(3),r(3)
    r(1)=a(2)*b(3)-a(3)*b(2)
    r(2)=a(3)*b(1)-a(1)*b(3)
    r(3)=a(1)*b(2)-a(2)*b(1)
  end subroutine Cross

  ! Signed distance point to plane/line
  subroutine Mix(a,b,c,m)
    implicit none
    real(8) :: a(3),b(3),c(3),r(3),m 
    r(1)=a(2)*b(3)-a(3)*b(2)
    r(2)=a(3)*b(1)-a(1)*b(3)
    r(3)=a(1)*b(2)-a(2)*b(1)
    m=(r(1)*c(1)+r(2)*c(2)+r(3)*c(3))/(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))
  end subroutine Mix

  ! Form local damping matrix for elements with viscous dampers
  subroutine FormLocalAbsC(el,side,m,indx)
    implicit none
    integer :: el,side,indx(:)
    real(8) :: m(:,:),matabs(dmn,dmn),vec1(dmn),vec2(dmn),vec3(dmn)
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    E=mat(1); nu=mat(2)
    dns=mat(3)
    select case(eltype) 
    case("tet")
       select case(side) 
       case(1)
          vec1=ecoords(2,:)-ecoords(1,:)
          call Cross(vec1,ecoords(4,:)-ecoords(1,:),vec3)
       case(2)
          vec1=ecoords(4,:)-ecoords(3,:)
          call Cross(vec1,ecoords(2,:)-ecoords(3,:),vec3)
       case(3)
          vec1=ecoords(4,:)-ecoords(1,:)
          call Cross(vec1,ecoords(3,:)-ecoords(1,:),vec3)
       case(4)
          vec1=ecoords(3,:)-ecoords(1,:)
          call Cross(vec1,ecoords(2,:)-ecoords(1,:),vec3)
       end select
       vec1=vec1/sqrt(sum(vec1*vec1))
       vec3=vec3/sqrt(sum(vec3*vec3))
       call Cross(vec1,vec3,vec2)
       matabs(:,1)=vec1; matabs(:,2)=vec2; matabs(:,3)=vec3
    case("hex")
       select case(side) 
       case(1)
          vec1=ecoords(2,:)-ecoords(1,:)
          call Cross(vec1,ecoords(5,:)-ecoords(1,:),vec3)
       case(2)
          vec1=ecoords(3,:)-ecoords(2,:)
          call Cross(vec1,ecoords(6,:)-ecoords(2,:),vec3)
       case(3)
          vec1=ecoords(4,:)-ecoords(3,:)
          call Cross(vec1,ecoords(7,:)-ecoords(3,:),vec3)
       case(4)
          vec1=ecoords(5,:)-ecoords(1,:)
          call Cross(vec1,ecoords(4,:)-ecoords(1,:),vec3)
       case(5)
          vec1=ecoords(4,:)-ecoords(1,:)
          call Cross(vec1,ecoords(2,:)-ecoords(1,:),vec3)
       case(6)
          vec1=ecoords(6,:)-ecoords(5,:)
          call Cross(vec1,ecoords(8,:)-ecoords(5,:),vec3)
       end select
       vec1=vec1/sqrt(sum(vec1*vec1))
       vec3=vec3/sqrt(sum(vec3*vec3))
       call Cross(vec1,vec3,vec2)
       matabs(:,1)=vec1; matabs(:,2)=vec2; matabs(:,3)=vec3
    end select
    call FormElAbsC(enodes,ecoords,side,matabs,E,nu,dns,m)
    call FormElIndx(enodes,indx)
  end subroutine FormLocalAbsC

  ! Translation matrix from global to face coordinate
  subroutine Glb2Face(ecoords,side,matrot,idface)
    implicit none  
    integer :: side,idface(:)
    real(8) :: ecoords(:,:),matrot(9),vec1(3),vec2(3),vec3(3),vec4(3),cntf(3),    &
       cntv(3)
    cntv=(/sum(ecoords(:,1)),sum(ecoords(:,2)),sum(ecoords(:,3))/)/dble(npel)
    select case(eltype)
    case("tet")
       select case(side) 
       case(1)
          vec1=ecoords(2,:)-ecoords(1,:)
          call Cross(vec1,ecoords(4,:)-ecoords(1,:),vec3)
          idface=(/1,2,4/)
       case(2)
          vec1=ecoords(4,:)-ecoords(3,:)
          call Cross(vec1,ecoords(2,:)-ecoords(3,:),vec3)
          idface=(/2,3,4/)
       case(3)
          vec1=ecoords(4,:)-ecoords(1,:)
          call Cross(vec1,ecoords(3,:)-ecoords(1,:),vec3)
          idface=(/1,3,4/)
       case(4)
          vec1=ecoords(3,:)-ecoords(1,:)
          call Cross(vec1,ecoords(2,:)-ecoords(1,:),vec3)
          idface=(/1,2,3/)
       end select
    case("hex")
       select case(side) 
       case(1)
          vec1=ecoords(2,:)-ecoords(1,:)
          call Cross(vec1,ecoords(5,:)-ecoords(1,:),vec3)
          idface=(/1,2,5,6/)
       case(2)
          vec1=ecoords(3,:)-ecoords(2,:)
          call Cross(vec1,ecoords(6,:)-ecoords(2,:),vec3)
          idface=(/2,3,6,7/)
       case(3)
          vec1=ecoords(4,:)-ecoords(3,:)
          call Cross(vec1,ecoords(7,:)-ecoords(3,:),vec3)
          idface=(/3,4,7,8/)
       case(4)
          vec1=ecoords(5,:)-ecoords(1,:)
          call Cross(vec1,ecoords(4,:)-ecoords(1,:),vec3)
          idface=(/1,4,5,8/)
       case(5)
          vec1=ecoords(4,:)-ecoords(1,:)
          call Cross(vec1,ecoords(2,:)-ecoords(1,:),vec3)
          idface=(/1,2,3,4/)
       case(6)
          vec1=ecoords(6,:)-ecoords(5,:)
          call Cross(vec1,ecoords(8,:)-ecoords(5,:),vec3)
          idface=(/5,6,7,8/)
       end select
    end select    
    vec1=vec1/sqrt(sum(vec1*vec1))
    vec3=vec3/sqrt(sum(vec3*vec3))
    cntf=(/sum(ecoords(idface,1)),sum(ecoords(idface,2)),sum(ecoords(idface,3))&
         /)/dble(size(idface))
    vec4=cntf-cntv; vec4=vec4/sqrt(sum(vec4*vec4))
    ! Outward positive
    vec3=vec3*sign(f1,vec3(1)*vec4(1)+vec3(2)*vec4(2)+vec3(3)*vec4(3))
    call Cross(vec1,vec3,vec2)
    matrot(1:3)=vec1; matrot(4:6)=vec2; matrot(7:9)=vec3
  end subroutine Glb2Face

  ! Impose roller boundary condition (Mat_K, Vec_F) 
  subroutine FixBnd(strng)
    implicit none  
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7)
#include "petsc.h"
#endif
    character(3) :: strng
    integer :: i,j,j1,idface(nps),dofface(3),idf(nps)
    real(8) :: kabs(eldof,eldof),valf(nps)
    do i=1,nabs_loc
       call FormLocalK(absel(i),k,indx); kabs=k
       indx=indxmap(indx,2)
       select case(eltype)
       case("tet") 
          select case(absid(i))
          case(1)
             idface=(/1,2,4/)
          case(2)
             idface=(/2,3,4/)
          case(3)
             idface=(/1,3,4/)
          case(4)
             idface=(/1,2,3/)
          end select
       case("hex")   
          select case(absid(i))  
          case(1)
             idface=(/1,2,5,6/)
          case(2)
             idface=(/2,3,6,7/)
          case(3)
             idface=(/3,4,7,8/)
          case(4)
             idface=(/1,4,5,8/)
          case(5)
             idface=(/1,2,3,4/)
          case(6)
             idface=(/5,6,7,8/)
          end select
       end select
       do j=1,nps
          dofface=(/((idface(j)-1)*3+j1,j1=1,3)/)
          kabs(dofface(absdir(i)),:)=-k(dofface(absdir(i)),:)
          kabs(:,dofface(absdir(i)))=-k(:,dofface(absdir(i)))
          kabs(dofface(absdir(i)),dofface(absdir(i)))=f0
          idf(j)=indx(dofface(absdir(i)))
       end do
       valf=f0
       if (strng=="all") call MatSetValues(Mat_K,eldof,indx,eldof,indx,kabs,   &
          Add_Values,ierr)
       call VecSetValues(Vec_F,nps,idf,valf,insert_Values,ierr)
    end do
    if (strng=="all") then
       call MatAssemblyBegin(Mat_K,Mat_Final_Assembly,ierr)
       call MatAssemblyend(Mat_K,Mat_Final_Assembly,ierr)
    end if
    call VecAssemblyBegin(Vec_F,ierr)
    call VecAssemblyEnd(Vec_F,ierr)
  end subroutine FixBnd

  ! Form RHS to free surface traction
  subroutine FreeSurfTrac
    implicit none  
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7)
#include "petsc.h"
#endif
    integer :: i,j,snodes(nps)
    real(8) :: mattmp(3,3),vectmp(3,1),matstr(3,3),vectrac(3),area
    call VecZeroEntries(Vec_F,ierr)
    do i=1,nsurf_loc
       do j=1,3  
          matstr(j,j)=surfdat(i,12+j)
       end do
       matstr(1,2)=surfdat(i,16); matstr(2,1)=surfdat(i,16)
       matstr(2,3)=surfdat(i,17); matstr(3,2)=surfdat(i,17)
       matstr(1,3)=surfdat(i,18); matstr(3,1)=surfdat(i,18)
       mattmp=reshape(surfmat(i,:),(/3,3/))
       matstr=matmul(matmul(transpose(mattmp),matstr),mattmp)
       vectmp(:,1)=matstr(:,3)
       ! Rotate to global coordinate and change sign 
       vectmp=-matmul(mattmp,vectmp)
       vectrac=vectmp(:,1) 
       enodes=nodes(surfel(i),:)
       ecoords=coords(enodes,:)
       call EdgeAreaNodes(enodes,ecoords,surfid(i),area,snodes)
       vectrac=vectrac*area/dble(nps)
       snodes=nl2g(snodes,2)
       do j=1,nps
          call ApplyNodalForce(snodes(j),vectrac)
       end do 
    end do
    call VecAssemblyBegin(Vec_F,ierr)
    call VecAssemblyEnd(Vec_F,ierr)
  end subroutine FreeSurfTrac 

  ! Observation/FD nodal base 
  subroutine GetObsNd
    implicit none
    integer :: neval,ob,el
    integer,allocatable :: nd_full(:,:),pick(:),oel_full(:)
    real(8) :: xmin,xmax,ymin,ymax,zmin,zmax,xmind,xmaxd,ymind,ymaxd,zmind,    &
       zmaxd,dd,du,dl,dr,df,db,d,eta,nu,psi,xob(dmn),N(npel),c,vec12(dmn),     &
       vec13(dmn),vec14(dmn),vec23(dmn),vec24(dmn),vec1o(dmn),vec2o(dmn),      &
       vec15(dmn),vec73(dmn),vec76(dmn),vec78(dmn),vec7o(dmn)
    real(8),allocatable :: N_full(:,:)
    logical :: p_in_dom, p_in_el
    c=0.125d0 
    ! Type of the evaluation Obs or FD
    neval=nobs
    allocate(pick(neval))
    pick=0
    select case(eltype) 
    case("tet"); allocate(nd_full(neval,4),N_full(neval,4))
    case("hex"); allocate(nd_full(neval,8),N_full(neval,8))
    allocate(oel_full(neval))
    end select
    xmind=minval(coords(:,1)); xmaxd=maxval(coords(:,1))
    ymind=minval(coords(:,2)); ymaxd=maxval(coords(:,2))
    zmind=minval(coords(:,3)); zmaxd=maxval(coords(:,3))
    do ob=1,neval ! Observation loop
       xob=ocoord(ob,:)
       p_in_dom=(xob(1)>=xmind .and. xob(1)<=xmaxd .and. xob(2)>=ymind .and.   &
          xob(2)<=ymaxd)
       if (dmn>2) p_in_dom=(p_in_dom .and. xob(3)>=zmind .and. xob(3)<=zmaxd)
       if (p_in_dom) then ! Point probably in domain
          do el=1,nels
             p_in_el=.false.
             enodes=nodes(el,:)
             ecoords=coords(enodes,:)
             xmin=minval(ecoords(:,1)); xmax=maxval(ecoords(:,1))
             ymin=minval(ecoords(:,2)); ymax=maxval(ecoords(:,2))
             zmin=minval(ecoords(:,3)); zmax=maxval(ecoords(:,3))
             ! Point probably in 3D cell
             if (xob(1)>=xmin .and. xob(1)<=xmax .and. xob(2)>=ymin .and.      &
                xob(2)<=ymax .and. xob(3)>=zmin .and. xob(3)<=zmax) then
                select case(eltype)
                case("tet")
                   vec12=ecoords(2,:)-ecoords(1,:)
                   vec13=ecoords(3,:)-ecoords(1,:)
                   vec14=ecoords(4,:)-ecoords(1,:)
                   vec23=ecoords(3,:)-ecoords(2,:)
                   vec24=ecoords(4,:)-ecoords(2,:)
                   vec1o=xob-ecoords(1,:)
                   vec2o=xob-ecoords(2,:)
                   call Mix(vec12,vec13,vec1o,dd)
                   call Mix(vec13,vec14,vec1o,dl)
                   call Mix(vec14,vec12,vec1o,df)
                   call Mix(vec24,vec23,vec2o,d)
                   ! Point in tet
                   if (dd>=f0 .and. dl>=f0 .and. df>=f0 .and. d>=f0) then 
                      call Mix(vec12,vec13,vec14,du)
                      call Mix(vec13,vec14,vec12,dr)
                      call Mix(vec14,vec12,vec13,db)
                      eta=dl/dr; nu=df/db; psi=dd/du
                      N(1)=f1-eta-nu-psi; N(2)=eta; N(3)=nu; N(4)=psi
                      p_in_el=.true.
                   end if
                case("hex")
                   vec12=ecoords(2,:)-ecoords(1,:)
                   vec14=ecoords(4,:)-ecoords(1,:)
                   vec15=ecoords(5,:)-ecoords(1,:)
                   vec73=ecoords(3,:)-ecoords(7,:)
                   vec76=ecoords(6,:)-ecoords(7,:)
                   vec78=ecoords(8,:)-ecoords(7,:)
                   vec1o=xob-ecoords(1,:)
                   vec7o=xob-ecoords(7,:)
                   call Mix(vec12,vec14,vec1o,dd)
                   call Mix(vec15,vec12,vec1o,df)
                   call Mix(vec14,vec15,vec1o,dl)
                   call Mix(vec76,vec78,vec7o,du)
                   call Mix(vec78,vec73,vec7o,db)
                   call Mix(vec73,vec76,vec7o,dr)
                   ! Point in hex
                   if (dd>=f0 .and. dl>=f0 .and. df>=f0 .and. du>=f0 .and.     &
                      dr>=f0 .and. db>=f0) then
                      eta=(dl-dr)/(dr+dl); nu=(df-db)/(df+db)
                      psi=(dd-du)/(dd+du)
                      N(1)=c*(f1-eta)*(f1-nu)*(f1-psi)
                      N(2)=c*(f1+eta)*(f1-nu)*(f1-psi)
                      N(3)=c*(f1+eta)*(f1+nu)*(f1-psi)
                      N(4)=c*(f1-eta)*(f1+nu)*(f1-psi)
                      N(5)=c*(f1-eta)*(f1-nu)*(f1+psi)
                      N(6)=c*(f1+eta)*(f1-nu)*(f1+psi)
                      N(7)=c*(f1+eta)*(f1+nu)*(f1+psi)
                      N(8)=c*(f1-eta)*(f1+nu)*(f1+psi)
                      p_in_el=.true.
                   end if
                end select
             end if ! Point probably in 3D cell
             ! Update local pick list
             if (p_in_el) then
                pick(ob)=ob
                N_full(ob,:)=N
                nd_full(ob,:)=enodes
                oel_full(ob)=el
                exit
             end if
          end do ! Element loop
       end if ! Point probably in domain
    end do ! Observation loop
    nobs_loc=size(pack(pick,pick/=0))
    allocate(ol2g(nobs_loc),ocoord_loc(nobs_loc,dmn),oel(nobs_loc))
    select case(eltype) 
    case("tet"); allocate(onlst(nobs_loc,4),oshape(nobs_loc,4))
    case("hex"); allocate(onlst(nobs_loc,8),oshape(nobs_loc,8))
    end select
    ol2g=pack(pick,pick/=0)
    ocoord_loc=ocoord(ol2g,:)
    onlst=nd_full(ol2g,:)
    oshape=N_full(ol2g,:)
    oel=oel_full(ol2g)
  end subroutine GetObsNd

  ! Form local index
  subroutine FormLocalIndx(enodes,indx)
    implicit none
    integer :: enodes(:),indx(:)
    call FormElIndx(enodes,indx)
  end subroutine FormLocalIndx

  ! Recover surface disp/stress
  subroutine RecoverSurf(i,el)
    implicit none
    integer :: i,el,j
    real(8) :: edisp(npel,3),estress(nip,6),matstr(3,3),mattmp(3,3),vecstr(6)
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    E=mat(1); nu=mat(2)
    call FormLocalIndx(enodes,indx)
    call CalcElStress(ecoords,uu(indx),E,nu,estress(:,:))
    edisp=reshape(uu(indx),(/npel,3/),(/f0,f0/),(/2,1/))
    vecstr=surfdat(i,13:18)+                                                   &
           (/sum(estress(:,1)),sum(estress(:,2)),sum(estress(:,3)),            &
             sum(estress(:,4)),sum(estress(:,5)),sum(estress(:,6))/)/dble(nip)
    surfdat(i,13:18)=vecstr                                         
    do j=1,3
       matstr(j,j)=vecstr(j) 
    end do
    matstr(1,2)=vecstr(4); matstr(2,1)=vecstr(4)
    matstr(2,3)=vecstr(5); matstr(3,2)=vecstr(5)
    matstr(1,3)=vecstr(6); matstr(3,1)=vecstr(6)
    mattmp=reshape(surfmat(i,:),(/3,3/))
    matstr=matmul(matmul(transpose(mattmp),matstr),mattmp)
    resid(i)=sqrt(sum(matstr(:,3)*matstr(:,3)))      
    select case(eltype)
    case("tet")
       select case(surfid(i)) 
       case(1)
          idface=(/1,2,4/)
       case(2)
          idface=(/2,3,4/)
       case(3)
          idface=(/1,3,4/)
       case(4)
          idface=(/1,2,3/)
       end select
    case("hex")
       select case(surfid(i)) 
       case(1)
          idface=(/1,2,5,6/)
       case(2)
          idface=(/2,3,6,7/)
       case(3)
          idface=(/3,4,7,8/)
       case(4)
          idface=(/1,4,5,8/)
       case(5)
          idface=(/1,2,3,4/)
       case(6)
          idface=(/5,6,7,8/)
       end select
    end select
    surfdat(i,10:12)=surfdat(i,10:12)+(/sum(edisp(idface,1)),                  &
                                        sum(edisp(idface,2)),                  &
                                        sum(edisp(idface,3))/)/dble(nps)
  end subroutine RecoverSurf
  
  ! Print message
  subroutine PrintMsg(msg)
    implicit none
    character(*) :: msg
    if (rank==0) print*,msg
  end subroutine PrintMsg

  ! Write results in ASCII VTK (legacy) format
  subroutine WriteOutput
    implicit none
    character(64) :: fmt
    character(256) :: name,name0,name1
    integer,save :: k=0
    integer :: i,j,j1,lnnds,lnels
    real(8),pointer :: field_val(:)
    field_val=>uu
    name0=output_file(:index(output_file,"/",BACK=.TRUE.))
    name1=output_file(index(output_file,"/",BACK=.TRUE.)+1:)
    write(name,'(A,I0,A,A,A)')trim(name0),rank,"_",trim(name1),".vtk"
    open(10,file=adjustl(name),status='replace')
    lnnds=size(coords,1)
    lnels=size(nodes,1)
    write(10,'(A)')"# vtk DataFile Version 2.0"
    write(10,'(A)')"File written by Esh3D"
    write(10,'(A)')"ASCII"
    write(10,'(A)')"DATASET UNSTRUCTURED_GRID"
    write(10,'(A,I0,A)')"POINTS ",lnnds," double"
    fmt="(3(F0.3,1X))"
    select case(dmn)
    case(2)
       do i=1,lnnds
          write(10,fmt)(/(coords(i,:)/km2m),f0/)
       end do
    case(3)
       do i=1,lnnds
          write(10,fmt)(/(coords(i,:)/km2m)/)
       end do
    end select
    write(10,'(A,I0,1X,I0)')"CELLS ",lnels,lnels*(npel+1)
    select case(npel)
    case(3); fmt="(I0,3(1X,I0))"
    case(4); fmt="(I0,4(1X,I0))"
    case(8); fmt="(I0,8(1X,I0))"
    end select
    do i=1,lnels
       write(10,fmt)npel,nodes(i,:)-1
    end do
    write(10,'(A,I0)')"CELL_TYPES ",lnels
    do i=1,lnels
       write(10,'(I0)')vtkid
    end do
    write(10,'(A,I0)')"POINT_DATA ",lnnds
    j=dmn
    write(10,'(A)')"VECTORS displacements double"
    fmt="(3(F0.6,1X))"
    select case(dmn)
    case(2)
       do i=1,lnnds
          j1=i*j
          write(10,fmt)(/field_val(j1-1),field_val(j1),f0/) ! 2D U
       end do
    case(3)
       do i=1,lnnds
          j1=i*j
          write(10,fmt)(/field_val(j1-2),field_val(j1-1),field_val(j1)/) ! 3D U
       end do
    end select
    close(10); k=k+1
  end subroutine WriteOutput

  ! Save Esh3D solution to H5 file
  subroutine EshSave
    implicit none  
    integer(hid_t) :: idfile,iddat,spc_dat
    integer(hsize_t) :: dim_dat(2)
    integer :: err
    character(256) :: name,name0,name1
    name0=output_file(:index(output_file,"/",BACK=.TRUE.))
    name1=output_file(index(output_file,"/",BACK=.TRUE.)+1:)
    write(name,'(A,A,A)')trim(name0),trim(name1),".h5"
    call h5open_f(err)
    call h5fcreate_f(trim(name),H5F_ACC_TRUNC_F,idfile,err)
    ! Obs grid (H5 inverse index order)
    dim_dat=(/3,nobs/)
    call h5screate_simple_f(2,dim_dat,spc_dat,err)
    call h5dcreate_f(idfile,"ocoord",h5t_native_double,spc_dat,iddat,err)
    call h5dwrite_f(iddat,h5t_native_double,transpose(ocoord)/km2m,dim_dat,err)
    call h5dclose_f(iddat,err)
    ! Obs data
    dim_dat=(/size(odat_glb,2),size(odat_glb,1)/)
    call h5screate_simple_f(2,dim_dat,spc_dat,err)
    call h5dcreate_f(idfile,"odat",h5t_native_double,spc_dat,iddat,err)
    call h5dwrite_f(iddat,h5t_native_double,transpose(odat_glb),dim_dat,err)
    call h5dclose_f(iddat,err)
    ! Ellip parameters
    dim_dat=(/size(ellip,2),size(ellip,1)/)
    call h5screate_simple_f(2,dim_dat,spc_dat,err)
    call h5dcreate_f(idfile,"ellip",h5t_native_double,spc_dat,iddat,err)
    call h5dwrite_f(iddat,h5t_native_double,transpose(ellip),dim_dat,err)
    call h5dclose_f(iddat,err)
    if (half) then
       ! Surface grid
       dim_dat=(/3,nsurf/)
       call h5screate_simple_f(2,dim_dat,spc_dat,err)
       call h5dcreate_f(idfile,"scoord",h5t_native_double,spc_dat,iddat,err)
       call h5dwrite_f(iddat,h5t_native_double,transpose(surfloc_glb)/km2m,    &
          dim_dat,err)
       call h5dclose_f(iddat,err)
       ! Surface data
       dim_dat=(/size(surfdat_glb,2),size(surfdat_glb,1)/)
       call h5screate_simple_f(2,dim_dat,spc_dat,err)
       call h5dcreate_f(idfile,"sdat",h5t_native_double,spc_dat,iddat,err)
       call h5dwrite_f(iddat,h5t_native_double,transpose(surfdat_glb),         &
          dim_dat,err)
       call h5dclose_f(iddat,err)
    end if
    call h5fclose_f(idfile,err)
    call h5close_f(err)
  end subroutine EshSave

end module global
