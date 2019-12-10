! Copyright (C) 2010-2015 Tabrez Ali 2015-present Chunfang Meng. All rights
! reserved. This file is part of Esh3D. See ../COPYING for license information.

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
  integer :: nnds,nels,ntrc,ntrc_loc,nobs,nobs_loc,nellip,nellip_loc,nrect,    &
     ntol,nrtol,nfix,nsolid,nfluid
  real(8) :: val,top,rstress(6),tol,rtol,mat(2)
  integer,allocatable :: nodes(:,:),work(:),onlst(:,:),surfel_glb(:),surfel(:),&
     surfside_glb(:),surfside(:),idface(:),oel(:),eel(:),bc(:,:),ndfix_glb(:), &
     ndfix(:),bcfix_glb(:,:),bcfix(:,:)
  real(8),allocatable :: coords(:,:),stress(:,:,:),vvec(:),ocoord(:,:),        &
     ocoord_loc(:,:),ellip(:,:),oshape(:,:),surftrc_glb(:,:),surftrc(:,:),     &
     surfloc_glb(:,:),surfloc(:,:),surfmat_glb(:,:),surfmat(:,:),surf_glb(:),  &
     surf(:),surfdat_glb(:,:),surfdat(:,:),resid(:),odat_glb(:,:),odat(:,:),   &
     rect(:,:),surfnrm_glb(:,:),instress(:,:),ellipeff(:,:),solfix(:,:),       &
     Keig(:,:),Feig(:),solok(:,:),Kvol(:,:),Fvol(:),Wsec(:,:),Esec(:),Vsec(:)
     !,KeigInv(:,:),KeigFldInv(:,:),EffEig(:,:)
  real(8),allocatable,target :: uu(:),uu0(:)
  character(12) :: stype
  character(256) :: output_file
  logical :: full,half,fini,incl,inho
  Vec :: Vec_F,Vec_U,Vec_Feig,Vec_Eig,Vec_incl,Vec_FixC,Vec_Fix,Vec_FixF,      &
     Vec_Fvol,Vec_Evol,Vec_dEig
  Mat :: Mat_K,Mat_Kfull,Mat_Keig,Mat_Kfld,Mat_Kvol
  KSP :: Krylov,KryInc,KryFld,KryVol
  PC :: PreCon
  ! Local element/side/node variables
  integer :: el,side,node
  real(8) :: E,nu
  integer,allocatable :: indx(:),enodes(:)
  real(8),allocatable :: k(:,:),f(:),ecoords(:,:)
  ! Variables for parallel code
  integer :: nprcs,rank,ierr
  integer,allocatable :: epart(:),npart(:) ! Partitioning
  ! L-G Mapping
  integer,allocatable :: nmap(:),emap(:),nl2g(:,:),indxmap(:,:),ol2g(:),el2g(:)
  Vec :: Seq_U
  IS :: From,To
  VecScatter :: Scatter
  real(8),pointer :: pntr(:)

contains

  ! Form local [K]
  subroutine FormLocalK(el,k,indx,kfix)
    implicit none
    integer :: el,indx(:)
    real(8) :: k(:,:)
    logical, optional :: kfix
    logical :: fix
    fix=.true.
    if (present(kfix)) fix=kfix
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    E=mat(1); nu=mat(2)
    call FormElK(ecoords,E,nu,k)
    if (fix) call FixBCinLocalK(el,k)
    call FormLocalIndx(enodes,indx)
  end subroutine FormLocalK

  ! Fix BCs (i.e., zero rows/columns) in local [K]
  subroutine FixBCinLocalK(el,k)
    implicit none
    integer :: el,j,j1,j2
    real(8) :: k(:,:)
    do j=1,npel
       do j1=1,dmn
          if (bc(nodes(el,j),j1)==0) then
             j2=dmn*j-dmn+j1
             val=k(j2,j2)
             k(j2,:)=f0; k(:,j2)=f0 ! Zero out rows and columns
             k(j2,j2)=val
          end if
       end do
    end do
  end subroutine FixBCinLocalK

  subroutine FormMatKfull
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7)
#include "petsc.h"
#endif
    integer :: i
    call MatDuplicate(Mat_K,MAT_DO_NOT_COPY_VALUES,Mat_Kfull,ierr)
    do i=1,nels
       call FormLocalK(i,k,indx,kfix=.false.)
       indx=indxmap(indx,2)
       call MatSetValues(Mat_Kfull,eldof,indx,eldof,indx,k,Add_Values,ierr)
    end do
    call MatAssemblyBegin(Mat_Kfull,Mat_Final_Assembly,ierr)
    call MatAssemblyEnd(Mat_Kfull,Mat_Final_Assembly,ierr)
  end subroutine FormMatKfull

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
    integer :: el,side,i,j,nd,snodes(nps),snode
    real(8) :: vvec(:),area
    real(8),allocatable :: vec(:)
    nd=size(vvec); allocate(vec(nd))
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    call EdgeAreaNodes(enodes,ecoords,side,area,snodes)
    vvec=vvec*area/dble(nps)
    do i=1,nps
       vec=vvec
       do j=1,nd
          if (bc(snodes(i),j)==0) vec(j)=f0
       end do
       snode=nl2g(snodes(i),2)
       call ApplyNodalForce(snode,vec)
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

  ! Translation matrix from global to face coordinate
  subroutine Glb2Face(ecoords,side,matrot,idface)
    implicit none
    integer :: side,idface(:)
    real(8) :: ecoords(:,:),matrot(9),vec1(3),vec2(3),vec3(3),vec4(3),cntf(3), &
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

  ! Form RHS to cancel residual traction/displacement surfdat, solfix -> Vec_F
  subroutine MatchSurf
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7)
#include "petsc.h"
#endif
    integer :: i,el,j
    real(8) :: mattmp(3,3),vectmp(3),vectmp0(3),matstr(3,3),vecstr(6),         &
       estress(nip,6)
    call VecZeroEntries(Vec_F,ierr) ! Incremental numerical RHS
    do i=1,ntrc_loc
       do j=1,3
          matstr(j,j)=surfdat(i,12+j)
       end do
       matstr(1,2)=surfdat(i,16); matstr(2,1)=surfdat(i,16)
       matstr(2,3)=surfdat(i,17); matstr(3,2)=surfdat(i,17)
       matstr(1,3)=surfdat(i,18); matstr(3,1)=surfdat(i,18)
       mattmp=reshape(surfmat(i,:),(/3,3/))
       matstr=matmul(matmul(transpose(mattmp),matstr),mattmp)
       vectmp=matstr(:,3)
       ! Subtract initial traction and change sign
       el=surfel(i)
       enodes=nodes(el,:)
       ecoords=coords(enodes,:)
       call FormLocalIndx(enodes,indx)
       call CalcElStress(ecoords,uu0(indx),mat(1),mat(2),estress(:,:))
       vecstr=(/sum(estress(:,1)),sum(estress(:,2)),sum(estress(:,3)),         &
                sum(estress(:,4)),sum(estress(:,5)),sum(estress(:,6))/)        &
              /dble(nip)
       do j=1,3
          matstr(j,j)=vecstr(j)
       end do
       matstr(1,2)=vecstr(4); matstr(2,1)=vecstr(4)
       matstr(2,3)=vecstr(5); matstr(3,2)=vecstr(5)
       matstr(1,3)=vecstr(6); matstr(3,1)=vecstr(6)
       matstr=matmul(matmul(transpose(mattmp),matstr),mattmp)
       vectmp0=matstr(:,3)
       vectmp=-(vectmp-vectmp0)
       ! Rotate to global coordinate
       vectmp=matmul(mattmp,vectmp)
       call ApplyTraction(surfel(i),surfside(i),vectmp)
    end do
    call VecAssemblyBegin(Vec_F,ierr)
    call VecAssemblyEnd(Vec_F,ierr)
    if (fini) call FixBndVecF
  end subroutine MatchSurf

  subroutine GetVecFixC
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7)
#include "petsc.h"
#endif
    integer :: i,j,j1
    call MatGetDiagonal(Mat_K,Vec_FixC,ierr)
    do i=1,size(coords,1)
       do j=1,3
          if (bc(i,j)/=0) then
             j1=(nl2g(i,2)-1)*3+j-1
             call VecSetValue(Vec_FixC,j1,f1,Insert_Values,ierr)
          end if
       end do
    end do
    call VecAssemblyBegin(Vec_FixC,ierr)
    call VecAssemblyEnd(Vec_FixC,ierr)
  end subroutine GetVecFixC

  ! Modify RHS for fixed value boundary condition: Mat_Kfull,solfix->Vec_F
  subroutine FixBndVecF
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7)
#include "petsc.h"
#endif
    integer :: i,j,j1
    real(8) :: val
    do i=1,nfix
       do j=1,3
          if (bcfix(i,j)==0) then
             val=solfix(i,j)-uu0((ndfix(i)-1)*3+j)
             j1=(nl2g(ndfix(i),2)-1)*3+j-1
             call VecSetValue(Vec_Fix,j1,val,Insert_Values,ierr)
             call VecSetValue(Vec_F,j1,-val,Insert_Values,ierr)
          end if
       end do
    end do
    call VecAssemblyBegin(Vec_Fix,ierr)
    call VecAssemblyEnd(Vec_Fix,ierr)
    call MatMult(Mat_Kfull,Vec_Fix,Vec_FixF,ierr)
    call VecAssemblyBegin(Vec_F,ierr)
    call VecAssemblyEnd(Vec_F,ierr)
    call VecPointWiseMult(Vec_F,Vec_F,Vec_FixC,ierr)
    do i=1,nfix
       do j=1,3
          if (bcfix(i,j)==0) then
             j1=(nl2g(ndfix(i),2)-1)*3+j-1
             call VecSetValue(Vec_FixF,j1,f0,Insert_Values,ierr)
          end if
       end do
    end do
    call VecAssemblyBegin(Vec_FixF,ierr)
    call VecAssemblyEnd(Vec_FixF,ierr)
    call VecAXPY(Vec_F,f1,Vec_FixF,ierr)
  end subroutine FixBndVecF

  ! Observation/inclusion nodal base
  subroutine GetObsNd(strng)
    implicit none
    character(2) :: strng
    integer :: neval,ob,el
    integer,allocatable :: nd_full(:,:),pick(:),oel_full(:)
    real(8) :: xmin,xmax,ymin,ymax,zmin,zmax,xmind,xmaxd,ymind,ymaxd,zmind,    &
       zmaxd,dd,du,dl,dr,df,db,d,eta,nu,psi,xob(dmn),N(npel),c,vec12(dmn),     &
       vec13(dmn),vec14(dmn),vec23(dmn),vec24(dmn),vec1o(dmn),vec2o(dmn),      &
       vec15(dmn),vec73(dmn),vec76(dmn),vec78(dmn),vec7o(dmn)
    real(8),allocatable :: N_full(:,:)
    logical :: p_in_dom, p_in_el
    c=0.125d0

    ! Type of the evaluation obs or inclusion
    if (strng=="ob") neval=nobs
    if (strng=="in") neval=nellip
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
       if (strng=="ob") then
          xob=ocoord(ob,:)
       elseif (strng=="in") then
          xob=ellip(ob,1:3)
       end if
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
    if (strng=="ob") then
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
    elseif (strng=="in") then
       nellip_loc=size(pack(pick,pick/=0))
       allocate(el2g(nellip_loc),eel(nellip_loc))
       el2g=pack(pick,pick/=0)
       eel=oel_full(el2g)
    end if
  end subroutine GetObsNd

  ! Form local index
  subroutine FormLocalIndx(enodes,indx)
    implicit none
    integer :: enodes(:),indx(:)
    call FormElIndx(enodes,indx)
  end subroutine FormLocalIndx

  ! Surface data and residual traction (uu, surfdat -> sufdata(:,10:18), resid)
  subroutine SurfSupResid
    implicit none
    integer :: i,el,j
    real(8) :: edisp(npel,3),estress(nip,6),matstr0(3,3),matstr(3,3),          &
       mattmp(3,3),vecstr(6)
    do i=1,ntrc_loc
       el=surfel(i)
       enodes=nodes(el,:)
       ecoords=coords(enodes,:)
       call FormLocalIndx(enodes,indx)
       call CalcElStress(ecoords,uu(indx),mat(1),mat(2),estress(:,:))
       edisp=reshape(uu(indx),(/npel,3/),(/f0,f0/),(/2,1/))
       ! Superpose with surfdat
       vecstr=surfdat(i,13:18)+                                                &
              (/sum(estress(:,1)),sum(estress(:,2)),sum(estress(:,3)),         &
                sum(estress(:,4)),sum(estress(:,5)),sum(estress(:,6))/)        &
              /dble(nip)
       surfdat(i,13:18)=vecstr
       do j=1,3
          matstr(j,j)=vecstr(j)
       end do
       matstr(1,2)=vecstr(4); matstr(2,1)=vecstr(4)
       matstr(2,3)=vecstr(5); matstr(3,2)=vecstr(5)
       matstr(1,3)=vecstr(6); matstr(3,1)=vecstr(6)
       mattmp=reshape(surfmat(i,:),(/3,3/))
       matstr=matmul(matmul(transpose(mattmp),matstr),mattmp)
       ! Subtract initial numerical traction to estimate residual
       call CalcElStress(ecoords,uu0(indx),mat(1),mat(2),estress(:,:))
       vecstr=(/sum(estress(:,1)),sum(estress(:,2)),sum(estress(:,3)),         &
                sum(estress(:,4)),sum(estress(:,5)),sum(estress(:,6))/)        &
              /dble(nip)
       do j=1,3
          matstr0(j,j)=vecstr(j)
       end do
       matstr0(1,2)=vecstr(4); matstr0(2,1)=vecstr(4)
       matstr0(2,3)=vecstr(5); matstr0(3,2)=vecstr(5)
       matstr0(1,3)=vecstr(6); matstr0(3,1)=vecstr(6)
       matstr=matstr-matmul(matmul(transpose(mattmp),matstr0),mattmp)
       resid(i)=sqrt(sum(matstr(:,3)*matstr(:,3)))
       select case(eltype)
       case("tet")
          select case(surfside(i))
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
          select case(surfside(i))
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
       surfdat(i,10:12)=surfdat(i,10:12)+(/sum(edisp(idface,1)),               &
                                           sum(edisp(idface,2)),               &
                                           sum(edisp(idface,3))/)/dble(nps)
    end do
  end subroutine SurfSupResid

  ! Superpose "fixed" boundary displacement
  subroutine FixSup
    implicit none
    integer :: i,j
    do i=1,nfix
       solfix(i,:3)=solfix(i,:3)+uu((/((ndfix(i)-1)*3+j,j=1,3)/))
    end do
  end subroutine FixSup

  ! Superpose observation disp/stress
  subroutine ObsSup
    implicit none
    integer :: ob,i,j,ind(npel),row(eldof)
    real(8) :: vectmp(3,1),strtmp(6),mattmp(3,npel),vecshp(npel,1),            &
       estress(nip,6)
    if ((nobs_loc)>0) then
       do ob=1,nobs_loc
          ind=onlst(ob,:)
          do i=1,npel
             row((/((i-1)*3+j,j=1,3)/))=(/((ind(i)-1)*3+j,j=1,3)/)
          end do
          mattmp=reshape(uu(row),(/3,npel/))
          vecshp=reshape(oshape(ob,:),(/npel,1/))
          vectmp=matmul(mattmp,vecshp)
          odat(ob,10:12)=odat(ob,10:12)+vectmp(:,1) ! Superpose
          enodes=nodes(oel(ob),:)
          ecoords=coords(enodes,:)
          call CalcElStress(ecoords,uu(row),mat(1),mat(2),estress(:,:))
          strtmp=(/sum(estress(:,1)),sum(estress(:,2)),sum(estress(:,3)),      &
                   sum(estress(:,4)),sum(estress(:,5)),sum(estress(:,6))/)     &
                   /dble(nip)
          odat(ob,13:18)=odat(ob,13:18)+strtmp
       end do
    end if
  end subroutine ObsSup

  ! Intrinsic fluid eigenstrains to RHS [Evol] -> [Feig]
  subroutine Evol2Feig(vm,ellip)
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7)
#include "petsc.h"
#endif
    integer :: i,j,j1,j2,idx(6)
    real(8) :: vm,ellip(:,:),vec(6),strain(6),Ch(6,6)
    call VecZeroEntries(Vec_Feig,ierr); Feig=f0
    call VecGetOwnershipRange(Vec_Evol,j1,j2,ierr)
    do i=1,nfluid
       idx=(/((i-1)*6+j-1,j=1,6)/)
       vec=f0
       if (idx(1)>=j1 .and. idx(6)<j2) call VecGetValues(Vec_Evol,6,idx,vec,   &
          ierr)
       call MPI_AllReduce(vec,strain,6,MPI_Real8,MPI_Sum,MPI_Comm_World,ierr)
       call Cmat(f3*ellip(nsolid+i,10)*(f1-f2*vm),vm,Ch)
       strain=matmul(Ch,strain)
       if (rank==nprcs-1) call VecSetValues(Vec_Feig,6,nsolid*6+idx,strain,    &
          Add_Values,ierr)
        Feig(nsolid*6+idx+1)=strain
    end do
    call VecAssemblyBegin(Vec_Feig,ierr)
    call VecAssemblyEnd(Vec_Feig,ierr)
  end subroutine Evol2Feig

  ! Update interacting eiginstrain Vec_Eig -> ellipeff(12:17)
  subroutine UpInhoEigen(eigen,fluid)
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7)
#include "petsc.h"
#endif
    logical,optional :: fluid
    logical :: fld
    integer :: i,n,j,j1,j2,idx(6)
    real(8) :: eigen(:,:),vec(6),strain(6)
    fld=.false.
    if (present(fluid)) fld=fluid
    call VecGetOwnershipRange(Vec_Eig,j1,j2,ierr)
    n=size(eigen,1)
    do i=1,n
       idx=(/((i-1)*6+j-1,j=1,6)/)
       vec=f0
       if (idx(1)>=j1 .and. idx(6)<j2) call VecGetValues(Vec_Eig,6,idx,vec,ierr)
       call MPI_AllReduce(vec,strain,6,MPI_Real8,MPI_Sum,MPI_Comm_World,ierr)
       if (fld) then
          eigen(i,:)=eigen(i,:)+strain
       else
         eigen(i,:)=strain
       end if
    end do
  end subroutine UpInhoEigen

  ! Vec_dEig -> Esec
  subroutine GetEigSec(Esec)
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7)
#include "petsc.h"
#endif
    integer :: i,j,j1,j2,idx(6)
    real(8) :: Esec(:),vec0(6),vec1(6)
    call VecGetOwnershipRange(Vec_dEig,j1,j2,ierr)
    do i=1,nsolid
       idx=(/((i-1)*6+j-1,j=1,6)/)
       vec0=f0
       if (idx(1)>=j1 .and. idx(6)<j2) call VecGetValues(Vec_dEig,6,idx,vec0,  &
          ierr)
       call MPI_AllReduce(vec0,vec1,6,MPI_Real8,MPI_Sum,MPI_Comm_World,ierr)
       Esec(idx+1)=vec1
    end do
  end subroutine GetEigSec

  ! L2 covergence of a n*nelem vector Vect = Vect + dVect
  subroutine ConvergeL2(Vect,dVect,n,nelem,conv)
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7)
#include "petsc.h"
#endif
    integer :: i,j,j1,j2,n,nelem,idx(nelem)
    real(8) :: vec0(nelem),vec1(nelem),tmp,resid,conv
    Vec :: Vect,dVect
    call VecGetOwnershipRange(Vect,j1,j2,ierr)
    resid=f0
    do i=1,(j2-j1)/nelem
       idx=j1+(/((i-1)*6+j-1,j=1,nelem)/)
       call VecGetValues(Vect,6,idx,vec0,ierr)
       call VecGetValues(dVect,6,idx,vec1,ierr)
       tmp=sqrt(sum(vec0*vec0))
       if (tmp>0) then
          tmp=sqrt(sum(vec1*vec1))/tmp
          if (tmp>resid) resid=tmp
       end if
    end do
    call MPI_AllReduce(resid,conv,1,MPI_Real8,MPI_Max,MPI_Comm_World,ierr)
  end subroutine ConvergeL2

  ! Evaluate inclusions stress changes uu (solok) -> instress(nellip,6)
  subroutine InStrEval(ok)
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7)
#include "petsc.h"
#endif
    integer :: ii,i,j,row(eldof)
    real(8) :: strtmp(6),estress(nip,6),ones(nellip),cnt
    logical :: okval
    logical,optional :: ok ! Include Okada source
    instress=f0 ! Reset inclusion stress
    ones=f0
    do ii=1,nellip_loc
       enodes=nodes(eel(ii),:)
       do i=1,npel
          row((/((i-1)*3+j,j=1,3)/))=(/((enodes(i)-1)*3+j,j=1,3)/)
       end do
       ecoords=coords(enodes,:)
       call CalcElStress(ecoords,uu(row),mat(1),mat(2),estress(:,:))
       strtmp=(/sum(estress(:,1)),sum(estress(:,2)),sum(estress(:,3)),         &
                sum(estress(:,4)),sum(estress(:,5)),sum(estress(:,6))/)        &
                /dble(nip)
       instress(el2g(ii),:)=strtmp
       ones(el2g(ii))=ones(el2g(ii))+f1
    end do
    if (present(ok)) then
       okval=ok
    else
       okval=.false.
    end if
    do ii=1,nellip
       call MPI_AllReduce(instress(ii,:),strtmp,6,MPI_Real8,MPI_Sum,           &
          MPI_Comm_World,ierr)
       call MPI_AllReduce(ones(ii),cnt,1,MPI_Real8,MPI_Sum,MPI_Comm_World,ierr)
       strtmp=strtmp/max(cnt,f1) ! Scale duplicates
       if (okval) strtmp=strtmp+solok(ii,4:9)
       instress(ii,:)=strtmp
    end do
  end subroutine InStrEval

  ! Print message
  subroutine PrintMsg(msg)
    implicit none
    character(*) :: msg
    if (rank==0) print*,msg
  end subroutine PrintMsg

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
    if (half .or. fini) then
       ! Surface grid
       dim_dat=(/3,ntrc/)
       call h5screate_simple_f(2,dim_dat,spc_dat,err)
       call h5dcreate_f(idfile,"scoord",h5t_native_double,spc_dat,iddat,err)
       call h5dwrite_f(iddat,h5t_native_double,transpose(surfloc_glb)/km2m,    &
          dim_dat,err)
       call h5dclose_f(iddat,err)
       ! Surface normal
       call h5screate_simple_f(2,dim_dat,spc_dat,err)
       call h5dcreate_f(idfile,"snorm",h5t_native_double,spc_dat,iddat,err)
       call h5dwrite_f(iddat,h5t_native_double,transpose(surfnrm_glb),         &
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
