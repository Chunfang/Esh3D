! Copyright (C) 2010-2015 Tabrez Ali 2015-2018 Chunfang Meng. All rights reserved.
! This file is part of Esh3D. See ../COPYING for license information.

program main

#include <petscversion.h>

  use esh3d  
  use global
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7)
  implicit none
#include "petsc.h"
#else
#include <petsc/finclude/petscksp.h>
  use petscksp
  implicit none
#endif
  character(256) :: input_file
  logical :: l
  integer,pointer :: null_i=>null()
  real(8),pointer :: null_r=>null()
  integer :: n,i,j,j1,j2,nodal_bw,status(MPI_STATUS_SIZE)
  integer,allocatable :: hit(:)

  call PetscInitialize(Petsc_Null_Character,ierr)

  call MPI_Comm_Rank(MPI_Comm_World,rank,ierr)
  call MPI_Comm_Size(MPI_Comm_World,nprcs,ierr)

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=6) 
  call PetscOptionsGetString(Petsc_Null_Character,'-f',input_file,l,ierr)
#elif (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==7)
  call PetscOptionsGetString(Petsc_Null_Object, Petsc_Null_Character,'-f',     &
     input_file,l,ierr)
#else
  call PetscOptionsGetString(Petsc_Null_Options,Petsc_Null_Character,'-f',     &
     input_file,l,ierr)
#endif
  if (.not. l) then
     call PrintMsg("Usage: [mpiexec -n <np>] defmod -f <input_filename>")
     go to 9
  end if

  ! Read input file parameters
  open(10,file=input_file,status='old')

  ! Output file name is the same as the input
  if (index(input_file,".inp")>0) then
     output_file=input_file(1:index(input_file,".inp")-1)
  else
     output_file=input_file
  end if

  call PrintMsg("Reading input ...")
  call ReadParameters

  if (half) then ! Half space problem
     ! Set element specific constants
     call InitializeElement
     ! Partition mesh using METIS, create mappings, and read on-rank mesh data
     allocate(npart(nnds),epart(nels)); epart=0; npart=0
     if (nprcs>1) then
        call PrintMsg("Partitioning mesh ...")
        if (rank==0) then
           allocate(nodes(1,npel*nels),work(nels+1)); work(1)=0
           do i=1,nels
              j=npel*(i-1)+1; n=npel*i; read(10,*)nodes(1,j:n); work(i+1)=n
           end do
           nodes=nodes-1
           call METIS_PartMeshNodal(nels,nnds,work,nodes,null_i,null_i,nprcs,     &
              null_r,null_i,n,epart,npart)
           deallocate(nodes,work)
           rewind(10); call ReadParameters
        end if
        call MPI_Bcast(npart,nnds,MPI_Integer,0,MPI_Comm_World,ierr)
        call MPI_Bcast(epart,nels,MPI_Integer,0,MPI_Comm_World,ierr)
     end if
     call PrintMsg("Reading mesh data ...")
     allocate(emap(nels),nmap(nnds)); emap=0; nmap=0
     ! Create original to local element mappings and read on-rank element data
     j=1
     do i=1,nels
        if (epart(i)==rank) then
           epart(i)=1; emap(i)=j; j=j+1
        else
           epart(i)=0
        end if
     end do
     n=sum(epart); allocate(nodes(n,npel)) ! id -> mtrl flag
     j=1
     do i=1,nels
        if (epart(i)==1) then
           read(10,*)nodes(j,:); j=j+1
        else
           read(10,*)val
        end if
     end do
     nels=n
     ! Create original to global nodal mappings and read on-rank + ghost nodes
     allocate(work(0:nprcs-1))
     j=0
     do i=1,nnds
        if (npart(i)==rank) j=j+1
     end do
     call MPI_AllGather(j,1,MPI_Integer,work,1,MPI_Integer,MPI_Comm_World,ierr)
     if (rank==0) n=1
     if (rank/=0) n=sum(work(0:rank-1))+1
     do i=1,nnds
        if (npart(i)==rank) then
           nmap(i)=n; n=n+1
        end if
     end do
     deallocate(work)
     allocate(work(nnds))
     call MPI_AllReduce(nmap,work,nnds,MPI_Integer,MPI_Sum,MPI_Comm_World,ierr)
     nmap=work
     npart=0
     do i=1,nels
        do j=1,npel
           npart(nodes(i,j))=1
        end do
     end do
     j=1; work=0
     do i=1,nnds
        if (npart(i)==1) then
           work(i)=j; j=j+1
        end if
     end do
     n=sum(npart); allocate(coords(n,dmn))
     j=1
     do i=1,nnds
        if (npart(i)==1) then
           read(10,*)coords(j,:); j=j+1
        else
           read(10,*)val
        end if
     end do
     coords=km2m*coords
     ! Re-number on-rank nodes and create local to global node and dof mappings
     do i=1,nels
        do j=1,npel
           nodes(i,j)=work(nodes(i,j))
        end do
     end do
     n=sum(npart); allocate(nl2g(n,2),indxmap(dmn*n,2))
     j=1
     do i=1,nnds
        if (work(i)==j) then
           nl2g(j,1)=j; nl2g(j,2)=nmap(i); j=j+1
        end if
     end do
     do i=1,n
        do j=1,dmn
           indxmap(dmn*i-j+1,:)=dmn*nl2g(i,:)-j ! 0 based index
        end do
     end do
     deallocate(work)

     deallocate(epart,npart)

     ! Initialize local element variables and global U
     allocate(ipoint(nip,dmn),weight(nip),k(eldof,eldof),m(eldof,eldof),       &
        f(eldof),indx(eldof),enodes(npel),ecoords(npel,dmn),vvec(dmn))
     call SamPts(ipoint,weight)
     n=dmn*nnds 
     call VecCreateMPI(Petsc_Comm_World,Petsc_Decide,n,Vec_U,ierr)
     allocate(stress(nels,nip,cdmn))

     ! Form stiffness matrix
     call PrintMsg("Forming [K] ...")
     nodal_bw=dmn*(nodal_bw+1)
     call MatCreateAIJ(Petsc_Comm_World,Petsc_Decide,Petsc_Decide,n,n,nodal_bw,&
        Petsc_Null_Integer,nodal_bw,Petsc_Null_Integer,Mat_K,ierr)
     call MatSetOption(Mat_K,Mat_New_Nonzero_Allocation_Err,Petsc_False,ierr)
     do i=1,nels
        call FormLocalK(i,k,indx)
        indx=indxmap(indx,2)
        call MatSetValues(Mat_K,eldof,indx,eldof,indx,k,Add_Values,ierr)
     end do
     call MatAssemblyBegin(Mat_K,Mat_Final_Assembly,ierr)
     call MatAssemblyEnd(Mat_K,Mat_Final_Assembly,ierr)

     ! Initialize and form mass matrix and its inverse
     call PrintMsg("Forming [M] & [M]^-1 ...")
     call MatCreateAIJ(Petsc_Comm_World,Petsc_Decide,Petsc_Decide,n,n,1,       &
        Petsc_Null_Integer,0,Petsc_Null_Integer,Mat_M,ierr)
     call MatSetOption(Mat_M,Mat_New_Nonzero_Allocation_Err,Petsc_False,ierr)
     do i=1,nels
        call FormLocalM(i,m,indx)
        indx=indxmap(indx,2)
        do j=1,eldof
           val=m(j,j)
           call MatSetValue(Mat_M,indx(j),indx(j),val,Add_Values,ierr)
        end do
     end do
     call MatAssemblyBegin(Mat_M,Mat_Final_Assembly,ierr)
     call MatAssemblyEnd(Mat_M,Mat_Final_Assembly,ierr)
     call MatDuplicate(Mat_M,Mat_Do_Not_Copy_Values,Mat_Minv,ierr)
     call MatGetDiagonal(Mat_M,Vec_U,ierr) ! Vec_U is used as a work vector
     call VecReciprocal(Vec_U,ierr)
     call MatDiagonalSet(Mat_Minv,Vec_U,Insert_Values,ierr)
     call VecZeroEntries(Vec_U,ierr)
     call MatScale(Mat_M,alpha,ierr)
     allocate(hit(nabs),absel_glb(nabs),absid_glb(nabs),absdir_glb(nabs)) 

     ! Read abs boundary
     hit=0
     do j=1,nabs 
        read(10,*)el,side,dir; el=emap(el)
        if (el/=0) then
           hit(j)=j; absel_glb(j)=el; absid_glb(j)=side; absdir_glb(j)=dir
           call FormLocalAbsC(el,side,m,indx)
           indx=indxmap(indx,2)
           do j1=1,eldof
              do j2=1,eldof 
                 val=m(j1,j2)
                 if (abs(val)>f0) call MatSetValue(Mat_M,indx(j1),indx(j2),val,&
                    Add_Values,ierr)
              end do
           end do
        end if
     end do
     call MatAssemblyBegin(Mat_M,Mat_Final_Assembly,ierr)
     call MatAssemblyEnd(Mat_M,Mat_Final_Assembly,ierr)
     nabs_loc=size(pack(hit,hit>0))
     allocate(absel(nabs_loc),absid(nabs_loc),absdir(nabs_loc))
     absel=absel_glb(pack(hit,hit>0))
     absid=absid_glb(pack(hit,hit>0))
     absdir=absdir_glb(pack(hit,hit>0))
     deallocate(hit,absel_glb,absid_glb,absdir_glb)

     ! Read free surface coordinate/orientation
     allocate(surfloc_glb(nsurf,3),hit(nsurf),surfid_glb(nsurf),               &
        surfmat_glb(nsurf,9),surfel_glb(nsurf),surf_glb(nsurf))
     hit=0; nsurf_loc=0
     allocate(idface(nps))
     do j=1,nsurf ! Read surf surface
        read(10,*)el,side,surf_glb(j); el=emap(el)
        if (el/=0) then ! Local surf element
           ecoords=coords(nodes(el,:),:)
           call Glb2Face(ecoords,side,surfmat_glb(j,:),idface)
           !surfloc_glb(j,:)=(/sum(ecoords(idface,1)),                          &
           !                   sum(ecoords(idface,2)),                          &
           !                   sum(ecoords(idface,3))/)/dble(nps)
           surfloc_glb(j,:)=(/sum(ecoords(:,1)),                               &
                              sum(ecoords(:,2)),                               &
                              sum(ecoords(:,3))/)/dble(npel)
           surfel_glb(j)=el 
           surfid_glb(j)=side
           hit(j)=j
        end if
     end do
     nsurf_loc=size(pack(hit,hit>0))
     allocate(surfloc(nsurf_loc,3)); surfloc=surfloc_glb(pack(hit,hit>0),:)
     allocate(surfid(nsurf_loc)); surfid=surfid_glb(pack(hit,hit>0))
     allocate(surfmat(nsurf_loc,9)); surfmat=surfmat_glb(pack(hit,hit>0),:)
     allocate(surfel(nsurf_loc)); surfel=surfel_glb(pack(hit,hit>0))
     allocate(surf(nsurf_loc)); surf=surf_glb(pack(hit,hit>0))
     deallocate(surfid_glb,surfmat_glb,surfel_glb,surf_glb,hit)
     deallocate(nmap,emap)

     ! Initialize arrays to communicate ghost node values
     call PrintMsg("Setting up solver ...")
     j=size(indxmap,1)
     call VecCreateSeq(Petsc_Comm_Self,j,Seq_U,ierr)
     call ISCreateGeneral(Petsc_Comm_Self,j,indxmap(:,2),Petsc_Copy_Values,    &
        From,ierr)
     call ISCreateGeneral(Petsc_Comm_Self,j,indxmap(:,1),Petsc_Copy_Values,To, &
        ierr)
     call VecScatterCreate(Vec_U,From,Seq_U,To,Scatter,ierr)
     allocate(uu(j),tot_uu(j)); uu=f0; tot_uu=f0

  end if ! Half space FE preparation  

  ! Full space Eshelby's solution
  read(10,*)rstress(:)
  do i=1,nellip 
     read(10,*)ellip(i,:) 
     ellip(i,1:6)=km2m*ellip(i,1:6) ! Centroids and semi-axises
  end do
  do i=1,nrect
      read(10,*)rect(i,:)
      rect(i,1:5)=km2m*rect(i,1:5)
  end do
  do i=1,nobs
     read(10,*)ocoord(i,:) 
     ocoord(i,:)=km2m*ocoord(i,:)
  end do
  close(10) ! End of input
  if (.not. half) then !rank==nprcs-1 .and. 
     call EshSol(mat(1),mat(2),rstress,ellip,ocoord,odat_glb(:,:9))
  end if

  if (half) then 
     call PrintMsg("Numerical half space correction ...")
     call VecGetOwnershipRange(Vec_U,j1,j2,ierr)
     if (rank==nprcs-1) print'(I0,A,I0,A)',j2," dofs on ",nprcs," processors."
     call VecDuplicate(Vec_U,Vec_F,ierr)
     allocate(surfdat(nsurf_loc,18)); surfdat=f0

     ! Find the topography thickness
     if (nsurf_loc>0) then
        val=maxval(surfloc(:,3))
     else
        val=maxval(coords(:,3))
     end if
     call MPI_AllReduce(val,top,1,MPI_Real8,MPI_Max,MPI_Comm_World,ierr)

     ! Full space solution at "free" surface 
     if (nsurf_loc>0) then 
        call EshSol(mat(1),mat(2),rstress,ellip,surfloc,surfdat(:,:9))
        if (nrect>0) then ! Add Okada fault solution
            call OkSol(mat(1),mat(2),rect,surfloc,top,surfdat(:,:9))
        end if
        surfdat(:,10:)=surfdat(:,:9)
     end if
     call FreeSurfTrac

!     allocate(edisp(eldof),estress(nip,6))
     ! Explicit solve with one-step traction
!     steps=int(ceiling(t/dt))
!     call VecDuplicate(Vec_U,Vec_Um,ierr)
!     call VecZeroEntries(Vec_Um,ierr)
!     call VecDuplicate(Vec_U,Vec_Up,ierr)
!     call VecDuplicateVecsF90(Vec_U,6,Vec_W,ierr)
!     do tstep=0,steps
!        call MatMult(Mat_K,Vec_U,Vec_W(1),ierr)
!        call VecAYPX(Vec_W(1),-f1,Vec_F,ierr)
!        call VecScale(Vec_W(1),dt**2,ierr)
!        call VecWAXPY(Vec_W(2),-f1,Vec_Um,Vec_U,ierr)
!        call MatMult(Mat_K,Vec_W(2),Vec_W(3),ierr)
!        call MatMult(Mat_M,Vec_W(2),Vec_W(4),ierr)
!        call VecAXPY(Vec_W(4),beta,Vec_W(3),ierr)
!        call VecScale(Vec_W(4),dt,ierr)
!        call VecWAXPY(Vec_W(5),-f1,Vec_W(4),Vec_W(1),ierr)
!        call MatMult(Mat_Minv,Vec_W(5),Vec_W(6),ierr)
!        call VecWAXPY(Vec_Up,f2,Vec_U,Vec_W(6),ierr)
!        call VecAXPY(Vec_Up,-f1,Vec_Um,ierr)
!        if (tstep==0) call VecZeroEntries(Vec_F,ierr) 
!        if (rank==0) print'(A12,I0,A,I0)'," Time Step ",tstep,'/',steps
!        call VecCopy(Vec_U,Vec_Um,ierr)
!        call VecCopy(Vec_Up,Vec_U,ierr)
!        call GetVec_U; tot_uu=tot_uu+uu
!        !if (frq>0) then
!        !   if (mod(tstep,frq)==0) call WriteOutput
!        !end if
!     end do

     ! Implicit solve
     call KSPCreate(Petsc_Comm_World,Krylov,ierr)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=4)
     call KSPSetOperators(Krylov,Mat_K,Mat_K,Different_Nonzero_Pattern,ierr)
#else
     call KSPSetOperators(Krylov,Mat_K,Mat_K,ierr)
#endif
     call SetupKSPSolver

     call FixBnd("all")
     allocate(resid(max(1,nsurf_loc))); resid=0
     do j=1,nsurf_loc
        call RecoverSurf(j,surfel(j))
     end do
     ! Initial Obs data
     call GetObsNd; allocate(odat(nobs_loc,18)); odat=f0
     if (nobs_loc>0) then
        call EshSol(mat(1),mat(2),rstress,ellip,ocoord_loc,odat(:,:9)) 
        if (nrect>0) then ! Add Okada fault solution
            call OkSol(mat(1),mat(2),rect,ocoord_loc,top,odat(:,:9))
        end if 
        odat(:,10:)=odat(:,:9)
     end if
     ! Initial residual traction
     call MPI_AllReduce(maxval(pack(resid,surf>0)),val,1,MPI_Real8,MPI_Max,    &
        MPI_Comm_World,ierr)
     if (val<tol) go to 8 
     i=0 
     do while(.true. .and. i<ntol+1)
        i=i+1
        call KSPSolve(Krylov,Vec_F,Vec_U,ierr)
        call GetVec_U; tot_uu=tot_uu+uu
        ! Superpose disp/stress and cancel residual traction
        do j=1,nsurf_loc
           call RecoverSurf(j,surfel(j))
        end do
        call ObsEval
        call MPI_AllReduce(maxval(pack(resid,surf>0)),val,1,MPI_Real8,MPI_Max, &
           MPI_Comm_World,ierr)
        if (val<tol) then
           if (rank==0) print('(A,X,I0,X,A,X,ES11.2E3,X,A,X,ES11.2E3)'),       &
               "Step",i,"converge ",val,"<",tol
           go to 8
        else
           if (rank==0) print('(A,X,I0,X,A,X,ES11.2E3,X,A,X,ES11.2E3)'),       &
               "Step",i,"residual traction",val,">",tol 
        end if    
        call FreeSurfTrac
        call FixBnd("rhs")
     end do 
8    call ObsGather
     allocate(surfdat_glb(nsurf,18),surfnrm_glb(nsurf,3))
     call EshGather(surfdat,surfdat_glb)
     call EshGather(surfloc,surfloc_glb)
     call EshGather(surfmat(:,7:9),surfnrm_glb)
  end if ! Half space case
  if (rank==nprcs-1) call EshSave 
9 call PetscFinalize(ierr)

contains

  ! Read simulation parameters
  subroutine ReadParameters
    implicit none
    integer,save :: k=0
    read(10,*)stype
    if (stype=="full") half=.false.
    if (stype=="half") half=.true.
    read(10,*)nellip,nrect,nobs
    if (k==0) allocate(ocoord(nobs,3),ellip(nellip,17),rect(nrect,9))
    if (half) then
       read(10,*)eltype,nodal_bw 
       read(10,*)nels,nnds,nabs,nsurf,frq
       read(10,*)t,dt,alpha,beta,rfrac,tol,ntol
       if (k==0) allocate(odat_glb(nobs,18),mat(3))
    else
       if (k==0) allocate(odat_glb(nobs,9),mat(2))
       nsurf=0
    end if
    odat_glb=f0
    read(10,*)mat(:) ! Material
    k=k+1
  end subroutine ReadParameters

  ! Scatter U and get all local values
  subroutine GetVec_U
    implicit none
    call VecScatterBegin(Scatter,Vec_U,Seq_U,Insert_Values,Scatter_Forward,ierr)
    call VecScatterEnd(Scatter,Vec_U,Seq_U,Insert_Values,Scatter_Forward,ierr)
    call VecGetArrayF90(Seq_U,pntr,ierr)
    uu=pntr
    call VecRestoreArrayF90(Seq_U,pntr,ierr)
  end subroutine GetVec_U

  ! Setup implicit solver
  subroutine SetupKSPSolver
    implicit none
    call KSPSetType(Krylov,"gmres",ierr)
    call KSPGetPC(Krylov,PreCon,ierr)
    call PCSetType(PreCon,"asm",ierr)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=4)
    call KSPSetTolerances(Krylov,1.0D-9,Petsc_Default_Double_Precision,        &
       Petsc_Default_Double_Precision,Petsc_Default_Integer,ierr)
#else
    call KSPSetTolerances(Krylov,1.0D-9,Petsc_Default_Real,Petsc_Default_Real, &
       Petsc_Default_Integer,ierr)
#endif
    call KSPSetFromOptions(Krylov,ierr)
  end subroutine SetupKSPSolver

  ! Evaluation observation disp/stress
  subroutine ObsEval 
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
  end subroutine ObsEval 

  ! Gather obsdat by laster rank
  subroutine ObsGather
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7)
#include "petsc.h"
#endif
    real(8) :: dattmp(nobs*18),buf(nobs*18),ones(nobs,18)
    odat_glb(ol2g,:)=odat
    dattmp=reshape(odat_glb,(/nobs*18/))
    call MPI_Reduce(dattmp,buf,nobs*18,MPI_Real8,MPI_Sum,nprcs-1,              &
       MPI_Comm_World,ierr)
    odat_glb=reshape(buf,(/nobs,18/))
    ones=f0; ones(ol2g,:)=f1
    dattmp=reshape(ones,(/nobs*18/))
    call MPI_Reduce(dattmp,buf,nobs*18,MPI_Real8,MPI_Sum,nprcs-1,              &
       MPI_Comm_World,ierr)
    ones=reshape(buf,(/nobs,18/))
    if (rank==nprcs-1) odat_glb=odat_glb/ones
  end subroutine ObsGather

  ! Gather 2D data by last rank
  subroutine EshGather(dat2D,dest2D)
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7)
#include "petsc.h"
#endif
    integer :: i,j,m,n,ncol
    real(8) :: dat2D(:,:),dest2D(:,:)
    real(8),allocatable :: buf(:)
    j=0
    do i=0,nprcs-2
       m=size(dat2D,1); n=size(dat2D,2); ncol=size(dest2D,2) 
       if (rank==i) call MPI_Send(m*n,1,MPI_Integer,nprcs-1,1234+i,            &
                                  MPI_Comm_World,ierr) 
       if (rank==nprcs-1) then 
          call MPI_Recv(n,1,MPI_Integer,i,1234+i,MPI_Comm_World,status,ierr)
          allocate(buf(n))
       end if
       if (rank==i) call MPI_Send(reshape(dat2D,(/m*n,1/)),m*n,MPI_Real8,      &
                                  nprcs-1,1235+i,MPI_Comm_World,ierr)
       if (rank==nprcs-1) then
          call MPI_Recv(buf,n,MPI_Real8,i,1235+i,MPI_Comm_World,status,ierr)
          n=n/ncol
          dest2D(j+1:j+n,:)=reshape(buf,(/n,ncol/))
          deallocate(buf)
          j=j+n
       end if
    end do
    if (rank==nprcs-1 .and. m>0) then
       n=size(dest2D,1) 
       dest2D(n-m+1:n,:)=dat2D
    end if
  end subroutine EshGather     

  subroutine TestEshInit  
    implicit none
    integer :: i
    real(8) :: vm,a(3),S4(6,6),PIvec(3)
    vm=0.25
    a=(/3.d0,2.d0,1.d0/)
    call EshS4(vm,a,S4,PIvec)
    do i=1,6
       print('(6(F0.4,1X))'), S4(i,:)
    end do
    print('(3(F0.4,1X))'), PIvec
  end subroutine TestEshInit     

  subroutine TestEshD4
    implicit none
    integer :: i,j,k
    real(8) :: vm,a(3),x(3),D4(3,3,3,3),fderphi(3),tderpsi(3,3,3),u(3),eigen(6)
    vm=0.25
    a=(/3.d0,2.d0,1.d0/)
    x=a+f1
    call EshD4(vm,a,x,D4,fderphi,tderpsi)
    do i=1,3
       do j=1,3 
          do k=1,3
             print('(3(F0.4,1X))'), D4(k,:,i,j)
          end do
       end do
    end do
    print('(3(F0.4,1X))'), fderphi
    do i=1,3
       do j=1,3
          print('(3(F0.4,1X))'), tderpsi(j,:,i)
       end do
    end do
    eigen=(/f0,f0,1.D-1,f0,f0,f0/)
    call EshDisp(vm,eigen,fderphi,tderpsi,u)
    print('(3(F0.8,1X))'), u
  end subroutine TestEshD4    

  subroutine TestEshSol
    implicit none
    integer :: i
    real(8) :: Em,vm,ellip(2,17),ocoord(2,3),sol(2,9),stress(6)
    Em=f8; vm=f1/f4; stress=f0
    ellip=reshape((/f0,f0,f0,f3,f2,f1,pi/f2,pi/f3,pi/f4,Em,vm,f0,f0,f1/f8,f0,  &
                    f0,f0,f0,f0,f0,f3,f2,f1,pi/f2,pi/f3,pi/f4,Em,vm,f0,f0,     &
                    f0*f1/f8,f0,f0,f0/),(/2,17/),(/f0,f0/),(/2,1/))                                                      
    ocoord=reshape((/f1,f0,f2,f2,f0,f0/),(/2,3/),(/f0,f0/),(/2,1/))
    call EshSol(Em,vm,stress,ellip,ocoord,sol)
    do i=1,size(sol,1) 
       print('(9(F0.8,1X))'),sol(i,:)
    end do
  end subroutine TestEshSol

end program main
