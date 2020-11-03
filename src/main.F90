! Copyright (C) 2010-2015 Tabrez Ali 2015-present Chunfang Meng. All rights
! reserved.  This file is part of Esh3D. See ../COPYING for license information.

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
  real(8) :: t1,t2
  integer :: n,i,j,j1,j2,n_incl,nodal_bw,incl_bw,status(MPI_STATUS_SIZE)
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

  t1=MPI_Wtime()
  call PrintMsg("Reading input ...")
  call ReadParameters

  if (half .or. fini) then ! Half or finite space problem
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
           call METIS_PartMeshNodal(nels,nnds,work,nodes,null_i,null_i,nprcs,  &
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
     n=sum(npart); allocate(coords(n,dmn),bc(n,dmn))
     j=1
     do i=1,nnds
        if (npart(i)==1) then
           read(10,*)coords(j,:),bc(j,:)
           j=j+1
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
     allocate(ipoint(nip,dmn),weight(nip),k(eldof,eldof),f(eldof),indx(eldof), &
        enodes(npel),ecoords(npel,dmn),vvec(dmn))
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
     if (fini) call FormMatKfull

     ! Read traction surface (surfel, surfside, surftrc, surf)
     allocate(surfel_glb(ntrc),surfside_glb(ntrc),surftrc_glb(ntrc,3),         &
        surf_glb(ntrc),surfloc_glb(ntrc,3),surfmat_glb(ntrc,9),hit(ntrc))
     hit=0; ntrc_loc=0
     allocate(idface(nps))
     do j=1,ntrc
        read(10,*)el,side,surftrc_glb(j,:),surf_glb(j); el=emap(el)
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
           surfside_glb(j)=side
           hit(j)=j
        end if
     end do
     ntrc_loc=size(pack(hit,hit>0))
     allocate(surfel(ntrc_loc)); surfel=surfel_glb(pack(hit,hit>0))
     allocate(surfside(ntrc_loc)); surfside=surfside_glb(pack(hit,hit>0))
     allocate(surftrc(ntrc_loc,3)); surftrc=surftrc_glb(pack(hit,hit>0),:)
     allocate(surf(ntrc_loc)); surf=surf_glb(pack(hit,hit>0))
     allocate(surfloc(ntrc_loc,3)); surfloc=surfloc_glb(pack(hit,hit>0),:)
     allocate(surfmat(ntrc_loc,9)); surfmat=surfmat_glb(pack(hit,hit>0),:)
     deallocate(surfel_glb,surfside_glb,surftrc_glb,surf_glb,surfmat_glb,hit)

     if (fini) then ! Nodes associated with fixed boundary bc = 0
        n=size(coords,1)
        allocate(ndfix_glb(n),bcfix_glb(n,3))
        ndfix_glb=0
        do i=1,n
           if (sum(bc(i,:))<3) then
              ndfix_glb(i)=i
              bcfix_glb(i,:)=bc(i,:)
           end if
        end do
        nfix=size(pack(ndfix_glb,ndfix_glb>0))
        allocate(ndfix(nfix),bcfix(nfix,3),solfix(nfix,9)); solfix=f0
        ndfix=ndfix_glb(pack(ndfix_glb,ndfix_glb>0))
        bcfix=bc(ndfix,:)
        deallocate(ndfix_glb,bcfix_glb)
     end if

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
     allocate(uu(j),uu0(j)); uu=f0; uu0=f0
  end if ! Half or finite space FE preparation

  if (inho) then
     read(10,*)rstress(:)
  elseif (incl) then
     rstress=f0
  end if
  do i=1,nellip
     if (incl) then
        if (i<=nsolid) then
           read(10,*)ellip(i,:9),ellip(i,12:)
           ellip(i,10:11)=mat
        else ! Inclusions can only be solid.
           read(10,*)val
        end if
     elseif (inho) then
        if (i<=nsolid) then
           read(10,*)ellip(i,:)
        else
           read(10,*)ellip(i,:10),ellip(i,12)
           ellip(i,11)=-f1 ! Negative "Poisson's ratio" indicates fluid.
           ellip(i,13:14)=f0; ellip(i,15:)=f0
        end if
     end if
     ellip(i,1:6)=km2m*ellip(i,1:6) ! Centroids and semi-axises
  end do

  ! Translate eigenstrain from inclusion coordinate to global coordinate
  call EigIncl2Glb(ellip)
  ! ellipeff with evolving eigen sources
  allocate(ellipeff(nellip,17)); ellipeff=ellip
  do i=1,nrect
      read(10,*)rect(i,:)
      rect(i,1:5)=km2m*rect(i,1:5)
  end do
  ! Okada solution at inclusion centroids
  if (inho .and. (half .or. fini) .and. nrect>0) allocate(solok(nellip,9))
  do i=1,nobs
     read(10,*)ocoord(i,:)
     ocoord(i,:)=km2m*ocoord(i,:)
  end do
  close(10) ! End of input

  ! Global linear system matrix for interacting inhomogeneities
  if (inho) then
     n=nellip*6
     incl_bw=n
     allocate(Keig(n,n),Feig(n))!,EffEig(nellip,6))
     call EshKeig(mat(1),mat(2),ellip,Keig)
     ! LAPACK inversion
     !allocate(KeigInv(n,n))
     !call MatInv(Keig,KeigInv)
     ! Entries of one stress tensor are not split by different ranks
     call VecCreateMPI(Petsc_Comm_World,Petsc_Decide,nellip,Vec_incl,ierr)
     call VecGetLocalSize(Vec_incl,n_incl,ierr)
     call VecDestroy(Vec_incl,ierr)
     call MatCreateAIJ(Petsc_Comm_World,n_incl*6,n_incl*6,n,n,incl_bw,         &
        Petsc_Null_Integer,incl_bw,Petsc_Null_Integer,Mat_Keig,ierr)
     call MatSetOption(Mat_Keig,Mat_New_Nonzero_Allocation_Err,Petsc_False,    &
        ierr)
     !call MatCreateDense(Petsc_Comm_World,n_incl*6,n_incl*6,n,n,               &
     !   Petsc_Null_Scalar,Mat_Keig,ierr)
     ! From [Keig]
     if (rank==nprcs-1) call MatSetValues(Mat_Keig,n,(/(i,i=0,n-1)/),n,        &
        (/(i,i=0,n-1)/),transpose(Keig),Insert_Values,ierr)
     call MatAssemblyBegin(Mat_Keig,Mat_Final_Assembly,ierr)
     call MatAssemblyEnd(Mat_Keig,Mat_Final_Assembly,ierr)
     call KSPCreate(Petsc_Comm_World,KryInc,ierr)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=4)
     call KSPSetOperators(KryInc,Mat_Keig,Mat_Keig,Different_Nonzero_Pattern,  &
        ierr)
#else
     call KSPSetOperators(KryInc,Mat_Keig,Mat_Keig,ierr)
#endif
     !call SetupKSPSolver(KryInc)
     call VecCreateMPI(Petsc_Comm_World,n_incl*6,n,Vec_Feig,ierr)
     call VecDuplicate(Vec_Feig,Vec_Eig,ierr)
     if (nellip-nsolid>0) then ! Has fluid inclusions
        call EshKeig(mat(1),mat(2),ellip,Keig,Kfluid=.true.)
        !allocate(KeigFldInv(n,n))
        !call MatInv(Keig,KeigFldInv)
        call MatDuplicate(Mat_Keig,Mat_Do_Not_Copy_Values,Mat_Kfld,ierr)
        if (rank==nprcs-1) call MatSetValues(Mat_Kfld,n,(/(i,i=0,n-1)/),n,     &
           (/(i,i=0,n-1)/),transpose(Keig),Insert_Values,ierr)
        call MatAssemblyBegin(Mat_Kfld,Mat_Final_Assembly,ierr)
        call MatAssemblyEnd(Mat_Kfld,Mat_Final_Assembly,ierr)
        call KSPCreate(Petsc_Comm_World,KryFld,ierr)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=4)
        call KSPSetOperators(KryFld,Mat_Kfld,Mat_Kfld,                         &
           Different_Nonzero_Pattern,ierr)
#else
        call KSPSetOperators(KryFld,Mat_Kfld,Mat_Kfld,ierr)
#endif
        !call SetupKSPSolver(KryFld)
        nfluid=(nellip-nsolid)
        n=6*nfluid; incl_bw=n
        allocate(Kvol(n,n),Fvol(n))
        call EshKvol(mat(1),mat(2),ellip(nsolid+1:,:),Kvol)
        ! Entries of one stress tensor are not split by different ranks
        call VecCreateMPI(Petsc_Comm_World,Petsc_Decide,nfluid,Vec_incl,ierr)
        call VecGetLocalSize(Vec_incl,n_incl,ierr)
        call VecDestroy(Vec_incl,ierr)
        call MatCreateAIJ(Petsc_Comm_World,n_incl*6,n_incl*6,n,n,incl_bw,      &
           Petsc_Null_Integer,incl_bw,Petsc_Null_Integer,Mat_Kvol,ierr)
        call MatSetOption(Mat_Kvol,Mat_New_Nonzero_Allocation_Err,             &
           Petsc_False,ierr)
        if (rank==nprcs-1) call MatSetValues(Mat_Kvol,n,(/(i,i=0,n-1)/),n,     &
           (/(i,i=0,n-1)/),transpose(Kvol),Insert_Values,ierr)
        call MatAssemblyBegin(Mat_Kvol,Mat_Final_Assembly,ierr)
        call MatAssemblyEnd(Mat_Kvol,Mat_Final_Assembly,ierr)
        call KSPCreate(Petsc_Comm_World,KryVol,ierr)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=4)
        call KSPSetOperators(KryVol,Mat_Kvol,Mat_Kvol,                         &
           Different_Nonzero_Pattern,ierr)
#else
        call KSPSetOperators(KryVol,Mat_Kvol,Mat_Kvol,ierr)
#endif
        !call SetupKSPSolver(KryVol)
        call VecCreateMPI(Petsc_Comm_World,n_incl*6,n,Vec_Fvol,ierr)
        call VecDuplicate(Vec_Fvol,Vec_Evol,ierr)
        if (nsolid>0) then ! Mixed solid fluid inclusions
           ! [V] = -[W][E]
           allocate(Wsec(n,6*nsolid),Esec(6*nsolid),Vsec(n))
           call EshWsec(mat(1),mat(2),ellip,Wsec)
           call VecDuplicate(Vec_Feig,Vec_dEig,ierr)
        end if
     end if
  end if

  ! Full space Eshelby's solution
  if (full) then
     if (inho) then ! Interacting eigenstrain at inhomogeneity centroids
        do i=1,nellip
           instress(i,:)=rstress
        end do
        ! Form [Feig]
        call EshFeig(mat(1),mat(2),instress,ellip,Feig,init=.true.)
        ! Non-Interacting effective eigenstrain
        !call EshEffEig(mat(1),mat(2),instress,ellip,EffEig,init=.true.)
        if (rank==nprcs-1) call VecSetValues(Vec_Feig,nellip*6,                &
           (/(i,i=0,nellip*6-1)/),Feig,Insert_Values,ierr)
        call VecAssemblyBegin(Vec_Feig,ierr)
        call VecAssemblyEnd(Vec_Feig,ierr)
        call KSPSolve(KryInc,Vec_Feig,Vec_Eig,ierr)
        call UpInhoEigen(ellipeff(:,12:17))
        if (nfluid>0) then ! Has fluid inclusion
           call GetFvol(ellipeff(nsolid+1:,:),Fvol,einit=ellip(nsolid+1:,12))
           if (rank==nprcs-1) call VecSetValues(Vec_Fvol,nfluid*6,             &
              (/(i,i=0,nfluid*6-1)/),Fvol,Insert_Values,ierr)
           call VecAssemblyBegin(Vec_Fvol,ierr)
           call VecAssemblyEnd(Vec_Fvol,ierr)
           call KSPSolve(KryVol,Vec_Fvol,Vec_Evol,ierr)
           call Evol2Feig(mat(2),ellip) ! Intrinsic fluid eigenstrains to RHS
           call KSPSolve(KryFld,Vec_Feig,Vec_Eig,ierr)
           if (nsolid>0) then ! Secodnary interaction
              j=0; call VecCopy(Vec_Eig,Vec_dEig,ierr)
              do while(j<nrtol)
                 call CoupleFSF(val) ! Fluid -> solid -> fluid coupling
                 if (val<rtol) then
                    if (rank==0) print('(A,X,I0,X,A,X,ES11.2E3,X,A,X,          &
                       &ES11.2E3)'),"Sub step",j+1,"converge",val,"<",rtol
                    exit
                 else
                    if (rank==0) print('(A,X,I0,X,A,X,ES11.2E3,X,A,X,          &
                       &ES11.2E3)'),"Sub step",j+1,"residual",val,">",rtol
                 end if
                 j=j+1
              end do
           end if
           call UpInhoEigen(ellipeff(:,12:17),fluid=.true.) ! Superpose
        end if
        call EshIncSol(mat(1),mat(2),ellipeff,ocoord,odat_glb(:,:9))
     else
        call EshIncSol(mat(1),mat(2),ellipeff(:nsolid,:),ocoord,odat_glb(:,:9))
     end if

  ! Tuncated space involving numerical boundary condition solutions.
  elseif (half .or. fini) then
     t2=MPI_Wtime(); if (rank==0) print'(F0.2,A)',t2-t1," seconds to assemble."
     call VecGetOwnershipRange(Vec_U,j1,j2,ierr)
     if (rank==nprcs-1) print'(I0,A,I0,A,I0,A)',j2-j1,"/",j2," DoFs across ",  &
        nprcs," processors."
     allocate(surfdat(ntrc_loc,18)); surfdat=f0

     ! Find the topography thickness
     if (ntrc_loc>0) then
        val=maxval(surfloc(:,3))
     else
        val=maxval(coords(:,3))
     end if
     call MPI_AllReduce(val,top,1,MPI_Real8,MPI_Max,MPI_Comm_World,ierr)

     ! Implicit solver
     call KSPCreate(Petsc_Comm_World,Krylov,ierr)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=4)
     call KSPSetOperators(Krylov,Mat_K,Mat_K,Different_Nonzero_Pattern,ierr)
#else
     call KSPSetOperators(Krylov,Mat_K,Mat_K,ierr)
#endif
     call SetupKSPSolver(Krylov)
     call VecDuplicate(Vec_U,Vec_F,ierr); call VecZeroEntries(Vec_F,ierr)
     if (fini) then ! Coefficient of fixed dofs
        call VecDuplicate(Vec_U,Vec_FixC,ierr)
        call GetVecFixC
        call VecDuplicate(Vec_U,Vec_Fix,ierr)
        call VecZeroEntries(Vec_Fix,ierr)
        call VecDuplicate(Vec_U,Vec_FixF,ierr)
     end if
     ! Solve initial boundary problem for nonzero surface loading
     allocate(resid(ntrc_loc)) ! Traction magnitudes
     resid=sqrt(surftrc(:,1)**2+surftrc(:,2)**2+surftrc(:,3)**2)
     call MPI_AllReduce(maxval(resid),val,1,MPI_Real8,MPI_Max,MPI_Comm_World,  &
        ierr)
     if (val>tol .and. inho) then
        call PrintMsg("Solving initial boundary condition problem ...")
        ! Apply background traction
        do i=1,ntrc_loc
           el=surfel(i); side=surfside(i); vvec=surftrc(i,:)
           call ApplyTraction(el,side,vvec)
        end do
        call VecAssemblyBegin(Vec_F,ierr)
        call VecAssemblyEnd(Vec_F,ierr)
        t1=MPI_Wtime()
        call KSPSolve(Krylov,Vec_F,Vec_U,ierr)
        t2=MPI_Wtime()
        if (rank==0) print'(F0.2,A)',t2-t1," seconds to converge."
        call GetVec_U; uu0=uu
     end if

     ! Interacting eigenstrain at inhomogeneity centroids, ellipeff(:,12:17)
     if (inho) then
        ! Stress perturbation at inclusion centroids
        call GetObsNd("in")
        if (nrect>0) then ! Okada stress at
           call OkSol(mat(1),mat(2),rect,ellip(:,:3),top,solok)
           call InStrEval(ok=.true.) ! uu (solok) -> instress
        else
           call InStrEval
        end if
        do i=1,nellip ! Add remote stress
           instress(i,:)=instress(i,:)+rstress
        end do
        call EshFeig(mat(1),mat(2),instress,ellip,Feig,init=.true.)
        !call EshEffEig(mat(1),mat(2),instress,ellip,EffEig,init=.true.)
        if (rank==nprcs-1) call VecSetValues(Vec_Feig,nellip*6,                &
           (/(i,i=0,nellip*6-1)/),Feig,Insert_Values,ierr)
        call VecAssemblyBegin(Vec_Feig,ierr)
        call VecAssemblyEnd(Vec_Feig,ierr)
        call KSPSolve(KryInc,Vec_Feig,Vec_Eig,ierr)
        call UpInhoEigen(ellipeff(:,12:17))
        if (nfluid>0) then ! Has fluid inclusion
           call GetFvol(ellipeff(nsolid+1:,:),Fvol,einit=ellip(nsolid+1:,12))
           if (rank==nprcs-1) call VecSetValues(Vec_Fvol,nfluid*6,             &
              (/(i,i=0,nfluid*6-1)/),Fvol,Insert_Values,ierr)
           call VecAssemblyBegin(Vec_Fvol,ierr)
           call VecAssemblyEnd(Vec_Fvol,ierr)
           call KSPSolve(KryVol,Vec_Fvol,Vec_Evol,ierr)
           call Evol2Feig(mat(2),ellip) ! Intrinsic fluid eigenstrains to RHS
           call KSPSolve(KryFld,Vec_Feig,Vec_Eig,ierr)
           if (nsolid>0) then ! Secodnary interaction
              j=0; call VecCopy(Vec_Eig,Vec_dEig,ierr)
              do while(j<nrtol)
                 call CoupleFSF(val) ! Fluid -> solid -> fluid coupling
                 if (val<rtol) then
                    if (rank==0) print('(A,X,I0,X,A,X,ES11.2E3,X,A,X,          &
                       &ES11.2E3)'),"Sub step",j+1,"converge",val,"<",rtol
                    exit
                 else
                    if (rank==0) print('(A,X,I0,X,A,X,ES11.2E3,X,A,X,          &
                       &ES11.2E3)'),"Sub step",j+1,"residual",val,">",rtol
                 end if
                 j=j+1
              end do
           end if
           call UpInhoEigen(ellipeff(:,12:17),fluid=.true.)
        end if
     end if

     ! Full space solution at traction surface
     if (ntrc_loc>0) then
        if (incl) then ! Inclusions have to be solid
           call EshIncSol(mat(1),mat(2),ellipeff,surfloc(:nsolid,:),           &
              surfdat(:,:9))
        else
           call EshIncSol(mat(1),mat(2),ellipeff,surfloc,surfdat(:,:9))
        end if
        ! Add Okada fault solution
        if (nrect>0) call OkSol(mat(1),mat(2),rect,surfloc,top,surfdat(:,:9))
        surfdat(:,10:)=surfdat(:,:9)
     end if
     ! Full space solution at fix boundary
     if (fini) then
        if (incl) then
           call EshIncSol(mat(1),mat(2),ellipeff(:nsolid,:),coords(ndfix,:),   &
              solfix)
        else
           call EshIncSol(mat(1),mat(2),ellipeff,coords(ndfix,:),solfix)
        end if
        if (nrect>0) call OkSol(mat(1),mat(2),rect,coords(ndfix,:),top,solfix)
     end if
     ! Full space solution at observation
     call GetObsNd("ob"); allocate(odat(nobs_loc,18)); odat=f0
     if (nobs_loc>0) then
        if (incl) then
           call EshIncSol(mat(1),mat(2),ellipeff(:nsolid,:),ocoord_loc,        &
              odat(:,:9))
        else
           call EshIncSol(mat(1),mat(2),ellipeff,ocoord_loc,odat(:,:9))
        end if
        if (nrect>0) call OkSol(mat(1),mat(2),rect,ocoord_loc,top,odat(:,:9))
        odat(:,10:)=odat(:,:9)
     end if

     ! Superpose, and estimate residual traction
     call SurfSupResid
     if (fini) call FixSup
     call ObsSup
     call MPI_AllReduce(maxval(pack(resid,surf>0)),val,1,MPI_Real8,MPI_Max,    &
        MPI_Comm_World,ierr)
     if (val<tol) go to 8
     if (rank==0) print('(A,X,ES11.2E3,X,A,X,ES11.2E3,A)'),                    &
        "Step 0 residual traction",val,">",tol,", run correction ..."
     call MatchSurf ! Cancel residual traction/displacement

     i=0; t1=MPI_Wtime() ! Half/finite space correction
     do while(i<ntol)
        call KSPSolve(Krylov,Vec_F,Vec_U,ierr)
        call GetVec_U
        if (inho) then ! Surface <-> inclusion interaction
           call InStrEval
           call EshFeig(mat(1),mat(2),instress,ellip,Feig)
           !call EshEffEig(mat(1),mat(2),instress,ellip,EffEig)
           if (rank==nprcs-1) call VecSetValues(Vec_Feig,nellip*6,             &
              (/(i,i=0,nellip*6-1)/),Feig,Insert_Values,ierr)
           call VecAssemblyBegin(Vec_Feig,ierr)
           call VecAssemblyEnd(Vec_Feig,ierr)
           call KSPSolve(KryInc,Vec_Feig,Vec_Eig,ierr)
           call UpInhoEigen(ellipeff(:,12:17))
           if (nfluid>0) then ! Has fluid inclusion
              call GetFvol(ellipeff(nsolid+1:,:),Fvol)
              if (rank==nprcs-1) call VecSetValues(Vec_Fvol,nfluid*6,          &
                 (/(i,i=0,nfluid*6-1)/),Fvol,Insert_Values,ierr)
              call VecAssemblyBegin(Vec_Fvol,ierr)
              call VecAssemblyEnd(Vec_Fvol,ierr)
              call KSPSolve(KryVol,Vec_Fvol,Vec_Evol,ierr)
              call Evol2Feig(mat(2),ellip) ! Intrinsic fluid eigenstrains to RHS
              call KSPSolve(KryFld,Vec_Feig,Vec_Eig,ierr)
              if (nsolid>0) then ! Secodnary interaction
                 j=0; call VecCopy(Vec_Eig,Vec_dEig,ierr)
                 do while(j<nrtol)
                    call CoupleFSF(val) ! Fluid -> solid -> fluid coupling
                    if (val<rtol) then
                       if (rank==0) print('(A,X,I0,X,A,X,ES11.2E3,X,A,X,       &
                          &ES11.2E3)'),"Sub step",j+1,"converge",val,"<",rtol
                       exit
                    else
                       if (rank==0) print('(A,X,I0,X,A,X,ES11.2E3,X,A,X,       &
                          &ES11.2E3)'),"Sub step",j+1,"residual",val,">",rtol
                    end if
                    j=j+1
                 end do
              end if
              call UpInhoEigen(ellipeff(:,12:17),fluid=.true.)
           end if
           call EshIncSol(mat(1),mat(2),ellipeff,surfloc,surfdat(:,10:))
           if (nobs_loc>0) call EshIncSol(mat(1),mat(2),ellipeff,ocoord_loc,   &
              odat(:,10:))
           if (fini) call EshIncSol(mat(1),mat(2),ellipeff,coords(ndfix,:),    &
              solfix)
        end if

        ! Superpose disp/stress to cancel residual traction
        call SurfSupResid
        if (fini) call FixSup
        call ObsSup
        call MPI_AllReduce(maxval(pack(resid,surf>0)),val,1,MPI_Real8,MPI_Max, &
           MPI_Comm_World,ierr)
        if (val<tol) then
           t2=MPI_Wtime()
           if (rank==0) print('(A,X,I0,X,A,X,ES11.2E3,X,A,X,ES11.2E3,X,A,X,    &
              &F0.2,X,A)'),"Step",i+1,"converge",val,"<",tol,"at",t2-t1,       &
              &"seconds."
           go to 8
        else
           if (rank==0) print('(A,X,I0,X,A,X,ES11.2E3,X,A,X,ES11.2E3)'),       &
               "Step",i+1,"residual traction",val,">",tol
        end if
        call MatchSurf
        i=i+1
     end do
     t2=MPI_Wtime()
     if (rank==0) print'(F0.2,A)',t2-t1," seconds to run correction."
8    call ObsGather
     allocate(surfdat_glb(ntrc,18),surfnrm_glb(ntrc,3))
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
    full=.false.;half=.false.;fini=.false.;incl=.false.;inho=.false.
    if (index(stype,"full")/=0) full=.true.
    if (index(stype,"half")/=0) half=.true.
    if (index(stype,"fini")/=0) fini=.true.
    if (index(stype,"incl")/=0) incl=.true.
    if (index(stype,"inho")/=0) inho=.true.
    if (full) then
       read(10,*)nellip,nsolid,nobs
    else
       read(10,*)nellip,nsolid,nrect,nobs
    end if
    ! Fluid-solid iteration
    if (nsolid>0 .and. nellip>nsolid) read(10,*)rtol,nrtol
    if (k==0) then
       allocate(ocoord(nobs,3))
       allocate(ellip(nellip,17),instress(nellip,6)); instress=f0
       if (half .or. fini) allocate(rect(nrect,9))
    end if
    if (half .or. fini) then
       read(10,*)eltype,nodal_bw
       read(10,*)nels,nnds,ntrc
       read(10,*)tol,ntol
       if (k==0) allocate(odat_glb(nobs,18))
    else
       if (k==0) allocate(odat_glb(nobs,9))
       ntrc=0; nrect=0
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
  subroutine SetupKSPSolver(Krylov)
    implicit none
    KSP :: Krylov
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

  ! Gather obsdat by last rank
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
    if (rank==nprcs-1) odat_glb=odat_glb/ones ! Scale duplicates, NaN allowed
    call MPI_Barrier(MPI_Comm_World,ierr) ! Prevent crash
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

  ! Fluid -> solid -> fluid coupling [dEig] => [Eig],resid
  subroutine CoupleFSF(resid)
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7)
#include "petsc.h"
#endif
    real(8) :: resid
    call GetEigSec(Esec) ! Vec_dEig => Esec
    if (rank==nprcs-1) then
       Vsec=-matmul(Wsec,Esec) ! nsolid -> nfluid
       call GetSecFvol(Vsec,Fvol) ! Vsec => Fvol
       call VecSetValues(Vec_Fvol,nfluid*6,(/(i,i=0,nfluid*6-1)/),Fvol,        &
          Insert_Values,ierr)
    end if
    call VecAssemblyBegin(Vec_Fvol,ierr)
    call VecAssemblyEnd(Vec_Fvol,ierr)
    call KSPSolve(KryVol,Vec_Fvol,Vec_Evol,ierr)
    call Evol2Feig(mat(2),ellip)
    call KSPSolve(KryFld,Vec_Feig,Vec_dEig,ierr)
    call VecAXPY(Vec_Eig,f1,Vec_dEig,ierr)
    call ConvergeL2(Vec_Eig,Vec_dEig,6,resid)
  end subroutine CoupleFSF

end program main
