! Copyright (C) 2010-2015 Tabrez Ali 2015-present Chunfang Meng. All rights
! reserved. This file is part of Esh3D. See ../COPYING for license information.

module local

  use elems
  implicit none
  real(8),allocatable :: ipoint(:,:),weight(:)

contains

  ! Form element [K]
  subroutine FormElK(ecoords,E,nu,k)
    implicit none
    real(8) :: ecoords(npel,dmn),E,nu,k(eldof,eldof),D(cdmn,cdmn),dN(dmn,npel),&
       detj,B(cdmn,eldof)
    integer :: i
    k=f0
    call DMat(D,E,nu)
    do i=1,nip
       call FormdNdetJ(ipoint(i,:),ecoords,dN,detj)
       call BMat(dN,B)
       k=k+matmul(transpose(B),matmul(D,B))*weight(i)*detj
    end do
  end subroutine FormElK

  ! Form element [M]
  subroutine FormElM(ecoords,rho,m)
    implicit none
    real(8) :: ecoords(npel,dmn),rho,m(eldof,eldof),ms(npel,npel),N(1,npel),   &
       detj,val
    integer :: i
    m=f0; ms=f0
    do i=1,nip
       call FormdetJ(ipoint(i,:),ecoords,detj)
       call ShapeFunc(N(1,:),ipoint(i,:))
       ms=ms+rho*matmul(transpose(N),N)*weight(i)*detj
    end do
    do i=1,dmn
       m(i::dmn,i::dmn)=ms
    end do
    do i=1,eldof
       val=sum(m(i,:)); m(i,:)=f0; m(i,i)=val ! Diagonalize [M]
    end do
  end subroutine FormElM

  ! Form element index
  subroutine FormElIndx(enodes,indx)
    implicit none
    integer :: enodes(npel),indx(eldof),j,j1
    do j=1,npel
       do j1=1,dmn
          indx(dmn*j-j1+1)=dmn*enodes(j)-j1+1
       end do
    end do
  end subroutine FormElIndx

  ! Calculate element stress
  subroutine CalcElStress(ecoords,edisp,E,nu,estress)
    implicit none
    real(8) :: ecoords(npel,dmn),edisp(eldof),E,nu,estrain(nip,cdmn),          &
       estress(nip,cdmn),D(cdmn,cdmn),dN(dmn,npel),detj,B(cdmn,eldof)
    integer :: i
    call DMat(D,E,nu)
    do i=1,nip
       call FormdNdetJ(ipoint(i,:),ecoords,dN,detj)
       call BMat(dN,B)
       estrain(i,:)=matmul(B,edisp); estress(i,:)=matmul(D,estrain(i,:))
    end do
  end subroutine CalcElStress

  ! Computes the strain displacement matrix 'B'
  subroutine BMat(dN,B)
    implicit none
    real(8) :: dN(dmn,npel),B(cdmn,eldof)
    integer :: j
    B=f0
    select case(dmn)
    case(2)
       do j=1,npel
          B(1,2*j-1)=dN(1,j); B(1,2*j)=f0
          B(2,2*j-1)=f0     ; B(2,2*j)=dN(2,j)
          B(3,2*j-1)=dN(2,j); B(3,2*j)=dN(1,j)
       end do
    case(3)
       do j=1,npel
          B(1,3*j-2)=dN(1,j); B(1,3*j-1)=f0     ; B(1,3*j)=f0
          B(2,3*j-2)=f0     ; B(2,3*j-1)=dN(2,j); B(2,3*j)=f0
          B(3,3*j-2)=f0     ; B(3,3*j-1)=f0     ; B(3,3*j)=dN(3,j)
          B(4,3*j-2)=dN(2,j); B(4,3*j-1)=dN(1,j); B(4,3*j)=f0
          B(5,3*j-2)=f0     ; B(5,3*j-1)=dN(3,j); B(5,3*j)=dN(2,j)
          B(6,3*j-2)=dN(3,j); B(6,3*j-1)=f0     ; B(6,3*j)=dN(1,j)
       end do
    end select
  end subroutine BMat

  ! Computes dN/dx(yz) and the determinant of Jacobian 'detj'
  subroutine FormdNdetJ(ipcoord,ecoords,dN,detj)
    implicit none
    real(8) :: ipcoord(dmn),detj,ecoords(npel,dmn),dN(dmn,npel),               &
       jacob(dmn,dmn),invj(dmn,dmn)
    call ShapeFuncd(dN,ipcoord) ! Form dN/de, dN/dn, (dN/ds)
    jacob=matmul(dN,ecoords)
    call MatDet(jacob,detj)
    call MatInv(jacob,invj)
    dN=matmul(invj,dN) ! Form dN/dx, dN/dy, (dN/dz)
  end subroutine FormdNdetJ

  ! Computes only the determinant of Jacobian 'detj' for use in Hs
  subroutine FormdetJ(ipcoord,ecoords,detj)
    implicit none
    real(8) :: ipcoord(dmn),detj,ecoords(npel,dmn),dN(dmn,npel),               &
       jacob(dmn,dmn)
    call ShapeFuncd(dN,ipcoord) ! Form dN/de, dN/dn, (dN/ds)
    jacob=matmul(dN,ecoords)
    call MatDet(jacob,detj)
  end subroutine FormdetJ

  ! Computes D, the matrix of elastic properties
  subroutine DMat(D,E,nu)
    implicit none
    real(8) :: D(:,:),E,nu
    if (size(D,1)==3) call DMat2d(D,E,nu)
    if (size(D,1)==6) call DMat3d(D,E,nu)
  end subroutine DMat

  ! Computes D, the matrix of elastic properties
  subroutine DMat2d(D,E,nu)
    implicit none
    real(8) :: D(3,3),E,nu,c
    c=E/((f1+nu)*(f1-f2*nu))
    D=c*reshape((/(f1-nu),nu,f0,nu,(f1-nu),f0,f0,f0,(f1-f2*nu)/f2/),           &
       shape=(/3,3/))
  end subroutine DMat2d

  ! Computes D, the matrix of elastic properties
  subroutine DMat3d(D,E,nu)
    implicit none
    real(8) :: D(6,6),E,nu,c
    c=E/((f1+nu)*(f1-f2*nu))
    D=c*reshape((/(f1-nu),nu,nu,f0,f0,f0,nu,(f1-nu),nu,f0,f0,f0,nu,nu,(f1-nu), &
       f0,f0,f0,f0,f0,f0,(f1-f2*nu)/f2,f0,f0,f0,f0,f0,f0,(f1-f2*nu)/f2,f0,f0,  &
       f0,f0,f0,f0,(f1-f2*nu)/f2/),shape=(/6,6/))
  end subroutine DMat3d

end module local
