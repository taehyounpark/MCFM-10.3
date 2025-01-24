!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function BoundaryConditionQQBARNS(x)
!     Interpolation routine for Boundary Condition data
!     taken from ancillary files at https://arxiv.org/abs/2207.07037
      implicit none
      include 'types.f'
      include 'first.f'
      integer,parameter::n=401
      real(dp)::BoundaryConditionQQBARNS,x,y
      real(dp),save::xa(n),ya(n),y2a(n)
      real(dp),parameter::yp1=1d30,ypn=1d30
!omp threadprivate(xa,ya,y2a)
      integer::j
      if (first) then
!$omp critical
      open(unit=10,file='../src/ptveto/BCxQQBARNS.dat',status='old') 
      do j=1,n
      read(10,*) XA(j),YA(j)
      enddo   
      close(10)
!$omp end critical
      first=.false.
      call spline(xa,ya,n,yp1,ypn,y2a)
      endif
      call splint(xa,ya,y2a,n,x,y)
      BoundaryConditionQQBARNS=y
      return
      end
      
