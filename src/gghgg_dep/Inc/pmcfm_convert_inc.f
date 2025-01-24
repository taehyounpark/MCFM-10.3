!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      implicit none
      include 'nlegborn.h'
      integer nlegs
      parameter (nlegs=nlegreal)
      real(dp):: pmcfm(mxpart,4)
      real(dp):: p(0:3,nlegs)
      integer nu

      pmcfm(:,:)=0

      do nu = 1,3
        pmcfm(1,nu)=-p(nu,1)
        pmcfm(2,nu)=-p(nu,2)
        pmcfm(3,nu)=p(nu,4)
        pmcfm(4,nu)=p(nu,5)
        pmcfm(5,nu)=p(nu,3)
      enddo

      pmcfm(1,4)=-p(0,1)
      pmcfm(2,4)=-p(0,2)
      pmcfm(3,4)=p(0,4)
      pmcfm(4,4)=p(0,5)
      pmcfm(5,4)=p(0,3)
      return


