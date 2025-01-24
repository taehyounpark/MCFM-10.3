!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      program run_gghgg
      use double_precision
      implicit none
      include 'nf.f'
      include 'masses.f'
      include 'nlegborn.h'
      integer nlegs
      parameter (nlegs=nlegreal)
      real * 8 p(0:3,nlegs)
      real(dp):: pmcfm(mxpart,4), msq(-nf:nf,-nf:nf),kappa
      integer rflav(nlegs)
      real * 8 amp2

      write(6,*)
      write(6,*)

      kappa=1d0/sqrt(94d0)

      p(0,1)=-15*kappa
      p(1,1)=-10*kappa
      p(2,1)=11*kappa
      p(3,1)=2*kappa

      p(0,2)=-9*kappa
      p(1,2)=8*kappa
      p(2,2)=kappa
      p(3,2)=-4*kappa

      p(0,4)=-21*kappa
      p(1,4)=4*kappa
      p(2,4)=-13*kappa
      p(3,4)=16*kappa

      p(0,5)=-7*kappa
      p(1,5)=2*kappa
      p(2,5)=-6*kappa
      p(3,5)=3*kappa


      p(0,1)=-15*kappa
      p(1,1)=-10*kappa
      p(2,1)=11*kappa
      p(3,1)=2*kappa

      p(0,2)=-9*kappa
      p(1,2)=8*kappa
      p(2,2)=kappa
      p(3,2)=-4*kappa

      p(0,4)=-21*kappa
      p(1,4)=4*kappa
      p(2,4)=-13*kappa
      p(3,4)=16*kappa

      p(0,5)=-7*kappa
      p(1,5)=2*kappa
      p(2,5)=-6*kappa
      p(3,5)=3*kappa

      !write(6,*) "p5",p(0,5)**2-p(1,5)**2-p(2,5)**2-p(3,5)**2

      hmass=5d0
      mt=1.5d0

      Gf=1.16639e-5_dp
      vevsq=1._dp/rt2/Gf
      as=0.118

      call testReal(p)

      !rflav(1)=0
      !rflav(2)=0
      !rflav(3)=25
      !rflav(4)=0
      !rflav(5)=0
      !call setreal(p,rflav,amp2)

c      write(6,*) 'p1 = ',pmcfm(1,:)
c      write(6,*) 'p2 = ',pmcfm(2,:)
c      write(6,*) 'p3 = ',pmcfm(3,:)
c      write(6,*) 'p4 = ',pmcfm(4,:)

        kappa=1d0/sqrt(94d0)

        pmcfm(1,4)=-15*kappa
        pmcfm(1,1)=-10*kappa
        pmcfm(1,2)=11*kappa
        pmcfm(1,3)=2*kappa

        pmcfm(2,4)=-9*kappa
        pmcfm(2,1)=8*kappa
        pmcfm(2,2)=kappa
        pmcfm(2,3)=-4*kappa

        pmcfm(5,4)=-21*kappa
        pmcfm(5,1)=4*kappa
        pmcfm(5,2)=-13*kappa
        pmcfm(5,3)=16*kappa

        pmcfm(6,4)=-7*kappa
        pmcfm(6,1)=2*kappa
        pmcfm(6,2)=-6*kappa
        pmcfm(6,3)=3*kappa

        hmass=5d0
        mt=1.5d0

c      write(6,*) 'p1 = ',pmcfm(1,:)
c      write(6,*) 'p2 = ',pmcfm(2,:)
c      write(6,*) 'p3 = ',pmcfm(5,:)
c      write(6,*) 'p4 = ',pmcfm(6,:)

      call gg_hgg_mass(pmcfm,msq)

      write(6,*)
      write(6,*)
      write(6,*) "gg initial state, MCFM:   ", msq(0,0)
      write(6,*) "gq initial state, MCFM:   ", msq(0,1)
      write(6,*) "qg initial state, MCFM:   ", msq(1,0)
      write(6,*) "qq initial state, MCFM:   ", msq(1,1)
      write(6,*) "qq' initial state, MCFM:   ", msq(1,2)
      write(6,*) "qqbar initial state, MCFM:   ", msq(1,-1)
      write(6,*) "qq'bar initial state, MCFM:   ", msq(1,-2)

      end
