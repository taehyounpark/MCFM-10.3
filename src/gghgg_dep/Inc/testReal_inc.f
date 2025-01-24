!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      use setreal_mcfm_generic
      implicit none
      include 'nlegborn.h'
      integer nlegs,i
      parameter (nlegs=nlegreal)
      real(dp):: p(0:3,nlegs)
      integer rflav(nlegs)
      real(dp):: amp2,amp2tmp

      amp2 = 0._dp
      amp2tmp = 0._dp

c gg initial state
        rflav(1)=0
        rflav(2)=0
        rflav(3)=25
        do i = 0,5
                rflav(4)=i
                rflav(5)=-i
                call setreal_mcfm(p,rflav,amp2tmp)
                amp2 = amp2+amp2tmp
        enddo

        write(6,*) "gg initial state, real.f: ",amp2

c gq initial state
      amp2 = 0._dp
      amp2tmp = 0._dp

        rflav(1)=0
        rflav(2)=1
        rflav(3)=25
        rflav(4)=1
        rflav(5)=0
        call setreal_mcfm(p,rflav,amp2tmp)
        amp2 = amp2+amp2tmp

      write(6,*) "gq initial state, real.f: ",amp2

c qg initial state
      amp2 = 0._dp
      amp2tmp = 0._dp

        rflav(1)=1
        rflav(2)=0
        rflav(3)=25
        rflav(4)=1
        rflav(5)=0
        call setreal_mcfm(p,rflav,amp2tmp)
        amp2 = amp2+amp2tmp

      write(6,*) "qg initial state, real.f: ",amp2

c qq initial state
      amp2 = 0._dp
      amp2tmp = 0._dp

        rflav(1)=1
        rflav(2)=1
        rflav(3)=25
        rflav(4)=1
        rflav(5)=1
        call setreal_mcfm(p,rflav,amp2tmp)
        amp2 = amp2+amp2tmp

        write(6,*) "qq initial state, real.f: ",amp2

c qq' initial state
      amp2 = 0._dp
      amp2tmp = 0._dp
        rflav(1)=1
        rflav(2)=2
        rflav(3)=25
        rflav(4)=1
        rflav(5)=2
        call setreal_mcfm(p,rflav,amp2tmp)
        amp2 = amp2+amp2tmp
      write(6,*) "qq' initial state, real.f: ",amp2

c q qbar initial state
      amp2 = 0._dp
      amp2tmp = 0._dp

        rflav(1)=1
        rflav(2)=-1
        rflav(3)=25
        do i = 0,5
                rflav(4)=i
                rflav(5)=-i
                call setreal_mcfm(p,rflav,amp2tmp)
                amp2 = amp2+amp2tmp
        enddo

        write(6,*) "qqbar initial state, real.f:",amp2

c qq'bar inital state
      amp2 = 0._dp
      amp2tmp = 0._dp
        rflav(1)=1
        rflav(2)=-2
        rflav(3)=25
        rflav(4)=1
        rflav(5)=-2
        call setreal_mcfm(p,rflav,amp2tmp)
        amp2 = amp2+amp2tmp
      write(6,*) "qq'bar initial state, real.f: ",amp2
      return
