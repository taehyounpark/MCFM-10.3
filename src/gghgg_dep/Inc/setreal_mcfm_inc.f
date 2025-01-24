!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      use spinor
      use pmcfm_convert_generic
      use haqgg_mass_generic
      use hgggg_mass_generic
      use haqaq_mass_generic
      use gghgg_dep_params
      implicit none
      include 'nlegborn.h'
      include 'zprods_decl.f'
      integer nlegs
      parameter (nlegs=nlegreal)
      integer rflav(nlegs)
      real(dp):: amp2,tmp
      real(dp):: p(0:3,5),pmcfm(mxpart,4)
      real(dp):: mtsq
      logical,save::first=.true.

c      if (first) then
c      first=.false.
c      write(6,*) 'vevsq',vevsq
c      write(6,*) 'mtsq',mtsq
c      write(6,*) 'as',as
c      endif
      amp2 = 0._dp


      call pmcfm_convert(p,pmcfm)

      !write(6,*) "p1mcfm^2",pmcfm(1,4)**2-pmcfm(1,1)**2-pmcfm(1,2)**2-pmcfm(1,3)**2
      !write(6,*) "p2mcfm^2",pmcfm(2,4)**2-pmcfm(2,1)**2-pmcfm(2,2)**2-pmcfm(2,3)**2
      !write(6,*) "p3mcfm^2",pmcfm(3,4)**2-pmcfm(3,1)**2-pmcfm(3,2)**2-pmcfm(3,3)**2
      !write(6,*) "p4mcfm^2",pmcfm(4,4)**2-pmcfm(4,1)**2-pmcfm(4,2)**2-pmcfm(4,3)**2
      !write(6,*) "p5mcfm^2",pmcfm(5,4)**2-pmcfm(5,1)**2-pmcfm(5,2)**2-pmcfm(5,3)**2

      !write(6,*) "p1mcfm",pmcfm(1,:)
      !write(6,*) "p2mcfm",pmcfm(2,:)
      !write(6,*) "p5mcfm",pmcfm(5,:)
      !write(6,*) "p3mcfm",pmcfm(3,:)
      !write(6,*) "p4mcfm",pmcfm(4,:)

c---fill spinor products
      call spinoru(4,pmcfm,za,zb)
c      call haqgg_mass(4,3,1,2,za,zb,amp2) ! Hggqa
c      write(6,*) amp2
c      call hgggg_mass(za,zb,amp2) ! Hgggg
c      write(6,*) amp2
c      call haqgg_mass(2,3,1,4,za,zb,amp2) ! Hgqqg
c      write(6,*) amp2
c      call haqgg_mass(4,2,3,1,za,zb,amp2) ! Hgqgq or Hgaga ! check
c      write(6,*) amp2
c      call haqgg_mass(1,3,2,4,za,zb,amp2) ! Hqgqg
c      write(6,*) amp2
c      call haqgg_mass(1,4,3,2,za,zb,amp2) ! Hqggq
c      write(6,*) amp2
c      call haqgg_mass(1,2,3,4,za,zb,amp2) ! Hqagg
c      write(6,*) amp2
c      call haqaq_mass(1,2,3,4,za,zb,tmp,amp2) ! Hqaaq
c      write(6,*) amp2
c      call haqaq_mass(1,2,4,3,za,zb,tmp,amp2) ! Hqaqa
c      write(6,*) amp2
c      call haqaq_mass(1,2,4,3,za,zb,amp2,tmp) ! Hqarb
c      write(6,*) amp2
c      call haqaq_mass(1,3,2,4,za,zb,tmp,amp2) ! Hqqqq
c      write(6,*) amp2
c      call  haqaq_mass(1,3,2,4,za,zb,amp2,tmp) ! Hqrqr
c      write(6,*) amp2
c      call haqaq_mass(1,4,2,3,za,zb,amp2,tmp) ! Hqrrq
c      write(6,*) amp2

      if ((rflav(1) == 0).and.(rflav(2) == 0)) then
               ! gg initial state
         if(rflav(4) /= 0) then
            ! g g -> H q qbar
            call haqgg_mass(4,3,1,2,za,zb,amp2) ! Hggqa
            amp2 = avegg*amp2
         else
            ! g g -> H g g
            call hgggg_mass(za,zb,amp2) ! Hgggg
            amp2 = avegg*amp2/2._dp
         endif

      else if ((rflav(1) == 0).and.(rflav(2) /= 0)) then
              ! gq initial state
         if(rflav(4) /= 0) then
            ! g q(bar) -> H q(bar) g
            call haqgg_mass(2,3,1,4,za,zb,amp2) ! Hgqqg
            amp2 = aveqg*amp2
         else
            ! g q(bar) -> H g q(bar)
            call haqgg_mass(4,2,3,1,za,zb,amp2) ! Hgqgq or Hgaga ! check
            amp2 = aveqg*amp2
         endif

      else if ((rflav(1) /= 0).and.(rflav(2) == 0)) then
              ! qg initial state
         if(rflav(4) /= 0) then
            ! q(bar) g -> H q(bar) g
            call haqgg_mass(1,3,2,4,za,zb,amp2) ! Hqgqg
            amp2 = aveqg*amp2
         else
            ! q(bar) g -> H g q(bar)
            call haqgg_mass(1,4,3,2,za,zb,amp2) ! Hqggq
            amp2 = aveqg*amp2
         endif

      else if ((rflav(1) == -rflav(2))) then
              ! qa initial state
       if(rflav(4) == 0) then
         ! q qbar -> H g g
         call haqgg_mass(1,2,3,4,za,zb,amp2) ! Hqagg
         amp2 = aveqq*amp2/2._dp
       else if (rflav(2) == rflav(4)) then
         ! q qbar -> H qbar q
         call haqaq_mass(1,2,3,4,za,zb,tmp,amp2) ! Hqaaq
         amp2 = aveqq*amp2
       else if (abs(rflav(1)) == abs(rflav(4))) then
         ! q qbar -> H q qbar
         call haqaq_mass(1,2,4,3,za,zb,tmp,amp2) ! Hqaqa
         amp2 = aveqq*amp2
       else
         ! q qbar -> H r rbar
         call haqaq_mass(1,2,4,3,za,zb,amp2,tmp) ! Hqarb
         amp2 = aveqq*amp2
       endif

      else if ((abs(rflav(1)) == abs(rflav(2))).and.
     &         (abs(rflav(1)) == abs(rflav(4)))) then
                ! qq initial state
         ! q(bar) + q(bar) -> H q(bar) q(bar)
         call haqaq_mass(1,3,2,4,za,zb,tmp,amp2) ! Hqqqq
         amp2 = aveqq*amp2/2._dp

      else if ((abs(rflav(1)) == abs(rflav(4))).and.
     &         (abs(rflav(2)) == abs(rflav(5)))) then
         ! q(bar) + q(bar)' -> H q(bar) q(bar)'
         call  haqaq_mass(1,3,2,4,za,zb,amp2,tmp) ! Hqrqr
         amp2 = aveqq*amp2

      else if ((abs(rflav(1)) == abs(rflav(5))).and.
     &         (abs(rflav(2)) == abs(rflav(4)))) then
         ! q(bar) + q(bar)' -> H q(bar)' q(bar)
         call haqaq_mass(1,4,2,3,za,zb,amp2,tmp) ! Hqrrq
         amp2 = aveqq*amp2
      else
              write(6,*) "rflav not matched with MCFM call"
      endif

c      amp2 = amp2/(st_alpha/(2d0*pi))

c      write(6,*) 'Hqrqr',aveqq*Hqrqr
c      write(6,*) 'Hqqqq',half*aveqq*Hqqqq
c      write(6,*) 'Hqarb',aveqq*Hqarb
c      write(6,*) 'Hqaqa',aveqq*Hqaqa
c      write(6,*)
c      write(6,*) 'Hqagg',aveqq*half*Hqagg
c      write(6,*) 'Hqgqg',aveqg*Hqgqg
c      write(6,*) 'Hgqqg',aveqg*Hgqqg
c      write(6,*) 'Hggqa',avegg*Hggqa
c      write(6,*)
c      end subroutine setreal_mcfm


