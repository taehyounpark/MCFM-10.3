!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module singletop_jet
      use ieee_arithmetic
      use types
      use singletop2_scale_m

#define ENABLE_NDOTP 0

      public :: singletop_jet_light_msqall

      public :: singletop_jet_light_virt_all

      public :: singletop_jet_light_z

      public :: singletop_jet_light_real_all
      public :: singletop_jet_light_gvec
      public :: singletop_jet_light_gs_all
      public :: singletop_jet_light_swap
      public :: singletop_jet_light_gvec_swap

      private :: singletop_jet_light ! for dipoles
      public :: ampsq_ugd_tdkb
      private :: ampsq_ugd_tdkb_v
      private :: amp_udqqb_tdkb
      private :: Augd_tdkb
      private :: subqcd_w_tdkb
      private :: W_tdkb_curr
      private :: ampsq_ugd_tdkb_n

      private

      contains

      subroutine singletop_jet_light_msqall(p,msqall)
      use singletop2_nnlo_vars
      implicit none
      include 'nf.f'
      include 'mxpart.f'
      include 'nwz.f'
      include 'zprods_com.f'
      include 'constants.f'
      include 'ewcouple.f'
      real(dp), intent(in) :: p(mxpart,4)
      real(dp), intent(out) :: msqall(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

      real(dp) :: fac,genfac
      real(dp) :: tmp(-nf:nf)

      integer :: i1,i2,m

      call spinoru(7,p,za,zb)

c factor representing sum over generations, i.e. (u,d~) and (c,s~) contributions
c set to 1._dp to reproduce Recola results for a specific final state
      genfac=2._dp

      msqall = 0._dp

      do m=1,maxbeams
          corr_on_beam = beams_enabled(m)
c loop over permutations of initial state:
c   cob=1    i1=1, i2=2,  light beam 1, heavy beam 2
c   cob=2    i1=2, i2=1,  light beam 2, heavy beam 1
          if (corr_on_beam == 1) then
            i1=1
            i2=2
            fac=2._dp*(fourpi*as_light_beam1)*cf*gw**8*xn**2
          else
            i1=2
            i2=1
            fac=2._dp*(fourpi*as_light_beam2)*cf*gw**8*xn**2
          endif

          tmp(:)=0._dp

          if (nwz == +1) then
            if (iand(partons_enabled, quarkChannel) > 0) then
                tmp([2,4]) = aveqq*fac*ampsq_ugd_tdkb(i1,i2,3,4,5,6,7,p)            ! u  b -> t d g
                tmp([-1,-3])= aveqq*fac*ampsq_ugd_tdkb(6,i2,3,4,5,i1,7,p)            ! d~ b -> t u g
            endif
            if (iand(partons_enabled, gluonChannel) > 0) then
                tmp(0) = genfac*aveqg*fac*ampsq_ugd_tdkb(7,i2,3,4,5,6,i1,p)       ! g  b -> t d u~  +  g  b -> t s c~
            endif
          else
            stop 'Abort in singletop_jet_light_msqall'
          endif

c update msq array
          if (corr_on_beam == 1) then
            msqall(:,5, 1,corr_on_beam)=tmp(:)
          else
            msqall(5,:, 1,corr_on_beam)=tmp(:)
          endif
      enddo

      end subroutine singletop_jet_light_msqall

      subroutine singletop_jet_light(p,msq)
      use singletop2_nnlo_vars
      implicit none
      include 'nf.f'
      include 'mxpart.f'
      include 'nwz.f'
      include 'zprods_com.f'
      include 'constants.f'
      include 'ewcouple.f'
      real(dp), intent(in) :: p(mxpart,4)
      real(dp), intent(out) :: msq(-nf:nf, -nf:nf)

      real(dp) :: fac,genfac
      real(dp) :: tmp(-nf:nf)

      integer :: i1,i2

      call spinoru(7,p,za,zb)

c factor representing sum over generations, i.e. (u,d~) and (c,s~) contributions
c set to 1._dp to reproduce Recola results for a specific final state
      genfac=2._dp

      tmp(:)=0._dp
      msq(:,:) = 0._dp

c loop over permutations of initial state:
c   iperm=1    i1=1, i2=2,  light beam 1, heavy beam 2
c   iperm=2    i1=2, i2=1,  light beam 2, heavy beam 1

      if (corr_on_beam == 1) then
        i1=1
        i2=2
        fac=2._dp*(fourpi*as_light_beam1)*cf*gw**8*xn**2
      else
        i1=2
        i2=1
        fac=2._dp*(fourpi*as_light_beam2)*cf*gw**8*xn**2
      endif

      if (nwz == +1) then
          tmp(2) = aveqq*fac*ampsq_ugd_tdkb(i1,i2,3,4,5,6,7,p)            ! u  b -> t d g
          tmp(-1)= aveqq*fac*ampsq_ugd_tdkb(6,i2,3,4,5,i1,7,p)            ! d~ b -> t u g
          tmp(4) = tmp(2)
          tmp(-3) = tmp(-1)
          tmp(0) = genfac*aveqg*fac*ampsq_ugd_tdkb(7,i2,3,4,5,6,i1,p)       ! g  b -> t d u~  +  g  b -> t s c~
      else
          write(6,*) 'Abort in'
      stop
      endif

c update msq array
      if (corr_on_beam == 1) then
        msq(:,5)=tmp(:)
      else
        msq(5,:)=tmp(:)
      endif

      end subroutine singletop_jet_light


      subroutine singletop_jet_light_swap(p,msq)
c This just calls the lowest order routine with p6 and p7 swapped, necessary for subtraction terms
      implicit none
      include 'nf.f'
      include 'mxpart.f'
      include 'nwz.f'
      include 'zprods_com.f'
      include 'constants.f'
      include 'ewcouple.f'
      real(dp), intent(in) :: p(mxpart,4)
      real(dp), intent(out) :: msq(-nf:nf, -nf:nf)

      real(dp) :: pswap(mxpart,4)

      pswap(:,:)=p(:,:)
      pswap(6,:)=p(7,:)
      pswap(7,:)=p(6,:)

      call singletop_jet_light(pswap,msq)

      end subroutine singletop_jet_light_swap

      subroutine singletop_jet_light_virt_all(p,msqall)
      use singletop2_nnlo_vars
      implicit none
      include 'nf.f'
      include 'mxpart.f'
      include 'nwz.f'
      include 'zprods_com.f'
      include 'constants.f'
      include 'ewcouple.f'
      include 'scheme.f'
      real(dp), intent(in) :: p(mxpart,4)
      real(dp), intent(out) :: msqall(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

      real(dp) :: fac,genfac
      real(dp) :: tmp(-nf:nf)

      integer :: i1,i2,m

      scheme = 'tH-V'
      call spinoru(7,p,za,zb)

c factor representing sum over generations, i.e. (u,d~) and (c,s~) contributions
c set to 1._dp to reproduce Recola results for a specific final state
      genfac=2._dp

      msqall(:,:,:,:) = 0._dp

c loop over permutations of initial state:
c   iperm=1    i1=1, i2=2,  light beam 1, heavy beam 2
c   iperm=2    i1=2, i2=1,  light beam 2, heavy beam 1
      do m=1,maxbeams
      corr_on_beam = beams_enabled(m)

      tmp(:) = 0._dp

      if (corr_on_beam == 1) then
        i1=1
        i2=2
        fac=2._dp*(fourpi*as_light_beam1)*cf*gw**8*xn**2
      else
        i1=2
        i2=1
        fac=2._dp*(fourpi*as_light_beam2)*cf*gw**8*xn**2
      endif

      if (nwz == +1) then
          tmp([2,4]) = aveqq*fac*ampsq_ugd_tdkb_v(i1,i2,3,4,5,6,7,p)      ! u  b -> t d g
          tmp([-1,-3]) = aveqq*fac*ampsq_ugd_tdkb_v(6,i2,3,4,5,i1,7,p)      ! d~ b -> t u g
          tmp(0) = genfac*aveqg*fac*ampsq_ugd_tdkb_v(7,i2,3,4,5,6,i1,p) ! g  b -> t d u~  +  g  b -> t s c~
      else
          stop 'Abort in singletop_jet_light'
      endif

c update msq array
      if (corr_on_beam == 1) then
        msqall(:,5, 1,1) = tmp(:)
      else
        msqall(5,:, 1,2) = tmp(:)
      endif

      enddo

      end subroutine singletop_jet_light_virt_all

      subroutine singletop_jet_light_real_all(p,msqall)
      use singletop2_nnlo_vars
      implicit none
      include 'nf.f'
      include 'mxpart.f'
      include 'nwz.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'constants.f'
      include 'ewcouple.f'
      include 'masses.f'

      real(dp), intent(in) :: p(mxpart,4)
      real(dp), intent(out) :: msqall(-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

      real(dp) :: facgg,facqq,propWsq,s1678,tmp(-nf:nf,max_bcontrib),
     & genfac
      integer :: iperm,i1,i2

      call spinoru(8,p,za,zb)

c factor representing sum over generations, i.e. (u,d~) and (c,s~) contributions
c set to 1._dp to reproduce Recola results for a specific final state
      genfac=2._dp

      msqall(:,:,:,:) = 0._dp

c loop over permutations of initial state:
c   iperm=1    i1=1, i2=2,  light beam 1, heavy beam 2
c   iperm=2    i1=2, i2=1,  light beam 2, heavy beam 1
      do iperm=1,maxbeams
          corr_on_beam = beams_enabled(iperm)

      tmp(:,:)=0._dp

      if (corr_on_beam == 1) then
        i1=1
        i2=2
        facqq=(fourpi*as_light_beam1)**2
        facgg=facqq
      else
        i1=2
        i2=1
        facqq=(fourpi*as_light_beam2)**2
        facgg=facqq
      endif
c overall propagator factor needs to change from 1/s1678 to 1/(s1678-wmass**2)
      s1678=s(i1,6)+s(i1,7)+s(i1,8)+s(6,7)+s(6,8)+s(7,8)
      propWsq=(s1678/(s1678-wmass**2))**2

      facqq=facqq*V*xn*gwsq**4*aveqq*propWsq
      facgg=facgg*V*xn**2*(gwsq/two)**4*propWsq

      if (iand(partons_enabled, quarkChannel) > 0) then

c note that labellings all correspond to heavy beam 2
      tmp([2,4], 1)=
c u b -> t d g g
     &   half*aveqq*facgg*w2jetsq_w_tdkb(i1,6,7,8,3,4,5,i2,za,zb)
c u b -> t d d d~
     &  +half*facqq*ampsq_udddb_tdkb(i1,i2,3,4,5,6,7,8,za,zb)
c u b -> t d u u~
     &  +facqq*ampsq_uduub_tdkb(i1,i2,3,4,5,6,7,8,za,zb)
c u b -> t s u c~
     &  +facqq*ampsq_udqqb_tdkb(8,i2,3,4,5,6,7,i1,za,zb)
c u b -> t d q q~, with q=s,c
     &  +facqq*ampsq_udqqb_tdkb(i1,i2,3,4,5,6,7,8,za,zb)*two
      tmp([2,4], 4) =
c u b -> t d b b~
     &  +facqq*ampsq_udqqb_tdkb(i1,i2,3,4,5,6,7,8,za,zb)

      tmp([1,3], 1)=
c d b -> t d u~ d
     & half*facqq*ampsq_udddb_tdkb(7,i2,3,4,5,6,8,i1,za,zb)
c d b -> t s c~ d
     & +facqq*ampsq_udqqb_tdkb(7,i2,3,4,5,6,8,i1,za,zb)      

c b b -> t d u~ b + b b -> t s c~ b
      tmp(5, 2)=genfac*facqq*ampsq_udqqb_tdkb(7,i2,3,4,5,6,8,i1,za,zb)

      tmp([-1,-3], 1)=
c d~ b -> t u~ g g
     &   half*aveqq*facgg*w2jetsq_w_tdkb(6,i1,7,8,3,4,5,i2,za,zb)
c d~ b -> t u~ d~ d
     &  +facqq*ampsq_udddb_tdkb(6,i2,3,4,5,i1,8,7,za,zb)          
c d~ b -> t u~ u~ u
     &  +half*facqq*ampsq_uduub_tdkb(6,i2,3,4,5,i1,8,7,za,zb)       
c d~ b -> t c~ d~ s
     &  +facqq*ampsq_udqqb_tdkb(6,i2,3,4,5,8,i1,7,za,zb)          
c d~ b -> t u~ q~ q, with q=s,c
     &  +facqq*ampsq_udqqb_tdkb(6,i2,3,4,5,i1,8,7,za,zb)*two
      tmp([-1,-3], 4) =
c d~ b -> t u~ b~ b
     &  +facqq*ampsq_udqqb_tdkb(6,i2,3,4,5,i1,8,7,za,zb)          

      tmp([-2,-4], 1)=
c u~ b -> t u~ d u~
     &  half*facqq*ampsq_uduub_tdkb(6,i2,3,4,5,7,i1,8,za,zb)       
c u~ b -> t c~ s u~
     &  +facqq*ampsq_udqqb_tdkb(6,i2,3,4,5,7,i1,8,za,zb)

      tmp(-5, 3)=genfac*facqq*ampsq_udqqb_tdkb(6,i2,3,4,5,7,i1,8,za,zb)       
c b~ b -> t u~ d b~ + b~ b -> t c~ s b~
      endif

      if (iand(partons_enabled, gluonChannel) > 0) then
c g b -> t d u~ g + g b -> t s c~ g
          tmp(0, 1)=genfac*aveqg*facgg*w2jetsq_w_tdkb(7,6,i1,8,3,4,5,i2,za,zb)    
      endif

c update msq array, taking care to include both b b contributions
      if (corr_on_beam == 1) then
          msqall(:,5,:,1) = tmp(:,:)
      else
          msqall(5,:,:,2) = tmp(:,:)
      endif

      enddo

      end subroutine singletop_jet_light_real_all

      subroutine singletop_jet_light_gvec(p,n,in,msq)
      use singletop2_nnlo_vars
      implicit none
      include 'nf.f'
      include 'mxpart.f'
      include 'nwz.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'constants.f'
      include 'ewcouple.f'
      real(dp), intent(in) :: p(mxpart,4), n(4)
      integer, intent(in) :: in
      real(dp), intent(out) :: msq(-nf:nf, -nf:nf)

      real(dp) :: fac,genfac
      complex(dp) :: zanb(mxpart,mxpart),zbna(mxpart,mxpart)

#if (ENABLE_NDOTP == 1)
      call checkndotp(p,n,in)
#endif

      call spinoru(7,p,za,zb)
      call spinork(7,p,zanb,zbna,n)

      if (corr_on_beam == 1) then
        fac=2._dp*(fourpi*as_light_beam1)*cf*gw**8*xn**2
      elseif (corr_on_beam == 2) then
        fac=2._dp*(fourpi*as_light_beam2)*cf*gw**8*xn**2
      endif

c factor representing sum over generations, i.e. (u,d~) and (c,s~) contributions
c set to 1._dp to reproduce Recola results for a specific final state
      genfac=2._dp

      msq(:,:) = 0._dp

      if (nwz == +1) then
          if     (in == 1) then
            if (corr_on_beam == 1) then
            msq(0,5) = genfac*aveqg*fac*ampsq_ugd_tdkb_n(7,2,3,4,5,6,1,zanb)
            endif
          elseif (in == 2) then
            if (corr_on_beam == 2) then
            msq(5,0) = genfac*aveqg*fac*ampsq_ugd_tdkb_n(7,1,3,4,5,6,2,zanb)
            endif
          elseif (in == 6) then
            if (corr_on_beam == 1) then
            msq([2,4],5) = aveqq*fac*ampsq_ugd_tdkb_n(1,2,3,4,5,7,6,zanb)
            msq([-1,-3],5)= aveqq*fac*ampsq_ugd_tdkb_n(7,2,3,4,5,1,6,zanb)
            endif

            if (corr_on_beam == 2) then
            msq(5,[2,4]) = aveqq*fac*ampsq_ugd_tdkb_n(2,1,3,4,5,7,6,zanb)
            msq(5,[-1,-3])= aveqq*fac*ampsq_ugd_tdkb_n(7,1,3,4,5,2,6,zanb)
            endif
          elseif (in == 7) then
            if (corr_on_beam == 1) then
            msq([2,4],5) = aveqq*fac*ampsq_ugd_tdkb_n(1,2,3,4,5,6,7,zanb)
            msq([-1,-3],5)= aveqq*fac*ampsq_ugd_tdkb_n(6,2,3,4,5,1,7,zanb)
            endif

            if (corr_on_beam == 2) then
            msq(5,[2,4]) = aveqq*fac*ampsq_ugd_tdkb_n(2,1,3,4,5,6,7,zanb)
            msq(5,[-1,-3])= aveqq*fac*ampsq_ugd_tdkb_n(6,1,3,4,5,2,7,zanb)
            endif
          else
            write(6,*) 'singletop_jet_light_gvec: unexpected value in =',in
            write(6,*) 'Abort in singletop_jet_light_gvec'
            stop
          endif
      else
          write(6,*) 'Abort in singletop_jet_light_gvec'
          stop
      endif

      end subroutine singletop_jet_light_gvec


      subroutine singletop_jet_light_gvec_swap(p,n,in,msq)
c This just calls the gvec routine with p6 and p7 swapped, necessary for subtraction terms
      implicit none
      include 'nf.f'
      include 'mxpart.f'
      include 'nwz.f'
      include 'zprods_com.f'
      include 'constants.f'
      include 'ewcouple.f'
      real(dp), intent(in) :: p(mxpart,4),n(4)
      real(dp), intent(out) :: msq(-nf:nf, -nf:nf)
      integer, intent(in) :: in

      real(dp) :: pswap(mxpart,4)

      pswap(:,:)=p(:,:)
      pswap(6,:)=p(7,:)
      pswap(7,:)=p(6,:)

      call singletop_jet_light_gvec(pswap,n,in,msq)

      end subroutine singletop_jet_light_gvec_swap


      function ampsq_ugd_tdkb(j1,j2,j3,j4,j5,j6,j7,p)
c--- Amplitude for 0 -> ubar(j1) + g(j7) + d(j6) [light line]
c---              + t (nu(j3) + e+(j4) + b(j5)) + bbar(j2)  [heavy line]
c---
      implicit none
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'masses.f'
      include 'scheme.f'
      include 'epinv.f'
      include 'b0.f'
      integer, intent(in) :: j1,j2,j3,j4,j5,j6,j7
      real(dp), intent(in) :: p(mxpart,4)
      real(dp):: ampsq_ugd_tdkb,s167,propW
      complex(dp):: amplop,amplom

c--- overall propagator factor needs to change from 1/s167 to 1/(s167-wmass**2)
      s167=s(j1,j6)+s(j1,j7)+s(j6,j7)
      propW=s167/(s167-wmass**2)

c Positive helicity gluon
c make sure to call these in this order to ensure correct sign of LO amplitude
      call Augd_tdkb(amplop,j1,j7,j6,j3,j4,j5,j2,za,zb,.false.)

c Negative helicity gluon
c make sure to call these in this order to ensure correct sign of LO amplitude
      call Augd_tdkb(amplom,j6,j7,j1,j3,j4,j5,j2,zb,za,.true.)

c Sum and include overall propagator factor
      ampsq_ugd_tdkb=cdabs(amplop)**2+cdabs(amplom)**2
      ampsq_ugd_tdkb=ampsq_ugd_tdkb*propW**2

      return
      end function ampsq_ugd_tdkb


      function ampsq_ugd_tdkb_v(j1,j2,j3,j4,j5,j6,j7,p)
c--- Amplitude for 0 -> ubar(j1) + g(j7) + d(j6) [light line]
c---              + t (nu(j3) + e+(j4) + b(j5)) + bbar(j2)  [heavy line]
c---
c--- with virtual corrections on the light line only
c---
      use singletop2_nnlo_vars
      implicit none
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'masses.f'
      include 'scheme.f'
      include 'epinv.f'
      include 'b0.f'
      integer, intent(in) :: j1,j2,j3,j4,j5,j6,j7
      real(dp), intent(in) :: p(mxpart,4)
      real(dp):: ampsq_ugd_tdkb_v,s167,propW
      complex(dp):: amplop,amplom,ampvp,ampvm,ampv_lc,ampv_slc

c--- overall propagator factor needs to change from 1/s167 to 1/(s167-wmass**2)
      s167=s(j1,j6)+s(j1,j7)+s(j6,j7)
      propW=s167/(s167-wmass**2)

c Positive helicity gluon
c make sure to call these in this order to ensure correct sign of LO amplitude
      call Augd_tdkb_v2(amplop,ampv_slc,j1,j6,j7,j3,j4,j5,j2,za,zb,.false.)
      call Augd_tdkb_v1(amplop,ampv_lc,j1,j7,j6,j3,j4,j5,j2,za,zb,.false.)

      ampvp=xn*ampv_lc-ampv_slc/xn
c--- strong coupling renormalization in dred scheme
      ampvp=ampvp-xn*(b0/xn*epinv-1._dp/6._dp)*amplop
      if (scheme  ==  'tH-V') then
c--- known translation rules
      ampvp=ampvp-(Cf+xn/6._dp)*amplop
      endif

c Negative helicity gluon
c make sure to call these in this order to ensure correct sign of LO amplitude
      call Augd_tdkb_v2(amplom,ampv_slc,j6,j1,j7,j3,j4,j5,j2,zb,za,.true.)
      call Augd_tdkb_v1(amplom,ampv_lc,j6,j7,j1,j3,j4,j5,j2,zb,za,.true.)

      ampvm=xn*ampv_lc-ampv_slc/xn
c--- strong coupling renormalization in dred scheme
      ampvm=ampvm-xn*(b0/xn*epinv-1._dp/6._dp)*amplom
      if (scheme  ==  'tH-V') then
c--- known translation rules
      ampvm=ampvm-(Cf+xn/6._dp)*amplom
      endif

c Sum and include overall propagator factor
      ampsq_ugd_tdkb_v=dble(ampvp*dconjg(amplop)+ampvm*dconjg(amplom))
      if     (corr_on_beam == 1) then
      ampsq_ugd_tdkb_v=ampsq_ugd_tdkb_v*(as_light_beam1/twopi)
      elseif (corr_on_beam == 2) then
      ampsq_ugd_tdkb_v=ampsq_ugd_tdkb_v*(as_light_beam2/twopi)
      endif
      ampsq_ugd_tdkb_v=ampsq_ugd_tdkb_v*propW**2

      return
      end function ampsq_ugd_tdkb_v

      subroutine Augd_tdkb(Alo,j1,j2,j3,j4,j5,j6,j7,za,zb,flip)
c--- Amplitude for 0 -> ubar(j1) + g(j2) + d(j3) [light line]
c---              + t (nu(j4) + e+(j5) + b(j6)) + bbar(j7)  [heavy line]
c---
c--- The amplitudes in this routine are obtained by extending the original
c--- BDK V+jet amplitudes
      implicit none
      include 'constants.f'
      include 'mxpart.f'
      include 'epinv.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      logical, intent(in) :: flip
      integer, intent(in) :: j1,j2,j3,j4,j5,j6,j7
      real(dp) :: s123
      complex(dp), intent(out) :: Alo

      s123=s(j1,j2)+s(j1,j3)+s(j2,j3)
c    -i * A5tree rewritten in a suitable way
c      A5lom =-za(j3,j4)**2/(za(j1,j2)*za(j2,j3)*za(j4,j5))
      Alo = -(za(j3,j1)*W_tdkb_curr(j3,j1,j4,j5,j6,j7,za,zb,flip)
     &       +za(j3,j2)*W_tdkb_curr(j3,j2,j4,j5,j6,j7,za,zb,flip))
     &       /(za(j1,j2)*za(j2,j3)*s123)

      return
      end subroutine Augd_tdkb


      function W_tdkb_curr(jx,jy,p3,p4,p5,p2,za,zb,flip)
      implicit none
c--- this function implements the replacement of a term
c--- of the form <jx j4> [jy j5] in the V(j4,j5)+1 jet amplitudes
c--- with the expression corresponding to replacing the
c--- current giving |jx> |jy] with the one for the line
c--- t( nu(p3) + e+(p4) + b(p5)) + bbar(p2)
c---
c--- if flip is true then we need to interchange the roles of jx and jy,
c--- and also interchange za and zb, which will already have been interchanged
c--- in the original call
      include 'cplx.h'
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'masses.f'
      logical, intent(in) :: flip
      integer, intent(in) :: jx,jy
      integer, intent(in) :: p2,p3,p4,p5
      real(dp):: s2345
      complex(dp):: W_tdkb_curr,zab2,propW

c--- statement function
      zab2(p2,p3,p4,p5)=za(p2,p3)*zb(p3,p5)+za(p2,p4)*zb(p4,p5)
      propw(s2345)=cone/cplx2((s2345-wmass**2),wmass*wwidth)

      s2345=s(p2,p3)+s(p2,p4)+s(p2,p5)+s(p3,p4)+s(p3,p5)+s(p4,p5)

      if (flip) then
      W_tdkb_curr=zb(p5,p3)/(im*twidth*mt)*propW(s(p3,p4))/2._dp
     &   *(2._dp*zab2(p4,p3,p5,jy)*za(jx,p2)
     &    -mt**2*(zab2(jx,p2,p3,jy)+zab2(jx,p4,p5,jy))/s2345*za(p4,p2))
      else
      W_tdkb_curr=za(p5,p3)/(im*twidth*mt)*propW(s(p3,p4))/2._dp
     &   *(2._dp*zab2(jx,p3,p5,p4)*zb(jy,p2)
     &    -mt**2*(zab2(jx,p2,p3,jy)+zab2(jx,p4,p5,jy))/s2345*zb(p4,p2))
      endif

      return
      end function W_tdkb_curr


      subroutine Augd_tdkb_v1(Alo,Av,j1,j2,j3,j4,j5,j6,j7,za,zb,flip)
c--- Amplitude for 0 -> ubar(j1) + g(j2) + d(j3) [light line]
c---              + t (nu(j4) + e+(j5) + b(j6)) + bbar(j7)  [heavy line]
c---
c--- with virtual corrections on the light line only
c---
c--- The amplitudes in this routine are obtained by extending the original
c--- BDK V+jet amplitudes
      use singletop2_nnlo_vars
      implicit none
      include 'constants.f'
      include 'mxpart.f'
      include 'epinv.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      logical, intent(in) :: flip
      integer, intent(in) :: j1,j2,j3,j4,j5,j6,j7
      complex(dp), intent(out) :: Alo, Av
      real(dp) :: s123,musq
      complex(dp) :: l12,l23,L0,L1,Lsm1,Anlo(-2:0)
      complex(dp) :: lnrat

      s123=s(j1,j2)+s(j1,j3)+s(j2,j3)
c    -i * A5tree rewritten in a suitable way
c      A5lom =-za(j3,j4)**2/(za(j1,j2)*za(j2,j3)*za(j4,j5))
      Alo =-(za(j3,j1)*W_tdkb_curr(j3,j1,j4,j5,j6,j7,za,zb,flip)
     &      +za(j3,j2)*W_tdkb_curr(j3,j2,j4,j5,j6,j7,za,zb,flip))
     &     /(za(j1,j2)*za(j2,j3)*s123)

c--- this is the other helicity, rewritten appropriately
c      Alo =-(zb(j1,j3)*W_tdkb_curr(j3,j1,j4,j5,j6,j7,za,zb,flip)
c     &      +zb(j1,j2)*W_tdkb_curr(j2,j1,j4,j5,j6,j7,za,zb,flip))
c     &     /(zb(j3,j2)*zb(j2,j1)*s123)

      if (corr_on_beam == 1) then
      musq=renscale_beam1_islight_onlight**2
      elseif (corr_on_beam == 2) then
      musq=renscale_beam2_islight_onlight**2
      endif
      l12=lnrat(musq,-s(j1,j2))
      l23=lnrat(musq,-s(j2,j3))

c--leading N
c      Vcc=
c     & -(epinv**2+epinv*l12+0.5_dp*l12**2)
c     & -(epinv**2+epinv*l23+0.5_dp*l23**2)
c     & -2._dp*(epinv+l23)-4._dp
      Anlo(-2)=-2._dp*Alo
      Anlo(-1)=-(l12+l23+2._dp)*Alo
      Anlo( 0)=(-0.5_dp*(l12**2+l23**2)-2._dp*l23-4._dp)*Alo

c      Fcc=za(j3,j4)**2/(za(j1,j2)*za(j2,j3)*za(j4,j5))
c     & *(Lsm1(-s(j1,j2),-s(j4,j5),-s(j2,j3),-s(j4,j5))
c     & -2._dp*za(j3,j1)*zb(j1,j5)*za(j5,j4)/za(j3,j4)
c     &   *L0(-s(j2,j3),-s(j4,j5))/s(j4,j5))

      Anlo(0)=Anlo(0)
     & +(za(j3,j1)*W_tdkb_curr(j3,j1,j4,j5,j6,j7,za,zb,flip)
     &  +za(j3,j2)*W_tdkb_curr(j3,j2,j4,j5,j6,j7,za,zb,flip))
     &  /(za(j1,j2)*za(j2,j3)*s123)
     &  *Lsm1(-s(j1,j2),-s123,-s(j2,j3),-s123)
     & -2._dp*za(j1,j3)/(za(j1,j2)*za(j2,j3))
     &  *W_tdkb_curr(j3,j1,j4,j5,j6,j7,za,zb,flip)
     &  *L0(-s(j2,j3),-s123)/s123

c      Vsc =0.5_dp*(epinv+l23)+1._dp
      Anlo(-1)=Anlo(-1)+0.5_dp*Alo
      Anlo( 0)=Anlo( 0)+(0.5_dp*l23+1._dp)*Alo

c      Fsc =za(j3,j4)*za(j3,j1)*zb(j1,j5)*za(j5,j4)
c     & /(za(j1,j2)*za(j2,j3)*za(j4,j5))*L0(-s(j2,j3),-s123)/s123
c     & +0.5_dp*(za(j3,j1)*zb(j1,j5))**2*za(j4,j5)
c     & /(za(j1,j2)*za(j2,j3))*L1(-s(j2,j3),-s123)/s123**2

      Anlo(0)=Anlo(0)
     & +za(j1,j3)/(za(j1,j2)*za(j2,j3))*W_tdkb_curr(j3,j1,j4,j5,j6,j7,za,zb,flip)
     &  *L0(-s(j2,j3),-s123)/s123
     & +0.5_dp*za(j1,j3)**2/(za(j1,j2)*za(j2,j3))
     &  *(zb(j1,j2)*W_tdkb_curr(j2,j1,j4,j5,j6,j7,za,zb,flip)
     &   +zb(j1,j3)*W_tdkb_curr(j3,j1,j4,j5,j6,j7,za,zb,flip))
     &  *L1(-s(j2,j3),-s123)/s123**2

c      A51=(Vcc+Vsc)*Alo+Fcc+Fsc

      Av=Anlo(-2)*epinv**2+Anlo(-1)*epinv+Anlo(0)

      return
      end subroutine Augd_tdkb_v1


      subroutine Augd_tdkb_v2(Alo,Av,j1,j2,j3,j4,j5,j6,j7,za,zb,flip)
      use singletop2_nnlo_vars
      implicit none
      include 'constants.f'
      include 'mxpart.f'
      include 'epinv.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      logical, intent(in) :: flip
      integer, intent(in) :: j1,j2,j3,j4,j5,j6,j7
      complex(dp), intent(out) :: Alo,Av
      real(dp) :: s123,musq
      complex(dp) :: l12,l123,L0,L1,Lsm1,Anlo(-2:0)
      complex(dp) :: lnrat

      s123=s(j1,j2)+s(j1,j3)+s(j2,j3)
c    -i * A5tree rewritten in a suitable way
c      Alom=za(j2,j4)**2/(za(j2,j3)*za(j3,j1)*za(j4,j5))
      Alo =-(za(j2,j1)*W_tdkb_curr(j2,j1,j4,j5,j6,j7,za,zb,flip)
     &      +za(j2,j3)*W_tdkb_curr(j2,j3,j4,j5,j6,j7,za,zb,flip))
     &     /(za(j1,j3)*za(j2,j3)*s123)

      if (corr_on_beam == 1) then
      musq=renscale_beam1_islight_onlight**2
      elseif (corr_on_beam == 2) then
      musq=renscale_beam2_islight_onlight**2
      endif
      l12=lnrat(musq,-s(j1,j2))
      l123=lnrat(musq,-s123)

c      Vcc=-(epinv**2+epinv*l12+0.5_dp*l12**2)
c     & -2._dp*(epinv+l123)-4._dp
      Anlo(-2)=-Alo
      Anlo(-1)=-(l12+2._dp)*Alo
      Anlo( 0)=-(0.5_dp*l12**2+2._dp*l123+4._dp)*Alo

c--subleading N
      Anlo(0)=Anlo(0)
     & -(za(j2,j1)*W_tdkb_curr(j2,j1,j4,j5,j6,j7,za,zb,flip)
     &  +za(j2,j3)*W_tdkb_curr(j2,j3,j4,j5,j6,j7,za,zb,flip))
     &  /(za(j2,j3)*za(j3,j1)*s123)
     &  *Lsm1(-s(j1,j2),-s123,-s(j1,j3), -s123)
     & +(za(j1,j2)*za(j3,j1)*W_tdkb_curr(j2,j1,j4,j5,j6,j7,za,zb,flip)
     &  +za(j1,j2)*za(j3,j2)*W_tdkb_curr(j2,j2,j4,j5,j6,j7,za,zb,flip)
     &  -za(j1,j2)*za(j2,j3)*W_tdkb_curr(j2,j2,j4,j5,j6,j7,za,zb,flip)
     &  -za(j1,j3)*za(j2,j3)*W_tdkb_curr(j2,j3,j4,j5,j6,j7,za,zb,flip))
     &  /(za(j2,j3)*za(j1,j3)**2*s123)
     &  *Lsm1(-s(j1,j2),-s123,-s(j2,j3),-s123)

     & +2._dp*zb(j1,j3)/(za(j1,j3)*s123)
     & *(za(j1,j2)*W_tdkb_curr(j2,j2,j4,j5,j6,j7,za,zb,flip)
     &  +za(j1,j3)*W_tdkb_curr(j2,j3,j4,j5,j6,j7,za,zb,flip))
     & *L0(-s(j2,j3),-s123)/s123

c     & +2._dp*zb(j1,j3)*za(j1,j4)*za(j2,j4)/(za(j1,j3)*za(j4,j5))
c     & *L0(-s(j2,j3),-s123)/s123

c      Vsc=0.5_dp*(epinv+l123)+0.5_dp
      Anlo(-1)=Anlo(-1)+0.5_dp*Alo
      Anlo( 0)=Anlo( 0)+0.5_dp*(l123+1._dp)*Alo

      Anlo(0)=Anlo(0)
     & +(za(j1,j2)*W_tdkb_curr(j1,j2,j4,j5,j6,j7,za,zb,flip)
     &  +za(j1,j3)*W_tdkb_curr(j1,j3,j4,j5,j6,j7,za,zb,flip))
     &  *za(j2,j3)/(za(j1,j3)**3*s123)
     &  *Lsm1(-s(j1,j2),-s123,-s(j2,j3),-s123)

     & -0.5_dp*zb(j1,j3)**2*za(j2,j3)/(za(j1,j3)*s123)
     & *(za(j1,j2)*W_tdkb_curr(j1,j2,j4,j5,j6,j7,za,zb,flip)
     &  +za(j1,j3)*W_tdkb_curr(j1,j3,j4,j5,j6,j7,za,zb,flip))
     & *L1(-s123,-s(j2,j3))/s(j2,j3)**2

     & +za(j2,j3)*zb(j3,j1)/(za(j1,j3)**2*s123)
     & *(za(j1,j2)*W_tdkb_curr(j1,j2,j4,j5,j6,j7,za,zb,flip)
     &  +za(j1,j3)*W_tdkb_curr(j1,j3,j4,j5,j6,j7,za,zb,flip))
     & *L0(-s123,-s(j2,j3))/s(j2,j3)

     & +za(j2,j1)*zb(j1,j3)*W_tdkb_curr(j3,j3,j4,j5,j6,j7,za,zb,flip)/za(j1,j3)
     & *L1(-s123,-s(j1,j2))/s(j1,j2)**2

     & -za(j2,j1)*zb(j1,j3)/(za(j1,j3)**2*s123)
     & *(za(j1,j2)*W_tdkb_curr(j3,j2,j4,j5,j6,j7,za,zb,flip)
     &  +za(j1,j3)*W_tdkb_curr(j3,j3,j4,j5,j6,j7,za,zb,flip))
     & *L0(-s123,-s(j1,j2))/s(j1,j2)

     & +0.5_dp*zb(j1,j3)
     & *(zb(j3,j1)*W_tdkb_curr(j1,j2,j4,j5,j6,j7,za,zb,flip)
     &  +zb(j3,j2)*W_tdkb_curr(j2,j2,j4,j5,j6,j7,za,zb,flip))
     & /(zb(j1,j2)*zb(j2,j3)*za(j1,j3)*s123)

     & +0.5_dp*zb(j2,j3)
     & *(zb(j3,j1)*W_tdkb_curr(j1,j1,j4,j5,j6,j7,za,zb,flip)
     &  +zb(j3,j2)*W_tdkb_curr(j2,j1,j4,j5,j6,j7,za,zb,flip))
     & /(zb(j1,j2)*zb(j2,j3)*za(j1,j3)*s123)

c     & -0.5_dp*(za(j4,j1)*zb(j1,j3))**2*za(j2,j3)/(za(j1,j3)*za(j4,j5))
c     & *L1(-s123,-s(j2,j3))/s(j2,j3)**2

c     & +za(j1,j4)**2*za(j2,j3)*zb(j3,j1)/(za(j1,j3)**2*za(j4,j5))
c     & *L0(-s123,-s(j2,j3))/s(j2,j3)

c     & -za(j2,j1)*zb(j1,j3)*za(j4,j3)*zb(j3,j5)/za(j1,j3)
c     & *L1(-s123,-s(j1,j2))/s(j1,j2)**2

c     & -za(j2,j1)*zb(j1,j3)*za(j3,j4)*za(j1,j4)/(za(j1,j3)**2*za(j4,j5))
c     & *L0(-s123,-s(j1,j2))/s(j1,j2)

c     & -0.5_dp*zb(j3,j5)*(zb(j1,j3)*zb(j2,j5)+zb(j2,j3)*zb(j1,j5))
c     & /(zb(j1,j2)*zb(j2,j3)*za(j1,j3)*zb(j4,j5))

c      A52=(Vcc+Vsc)*A5lom+Fcc+Fsc

      Av=-(Anlo(-2)*epinv**2+Anlo(-1)*epinv+Anlo(0))

      return
      end subroutine Augd_tdkb_v2


      function w2jetsq_w_tdkb(i1,i2,i5,i6,in,ie,ib,ic,za,zb)
      implicit none
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'lc.f'
      include 'mmsq_cs.f'
      complex(dp):: qcd1(-1:1,-1:1),qcd2(-1:1,-1:1),qed(-1:1,-1:1)
      real(dp):: w2jetsq_w_tdkb,msq1,msq2,msqq
      integer, intent(in) :: i1,i2,i5,i6,in,ie,ib,ic

      call subqcd_w_tdkb(i1,i2,i5,i6,in,ie,ib,ic,za,zb,qcd1)
      call subqcd_w_tdkb(i1,i2,i6,i5,in,ie,ib,ic,za,zb,qcd2)

      qed(+1,+1)=qcd1(+1,+1)+qcd2(+1,+1)
      qed(+1,-1)=qcd1(+1,-1)+qcd2(-1,+1)
      qed(-1,+1)=qcd1(-1,+1)+qcd2(+1,-1)
      qed(-1,-1)=qcd1(-1,-1)+qcd2(-1,-1)

      msq1= abs(qcd1(+1,+1))**2+abs(qcd1(+1,-1))**2
     &     +abs(qcd1(-1,+1))**2+abs(qcd1(-1,-1))**2

      msq2= abs(qcd2(+1,+1))**2+abs(qcd2(+1,-1))**2
     &     +abs(qcd2(-1,+1))**2+abs(qcd2(-1,-1))**2

      msqq= abs( qed(+1,+1))**2+abs( qed(+1,-1))**2
     &     +abs( qed(-1,+1))**2+abs( qed(-1,-1))**2

      w2jetsq_w_tdkb=msq1+msq2-msqq/xnsq

c      write(6,*) 'w2jetsq_w_tdkb',w2jetsq_w_tdkb

      return
      end function w2jetsq_w_tdkb

      subroutine subqcd_w_tdkb(i1,i2,i5,i6,in,ie,ib,ic,za,zb,amp)
      implicit none
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
c*******************************************************************
c     the matrix elements of the
c     helicity amplitudes for the QCD process
c     q(-p1)+qbar(-p2) --> l(p3)+abar(p4)+g(p5)+g(p6)
c     multiplied by ((a+l)^2-M**2)/(a+l)^2/g^4
c     one colour ordering only
c     left-on quark line only
c*******************************************************************
c--- JC, Oct 8 2019
c---   generalized so that a current can be attached in place of the
c---   usual leptonic W decay
c---   in this case: l(p3)+abar(p4) is replaced by
c---             t (nu(in) + e+(ie) + b(ib)) + bbar(ic)
      logical:: flip
      integer, intent(in) :: i1,i2,i5,i6,in,ie,ib,ic
      integer :: i3,i4
      real(dp):: s156,s56,s256,s1256
      complex(dp), intent(out) :: amp(-1:1,-1:1)
      complex(dp):: zba2,zba3
c begin statement functions
      zba2(i1,i2,i3,i4)=zb(i1,i2)*za(i2,i4)+zb(i1,i3)*za(i3,i4)
      zba3(i1,i2,i3,i4,i5)=
     & zb(i1,i2)*za(i2,i5)+zb(i1,i3)*za(i3,i5)+zb(i1,i4)*za(i4,i5)
c end statement functions

      s156=s(i1,i5)+s(i1,i6)+s(i5,i6)
      s256=s(i2,i5)+s(i2,i6)+s(i5,i6)
      s56=s(i5,i6)
      s1256=s156+s256+s(i1,i2)-s56

      flip=.false.

c      amp(-1,-1)=-four/s1256
c     & *zb(i1,i4)*zba3(i1,i2,i6,i5,i3)/zb(i1,i5)/zb(i2,i6)/zb(i6,i5)
      amp(-1,-1)=-four/s1256*(
     & +zb(i1,i2)*W_tdkb_curr(i2,i1,in,ie,ib,ic,za,zb,flip)
     & +zb(i1,i6)*W_tdkb_curr(i6,i1,in,ie,ib,ic,za,zb,flip)
     & +zb(i1,i5)*W_tdkb_curr(i5,i1,in,ie,ib,ic,za,zb,flip)
     & )/zb(i1,i5)/zb(i2,i6)/zb(i6,i5)

c      amp(+1,+1)=+four/s1256
c     & *za(i2,i3)*zba3(i4,i1,i6,i5,i2)/za(i1,i5)/za(i2,i6)/za(i6,i5)
      amp(+1,+1)=+four/s1256*(
     & -za(i1,i2)*W_tdkb_curr(i2,i1,in,ie,ib,ic,za,zb,flip)
     & -za(i6,i2)*W_tdkb_curr(i2,i6,in,ie,ib,ic,za,zb,flip)
     & -za(i5,i2)*W_tdkb_curr(i2,i5,in,ie,ib,ic,za,zb,flip)
     & )/za(i1,i5)/za(i2,i6)/za(i6,i5)

c      amp(+1,-1)=four/(s1256*s(i5,i6))
c     & *(zba2(i4,i1,i5,i6)*zba2(i5,i2,i6,i3)/za(i1,i5)/zb(i2,i6)
c     & -za(i1,i6)*za(i2,i3)*zb(i1,i5)*zba2(i4,i1,i5,i6)/za(i1,i5)/s156
c     & -za(i2,i6)*zb(i1,i4)*zb(i2,i5)*zba2(i5,i2,i6,i3)/zb(i2,i6)/s256)
      amp(+1,-1)=four/(s1256*s(i5,i6))*((
     &  -za(i1,i6)*zb(i5,i2)*W_tdkb_curr(i2,i1,in,ie,ib,ic,za,zb,flip)
     &  -za(i1,i6)*zb(i5,i6)*W_tdkb_curr(i6,i1,in,ie,ib,ic,za,zb,flip)
     &  -za(i5,i6)*zb(i5,i2)*W_tdkb_curr(i2,i5,in,ie,ib,ic,za,zb,flip)
     &  -za(i5,i6)*zb(i5,i6)*W_tdkb_curr(i6,i5,in,ie,ib,ic,za,zb,flip)
     &    )/za(i1,i5)/zb(i2,i6)
     & -za(i1,i6)*zb(i1,i5)*(
     &  -za(i1,i6)*W_tdkb_curr(i2,i1,in,ie,ib,ic,za,zb,flip)
     &  -za(i5,i6)*W_tdkb_curr(i2,i5,in,ie,ib,ic,za,zb,flip)
     &    )/za(i1,i5)/s156
     & -za(i2,i6)*zb(i2,i5)*(
     &  +zb(i5,i2)*W_tdkb_curr(i2,i1,in,ie,ib,ic,za,zb,flip)
     &  +zb(i5,i6)*W_tdkb_curr(i6,i1,in,ie,ib,ic,za,zb,flip)
     &  )/zb(i2,i6)/s256)

c      amp(-1,+1)=four/(s1256*s(i5,i6))
c     & *(za(i2,i3)*za(i2,i5)*zb(i1,i4)*zb(i1,i6)/za(i2,i6)/zb(i1,i5)
c     & -za(i2,i3)*zb(i1,i6)**2*zba2(i4,i1,i6,i5)/zb(i1,i5)/s156
c     & -za(i2,i5)**2*zb(i1,i4)*zba2(i6,i2,i5,i3)/za(i2,i6)/s256)
      amp(-1,+1)=four/(s1256*s(i5,i6))
     & *(za(i2,i5)*zb(i1,i6)/za(i2,i6)/zb(i1,i5)*W_tdkb_curr(i2,i1,in,ie,ib,ic,za,zb,flip)
     & -zb(i1,i6)**2*(
     &   -za(i1,i5)*W_tdkb_curr(i2,i1,in,ie,ib,ic,za,zb,flip)
     &   -za(i6,i5)*W_tdkb_curr(i2,i6,in,ie,ib,ic,za,zb,flip)
     &  )/zb(i1,i5)/s156
     & -za(i2,i5)**2*(
     &   +zb(i6,i2)*W_tdkb_curr(i2,i1,in,ie,ib,ic,za,zb,flip)
     &   +zb(i6,i5)*W_tdkb_curr(i5,i1,in,ie,ib,ic,za,zb,flip)
     &  )/za(i2,i6)/s256)


      return
      end subroutine subqcd_w_tdkb


      function ampsq_udqqb_tdkb(j1,j2,j3,j4,j5,j6,j7,j8,za,zb)
c matrix elements squared for the process
c 0 -> ubar(i1) + d(i6) + q(i7) + qbar(i8) + t (nu(3) + e+(i4) + b(i5)) + bbar(i2)

c where q(i7) is not a down quark, i.e. the non-identical process

c where a W is radiated from the (d,ubar) line and connects to the (t, bbar) line,
c and quark helicity(i7) is summed over, all other helicities are fixed (left-handed i5, i6)

c remaining overall factor is [(s1678/(s1678-MW^2))*gw^4*gs^2]^2*V*xn
      implicit none
      include 'cplx.h'
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer, intent(in) :: j1,j2,j3,j4,j5,j6,j7,j8
      real(dp):: ampsq_udqqb_tdkb

      ampsq_udqqb_tdkb=cdabs(amp_udqqb_tdkb(j1,j2,j3,j4,j5,j6,j7,j8,za,zb))**2
     &            +cdabs(amp_udqqb_tdkb(j1,j2,j3,j4,j5,j6,j8,j7,za,zb))**2

      return
      end function ampsq_udqqb_tdkb


      function ampsq_udddb_tdkb(j1,j2,j3,j4,j5,j6,j7,j8,za,zb)
c matrix elements squared for the process
c 0 -> ubar(i1) + d(i6) + d(i7) + dbar(i8) + t (nu(3) + e+(i4) + b(i5)) + bbar(i2)

c i.e. an identical-quark case of the amplitude-squared above

c where a W is radiated from the (d,ubar) line and connects to the (t, bbar) line,
c and quark helicity(i7) is summed over, all other helicities are fixed (left-handed i5, i6)

c remaining overall factor is [(s1678/(s1678-MW^2))*gw^4*gs^2]^2*V*xn
      implicit none
      include 'cplx.h'
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer, intent(in) :: j1,j2,j3,j4,j5,j6,j7,j8
      real(dp):: ampsq_udddb_tdkb
      complex(dp):: amporig,ampswap

c left-handed amplitude must also include interference between the two diagrams
      amporig=amp_udqqb_tdkb(j1,j2,j3,j4,j5,j6,j7,j8,za,zb)
      ampswap=amp_udqqb_tdkb(j1,j2,j3,j4,j5,j7,j6,j8,za,zb)
      ampsq_udddb_tdkb=cdabs(amporig)**2+cdabs(ampswap)**2+two/xn*real(amporig*conjg(ampswap),dp)

c right-handed amplitude cannot interfere
      amporig=amp_udqqb_tdkb(j1,j2,j3,j4,j5,j6,j8,j7,za,zb)
      ampswap=amp_udqqb_tdkb(j1,j2,j3,j4,j5,j7,j8,j6,za,zb)
      ampsq_udddb_tdkb=ampsq_udddb_tdkb+cdabs(amporig)**2+cdabs(ampswap)**2

      return
      end function ampsq_udddb_tdkb


      function ampsq_uduub_tdkb(j1,j2,j3,j4,j5,j6,j7,j8,za,zb)
c matrix elements squared for the process
c 0 -> ubar(i1) + d(i6) + u(i7) + ubar(i8) + t (nu(3) + e+(i4) + b(i5)) + bbar(i2)

c i.e. an identical-quark case of the amplitude-squared above

c where a W is radiated from the (d,ubar) line and connects to the (t, bbar) line,
c and quark helicity(i7) is summed over, all other helicities are fixed (left-handed i5, i6)

c remaining overall factor is [(s1678/(s1678-MW^2))*gw^4*gs^2]^2*V*xn
      implicit none
      include 'cplx.h'
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer, intent(in) :: j1,j2,j3,j4,j5,j6,j7,j8
      real(dp):: ampsq_uduub_tdkb
      complex(dp):: amporig,ampswap

c left-handed amplitude must also include interference between the two diagrams
      amporig=amp_udqqb_tdkb(j1,j2,j3,j4,j5,j6,j7,j8,za,zb)
      ampswap=amp_udqqb_tdkb(j8,j2,j3,j4,j5,j6,j7,j1,za,zb)
      ampsq_uduub_tdkb=cdabs(amporig)**2+cdabs(ampswap)**2+two/xn*real(amporig*conjg(ampswap),dp)

c right-handed amplitude cannot interfere
      amporig=amp_udqqb_tdkb(j1,j2,j3,j4,j5,j6,j8,j7,za,zb)
      ampswap=amp_udqqb_tdkb(j8,j2,j3,j4,j5,j6,j1,j7,za,zb)
      ampsq_uduub_tdkb=ampsq_uduub_tdkb+cdabs(amporig)**2+cdabs(ampswap)**2

      return
      end function ampsq_uduub_tdkb


      function amp_udqqb_tdkb(j1,j2,j3,j4,j5,j6,j7,j8,za,zb)
c amplitude for the process
c 0 -> ubar(i1) + d(i6) + q(i7) + qbar(i8) + t (nu(3) + e+(i4) + b(i5)) + bbar(i2)

c where a W is radiated from the (d,ubar) line and connects to the (t, bbar) line,
c and quark(i7) is left-handed and all the other helicities are fixed (left-handed i5, i6)

c remaining overall factor is (s1678/(s1678-MW^2))*gw^4*gs^2
      implicit none
      include 'cplx.h'
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'masses.f'
      integer, intent(in) ::  j1,j2,j3,j4,j5,j6,j7,j8
      real(dp):: s78,s178,s678,s1678
      complex(dp):: amp_udqqb_tdkb,zab2,propw

c--- statement function
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      propw(s78)=cone/cplx2((s78-wmass**2),wmass*wwidth)

      s1678=s(j1,j6)+s(j1,j7)+s(j1,j8)+s(j6,j7)+s(j6,j8)+s(j7,j8)
      s178=s(j1,j7)+s(j1,j8)+s(j7,j8)
      s678=s(j6,j7)+s(j6,j8)+s(j7,j8)
      s78=s(j7,j8)

      amp_udqqb_tdkb=propw(s(j3,j4))/(im*mt*twidth)*za(j5,j3)/(s78*s1678)*(
     &  za(j6,j7)/s678*(zb(j8,j6)*zab2(j6,j3,j5,j4)+zb(j8,j7)*zab2(j7,j3,j5,j4))*zb(j2,j1)
     & -zb(j8,j1)/s178*zab2(j6,j3,j5,j4)*zab2(j7,j1,j8,j2))

      return
      end function amp_udqqb_tdkb

      function ampsq_ugd_tdkb_n(j1,j2,j3,j4,j5,j6,j7,zanb)
c--- Amplitude squared for 0 -> ubar(j1) + g(j7) + d(j6) [light line]
c---                     + t (nu(j3) + e+(j4) + b(j5)) + bbar(j2)  [heavy line]
c---
c--- but with the polarization vector of g(j7) replaced by an arbitrary vector n
c--- that satisfies n.g = 0; spinor products <a|n|b] are passed in via zanb
c    overall factor is [(s167/(s167-MW^2))*gw^4*gs]^2*V*xn^2
      implicit none
      include 'cplx.h'
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'masses.f'
      real(dp) :: ampsq_ugd_tdkb_n
      integer, intent(in) ::  j1,j2,j3,j4,j5,j6,j7
      complex(dp), intent(in) :: zanb(mxpart,mxpart)

      real(dp) :: s167
      complex(dp) :: amp,zab2,propW

c--- statement function
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      propw(s167)=cone/cplx2((s167-wmass**2),wmass*wwidth)

      s167=s(j1,j6)+s(j1,j7)+s(j6,j7)
      amp=propW(s(j3,j4))/(im*mt*twidth)*za(j5,j3)/rt2*(
     & (zanb(j6,j6)*zab2(j6,j3,j5,j4)+zanb(j6,j7)*zab2(j7,j3,j5,j4))*zb(j2,j1)/s(j6,j7)
     & -zab2(j6,j3,j5,j4)*(zb(j2,j1)*zanb(j1,j1)+zb(j2,j7)*zanb(j7,j1))/s(j1,j7))

      ampsq_ugd_tdkb_n=cdabs(amp/(s167-wmass**2))**2

      return
      end function ampsq_ugd_tdkb_n

      subroutine singletop_jet_light_gs_all(p,ndmx,msqall)
        use singletop2_nnlo_vars
        use dipoles, only: dips_new
        implicit none
        include 'constants.f'
        include 'nf.f'
        include 'mxpart.f'
        include 'ptilde.f'
        include 'qqgg.f'
        real(dp), intent(in) :: p(mxpart,4)
        integer, intent(in) :: ndmx
        real(dp), intent(out) :: msqall(ndmx,-nf:nf,-nf:nf,max_bcontrib,max_corr_on_beam)

        real(dp) :: msq(-nf:nf,-nf:nf), msqv(-nf:nf,-nf:nf)
        real(dp) :: sub(4), subv

        real(dp), parameter :: nfmhalf = real(nf,dp) - 0.5_dp
        logical :: quark, gluon

#define noglue [-3,-1,2,4]

      ndmax=32

      corr_islight = .true.
      msqall = 0._dp
      quark = iand(partons_enabled, quarkChannel) > 0
      gluon = iand(partons_enabled, gluonChannel) > 0

c Note that swaps are required to get the outgoing quark and antiquark in the canonical order

      if (any(beams_enabled == 1)) then

      corr_beam1=.true.
      corr_on_beam=1

      if (quark) then
        call dips_new(1,p,1,7,6,sub,msq,singletop_jet_light)
        msqall(1,noglue,5, 1,1) = -msq(noglue,5)*(sub(qq))/xn*half
      endif

      call dips_new(2,p,1,8,6,sub,msq,singletop_jet_light,subv,msqv,singletop_jet_light_gvec)
      if (quark) msqall(2,noglue,5, 1,1)=-msq(noglue,5)*sub(qq)/xn*half
      if (gluon) msqall(2,0,5, 1,1)=xn*(msq(0,5)*sub(gg)+msqv(0,5)*subv)

      call dips(3,p,1,7,8,sub,subv,msq,msqv,singletop_jet_light,singletop_jet_light_gvec)
      if (quark) msqall(3,noglue,5, 1,1)=xn*msq(noglue,5)*sub(qq)*half
      if (gluon) msqall(3,0,5, 1,1)=two*tr*(msq(2,5)+msq(4,5))*sub(qg)
      if (quark) msqall(3,[2,4],5, 1,1)=msqall(3,[2,4],5, 1,1)+two*Cf*(msq(0,5)*sub(gq)+msqv(0,5)*subv)

      call dips(25,p,1,8,7,sub,subv,msq,msqv,singletop_jet_light,singletop_jet_light_gvec) !3, finished
      if (quark) msqall(25,noglue,5, 1,1)=xn*msq(noglue,5)*sub(qq)*half
      if (gluon) msqall(25,0,5, 1,1)=xn*(msq(0,5)*sub(gg)+msqv(0,5)*subv)
      if (quark) msqall(25,[1,3],5, 1,1)=msqall(25,[1,3],5, 1,1)+1.5_dp*Cf*(msq(0,5)*sub(gq)+msqv(0,5)*subv)
      if (quark) msqall(25,5,5, 1,1)=two*Cf*(msq(0,5)*sub(gq)+msqv(0,5)*subv)

      if (quark) then
        call dips_new(27,p,6,7,1,sub,msq,singletop_jet_light) ! 1, finished
        msqall(27,noglue,5, 1,1)=-msq(noglue,5)*(sub(qq))/xn*half
      endif

      call dips_new(28,p,6,8,1,sub,msq,singletop_jet_light) ! 2, finished
      if (quark) msqall(28,noglue,5, 1,1)=-msq(noglue,5)*(sub(qq))/xn*half
      if (gluon) msqall(28,0,5, 1,1)=xn*(msq(0,5)*(sub(qq)))

      call dips(29,p,7,8,1,sub,subv,msq,msqv,singletop_jet_light,singletop_jet_light_gvec) ! 3, finished
      if (quark) msqall(29,noglue,5, 1,1)=xn*(msq(noglue,5)*(half*sub(gg))+half*msqv(noglue,5)*subv
     &            +msq(noglue,5)*half*sub(gg)+half*msqv(noglue,5)*subv)*half
     &            +tr*nfmhalf*(msq(noglue,5)*sub(gq)-msqv(noglue,5)*subv)
      if (gluon) msqall(29,0,5, 1,1)=xn*msq(0,5)*sub(qq)

      if (quark) then
        call dips(10,p,8,6,1,sub,subv,msq,msqv,singletop_jet_light,singletop_jet_light_gvec)
        msqall(10,noglue,5, 1,1)=tr*half*(msq(noglue,5)*sub(gq)-msqv(noglue,5)*subv)
      endif

      if (quark) then
        call dips(12,p,1,7,8,sub,subv,msq,msqv,singletop_jet_light_swap,singletop_jet_light_gvec_swap)
        msqall(12,[-3,-1],5, 1,1)=two*Cf*(msq(0,5)*sub(gq)+msqv(0,5)*subv)
      endif

      if (quark) then
        call dips(13,p,1,6,8,sub,subv,msq,msqv,singletop_jet_light_swap,singletop_jet_light_gvec_swap)
        msqall(13,[1,3],5, 1,1)=half*Cf*(msq(0,5)*sub(gq)+msqv(0,5)*subv)
      endif

      if (quark) then
        call dips(14,p,1,8,7,sub,subv,msq,msqv,singletop_jet_light_swap,singletop_jet_light_gvec_swap)
        msqall(14,[-4,-2],5, 1,1)=1.5_dp*Cf*(msq(0,5)*sub(gq)+msqv(0,5)*subv)
        msqall(14,-5,5, 1,1)=two*Cf*(msq(0,5)*sub(gq)+msqv(0,5)*subv)
      endif

      call dips(15,p,1,6,8,sub,subv,msq,msqv,singletop_jet_light,singletop_jet_light_gvec)
      if (gluon) msqall(15,0,5, 1,1)=two*tr*(msq(-1,5)+msq(-3,5))*sub(qg)
      if (quark) msqall(15,[-4,-2],5, 1,1)=half*Cf*(msq(0,5)*sub(gq)+msqv(0,5)*subv)

      if (quark) then
        call dips_new(7,p,6,7,8,sub,msq,singletop_jet_light)
        msqall(7,noglue,5, 1,1)=xn*msq(noglue,5)*sub(qq)*half
      endif

      call dips_new(8,p,6,8,7,sub,msq,singletop_jet_light)
      if (quark) msqall(8,noglue,5, 1,1)=xn*msq(noglue,5)*sub(qq)*half
      if (gluon) msqall(8,0,5, 1,1)=-msq(0,5)*sub(qq)/xn

      call dips(9,p,7,8,6,sub,subv,msq,msqv,singletop_jet_light,singletop_jet_light_gvec)
      if (quark) msqall(9,noglue,5, 1,1)=xn*(msq(noglue,5)*sub(gg)+msqv(noglue,5)*subv)*half
     &                   +tr*nfmhalf*(msq(noglue,5)*sub(gq)-msqv(noglue,5)*subv)
      if (gluon) msqall(9,0,5, 1,1)=-msq(0,5)*sub(qq)/xn

      if (quark) then
        call dips(11,p,8,6,7,sub,subv,msq,msqv,singletop_jet_light,singletop_jet_light_gvec)
        msqall(11,noglue,5, 1,1)=tr*half*(msq(noglue,5)*sub(gq)-msqv(noglue,5)*subv)
      endif

      endif

      if (any(beams_enabled == 2)) then

      corr_beam1=.false.
      corr_on_beam=2

      if (quark) then
        call dips_new(4,p,2,7,6,sub,msq,singletop_jet_light)
        msqall(4,5,noglue, 1,2)=-msq(5,noglue)*(sub(qq))/xn*half
      endif

      call dips(5,p,2,8,6,sub,subv,msq,msqv,singletop_jet_light,singletop_jet_light_gvec)
      if (quark) msqall(5,5,noglue, 1,2)=-msq(5,noglue)*(sub(qq))/xn*half
      if (gluon) msqall(5,5,0, 1,2)=xn*(msq(5,0)*(sub(gg))+msqv(5,0)*subv)

      call dips(6,p,2,7,8,sub,subv,msq,msqv,singletop_jet_light,singletop_jet_light_gvec)
      if (quark) msqall(6,5,noglue, 1,2)=xn*(msq(5,noglue)*(sub(qq)))*half
      if (quark) msqall(6,5,[2,4], 1,2)=msqall(6,5,[2,4], 1,2)+two*Cf*(msq(5,0)*sub(gq)+msqv(5,0)*subv)
      if (gluon) msqall(6,5,0, 1,2)=two*tr*(msq(5,2)+msq(5,4))*sub(qg)

      call dips(26,p,2,8,7,sub,subv,msq,msqv,singletop_jet_light,singletop_jet_light_gvec) !6, finished
      if (quark) msqall(26,5,noglue, 1,2)=xn*msq(5,noglue)*sub(qq)*half
      if (quark) msqall(26,5,[1,3], 1,2)=msqall(26,5,[1,3], 1,2)+1.5_dp*Cf*(msq(5,0)*sub(gq)+msqv(5,0)*subv)
      if (quark) msqall(26,5,5, 1,2)=two*Cf*(msq(5,0)*sub(gq)+msqv(5,0)*subv)
      if (gluon) msqall(26,5,0, 1,2)=xn*(msq(5,0)*(sub(gg))+msqv(5,0)*subv)

      if (quark) then
        call dips_new(30,p,6,7,2,sub,msq,singletop_jet_light) ! 4, finished
        msqall(30,5,noglue, 1,2)=-msq(5,noglue)*(sub(qq))/xn*half
      endif

      call dips_new(31,p,6,8,2,sub,msq,singletop_jet_light) ! 5, finished
      if (quark) msqall(31,5,noglue, 1,2)=-msq(5,noglue)*(sub(qq))/xn*half
      if (gluon) msqall(31,5,0, 1,2)=xn*(msq(5,0)*(sub(qq)))

      call dips(32,p,7,8,2,sub,subv,msq,msqv,singletop_jet_light,singletop_jet_light_gvec) ! 6, finished
      if (quark) msqall(32,5,noglue, 1,2)=xn*(msq(5,noglue)*(half*sub(gg))+half*msqv(5,noglue)*subv
     &            +msq(5,noglue)*half*sub(gg)+half*msqv(5,noglue)*subv)*half
     &            +tr*nfmhalf*(msq(5,noglue)*sub(gq)-msqv(5,noglue)*subv)
      if (gluon) msqall(32,5,0, 1,2)=xn*(msq(5,0)*(sub(qq)))

      call dips(16,p,2,6,8,sub,subv,msq,msqv,singletop_jet_light,singletop_jet_light_gvec)
      if (quark) msqall(16,5,[-4,-2], 1,2)=half*Cf*(msq(5,0)*sub(gq)+msqv(5,0)*subv)
      if (gluon) msqall(16,5,0, 1,2)=two*tr*(msq(5,-1)+msq(5,-3))*sub(qg)

      if (quark) then
        call dips(17,p,8,6,2,sub,subv,msq,msqv,singletop_jet_light,singletop_jet_light_gvec)
        msqall(17,5,noglue, 1,2)=tr*half*(msq(5,noglue)*sub(gq)-msqv(5,noglue)*subv)
      endif

      if (quark) then
        call dips(18,p,2,7,8,sub,subv,msq,msqv,singletop_jet_light_swap,singletop_jet_light_gvec_swap)
        msqall(18,5,[-3,-1], 1,2)=two*Cf*(msq(5,0)*sub(gq)+msqv(5,0)*subv)
      endif

      if (quark) then
        call dips(19,p,2,6,8,sub,subv,msq,msqv,singletop_jet_light_swap,singletop_jet_light_gvec_swap)
        msqall(19,5,[1,3], 1,2)=half*Cf*(msq(5,0)*sub(gq)+msqv(5,0)*subv)
      endif

      if (quark) then
        call dips(20,p,2,8,7,sub,subv,msq,msqv,singletop_jet_light_swap,singletop_jet_light_gvec_swap)
        msqall(20,5,[-4,-2], 1,2)=1.5_dp*Cf*(msq(5,0)*sub(gq)+msqv(5,0)*subv)
        msqall(20,5,-5, 1,2)=two*Cf*(msq(5,0)*sub(gq)+msqv(5,0)*subv)
      endif

      if (quark) then
        call dips_new(21,p,6,7,8,sub,msq,singletop_jet_light) ! 7, finished
        msqall(21,5,noglue, 1,2)=xn*msq(5,noglue)*sub(qq)*half
      endif

      call dips_new(22,p,6,8,7,sub,msq,singletop_jet_light) ! 8, finished
      if (quark) msqall(22,5,noglue, 1,2)=xn*msq(5,noglue)*sub(qq)*half
      if (gluon) msqall(22,5,0, 1,2)=-msq(5,0)*sub(qq)/xn

      call dips(23,p,7,8,6,sub,subv,msq,msqv,singletop_jet_light,singletop_jet_light_gvec) ! 9, finished
      if (quark) msqall(23,5,noglue, 1,2)=xn*(msq(5,noglue)*sub(gg)+msqv(5,noglue)*subv)*half
     &                           +tr*nfmhalf*(msq(5,noglue)*sub(gq)-msqv(5,noglue)*subv)
      if (gluon) msqall(23,5,0, 1,2)=-msq(5,0)*sub(qq)/xn

      if (quark) then
        call dips(24,p,8,6,7,sub,subv,msq,msqv,singletop_jet_light,singletop_jet_light_gvec) ! 11, finished
        msqall(24,5,noglue, 1,2)=tr*half*(msq(5,noglue)*sub(gq)-msqv(5,noglue)*subv)
      endif

      endif

      end subroutine singletop_jet_light_gs_all

      subroutine singletop_jet_light_z(p,z)
      use singletop2_nnlo_vars
      implicit none
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'PR_new.f'
      include 'PR_stop.f'
      include 'agq.f'
      integer:: is
      real(dp), intent(in) :: p(mxpart,4),z
      real(dp):: dot,xl16,xl17,xl26,xl27,ff_qq,ff_gg,fi_qq,fi_gg,if_qq,if_gg,if_gq,if_qg
      real(dp) :: xl67_beam1, xl67_beam2

      xl16=log(-two*dot(p,1,6)/renscale_beam1_islight_onlight**2)
      xl17=log(-two*dot(p,1,7)/renscale_beam1_islight_onlight**2)

      xl26=log(-two*dot(p,2,6)/renscale_beam2_islight_onlight**2)
      xl27=log(-two*dot(p,2,7)/renscale_beam2_islight_onlight**2)

      !if (corr_on_beam == 1) then
      xl67_beam1 = log(+two*dot(p,6,7)/renscale_beam1_islight_onlight**2)
      !elseif (corr_on_beam == 2) then
      xl67_beam2 = log(+two*dot(p,6,7)/renscale_beam2_islight_onlight**2)
      !endif

c 2-quark, 2-gluon contribution and final-state g->qqb splittings from 4-quark contribution
      do is=1,3
      !if (corr_on_beam == 1) then
      B1(q,q,b,is) = (as_light_beam1/4._dp/pi)*xn*(if_qq(z,xl17,is)+half*fi_gg(z,xl17,is)
     &                    +ff_qq(z,xl67_beam1,is)+half*ff_gg(z,xl67_beam1,is))
     &            -(as_light_beam1/4._dp/pi)/xn*(if_qq(z,xl16,is)+fi_qq(z,xl16,is))
      B1(g,g,b,is) = (as_light_beam1/4._dp/pi)*xn*(if_gg(z,xl16,is)+fi_qq(z,xl16,is)
     &                    +if_gg(z,xl17,is)+fi_qq(z,xl17,is))
     &            -(as_light_beam1/4._dp/pi)/xn*(two*ff_qq(z,xl67_beam1,is))
      !endif

      !if (corr_on_beam == 2) then
      B2(q,q,b,is) = (as_light_beam2/4._dp/pi)*xn*(if_qq(z,xl27,is)+half*fi_gg(z,xl27,is)
     &                    +ff_qq(z,xl67_beam2,is)+half*ff_gg(z,xl67_beam2,is))
     &            -(as_light_beam2/4._dp/pi)/xn*(if_qq(z,xl26,is)+fi_qq(z,xl26,is))
      B2(g,g,b,is) = (as_light_beam2/4._dp/pi)*xn*(if_gg(z,xl26,is)+fi_qq(z,xl26,is)
     &                    +if_gg(z,xl27,is)+fi_qq(z,xl27,is))
     &            -(as_light_beam2/4._dp/pi)/xn*(two*ff_qq(z,xl67_beam2,is))
      !endif
      enddo

c off diagonal terms from 4-quark contribution
      do is=1,3
      !if (corr_on_beam == 1) then
      B1(g,q,b,is) = (as_light_beam1/4._dp/pi)*two*Cf*if_gq(z,xl17,is)
      B1(q,g,b,is) = (as_light_beam1/4._dp/pi)*two*Tr*if_qg(z,xl17,is)
      !endif

      !if (corr_on_beam == 2) then
      B2(g,q,b,is) = (as_light_beam2/4._dp/pi)*two*Cf*if_gq(z,xl27,is)
      B2(q,g,b,is) = (as_light_beam2/4._dp/pi)*two*Tr*if_qg(z,xl27,is)
      !endif
      enddo

      return
      end subroutine singletop_jet_light_z



      end module singletop_jet


