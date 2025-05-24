!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

c--- Author: J. M. Campbell, September 2013
c--- Effect of interference between gg -> H -> ZZ signal process
c--- and gg -> ZZ NNLO contribution to continuum background
c--- The effect of massive bottom and top quark loops is included
      subroutine gg_ZZ_int(p,msq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'noglue.f'
      include 'interference.f'
      include 'zcouple_cms.f'
      include 'bsm_higgs.f'
      integer:: h1,h2,h34,h56,j,k
      integer,parameter:: i4(2)=(/4,5/),i5(2)=(/5,4/),
     & jkswitch(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
      real(dp):: p(mxpart,4),wt,wt2,msq(fn:nf,fn:nf),msqgg,fac,
     & pswap(mxpart,4),oprat
      complex(dp)::
     & Mloop_uptype(2,2,2,2),Mloop_dntype(2,2,2,2),
     & Mloop_bquark(2,2,2,2),Mloop_tquark(2,2,2,2),
     & Sloop_uptype(2,2,2,2),Sloop_dntype(2,2,2,2),
     & Sloop_bquark(2,2,2,2),Sloop_tquark(2,2,2,2),
     & ggH_bquark(2,2,2,2),ggH_tquark(2,2,2,2),Acont,AHiggs,
     & ggH_bquark_swap(2,2,2,2),ggH_tquark_swap(2,2,2,2),AHiggs_swap,
     & Acont_swap,
     & Mamp,Samp,
     & Mloop_c6_propagator(2,2,2,2),
     & Mloop_c6_propagator_swap(2,2,2,2),
     & AHiggs_c6,AHiggs_c6_swap,
     & Mloop_c6_decay(2,2,2,2),
     & Mloop_c6_decay_swap(2,2,2,2),
     & Mloop_c6_production(2,2,2,2),
     & Mloop_c6_production_swap(2,2,2,2),
     & Mloop_c6_width(2,2,2,2),
     & Mloop_c6_width_swap(2,2,2,2),
     & Mloop_SMEFT(2,2,2,2),
     & Mloop_SMEFT_swap(2,2,2,2),
     & AHiggs_ct,AHiggs_ct_swap,
     & AHiggs_cg,AHiggs_cg_swap,
     & AHiggs_SMEFT,AHiggs_SMEFT_swap
      logical:: includegens1and2,includebottom,includetop

c--- set this to true to include generations 1 and 2 of (light) quarks
      includegens1and2=.true.
c--- set this to true to include massive bottom quark
      includebottom=.true.
c--- set this to true to include massive top quark
      includetop=.true.

c--- if neither contribution is included print warning message and stop
      if ((includegens1and2 .eqv. .false.) .and.
     &    (includebottom    .eqv. .false.) .and.
     &    (includetop       .eqv. .false.)) then
         write(6,*) 'Box loop is set to zero, please edit gg_ZZ_int.f'
         stop
      endif
c--- if noglue print warning message and stop
      if (noglue) then
         write(6,*) 'Please set noglue .false. in input file'
         stop
      endif

      msq(:,:)=0._dp

c      if (pttwo(3,4,p) < 7._dp) return ! Kauer gg2VV cut on |H+C|^2

      call getggZZamps(p,includegens1and2,includebottom,includetop,
     & Mloop_uptype,Mloop_dntype,Mloop_bquark,Mloop_tquark)

      call getggHZZamps(p,ggH_bquark,ggH_tquark,Mloop_c6_propagator,Mloop_c6_decay,Mloop_c6_production,Mloop_c6_width,Mloop_SMEFT)

      if (interference) then
c--- for interference, compute amplitudes after 4<->6 swap
       pswap(1,:)=p(1,:)
       pswap(2,:)=p(2,:)
       pswap(3,:)=p(3,:)
       pswap(4,:)=p(6,:)
       pswap(5,:)=p(5,:)
       pswap(6,:)=p(4,:)
       call getggZZamps(pswap,includegens1and2,includebottom,includetop,
     &  Sloop_uptype,Sloop_dntype,Sloop_bquark,Sloop_tquark)
       call getggHZZamps(pswap,ggH_bquark_swap,ggH_tquark_swap,Mloop_c6_propagator_swap,Mloop_c6_decay_swap,Mloop_c6_production_swap,Mloop_c6_width_swap,Mloop_SMEFT_swap)
      endif

      msqgg=0._dp
      do h1=1,2
      do h2=1,2
      do h34=1,2
      do h56=1,2

c--- compute total continuum amplitude
      Acont=
     &  +two*Mloop_uptype(h1,h2,h34,h56)
     &  +two*Mloop_dntype(h1,h2,h34,h56)
     &      +Mloop_bquark(h1,h2,h34,h56)
     &      +Mloop_tquark(h1,h2,h34,h56)
c--- compute total Higgs amplitude
      AHiggs=
     &  +ggH_bquark(h1,h2,h34,h56)
     &  +ggH_tquark(h1,h2,h34,h56)
c--- compute total c6 correction to Higgs amplitude        
      AHiggs_c6=
     &  +Mloop_c6_propagator(h1,h2,h34,h56)
     &  +Mloop_c6_decay(h1,h2,h34,h56)
     &  +Mloop_c6_production(h1,h2,h34,h56)
     &  +Mloop_c6_width(h1,h2,h34,h56)
c--- compute total ct correction to Higgs amplitude
      AHiggs_ct=ct*ggH_tquark(h1,h2,h34,h56)
c--- compute total cg correction to Higgs amplitude
      AHiggs_cg=cg*Mloop_SMEFT(h1,h2,h34,h56)
c--- compute total SMEFT correction to Higgs amplitude
      AHiggs_SMEFT=AHiggs_c6+AHiggs_ct+AHiggs_cg
c-------------------------------      

c---- This only accumulates contributions from the interference
      
      if (interference .eqv. .false.) then
c---  normal case
        msqgg=msqgg
     &        +abs(Acont+AHiggs+AHiggs_SMEFT)**2
     &        -abs(Acont)**2
     &        -abs(AHiggs+AHiggs_SMEFT)**2
!     &        +two*real(conjg(Acont)*AHiggs)
!     &        +two*real(conjg(AHiggs)*AHiggs_c6)
!     &        +two*real(conjg(Acont)*AHiggs_c6)
      else
c--- with interference
        Acont_swap=
     &  +two*Sloop_uptype(h1,h2,h34,h56)
     &  +two*Sloop_dntype(h1,h2,h34,h56)
     &      +Sloop_bquark(h1,h2,h34,h56)
     &      +Sloop_tquark(h1,h2,h34,h56)
        AHiggs_swap=
     &  +ggH_bquark_swap(h1,h2,h34,h56)
     &  +ggH_tquark_swap(h1,h2,h34,h56)
c--- compute total c6 correction to Higgs amplitude
        AHiggs_c6_swap=
     &  +Mloop_c6_propagator_swap(h1,h2,h34,h56)
     &  +Mloop_c6_decay_swap(h1,h2,h34,h56)
     &  +Mloop_c6_production_swap(h1,h2,h34,h56)
     &  +Mloop_c6_width_swap(h1,h2,h34,h56)  
c--- compute total ct Higgs amplitude
        AHiggs_ct_swap=ct*ggH_tquark_swap(h1,h2,h34,h56)
c--- compute total cg Higgs amplitude
        AHiggs_cg_swap=cg*Mloop_SMEFT_swap(h1,h2,h34,h56)
c--- compute total SMEFT Higgs amplitude
        AHiggs_SMEFT_swap=AHiggs_c6_swap+AHiggs_ct_swap+AHiggs_cg_swap

!        Mamp=Acont+AHiggs+AHiggs_c6
!        Samp=Acont_swap+AHiggs_swap+AHiggs_c6_swap

        Mamp=Acont+AHiggs+AHiggs_SMEFT
        Samp=Acont_swap+AHiggs_swap+AHiggs_SMEFT_swap

        if (h34 == h56) then
          oprat=1._dp-two*real(conjg(Mamp)*Samp)
     &                 /(abs(Mamp)**2+abs(Samp)**2)
        else
          oprat=1._dp
        endif
        if (bw34_56) then
          msqgg=msqgg+two*abs(Mamp)**2*oprat
        else
          msqgg=msqgg+two*abs(Samp)**2*oprat
        endif
c--- subtract |Acont|^2
        Mamp=Acont
        Samp=Acont_swap
        if (h34 == h56) then
          oprat=1._dp-two*real(conjg(Mamp)*Samp)
     &                 /(abs(Mamp)**2+abs(Samp)**2)
        else
          oprat=1._dp
        endif
        if (bw34_56) then
          msqgg=msqgg-two*abs(Mamp)**2*oprat
        else
          msqgg=msqgg-two*abs(Samp)**2*oprat
        endif
c--- subtract |AHiggs|^2
        if (h1 == h2) then
c------ Higgs amplitudes only non-zero for h1=h2
        Mamp=AHiggs+AHiggs_SMEFT
        Samp=AHiggs_swap+AHiggs_SMEFT_swap
        if (h34 == h56) then
          oprat=1._dp-two*real(conjg(Mamp)*Samp)
     &                 /(abs(Mamp)**2+abs(Samp)**2)
        else
          oprat=1._dp
        endif
        if (bw34_56) then
          msqgg=msqgg-two*abs(Mamp)**2*oprat
        else
          msqgg=msqgg-two*abs(Samp)**2*oprat
        endif
        endif
      endif

      enddo
      enddo
      enddo
      enddo

c--- overall factor extracted (c.f. getggZZamps.f and getggHZZamps.f )
      fac=avegg*V*(four*abs(zesq)*gsq/(16._dp*pisq)*abs(zesq))**2

      msq(0,0)=msqgg*fac*vsymfact

      return
      end


