!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

c--- Author: J. M. Campbell, September 2013
c--- Matrix element squared for gg -> H -> ZZ signal process
c--- The exact result for massive bottom and top quark loops is included
      subroutine gg_hZZ_tb(p,msq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'interference.f'
      include 'zcouple_cms.f'
      integer:: h1,h2,h34,h56
      real(dp):: p(mxpart,4),wt,wt2,msq(fn:nf,fn:nf),msqgg,fac,
     & pswap(mxpart,4),oprat
      complex(dp)::
     & Mloop_uptype(2,2,2,2),Mloop_dntype(2,2,2,2),
     & Mloop_bquark(2,2,2,2),Mloop_tquark(2,2,2,2),
     & Sloop_uptype(2,2,2,2),Sloop_dntype(2,2,2,2),
     & Sloop_bquark(2,2,2,2),Sloop_tquark(2,2,2,2),
     & ggH_bquark(2,2,2,2),ggH_tquark(2,2,2,2),Acont,AHiggs,
     & ggH_bquark_swap(2,2,2,2),ggH_tquark_swap(2,2,2,2),AHiggs_swap,
     & Acont_swap,Mamp,Samp,
     & Mloop_c6_propagator(2,2,2,2),
     & Mloop_c6_propagator_swap(2,2,2,2),
     & AHiggs_c6,AHiggs_c6_swap,
     & Mloop_c6_decay(2,2,2,2),
     & Mloop_c6_decay_swap(2,2,2,2),
     & Mloop_c6_production(2,2,2,2),
     & Mloop_c6_production_swap(2,2,2,2),
     & Mloop_c6_width(2,2,2,2),
     & Mloop_c6_width_swap(2,2,2,2)

      msq(:,:)=0._dp
      call getggHZZamps(p,ggH_bquark,ggH_tquark,Mloop_c6_propagator,Mloop_c6_decay,Mloop_c6_production,Mloop_c6_width)

      if (interference) then
c--- for interference, compute amplitudes after 4<->6 swap
       pswap(1,:)=p(1,:)
       pswap(2,:)=p(2,:)
       pswap(3,:)=p(3,:)
       pswap(4,:)=p(6,:)
       pswap(5,:)=p(5,:)
       pswap(6,:)=p(4,:)
       call getggHZZamps(pswap,ggH_bquark_swap,ggH_tquark_swap,Mloop_c6_propagator_swap,Mloop_c6_decay_swap,Mloop_c6_production_swap,Mloop_c6_width_swap)
      endif

      msqgg=0._dp
      do h1=1,2
      h2=h1
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

c-------------------------------
      if (interference .eqv. .false.) then
c--- normal case
      msqgg=msqgg+abs(AHiggs+AHiggs_c6)**2
      else
c--- with interference
        AHiggs_swap=
     &  +ggH_bquark_swap(h1,h2,h34,h56)
     &  +ggH_tquark_swap(h1,h2,h34,h56)

        AHiggs_c6_swap=
     &  +Mloop_c6_propagator_swap(h1,h2,h34,h56)
     &  +Mloop_c6_decay_swap(h1,h2,h34,h56)
     &  +Mloop_c6_production_swap(h1,h2,h34,h56)
     &  +Mloop_c6_width_swap(h1,h2,h34,h56) 
     

        if (h34 == h56) then   
          oprat=1._dp-two*real(conjg(AHiggs+AHiggs_c6)*(AHiggs_swap+AHiggs_c6_swap))
     &  /(abs(AHiggs+AHiggs_c6)**2+abs(AHiggs_swap+AHiggs_c6_swap)**2)
        else
          oprat=1._dp
        endif
        if (bw34_56) then
          msqgg=msqgg+two*abs(AHiggs+AHiggs_c6)**2*oprat
        else
          msqgg=msqgg+two*abs(AHiggs_swap+AHiggs_c6_swap)**2*oprat
        endif
      endif
      enddo
      enddo
      enddo

c--- overall factor extracted (c.f. getggHZZamps.f)
      fac=avegg*V*(four*abs(zesq)*gsq/(16._dp*pisq)*abs(zesq))**2

      msq(0,0)=msqgg*fac*vsymfact

      ! print*,"|||||||||||||||| phase space point: ||||||||||||||||||"
      ! print*,p(1,:)
      ! print*,p(2,:)
      ! print*,p(3,:)
      ! print*,p(4,:)
      ! print*,p(5,:)
      ! print*,p(6,:)
      ! print*,""
      ! print*,"msq(c6) gg_hzz_tb: ", msq(0,0)
      ! call writecsv(p,msq,wt)

      return
      end


