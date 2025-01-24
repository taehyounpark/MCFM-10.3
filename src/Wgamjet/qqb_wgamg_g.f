!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_wgamg_g(p,msq)
      implicit none
c-----Author:Keith Ellis, May 2020
c---- Matrix element for W gam production
c----averaged over initial colours and spins
c For nwz=+1
c     u(-p1)+dbar(-p2)-->W^+(nu(p3)+e^+(p4)) + gamma(p5) + g(p6)
c For nwz=-1
c     ubar(-p1)+d(-p2)-->W^-(e^-(p3)+nbar(p4)) + gamma(p5) + g(p6)
c---
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'ckm.f'
      include 'zprods_com.f'
      include 'zcouple_cms.f'
      include 'nwz.f'
      include 'qqxqq.f'
      include 'nflav.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac
      real(dp):: qbq,qqb,qg,gq,qbg,gqb,gg,ubdggmsq,qqxqqsq(12)

      call spinoru(7,p,za,zb)

      msq(:,:)=0._dp

c DEBUG: 4-quark diagrams only
c      goto 33

      if (nwz == -1) then
      qqb=half*aveqq*xn*ubdggmsq(1,2,3,4,5,6,7,za,zb)
      qbq=half*aveqq*xn*ubdggmsq(2,1,3,4,5,6,7,za,zb)
      qg =aveqg*xn*ubdggmsq(1,6,3,4,5,2,7,za,zb)
      gqb=aveqg*xn*ubdggmsq(6,2,3,4,5,1,7,za,zb)
      gq =aveqg*xn*ubdggmsq(2,6,3,4,5,1,7,za,zb)
      qbg=aveqg*xn*ubdggmsq(6,1,3,4,5,2,7,za,zb)
      gg =avegg*xn*ubdggmsq(7,6,3,4,5,1,2,za,zb)
      elseif (nwz == +1) then
      qqb=half*aveqq*xn*ubdggmsq(2,1,4,3,5,6,7,zb,za)
      qbq=half*aveqq*xn*ubdggmsq(1,2,4,3,5,6,7,zb,za)
      qg =aveqg*xn*ubdggmsq(6,1,4,3,5,2,7,zb,za)
      gqb=aveqg*xn*ubdggmsq(2,6,4,3,5,1,7,zb,za)
      gq =aveqg*xn*ubdggmsq(6,2,4,3,5,1,7,zb,za)
      qbg=aveqg*xn*ubdggmsq(1,6,4,3,5,2,7,zb,za)
      gg =avegg*xn*ubdggmsq(6,7,4,3,5,1,2,zb,za)
      endif


      do j=-4,4
      do k=-4,4
          if ((j > 0) .and. (k < 0)) then
            msq(j,k)=Vsq(j,k)*qqb
          elseif ((j == 0) .and. (k < 0)) then
            msq(j,k)=Vsum(k)*gqb
          elseif ((j == 0) .and. (k > 0)) then
            msq(j,k)=Vsum(k)*gq
          elseif ((j < 0) .and. (k > 0)) then
            msq(j,k)=Vsq(j,k)*qbq
          elseif ((j > 0) .and. (k == 0)) then
            msq(j,k)=Vsum(j)*qg
         elseif ((j < 0) .and. (k == 0)) then
            msq(j,k)=Vsum(j)*qbg
         elseif ((j == 0) .and. (k == 0)) then
            msq(j,k)=2._dp*gg
          endif
      enddo
      enddo

c      goto 99

c   33 continue

c---*****************4 quark diagrams********************

c     dc_uc1=1,ds_us1=2,ds_dc1=3,us_uc1=4,
c     dc_uc2=5,ds_us2=6,ds_dc2=7,us_uc2=8,
c     du_uu= 9,dd_ud=10,dd_du=11,ud_uu=12

c     q q --> q q
      if (nwz == -1) then
      call a7Wgamfourqsq(1,2,3,4,5,6,7,za,zb,qqxqqsq)
      msq(1,1)=qqxqqsq(dd_ud)
      msq(1,2)=half*qqxqqsq(du_uu)
      msq(1,3)=qqxqqsq(ds_us1)+qqxqqsq(ds_dc1)
      msq(1,4)=qqxqqsq(dc_uc1)
      msq(2,1)=half*qqxqqsq(ud_uu)
      msq(2,3)=qqxqqsq(us_uc1)

      msq(3,1)=msq(1,3)
      msq(3,2)=msq(1,4)
      msq(3,3)=msq(1,1)
      msq(3,4)=msq(1,2)
      msq(4,1)=msq(2,3)
      msq(4,3)=msq(2,1)

      if (nflav == 5) then
      msq(5,1)=qqxqqsq(ds_dc1)
      msq(5,3)=msq(5,1)
      msq(1,5)=qqxqqsq(ds_us1)
      msq(3,5)=msq(1,5)
      endif

      elseif (nwz == +1) then
      call a7Wgamfourqsq(6,7,4,3,5,1,2,zb,za,qqxqqsq)
      msq(2,2)=qqxqqsq(dd_ud)
      msq(2,1)=half*qqxqqsq(du_uu)
      msq(2,4)=qqxqqsq(ds_us1)+qqxqqsq(ds_dc1)
      msq(2,3)=qqxqqsq(dc_uc1)
      msq(1,2)=half*qqxqqsq(ud_uu)
      msq(1,4)=qqxqqsq(us_uc1)

      msq(4,2)=msq(2,4)
      msq(4,1)=msq(2,3)
      msq(4,4)=msq(2,2)
      msq(4,3)=msq(2,1)
      msq(3,2)=msq(1,4)
      msq(3,4)=msq(1,2)

      if (nflav==5) then
      msq(2,5)=msq(2,3)
      msq(4,5)=msq(2,3)
      msq(5,2)=msq(3,2)
      msq(5,4)=msq(3,2)
      endif
      endif


c     qb q --> qb q
      if (nwz == -1) then
      call a7Wgamfourqsq(6,2,3,4,5,1,7,za,zb,qqxqqsq)
      msq(-1,1)=qqxqqsq(dd_du)
      msq(-1,3)=qqxqqsq(ds_dc1)
      msq(-2,4)=qqxqqsq(dc_uc1)
      msq(-2,1)=msq(-2,1)
     & +qqxqqsq(ud_uu)+qqxqqsq(dd_ud)+qqxqqsq(ds_us2)+qqxqqsq(us_uc2)
      msq(-2,2)=qqxqqsq(du_uu)
      msq(-2,3)=qqxqqsq(us_uc1)+qqxqqsq(ds_us1)

      msq(-4,1)=msq(-2,3)
      msq(-4,2)=msq(-2,4)
      msq(-4,3)=msq(-2,1)
      msq(-4,4)=msq(-2,2)
      msq(-3,1)=msq(-1,3)
      msq(-3,3)=msq(-1,1)

      if (nflav == 5) then
      msq(-5,1)=msq(-3,1)
      msq(-5,3)=msq(-3,1)
      msq(-2,5)=qqxqqsq(ds_us1)
      msq(-4,5)=msq(-2,5)
      msq(-2,1)=msq(-2,1)+qqxqqsq(ds_us2)
      msq(-4,3)=msq(-2,1)
      endif

      elseif (nwz == +1) then
      call a7Wgamfourqsq(1,7,4,3,5,6,2,zb,za,qqxqqsq)
      msq(-2,2)=qqxqqsq(dd_du)
      msq(-2,4)=qqxqqsq(ds_dc1)
      msq(-1,3)=qqxqqsq(dc_uc1)
      msq(-1,2)=msq(-1,2)
     & +qqxqqsq(ud_uu)+qqxqqsq(dd_ud)+qqxqqsq(ds_us2)+qqxqqsq(us_uc2)
      msq(-1,1)=qqxqqsq(du_uu)
      msq(-1,4)=qqxqqsq(us_uc1)+qqxqqsq(ds_us1)

      msq(-3,2)=msq(-1,4)
      msq(-3,1)=msq(-1,3)
      msq(-3,4)=msq(-1,2)
      msq(-3,3)=msq(-1,1)
      msq(-4,2)=msq(-2,4)
      msq(-4,4)=msq(-2,2)

      if (nflav == 5) then
      msq(-1,5)=msq(-1,3)
      msq(-3,5)=msq(-1,3)
      msq(-5,2)=qqxqqsq(us_uc1)
      msq(-5,4)=msq(-5,2)
      msq(-1,2)=msq(-1,2)+qqxqqsq(us_uc2)
      msq(-3,4)=msq(-1,2)
      endif
      endif

c     q qb --> q qb
      if (nwz == -1) then
      call a7Wgamfourqsq(1,7,3,4,5,6,2,za,zb,qqxqqsq)
      msq(1,-4)=qqxqqsq(dc_uc1)+qqxqqsq(ds_dc1)
      msq(1,-3)=qqxqqsq(ds_us1)
      msq(1,-2)=msq(1,-2)
     & +qqxqqsq(du_uu)+qqxqqsq(dd_du)+qqxqqsq(dc_uc2)+qqxqqsq(ds_dc2)
      msq(1,-1)=qqxqqsq(dd_ud)
      msq(2,-2)=qqxqqsq(ud_uu)
      msq(2,-4)=qqxqqsq(us_uc1)

      msq(3,-4)=msq(1,-2)
      msq(3,-3)=msq(1,-1)
      msq(3,-2)=msq(1,-4)
      msq(3,-1)=msq(1,-3)
      msq(4,-2)=msq(2,-4)
      msq(4,-4)=msq(2,-2)

      if (nflav == 5) then
      msq(5,-2)=qqxqqsq(ds_dc1)
      msq(5,-4)=msq(5,-2)
      msq(1,-5)=msq(1,-3)
      msq(3,-5)=msq(1,-3)
      msq(1,-2)=msq(1,-2)+qqxqqsq(ds_dc2)
      msq(3,-4)=msq(1,-2)
      endif

      elseif (nwz == +1) then
      call a7Wgamfourqsq(6,2,4,3,5,1,7,zb,za,qqxqqsq)
      msq(2,-3)=qqxqqsq(dc_uc1)+qqxqqsq(ds_dc1)
      msq(2,-4)=qqxqqsq(ds_us1)
      msq(2,-1)=msq(2,-1)
     & +qqxqqsq(du_uu)+qqxqqsq(dd_du)+qqxqqsq(dc_uc2)+qqxqqsq(ds_dc2)
      msq(2,-2)=qqxqqsq(dd_ud)
      msq(1,-1)=qqxqqsq(ud_uu)
      msq(1,-3)=qqxqqsq(us_uc1)

      msq(4,-3)=msq(2,-1)
      msq(4,-4)=msq(2,-2)
      msq(4,-1)=msq(2,-3)
      msq(4,-2)=msq(2,-4)
      msq(3,-1)=msq(1,-3)
      msq(3,-3)=msq(1,-1)

      if (nflav == 5) then
      msq(5,-1)=msq(3,-1)
      msq(5,-3)=msq(3,-1)
      msq(2,-5)=qqxqqsq(dc_uc1)
      msq(4,-5)=msq(2,-5)
      msq(2,-1)=msq(2,-1)+qqxqqsq(dc_uc2)
      msq(4,-3)=msq(2,-1)
      endif
      endif

c     qb qb --> qb qb
      if (nwz == -1) then
      call a7Wgamfourqsq(6,7,3,4,5,1,2,za,zb,qqxqqsq)
      msq(-4,-3)=qqxqqsq(ds_dc1)
      msq(-2,-4)=qqxqqsq(dc_uc1)+qqxqqsq(us_uc1)
      msq(-2,-3)=qqxqqsq(ds_us1)
      msq(-2,-2)=qqxqqsq(du_uu)
      msq(-2,-1)=half*qqxqqsq(dd_ud)
      msq(-1,-4)=qqxqqsq(ds_dc1)
      msq(-1,-2)=half*qqxqqsq(dd_du)

      msq(-4,-4)=msq(-2,-2)
      msq(-4,-3)=msq(-2,-1)
      msq(-4,-2)=msq(-2,-4)
      msq(-4,-1)=msq(-2,-3)
      msq(-3,-4)=msq(-1,-2)
      msq(-3,-2)=msq(-1,-4)

      if (nflav == 5) then
      msq(-5,-2)=msq(-3,-2)
      msq(-5,-4)=msq(-3,-2)
      msq(-2,-5)=msq(-2,-3)
      msq(-4,-5)=msq(-2,-3)
      endif

      elseif (nwz == +1) then
      call a7Wgamfourqsq(1,2,4,3,5,6,7,zb,za,qqxqqsq)
      msq(-3,-4)=qqxqqsq(ds_dc1)
      msq(-1,-3)=qqxqqsq(dc_uc1)+qqxqqsq(us_uc1)
      msq(-1,-4)=qqxqqsq(ds_us1)
      msq(-1,-1)=qqxqqsq(du_uu)
      msq(-1,-2)=half*qqxqqsq(dd_ud)
      msq(-2,-3)=qqxqqsq(ds_dc1)
      msq(-2,-1)=half*qqxqqsq(dd_du)

      msq(-3,-3)=msq(-1,-1)
      msq(-3,-4)=msq(-1,-2)
      msq(-3,-1)=msq(-1,-3)
      msq(-3,-2)=msq(-1,-4)
      msq(-4,-3)=msq(-2,-1)
      msq(-4,-1)=msq(-2,-3)

      if (nflav == 5) then
      msq(-5,-1)=qqxqqsq(us_uc1)
      msq(-5,-3)=msq(-5,-1)
      msq(-1,-5)=qqxqqsq(dc_uc1)
      msq(-3,-5)=msq(-1,-5)
      endif
      endif

c   99 continue

c     Normalize both qqgg and qqqq matrix element
      fac=V*gsq**2*abs((zesq/zxw)**2*zesq)/2._dp
      msq(:,:)=fac*msq(:,:)
      return
      end

