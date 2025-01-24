!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module singletop_jetdeps
      use types
      use singletop2_nnlo_vars
      use singletop2_scale_m
      implicit none

      private

      public :: qb_wtq
      public :: qb_wtq_intf

      public :: wqq_sc_intf

      public :: interdk_gen
      public :: interdk_qq_gen
      public :: interdk_qqid_gen

      public :: inter_gen
      public :: inter_qq_gen
      public :: inter_qqid_gen

      contains

      subroutine wqq_sc_intf(is,ic,ie,in,iq,iqb,p,msq)
          use types
          use singletop2_nnlo_vars
          use singletop2_scale_m
      implicit none

c Updated of routine wqq_sc to account for identical-particle interference
c between is and iqb (J. Campbell, November 2019)

c Author: J. Campbell, April 2004.

c     s(-p1)+cbar(-p2) --> l(p3)+abar(p4)+s(p5)+sb(p6)
c     with c massive
c     multiplied by (((a+l)^2-M**2)/(a+l)^2)^2*g^4/gwsq^2/2
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'masses.f'
      include 'zprods_decl.f'
      real(dp):: msq,p(mxpart,4),q(mxpart,4),mch
      real(dp):: s134,s156,se,sn,en,sq,sqb,qqb,cq,afac,prop,facqq
      complex(dp):: ampmm(2),ampmp(2),amppm(2),amppp(2)
      integer:: is,ic,ie,in,iq,iqb,nu,j,js,jqb

      real(dp) :: gsq

      en=p(ie,4)*p(in,4)-p(ie,1)*p(in,1)-p(ie,2)*p(in,2)-p(ie,3)*p(in,3)

      cq=p(iq,4)*p(ic,4)-p(iq,1)*p(ic,1)-p(iq,2)*p(ic,2)-p(iq,3)*p(ic,3)
      mch=sqrt(abs(p(ic,4)**2-p(ic,1)**2-p(ic,2)**2-p(ic,3)**2))

c----define massless momenta
      afac=0.5_dp*mch**2/cq

      do nu=1,4
      do j=1,6
        if (j==ic) then
          q(j,nu)=p(ic,nu)-p(iq,nu)*afac
        else
          q(j,nu)=p(j,nu)
        endif
      enddo
      enddo

      call spinoru(6,q,za,zb)

c loop over identical-particle exchange
      do j=1,2
      if (j == 1) then
        js=is
        jqb=iqb
      else
        js=iqb
        jqb=is
      endif

      se=p(js,4)*p(ie,4)-p(js,1)*p(ie,1)-p(js,2)*p(ie,2)-p(js,3)*p(ie,3)
      sn=p(js,4)*p(in,4)-p(js,1)*p(in,1)-p(js,2)*p(in,2)-p(js,3)*p(in,3)
      sq=p(js,4)*p(iq,4)-p(js,1)*p(iq,1)-p(js,2)*p(iq,2)-p(js,3)*p(iq,3)
      sqb=p(js,4)*p(jqb,4)-p(js,1)*p(jqb,1)
     &   -p(js,2)*p(jqb,2)-p(js,3)*p(jqb,3)
      qqb=p(iq,4)*p(jqb,4)-p(iq,1)*p(jqb,1)
     &   -p(iq,2)*p(jqb,2)-p(iq,3)*p(jqb,3)

      s134=2._dp*(se+sn+en)-mch**2
      s156=2._dp*(sq+sqb+qqb)

      ampmm(j) =
     & 2*za(ie,ic)*za(iq,ic)*zb(js,in)*zb(jqb,ic)*s134**(-1) - 2*za(ie,
     & ic)*za(iq,jqb)*zb(jqb,js)*zb(jqb,in)*s156**(-1) - 2*za(ie,ic)*
     & za(iq,js)*zb(js,in)*zb(jqb,js)*s156**(-1) + 2*za(iq,ie)*za(iq,ic
     & )*zb(js,in)*zb(iq,jqb)*s134**(-1)

      amppm(j) =
     & 2*za(iq,ie)*za(iq,jqb)*zb(jqb,js)*zb(jqb,in)/za(iq,ic)*
     & s156**(-1)*mch - 2*za(iq,ie)*zb(js,in)*zb(jqb,ic)*s134**(-1)*mch
     &  + 2*za(iq,js)*za(iq,ie)*zb(js,in)*zb(jqb,js)/za(iq,ic)*
     & s156**(-1)*mch

      ampmp(j) =
     & 2*za(ie,ic)*za(iq,jqb)*zb(iq,js)*zb(iq,in)*s156**(-1) + 2*za(ie,
     & ic)*za(jqb,ic)*zb(js,in)*zb(iq,ic)*s134**(-1) - 2*za(ie,ic)*za(
     & jqb,js)*zb(js,in)*zb(iq,js)*s156**(-1) - 2*za(jqb,ie)*za(jqb,ic)
     & *zb(js,in)*zb(iq,jqb)*s134**(-1)

      amppp(j) =
     &  - 2*za(ie,ic)*za(iq,jqb)*zb(js,in)*zb(iq,ic)/za(iq,ic)*
     & s134**(-1)*mch - 2*za(iq,ie)*za(iq,jqb)*zb(iq,js)*zb(iq,in)/za(
     & iq,ic)*s156**(-1)*mch + 2*za(iq,ie)*za(jqb,js)*zb(js,in)*zb(iq,js
     & )/za(iq,ic)*s156**(-1)*mch + 2*za(iq,jqb)*za(jqb,ie)*zb(js,in)*
     & zb(iq,jqb)/za(iq,ic)*s134**(-1)*mch - 2*za(jqb,ie)*zb(js,in)*zb(
     & iq,ic)*s134**(-1)*mch

c-- insert the gluon propagator
      ampmm(j)=ampmm(j)/2._dp/qqb
      ampmp(j)=ampmp(j)/2._dp/qqb
      amppm(j)=amppm(j)/2._dp/qqb
      amppp(j)=amppp(j)/2._dp/qqb

      enddo

      if (corr_on_beam == 1) then
          gsq = 4*pi*as_heavy_beam1
      else
          gsq = 4*pi*as_heavy_beam2
      endif

c--- overall factor from sc_Wqq.frm
      prop=1._dp/((two*en-wmass**2)**2+wmass**2*wwidth**2)
      facqq=aveqq*xn*cf*(gsq*gwsq)**2/2._dp*prop

c Note: no identical-particle factor of 1/2 applied yet
      msq=facqq*real(ampmm(1)*conjg(ampmm(1))+ampmp(1)*conjg(ampmp(1))
     &              +amppm(1)*conjg(amppm(1))+amppp(1)*conjg(amppp(1))
     &              +ampmm(2)*conjg(ampmm(2))+ampmp(2)*conjg(ampmp(2))
     &              +amppm(2)*conjg(amppm(2))+amppp(2)*conjg(amppp(2))
     &      +two/xn*(ampmm(1)*conjg(ampmm(2))+amppm(1)*conjg(amppm(2))))

      return
      end

      subroutine interdk_gen(p,u,g1,i3,i4,i5,b,d,g2,mq,ampsq)
c--- Matrix-element squared for the process:
c---
c--- 0  ->  u~ + g1 + d + t(i3,i4,i5) + b~ + g2
c---
c---        u -------------- d
c---                 $
c---                 $
c---        b -------------- t
c---             \       \
c---              g1      g2
c---
c--- where gluon radiation is *only* attached to the (t,b) line
c--- and without any initial state averaging over helicites and colors
c---
c--- Based on earlier routines for, alternatively,
c---  (1) t-channel single top production with a massive b-quark
c---  (2) W+t associated single top production (more efficient)
c---
      use types
      use singletop2_nnlo_vars
      use singletop2_scale_m
      implicit none
      include 'constants.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'mxpart.f'
      include 'nf.f'
      integer:: u,g1,i3,i4,i5,b,d,g2,ia,ig
      real(dp) :: p(mxpart,4),mq,ampsq,ampsq_tb,ampsq_Wt,sud,s34,dot,prop,fac
      complex(dp) :: amp_ga(2,2),amp_ag(2,2)
      integer, parameter :: iflag = 2 ! Choice of routine above, or 3 for comparison

      real(dp) :: gsq

      if (corr_on_beam == 1) then
          gsq = 4*pi*as_heavy_beam1
      else
          gsq = 4*pi*as_heavy_beam2
      endif

      if ((iflag ==1) .or. (iflag == 3)) then
        ! please reinstate only for debugging
        !as_L=zip
        error stop __FILE__//" interdk_gen bad iflag: not to be used in production, interdk depends on uncontrolled gsq"
        call interdk(p,u,g1,i3,i4,b,d,g2,mq,zip,ampsq_tb)
c-- only have to remove initial-state averaging
        ampsq_tb=ampsq_tb/aveqg
      endif
      if ((iflag ==2) .or. (iflag == 3)) then
        call gs_wt_prog(mt,zip,p,g1,b,d,u,i3,i4,i5,g2,amp_ga)
        call gs_wt_prog(mt,zip,p,g2,b,d,u,i3,i4,i5,g1,amp_ag)
        ampsq_Wt=zip
c-- must square and sum, then apply overall factors
        do ia=1,2
        do ig=1,2
          ampsq_Wt=ampsq_Wt+xn*cf**2*(
     &          + abs(amp_ag(ia,ig))**2 + abs(amp_ga(ig,ia))**2
     &          - one/xn/cf*real(amp_ag(ia,ig)*conjg(amp_ga(ig,ia))))
        enddo
        enddo
        fac=gsq**2*gwsq**4*xn
        sud=two*dot(p,u,d)
        if (sud < zip) then
          prop=(sud-wmass**2)**2
        else
          prop=(sud-wmass**2)**2+(wmass*wwidth)**2
        endif
        s34=two*dot(p,i3,i4)
        if (s34 < zip) then
          prop=prop*(s34-wmass**2)**2
        else
          prop=prop*((s34-wmass**2)**2+(wmass*wwidth)**2)
        endif
        prop=prop*(mt*twidth)**2
        ampsq_Wt=fac/prop*ampsq_Wt
      endif

      if     (iflag == 1) then
        ampsq=ampsq_tb
      elseif (iflag == 2) then
        ampsq=ampsq_Wt
      elseif (iflag == 3) then
        if (i5 /= 5) then
          write(6,*) 'Comparison not possible for i5 /= 5'
          stop
        endif
        write(6,*) 'interdk_gen comparison for: ',u,g1,i3,i4,i5,b,d,g2
        write(6,*) 'ampsq_tb,ampsq_Wt,ampsq_tb/ampsq_Wt',ampsq_tb,ampsq_Wt,ampsq_tb/ampsq_Wt
      endif

      return
      end

      subroutine inter_qqid_gen(p,u,q1,t,b,d,q2,mq,ampsq)
c--- Matrix-element squared for the process:
c---
c--- 0  ->  u~ + b~ + d + t + b~ + b
c---
c---        u -------------- d
c---                 $
c---                 $
c---        b -------------- t
c---              \
c---               \---- b
c---                \
c---                 \
c---                  b~
c---
c--- where gluon radiation is *only* attached to the (t,b) line
c--- and without any initial state averaging over helicites and colors
c---
      use types
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'mxpart.f'
      include 'nf.f'
      integer:: u,q1,q2,d,b,t
      real(dp) :: p(mxpart,4),mq,ampsq,sud,dot,prop

      call wqq_sc_intf(b,t,d,u,q2,q1,p,ampsq)
c-- must apply overall factors and also remove initial-state averaging
      sud=two*dot(p,u,d)
      prop=((sud-wmass**2)**2+(wmass*wwidth)**2)/(sud-wmass**2)**2
      ampsq=xn*prop*ampsq/aveqq

      return
      end

      subroutine inter_qq_gen(p,u,q1,t,b,d,q2,mq,ampsq)
c--- Matrix-element squared for the process:
c---
c--- 0  ->  u~ + q1~ + d + t + b~ + q2
c---
c---        u -------------- d
c---                 $
c---                 $
c---        b -------------- t
c---              \
c---               \---- q2
c---                \
c---                 \
c---                  q1~
c---
c--- where gluon radiation is *only* attached to the (t,b) line
c--- and without any initial state averaging over helicites and colors
c---
c--- Based on earlier routines for, alternatively,
c---  (1) t-channel single top production with a massive b-quark
c---  (2) W+t associated single top production (more efficient)
c---
      use types
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'mxpart.f'
      include 'nf.f'
      integer:: u,q1,q2,d,b,t
      real(dp) :: p(mxpart,4),mq,ampsq,ampsq_tb,ampsq_Wt,sud,dot,prop,
     & ampsqarray(4),intfarray(4)
      integer, parameter :: iflag = 2 ! Choice of routine above, or 3 for comparison

      if ((iflag ==1) .or. (iflag == 3)) then
        call inter_qq(p,u,q1,t,b,d,q2,mq,zip,ampsqarray,intfarray)
c-- only have to remove initial-state averaging
        ampsq_tb=ampsqarray(1)/aveqq
      endif
      if ((iflag ==2) .or. (iflag == 3)) then
        call wqq_sc(b,t,d,u,q2,q1,p,ampsq_Wt)
c-- must apply overall factors and also remove initial-state averaging
        sud=two*dot(p,u,d)
        prop=((sud-wmass**2)**2+(wmass*wwidth)**2)/(sud-wmass**2)**2
        ampsq_Wt=xn*prop*ampsq_Wt/aveqq
      endif

      if     (iflag == 1) then
        ampsq=ampsq_tb
      elseif (iflag == 2) then
        ampsq=ampsq_Wt
      elseif (iflag == 3) then
        write(6,*) 'inter_qq_gen comparison for: ',u,q1,t,b,d,q2
        write(6,*) 'ampsq_tb,ampsq_Wt,ampsq_tb/ampsq_Wt',ampsq_tb,ampsq_Wt,ampsq_tb/ampsq_Wt
      endif

      return
      end

      subroutine inter_gen(p,u,g1,t,b,d,g2,mq,ampsq)
c--- Matrix-element squared for the process:
c---
c--- 0  ->  u~ + g1 + d + t + b~ + g2
c---
c---        u -------------- d
c---                 $
c---                 $
c---        b -------------- t
c---             \       \
c---              g1      g2
c---
c--- where gluon radiation is *only* attached to the (t,b) line
c--- and without any initial state averaging over helicites and colors
c---
c--- Based on earlier routines for, alternatively,
c---  (1) t-channel single top production with a massive b-quark
c---  (2) W+t associated single top production (more efficient)
c---
      use types
      use singletop2_nnlo_vars
      use singletop2_scale_m
      implicit none
      include 'constants.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'mxpart.f'
      include 'nf.f'
      integer:: u,g1,g2,d,b,t
      real(dp) :: p(mxpart,4),mq,ampsq,ampsq_tb,ampsq_Wt,sud,dot,prop,fac
      integer, parameter :: iflag = 2 ! Choice of routine above, or 3 for comparison
      real(dp) :: gsq

      if (corr_on_beam == 1) then
          gsq = 4*pi*as_heavy_beam1
      else
          gsq = 4*pi*as_heavy_beam2
      endif

      if ((iflag ==1) .or. (iflag == 3)) then
        ! please reinstate only for debugging
        !as_L=zip
        call inter(p,u,g1,t,b,d,g2,mq,zip,ampsq_tb)
c-- only have to remove initial-state averaging
        ampsq_tb=ampsq_tb/aveqg
      endif
      if ((iflag ==2) .or. (iflag == 3)) then
        call w2jetsq_mass(b,t,d,u,g1,g2,p,ampsq_Wt)
c-- must apply overall factors
        sud=two*dot(p,u,d)
        prop=sud**2/(sud-wmass**2)**2
        fac=V*xn/four*(gwsq/two*gsq)**2*prop*xn
        ampsq_Wt=fac*ampsq_Wt
      endif

      if     (iflag == 1) then
        ampsq=ampsq_tb
      elseif (iflag == 2) then
        ampsq=ampsq_Wt
      elseif (iflag == 3) then
        write(6,*) 'inter_gen comparison for: ',u,g1,t,b,d,g2
        write(6,*) 'ampsq_tb,ampsq_Wt,ampsq_tb/ampsq_Wt',ampsq_tb,ampsq_Wt,ampsq_tb/ampsq_Wt
      endif

      return
      end

      subroutine interdk_qq_gen(p,u,q1,i3,i4,i5,b,d,q2,mq,ampsq)
c--- Matrix-element squared for the process:
c---
c--- 0  ->  u~ + q1~ + d + t(i3,i4,i5) + b~ + q2
c---
c---        u -------------- d
c---                 $
c---                 $
c---        b -------------- t
c---              \
c---               \---- q2
c---                \
c---                 \
c---                  q1~
c---
c--- where gluon radiation is *only* attached to the (t,b) line
c--- and without any initial state averaging over helicites and colors
c---
c--- Based on earlier routines for, alternatively,
c---  (1) t-channel single top production with a massive b-quark
c---  (2) W+t associated single top production (more efficient)
c---
      use types
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'mxpart.f'
      include 'nf.f'
      integer:: u,q1,i3,i4,i5,b,d,q2
      real(dp) :: p(mxpart,4),mq,ampsq,ampsq_tb,ampsq_Wt,sud,s34,dot,prop,
     & ampsqarray(4),intfarray(4)
      integer, parameter :: iflag = 2 ! Choice of routine above, or 3 for comparison

      if ((iflag ==1) .or. (iflag == 3)) then
        call interdk_qq(p,u,q1,i3,i4,b,d,q2,mq,zip,ampsqarray,intfarray)
c-- only have to remove initial-state averaging
        ampsq_tb=ampsqarray(1)/aveqq
      endif
      if ((iflag ==2) .or. (iflag == 3)) then
        call qb_wtq(mt,twidth,p,q1,b,d,u,i3,i4,i5,q2,ampsq_Wt)
c-- must apply overall factors and also remove initial-state averaging,
c-- also adjusting W propagators if no W width is necessary
        sud=two*dot(p,u,d)
        if (sud < zip) then
          prop=((sud-wmass**2)**2+(wmass*wwidth)**2)/(sud-wmass**2)**2
        else
          prop=1._dp
        endif
        s34=two*dot(p,i3,i4)
        if (s34 < zip) then
          prop=prop*((s34-wmass**2)**2+(wmass*wwidth)**2)/(s34-wmass**2)**2
        endif
        ampsq_Wt=xn*prop*ampsq_Wt/aveqq
      endif

      if     (iflag == 1) then
        ampsq=ampsq_tb
      elseif (iflag == 2) then
        ampsq=ampsq_Wt
      elseif (iflag == 3) then
        if (i5 /= 5) then
          write(6,*) 'Comparison not possible for i5 /= 5'
          stop
        endif
        write(6,*) 'interdk_qq_gen comparison for: ',u,q1,i3,i4,i5,b,d,q2
        write(6,*) 'ampsq_tb,ampsq_Wt,ampsq_tb/ampsq_Wt',ampsq_tb,ampsq_Wt,ampsq_tb/ampsq_Wt
      endif

      return
      end

      subroutine interdk_qqid_gen(p,u,q1,i3,i4,i5,b,d,q2,mq,ampsq)
c--- Matrix-element squared for the process:
c---
c--- 0  ->  u~ + b~ + d + t(i3,i4,i5) + b~ + b
c---
c---        u -------------- d
c---                 $
c---                 $
c---        b -------------- t
c---              \
c---               \---- q2
c---                \
c---                 \
c---                  q1~
c---
c--- where gluon radiation is *only* attached to the (t,b) line
c--- and without any initial state averaging over helicites and colors
c---
      use types
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'mxpart.f'
      include 'nf.f'
      integer:: u,q1,i3,i4,i5,b,d,q2
      real(dp) :: p(mxpart,4),mq,ampsq,sud,s34,dot,prop

      call qb_wtq_intf(mt,twidth,p,q1,b,d,u,i3,i4,i5,q2,ampsq)
c-- must apply overall factors and also remove initial-state averaging,
c-- also adjusting W propagators if no W width is necessary
      sud=two*dot(p,u,d)
      if (sud < zip) then
        prop=((sud-wmass**2)**2+(wmass*wwidth)**2)/(sud-wmass**2)**2
      else
        prop=1._dp
      endif
      s34=two*dot(p,i3,i4)
      if (s34 < zip) then
        prop=prop*((s34-wmass**2)**2+(wmass*wwidth)**2)/(s34-wmass**2)**2
      endif
      ampsq=xn*prop*ampsq/aveqq

      return
      end

      subroutine qb_wtq_intf(mq,qwidth,p,ix,is,ie,in,jn,je,jb,iy,msq)
      use types
      implicit none

c Updated of routine qb_wtq to account for identical-particle interference
c between is and ix (J. Campbell, November 2019)

c***********************************************************************
c     Author: Francesco Tramontano                                     *
c     February, 2005.                                                  *
c     Real correction to W+t, radiation in production                  *
c                                                                      *
c    Matrix element squared and averaged over initial colours and spins*
c--- W+t production (nwz=-1)                                           *
c     q(-ix) + b(-is) --> W + t(pneb) + q(iy)                          *
c                         |   |                                        *
c                         |   --> nu(jn) + e^+(jn) + b(jb)             *
c                         |                                            *
c                         --> e^-(ie) + nubar(in)                      *
c--- W+tbar production (nwz=+1)                                        *
c     q(-p1) + qbar(-p2) --> W + t(p567) + f(p8)                       *
c                            |   |                                     *
c                            |   --> e^-(p5) + nubar(p6) + b(p7)       *
c                            |                                         *
c                            --> nu(p3) + e^+(p4)                      *
c***********************************************************************
c---- helicities: 1=minus 2=plus

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'masses.f'
      integer:: ix,is,ie,in,jn,je,jb,iy
      real(dp):: p(mxpart,4),fac,prop,dot,msq,mq,qwidth
      complex(dp):: ampa(2),ampb(2)
      real(dp) :: gsq

      if (corr_on_beam == 1) then
          gsq = 4*pi*as_heavy_beam1
      else
          gsq = 4*pi*as_heavy_beam2
      endif


      fac=aveqq*gsq**2*gw**8*xn*cf*half
      prop=(two*dot(p,ie,in)-wmass**2)**2+(wmass*wwidth)**2
      prop=prop*((two*dot(p,jn,je)-wmass**2)**2+(wmass*wwidth)**2)
      prop=prop*((two*(dot(p,jn,je)+dot(p,jn,jb)+dot(p,je,jb))-mq**2)**2
     & +(mq*qwidth)**2)
c Note: no identical-particle factor of 1/2 applied yet
      call amps_4quark(mq,qwidth,p,ix,is,ie,in,jn,je,jb,iy,ampa)
      call amps_4quark(mq,qwidth,p,is,ix,ie,in,jn,je,jb,iy,ampb)
      msq=abs(ampa(1))**2+abs(ampa(2))**2
     &   +abs(ampb(1))**2+abs(ampb(2))**2
     &   +real(ampa(1)*conjg(ampb(1)))*two/xn
      msq=fac*msq/prop
      return
      end

      subroutine qb_wtq(mq,qwidth,p,ix,is,ie,in,jn,je,jb,iy,msq)
      implicit none
      include 'types.f'
c***********************************************************************
c     Author: Francesco Tramontano                                     *
c     February, 2005.                                                  *
c     Real correction to W+t, radiation in production                  *
c                                                                      *
c    Matrix element squared and averaged over initial colours and spins*
c--- W+t production (nwz=-1)                                           *
c     q(-ix) + b(-is) --> W + t(pneb) + q(iy)                          *
c                         |   |                                        *
c                         |   --> nu(jn) + e^+(jn) + b(jb)             *
c                         |                                            *
c                         --> e^-(ie) + nubar(in)                      *
c--- W+tbar production (nwz=+1)                                        *
c     q(-p1) + qbar(-p2) --> W + t(p567) + f(p8)                       *
c                            |   |                                     *
c                            |   --> e^-(p5) + nubar(p6) + b(p7)       *
c                            |                                         *
c                            --> nu(p3) + e^+(p4)                      *
c***********************************************************************
c---- helicities: 1=minus 2=plus

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'masses.f'
      integer:: ix,is,ie,in,jn,je,jb,iy
      real(dp):: p(mxpart,4),fac,prop,dot,msq,mq,qwidth
      complex(dp):: amp(2)
      real(dp) :: gsq

      if (corr_on_beam == 1) then
          gsq = 4._dp*pi*as_heavy_beam1
      else
          gsq = 4._dp*pi*as_heavy_beam2
      endif

      fac=aveqq*gsq**2*gw**8*xn*cf*half
      prop=(two*dot(p,ie,in)-wmass**2)**2+(wmass*wwidth)**2
      prop=prop*((two*dot(p,jn,je)-wmass**2)**2+(wmass*wwidth)**2)
      prop=prop*((two*(dot(p,jn,je)+dot(p,jn,jb)+dot(p,je,jb))-mq**2)**2
     & +(mq*qwidth)**2)
      call amps_4quark(mq,qwidth,p,ix,is,ie,in,jn,je,jb,iy,amp)
      msq=abs(amp(1))**2+abs(amp(2))**2
      msq=fac*msq/prop
      return
      end


      subroutine amps_4quark(mq,qwidth,p,ix,is,ie,in,jn,je,jb,iy,amp)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j,ix,is,ie,in,jn,je,jb,iy
      real(dp):: p(mxpart,4),dot,xDy,txy2,sxy2,t(4),q(mxpart,4),
     & mq,qwidth
      complex(dp):: amp(2),zab(mxpart,mxpart),zba(mxpart,mxpart)
      do j=1,4
      t(j)=p(jn,j)+p(je,j)+p(jb,j)
      q(1,j)=t(j)
      q(2,j)=p(ix,j)
      q(3,j)=p(iy,j)
      enddo
      xDy=dot(p,ix,iy)
      sxy2=two*(dot(p,is,ix)+dot(p,is,iy)+dot(p,ix,iy))
      txy2=two*(dot(q,1,2)+dot(q,1,3)+dot(q,2,3))
      call spinoru(8,p,za,zb)
      call spinork(8,p,zab,zba,t)
      amp(1)=  + zab(iy,je)*xDy**(-1) * ( za(jb,jn)*za(iy,ie
     &    )*zb(in,is)*zb(ix,iy) )
      amp(1) = amp(1) + zab(ie,je)*sxy2**(-1)*xDy**(-1)*txy2 * (
     &     - za(jb,jn
     &    )*za(ix,iy)*zb(in,ix)*zb(ix,is) - za(jb,jn)*za(iy,iy)*zb(in,
     &    iy)*zb(ix,is) - za(jb,jn)*za(is,iy)*zb(in,is)*zb(ix,is) )
      amp(1) = amp(1) + zab(ie,je)*zab(iy,ix)*xDy**(-1) * (
     &    za(jb,jn)*zb(in,is) )

      amp(2)=  + zab(ix,je)*xDy**(-1) * ( za(jb,jn)*za(ix,ie
     &    )*zb(in,is)*zb(iy,ix) )
      amp(2) = amp(2) + zab(ie,je)*sxy2**(-1)*xDy**(-1)*txy2 * (
     &     - za(jb,jn
     &    )*za(ix,ix)*zb(in,ix)*zb(iy,is) - za(jb,jn)*za(iy,ix)*zb(in,
     &    iy)*zb(iy,is) - za(jb,jn)*za(is,ix)*zb(in,is)*zb(iy,is) )
      amp(2) = amp(2) + zab(ie,je)*zab(ix,iy)*xDy**(-1) * (
     &    za(jb,jn)*zb(in,is) )

c--- Use the "overall scheme" to include the top width when necessary
c      if (txy2+mq**2 > zero) then
c        txy2=sqrt(txy2**2+0*(mq*qwidth)**2)
c      endif
      amp(1)=amp(1)/txy2
      amp(2)=amp(2)/txy2

      return
      end


      end module
