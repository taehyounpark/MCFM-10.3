!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine a6Wpgam_v(i1,i2,i3,i4,i5,i6,za,zb,amplo,amp)
c Matrix element squared for the procoess
c  q(i1) + qb(i2) -> l(p3)+lbar(p4)+ gamma(i5) + g(i6)
c with full top mass dependence retained
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'nf.f'
      include 'zprods_decl.f'
      include 'ewcharge.f'
      include 'epinv.f'
      include 'scheme.f'
      include 'b0.f'
      integer i1,i2,i3,i4,i5,i6
      complex(dp) :: amp(2,2),amplo(2,2),ampLC,ampSL,
     & ampdklc(2,2),ampdksl(2,2),ampQu,ampQd !      amp(h5,h6)

c DEBUG: just Qu for now
c      Q(1)=zip
c      Q(2)=1d0
c DEBUG: just Qd for now
c      Q(1)=1d0
c      Q(2)=zip

c DEBUG: to produce KCheck numbers
c      epinv=0d0
c Terms proportional to Qu or Qd

c (+, +)
c      call WpampLC_uppe_pp_full(i1,i2,i3,i4,i5,i6,za,zb,ampQux)
c      call WpampLC_down_pp_full(i1,i2,i3,i4,i5,i6,za,zb,ampQdx)
c      write(6,*) 'LC_uppe_pp',ampQu
c      write(6,*) 'LC_down_pp',ampQd

      call WpampLC_uppe_pp_paper(i1,i2,i3,i4,i5,i6,za,zb,ampQu)
      call WpampLC_down_pp_paper(i1,i2,i3,i4,i5,i6,za,zb,ampQd)
      ampLC=Q(2)*ampQu+Q(1)*ampQd
c      write(6,*) 'LC pp Qu',ampQux/ampQu
c      write(6,*) 'LC pp Qd',ampQdx/ampQd

c      call WpampSL_uppe_pp_full(i1,i2,i3,i4,i5,i6,za,zb,ampQux)
c      call WpampSL_down_pp_full(i1,i2,i3,i4,i5,i6,za,zb,ampQdx)
c      write(6,*) 'SL_uppe_pp',ampQux
c      write(6,*) 'SL_down_pp',ampQdx

      call WpampSL_uppe_pp_paper(i1,i2,i3,i4,i5,i6,za,zb,ampQu)
      call WpampSL_down_pp_paper(i1,i2,i3,i4,i5,i6,za,zb,ampQd)
      ampSL=Q(2)*ampQu+Q(1)*ampQd
c      write(6,*) 'SL pp Qu',ampQux/ampQu
c      write(6,*) 'SL pp Qd',ampQdx/ampQd

      amp(2,2)=ampLC*xn+ampSL/xn

c (-, -)
      call WpampLC_uppe_pp_paper(i2,i1,i4,i3,i5,i6,zb,za,ampQd)
      call WpampLC_down_pp_paper(i2,i1,i4,i3,i5,i6,zb,za,ampQu)
      ampLC=Q(2)*ampQu+Q(1)*ampQd

      call WpampSL_uppe_pp_paper(i2,i1,i4,i3,i5,i6,zb,za,ampQd)
      call WpampSL_down_pp_paper(i2,i1,i4,i3,i5,i6,zb,za,ampQu)
      ampSL=Q(2)*ampQu+Q(1)*ampQd

      amp(1,1)=ampLC*xn+ampSL/xn

c (-, +)
c      call WpampLC_uppe_mp_full(i1,i2,i3,i4,i5,i6,za,zb,ampQux)
c      call WpampLC_down_mp_full(i1,i2,i3,i4,i5,i6,za,zb,ampQdx)
c      write(6,*) 'LC_uppe_mp',ampQux
c      write(6,*) 'LC_down_mp',ampQdx

      call WpampLC_uppe_mp_paper(i1,i2,i3,i4,i5,i6,za,zb,ampQu)
      call WpampLC_down_mp_paper(i1,i2,i3,i4,i5,i6,za,zb,ampQd)
      ampLC=Q(2)*ampQu+Q(1)*ampQd
c      write(6,*) 'LC mp Qu',ampQux/ampQu
c      write(6,*) 'LC mp Qd',ampQdx/ampQd

c      call WpampSL_uppe_mp_full(i1,i2,i3,i4,i5,i6,za,zb,ampQux)
c      call WpampSL_down_mp_full(i1,i2,i3,i4,i5,i6,za,zb,ampQdx)
c      write(6,*) 'SL uppe mp full',ampQux
c      write(6,*) 'SL down mp full',ampQdx
      call WpampSL_uppe_mp_paper(i1,i2,i3,i4,i5,i6,za,zb,ampQu)
      call WpampSL_down_mp_paper(i1,i2,i3,i4,i5,i6,za,zb,ampQd)
      ampSL=Q(2)*ampQu+Q(1)*ampQd
c      write(6,*) 'SL mp Qu',ampQux/AmpQu
c      write(6,*) 'SL mp Qd',ampQdx/AmpQd
c      pause

      amp(1,2)=ampLC*xn+ampSL/xn

c (+, -)
      call WpampLC_uppe_mp_paper(i2,i1,i4,i3,i5,i6,zb,za,ampQd)
      call WpampLC_down_mp_paper(i2,i1,i4,i3,i5,i6,zb,za,ampQu)
      ampLC=Q(2)*ampQu+Q(1)*ampQd

      call WpampSL_uppe_mp_paper(i2,i1,i4,i3,i5,i6,zb,za,ampQd)
      call WpampSL_down_mp_paper(i2,i1,i4,i3,i5,i6,zb,za,ampQu)
      ampSL=Q(2)*ampQu+Q(1)*ampQd

      amp(2,1)=ampLC*xn+ampSL/xn

c Comparison with calculation by coefficients
c      call a6Wpgam_v_coeffs(i1,i2,i3,i4,i5,i6,za,zb,ampold)
c      write(6,*) 'amp(1,1)/ampold(1,1)',amp(1,1)/ampold(1,1)
c      write(6,*) 'amp(1,2)/ampold(1,2)',amp(1,2)/ampold(1,2)
c      write(6,*) 'amp(2,1)/ampold(2,1)',amp(2,1)/ampold(2,1)
c      write(6,*) 'amp(2,2)/ampold(2,2)',amp(2,2)/ampold(2,2)
c      pause


c Terms proportional to (Qu-Qd)

c      call a6Wgamxndk(i2,i1,i4,i3,i5,i6,zb,za,ampdklc)
c      call a6Wgamxninvdk(i2,i1,i4,i3,i5,i6,zb,za,ampdksl)
      call a6Wpgamxndk(i1,i2,i3,i4,i5,i6,za,zb,ampdklc)
      call a6Wpgamxninvdk(i1,i2,i3,i4,i5,i6,za,zb,ampdksl)

c      write(6,*) 'LC_dk_pp',ampdklc(2,2)
c      write(6,*) 'LC_dk_mm',ampdklc(1,1)
c      write(6,*) 'LC_dk_pm',ampdklc(2,1)
c      write(6,*) 'LC_dk_mp',ampdklc(1,2)

c      write(6,*) 'SL_dk_pp',ampdksl(2,2)
c      write(6,*) 'SL_dk_mm',ampdksl(1,1)
c      write(6,*) 'SL_dk_pm',ampdksl(2,1)
c      write(6,*) 'SL_dk_mp',ampdksl(1,2)
c      stop

      amp(:,:)=amp(:,:)+ampdklc(:,:)*xn+ampdksl(:,:)/xn

c      write(6,*) 'amp(2,2)',amp(2,2)
c      write(6,*) 'amp(1,1)',amp(1,1)
c      write(6,*) 'amp(1,2)',amp(1,2)
c      write(6,*) 'amp(2,1)',amp(2,1)

      call a6Wpgam(i1,i2,i3,i4,i5,i6,za,zb,amplo)

c--- strong coupling renormalization in dred scheme
      amp(:,:)=amp(:,:)-xn*(b0/xn*epinv-1._dp/6._dp)*amplo(:,:)
      if (scheme  ==  'tH-V') then
c--- known translation rules
        amp(:,:)=amp(:,:)-(Cf+xn/6d0)*amplo(:,:)
      endif

c      stop

      return
      end

