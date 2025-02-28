!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine epem3j_v(p,msq)
      implicit none
      include 'types.f'

c--- simple modification of qqb_w_g.f: permuted 1 and 4, 2 and 3
c--- to switch leptons with quarks and added a factor of Nc

c----Matrix element for W + jet production
c----in order alpha_s^2
c----averaged over initial colours and spins
c     q(-p1)+qbar(-p2)-->W^+(nu(p3)+e^+(p4))+g(p5)
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'epinv.f'
      include 'scheme.f'
      include 'nflav.f'

      real(dp):: msq(-nf:nf,-nf:nf),msq0(-nf:nf,-nf:nf),
     & p(mxpart,4),fac,sw,prop,subuv_lc,subuv_tr,
     & qqbWg_lc,qqbWg_slc

      integer,parameter::iqqbg(5)=(/4,3,2,1,5/)

c--set msq=0 to initialize
      msq(:,:)=0._dp
      scheme='dred'

c--calculate spinor and dot-products (using BDK type notation)
      call spinoru(5,p,za,zb)

c--- calculate lowest order
      call epem3j(p,msq0)

c--- UV counterterm contains the finite renormalization to arrive
c--- at MS bar scheme.
c      subuv=ason2pi*ca*(epinv*(11._dp-4._dp*tr*real(nflav)/ca)-1._dp)/6._dp
      subuv_lc=ason2pi*ca*(epinv*11._dp-1._dp)/6._dp
      subuv_tr=ason2pi*ca
     & *(epinv*(-4._dp*tr*real(nflav,kind=dp)/ca))/6._dp

c--- calculate propagator
      sw=s(1,2)
      prop=sw**2/((sw-wmass**2)**2+(wmass*wwidth)**2)

      fac=2._dp*cf*xnsq*xn*gwsq**2*gsq*prop

      call virt5colsep(iqqbg,za,zb,qqbWg_lc,qqbWg_slc)

c--- leading colour (Ca) goes in (0,0)
c--- subleading colour (2*Cf-Ca) goes in (0,1)
c--- Tr terms go in (1,0)
      msq(0,0)=aveqq*fac*qqbWg_lc-subuv_lc*msq0(0,0)
      msq(0,1)=aveqq*fac*qqbWg_slc
      msq(1,0)=-subuv_tr*msq0(0,0)

      return
      end


c--- this is a copy of the routine virt5.f, but changed from a function
c---  to a subroutine and modified to extract different colour pieces
      subroutine virt5colsep(ip,za,zb,virt5lc,virt5slc)
      implicit none
      include 'types.f'

c***********************************************************************
c     Author: R.K. Ellis                                               *
c     July, 1999.                                                      *
c   Given za and zb calculate the                                      *
c   the interference of the amplitude for the process                  *
c   0--> qb_R(1)+q_L(2)+l_L(3)+a_R(4)+g_L/R(5)                         *
c   at one loop with the corresponding lowest order amplitude          *
c   summed over the polarizations of the emitted gluon                 *
c   Virtual terms are in units of
c   (as/4/pi) (4 pi)^ep Gamma(1+ep)*Gamma(1-ep)^2/Gamma(1-2*ep)
c***********************************************************************
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'zprods_decl.f'
      integer:: ip(5)
      complex(dp):: A5LOm,A5NLOm_lc,A5NLOm_slc,
     &               A5LOp,A5NLOp_lc,A5NLOp_slc
      real(dp):: virt5lc,virt5slc

c   0--> qb_R(1)+q_L(2)+l_L(3)+a_R(4)+g_L(5)
      call A5NLOcolsep(ip(1),ip(2),ip(3),ip(4),ip(5),za,zb,
     &                 A5LOm,A5NLOm_lc,A5NLOm_slc)
c   0--> qb_R(1)+q_L(2)+l_L(3)+a_R(4)+g_R(5)
      call A5NLOcolsep(ip(2),ip(1),ip(4),ip(3),ip(5),zb,za,
     &                 A5LOp,A5NLOp_lc,A5NLOp_slc)

      virt5lc=ason2pi*(
     & Dble(conjg(A5LOp)*A5NLOp_lc)+Dble(conjg(A5LOm)*A5NLOm_lc))
      virt5slc=ason2pi*(
     & Dble(conjg(A5LOp)*A5NLOp_slc)+Dble(conjg(A5LOm)*A5NLOm_slc))

      return
      end



c--- this is a copy of the routine a5nlo.f, but modified to return
c---  different colour pieces separately
      subroutine A5NLOcolsep(j1,j2,j3,j4,j5,za,zb,
     &                       A5LOm,A5NLOm_lc,A5NLOm_slc)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5
      complex(dp):: A51,A52,A5NLOm_lc,A5NLOm_slc,A5LOm

c As originally written, the functions A51, A52 correspond to
c 0 --> q_R(1)+qb_L(3)+g_R(2)+ebar_L(4)+e_R(5)
c with all RH couplings
c However we want it in our
c standard form
c       0--> qb_R(1)+q_L(2)++e_L(3)+ebar_R(4)+g_L(5)
c with all LH couplings

c so we have made the changes

c                    'q+g+qb-'   (A51)
c                   (1 ---> 2)
c                   (2 ---> 5)
c                   (3 ---> 1)
c                   (4 ---> 4)
c                   (5 ---> 3)

c                    'q+qb-g+'   (A52)
c                   (1 ---> 2)
c                   (2 ---> 1)
c                   (3 ---> 5)
c                   (4 ---> 4)
c                   (5 ---> 3)

c  and also exchanged za and zb.

c--- corresponds to (1V.1) times minus i, with the (A51) change
      A5LOm=-zb(j1,j4)**2/(zb(j2,j5)*zb(j5,j1)*zb(j4,j3))
      A5NLOm_lc=A51(j2,j5,j1,j4,j3,zb,za)
      A5NLOm_slc=A52(j2,j1,j5,j4,j3,zb,za)/xnsq

      return
      end
