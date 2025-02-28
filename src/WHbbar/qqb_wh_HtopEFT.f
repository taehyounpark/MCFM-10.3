!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c===== C. Williams Sept 2015
c===== routine which calcualtes the process
c====   q(i1)^++qb(i2)^-=>ell(i3)^-+ell(i4)^+
c     + Higgs where superscripts denote helicity.
c===== this is the piece which goes like the Higgs-top coupling in the EFT (i.e. O(alpha_s^2) in QCD)

c===== Note Higgs decay not handeled in this routine
c===== routine returns msq with no color averaging
      function qqb_wh_HtopEFT(i1,i2,i3,i4,p)
      implicit none
      include 'types.f'
      real(dp) qqb_wh_HtopEFT
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'cutoff.f'
      integer:: i1,i2,i3,i4
      real(dp) :: p(mxpart,4)
      real(dp) :: cutoff_orig
      complex(dp) :: qqb_WH_HtopEFT_loop,EFT_loop
      complex(dp) :: qqb_WH_treeamp,tree
      complex(dp) :: fac_eft,fac_tree
      complex(dp) :: prop_12,prop_34,Ch
      include 'cplx.h'

      qqb_wh_HtopEFT=zip
      EFT_loop=czip
      tree=czip
      cutoff_orig=cutoff

      Ch=as/(3._dp*pi*sqrt(vevsq))

c==== factors
      fac_eft=-ason4pi*gw*(2._dp*cf)*Ch*cone*gw
      fac_tree=gw**2*wmass*gw


c==== just need spinors for 1-4
      call spinoru(4,p,za,zb)

c=====this routine runs in NNLO mode, and can have very small
c=====values of s, since its non-singular lets ensure stability
c==== by doing a smalls cut.
      cutoff=1.E-3_dp
      call smalls(s,4,*999)
      cutoff=cutoff_orig

c======= props
      prop_34=s(i3,i4)/cplx2(s(i3,i4)-wmass**2,wmass*wwidth)
      prop_12=s(i1,i2)/cplx2(s(i1,i2)-wmass**2,wmass*wwidth)

c======== call amplitude for the loop part
      EFT_loop=qqb_WH_HtopEFT_loop(i1,i2,i3,i4,za,zb)


c======= call amplitude for the LO intf.
      tree=qqb_WH_treeamp(i1,i2,i3,i4,za,zb)

c===== dress with factors and propagators
      tree=tree*fac_tree*prop_12*prop_34

c      EFT_loop=tree
      EFT_loop=EFT_loop*fac_EFT*prop_34

c====== perform interference and include overall color charge
      qqb_WH_HtopEFT=real(xn*(conjg(tree)*EFT_loop+conjg(EFT_loop)*tree)
     &     ,dp)

      return

 999  continue
      cutoff=cutoff_orig
      return

      end

      function qqb_WH_treeamp(i1,i2,i3,i4,za,zb)
c==== amplitude for intf. with EFT loop
c     q(i1)^++qb(i2)^-=>ell(i3)^-+ell(i4)^+
      implicit none
      include 'types.f'
      complex(dp) qqb_WH_treeamp
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer i1,i2,i3,i4


      qqb_WH_treeamp=zb(i1,i4)*za(i3,i2)/(s(i1,i2)*s(i3,i4))


      return
      end


      function qqb_WH_HtopEFT_loop(i1,i2,i3,i4,za,zb)
      implicit none
      include 'types.f'
      complex(dp) qqb_WH_HtopEFT_loop
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer i1,i2,i3,i4
      complex(dp) Lsm1_2me,lnrat
      complex(dp) Box,bub134,bub234,rat
      complex(dp) zab2
      real(dp) t,mhsq

      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i1,i3)

      mhsq=t(i1,i3,i4)+t(i2,i3,i4)-s(i3,i4)+s(i1,i2)


c==== box coefficient
      box= (za(i2,i3)**2/(za(i1,i2)*za(i3,i4)) +
     &    zb(i4,i1)**2/(zb(i2,i1)*zb(i4,i3)))/2._dp

      box=box*Lsm1_2me(t(i1,i3,i4),t(i2,i3,i4),s(i3,i4),mhsq)

c=====completed bubble coefficients

      bub134= ((za(i2,i3)*zb(i4,i1))/zab2(i1,i3,i4,i1) -
     &    (za(i1,i2)*za(i3,i4)*zb(i4,i1)**2)/
     &     (2.*zab2(i1,i3,i4,i1)**2))/2._dp

      bub234= ((za(i2,i3)*zb(i4,i1))/zab2(i2,i3,i4,i2) -
     &    (za(i2,i3)**2*zb(i2,i1)*zb(i4,i3))/
     &     (2.*zab2(i2,i3,i4,i2)**2))/2._dp

c      write(6,*) 'bub 134',bub134
c      write(6,*) 'bub 234',bub234


      bub234=bub234*lnrat(-s(i3,i4),-t(i2,i3,i4))
      bub134=bub134*lnrat(-s(i3,i4),-t(i1,i3,i4))
c==== rational
      rat=   (-(za(i2,i3)**2*zb(i2,i1))/
     &     (2.*za(i3,i4)*zab2(i2,i3,i4,i2)) -
     &    (za(i1,i2)*zb(i4,i1)**2)/
     &     (2.*zab2(i1,i3,i4,i1)*zb(i4,i3)))/2.

c      write(6,*) 'rational ',rat

      qqb_WH_HtopEFT_loop=box+bub134+bub234+rat
c      write(6,*) 'total ',qqb_WH_HtopEFT_loop
c      write(6,*) 'CC ',box+(bub134+bub234)
c      write(6,*) 'bubs',(bub134+bub234)
c      write(6,*) 'box',box
c      write(6,*) 'rat ',rat
      return
      end


