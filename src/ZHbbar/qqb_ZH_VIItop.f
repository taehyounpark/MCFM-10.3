!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function qqb_ZH_VIItop(i1,i2,i3,i4,p,j)
c===== C. Williams Sept 2015
c===== routine which calculates the VII style contributions to
c==== qqb => ZH (see e.g. Brein, Harlander Wisemann and Zirke)
      implicit none
      include 'types.f'
      real(dp) :: qqb_ZH_VIItop
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'nf.f'
      include 'zcouple.f'
      include 'cutoff.f'
      integer:: i1,i2,i3,i4,j
      integer h1,h2
      real(dp) :: p(mxpart,4),cl(2),cq(2),cutoff_orig
      complex(dp) :: qqb_ZH_asymt,asympt(2,2)
      complex(dp) :: qqb_WH_treeamp,tree(2,2)
      complex(dp) :: fac_asym,fac_tree
      complex(dp) :: prop_12,prop_34
      real(dp) :: ccou(2,2)

      include 'cplx.h'

      qqb_ZH_VIItop=zip


c==== debug, routine is written for leptons, needs generalization...
      cl(1)=l1
      cl(2)=r1
c      if((nproc /= 101)
c     &     .or.(nproc /= 104).or.(nproc /= 106).or.(nproc /= 109)) then
c         write(6,*) 'Warning, only leptonic decays implemented in '
c         write(6,*) 'file qqb_ZH_VIItop.f'
c         write(6,*) 'nproc = 'nproc,' not yet supported'
c         stop
c      endif

      cq(1)=-l(j)
      cq(2)=r(j)


c==== factors
      fac_asym=-(ason4pi)**2*Gf*wmass*(2._dp*cf)*sqrt(esq)*(l(4)-r(4))*8._dp
      fac_tree=esq*wmass*gw/(one-xw)*2._dp

c==== just need spinors for 1-4
      call spinoru(4,p,za,zb)

c=====this routine runs in NNLO mode, and can have very small
c=====values of s, since its non-singular lets ensure stability
c==== by doing a smalls cut.
      cutoff_orig=cutoff
      cutoff=1.E-3_dp
      call smalls(s,4,*999)
      cutoff=cutoff_orig

c======= props
      prop_34=s(i3,i4)/cplx2(s(i3,i4)-zmass**2,zmass*zwidth)
      prop_12=s(i1,i2)/cplx2(s(i1,i2)-zmass**2,zmass*zwidth)

c======== cal
c==== fill tree array and assympt

      tree(1,1)=qqb_WH_treeamp(i2,i1,i3,i4,za,zb)
      tree(2,1)=qqb_WH_treeamp(i1,i2,i3,i4,za,zb)
      tree(1,2)=qqb_WH_treeamp(i2,i1,i4,i3,za,zb)
      tree(2,2)=qqb_WH_treeamp(i1,i2,i4,i3,za,zb)

      asympt(1,1)=qqb_ZH_asymt(i2,i1,i3,i4,za,zb)
      asympt(2,1)=qqb_ZH_asymt(i1,i2,i3,i4,za,zb)
      asympt(1,2)=qqb_ZH_asymt(i2,i1,i4,i3,za,zb)
      asympt(2,2)=qqb_ZH_asymt(i1,i2,i4,i3,za,zb)


      asympt(:,:)=asympt(:,:)*prop_34*fac_asym
      tree(:,:)=tree(:,:)*prop_12*prop_34*fac_tree

c===== chiral couplings
      qqb_ZH_VIItop=zip
      do h1=1,2
         do h2=1,2
            ccou(h1,h2)=cq(h1)*cl(h2)**2
            qqb_ZH_VIItop=qqb_ZH_VIItop+
     &    xn*ccou(h1,h2)*real((conjg(tree(h1,h2))*asympt(h1,h2)
     &           +conjg(asympt(h1,h2))*tree(h1,h2)),dp)

         enddo
      enddo

      return

 999  continue
      cutoff=cutoff_orig
      return

      end


      function qqb_ZH_asymt(i1,i2,i3,i4,za,zb)
c==== amplitude for intf. with EFT loop
c     q(i1)^++qb(i2)^-=>ell(i3)^-+ell(i4)^+
      implicit none
      include 'types.f'
      complex(dp)::  qqb_ZH_asymt
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: i1,i2,i3,i4


      qqb_ZH_asymt=zb(i1,i4)*za(i3,i2)/(s(i3,i4))

      return
      end

