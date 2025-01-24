!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function A53_VH(j1,j2,j3,j4,j5,za,zb,mtsq)
c  Amplitudes taken from Appendix IV of
c  Z.~Bern, L.~J.~Dixon and D.~A.~Kosower,
c  %``One loop amplitudes for e+ e- to four partons,''
c  Nucl.\ Phys.\ B {\bf 513}, 3 (1998)
c  [hep-ph/9708239].
c  Eq.(IV.11)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      complex(dp):: A53_VH,L1,kinstructure,kinstructure2,F1anom,F3anom,G1,G3,MZsq
      integer:: j1,j2,j3,j4,j5
      real(dp):: s45,s12,mtsq,musq
c      logical,parameter::largemassexp=.false.
      logical,parameter::tandbtogether=.true.
      s12=s(j1,j2)
      s45=s(j1,j2)+s(j2,j3)+s(j3,j1)
      kinstructure=two*zb(j5,j3)*zb(j3,j1)*za(j2,j4)/s45
c      kinstructure2=
c     & -zb(j1,j3)*(za(j4,j1)*zb(j1,j5)+za(j4,j2)*zb(j2,j5)+za(j4,j3)*zb(j3,j5))/s45
      kinstructure2=(za(j4,j1)*zb(j1,j5)+za(j4,j2)*zb(j2,j5)+za(j4,j3)*zb(j3,j5))
     & *zb(j1,j3)**2/zb(j1,j2)/s45
      MZsq=zmass**2-im*zmass*zwidth
      A53_VH=czip
c     musq is irrelevant, set to some value
      musq=mtsq
      if (tandbtogether) then
c        A53_VH=-kinstructure*L1(-s12,-s45)/two/s45
c        A53_VH=A53_VH-kinstructure2*F3anom(s12,s45,0._dp,musq)
        G1=F1anom(s12,s45,0._dp,musq)
        G3=F3anom(s12,s45,0._dp,musq)
        A53_VH=-(kinstructure*G1+kinstructure2*(s12/MZsq*G1+(s45/MZsq-one)*G3))
      endif
      if (mtsq == 0._dp) then
         A53_VH=+kinstructure*L1(-s12,-s45)/two/s45
         A53_VH=A53_VH+kinstructure2*F3anom(s12,s45,0._dp,musq)
      else
c         if (largemassexp) then
c            A53_VH=A53_VH+kinstructure/(24._dp*mtsq)
c     &      *(1._dp+(2._dp*s12+s45)/15._dp/mtsq*zip
c     &      +(2._dp*s12*s45+3._dp*s12**2+s45**2)/140._dp/mtsq**2*zip)
c         else
c            A53_VH=A53_VH+kinstructure*F1anom(s12,s45,mtsq,musq)
c            A53_VH=A53_VH+kinstructure2*F3anom(s12,s45,mtsq,musq)
        G1=F1anom(s12,s45,mtsq,musq)
        G3=F3anom(s12,s45,mtsq,musq)
        A53_VH=A53_VH+(kinstructure*G1+kinstructure2*(s12/MZsq*G1+(s45/MZsq-one)*G3))
c         endif
      endif
      return
      end
