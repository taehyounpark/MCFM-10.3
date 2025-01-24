!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine tildep0(z,tp0)
      implicit none
      include 'types.f'
      include 'singletlabels.f'
c     gg=0,qqV=1,qbV=2,qqS=3,qg=4,gq=5,qqDS=6
      include 'distributions.f'
c     dmin=-2,dmax=1,rglr=-2,delt=-1,plus=0,lpls=1
      real(dp):: z,omz,opzsq,tp0(0:6,dmin:dmax)

c     initialize to zero
      tp0(:,:)=0

      omz=1._dp-z
      opzsq=1._dp+z**2
c     1403.6451v2,Eq C7
c     function called for both positive and negative z
      if (z > 0._dp) then
      tp0(qqV,plus)=2._dp;tp0(qqV,rglr)=-1._dp-z
      tp0(gg,plus)=1._dp; tp0(gg,rglr)=-1._dp+omz*opzsq/z
      else
      tp0(qqV,rglr)=opzsq/omz
      tp0(gg,rglr)=(opzsq-z)**2/z/omz
      endif
      tp0(qg,rglr)=z**2+omz**2
      tp0(gq,rglr)=(1._dp+omz**2)/z

      return
      end
