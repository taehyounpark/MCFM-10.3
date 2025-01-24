!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine CggHsq_onestep(Mhsq,musq,coeff0,coeff)
      use LHAPDF, only: getalphas
      implicit none
c   Wilson Coefficients for gg->H
c   Taken from Ahrens,Becher,Neubert,Yang, 0809.4283v3
c   Also allows for timelike value of musq and corresponding complex alphas
      include 'types.f'
      include 'constants.f'
      include 'masses.f'
      include 'nflav.f'
      include 'scet_const.f'
      real(dp),intent(in)::Mhsq,musq
      real(dp),intent(out)::coeff0,coeff(2)
      real(dp):: z
      complex(dp)::L,lnrat
      complex(dp)::Cggh(0:2)
      complex(dp)::ason4pi
      complex(dp):: CggH1,CggH2,ggHf1,ggHf2,ggHf0
      complex(dp)::getason2pic

      z=Mhsq/(4._dp*mt**2)
      Cggh(1)=CggH1(MHsq,musq)+ggHf1(z)
      Cggh(2)=CggH2(MHsq,musq,z)+ggHf2(z)

! allow for complex alpha-s (timelike musq)
      ason4pi=getason2pic(musq)/2._dp

! apply factors of ason4pi
      Cggh(1)=Cggh(1)*ason4pi
      Cggh(2)=Cggh(2)*ason4pi**2

! overall factor from leading order, including reweighting from
! (always real) alphas in EFT to possibly-complex alphas for musq < 0
      CggH(0)=ggHf0(z)*ason4pi/(getalphas(sqrt(abs(musq)))/fourpi)

      coeff0=abs(CggH(0))**2
      coeff(1)=real(2*Cggh(1),dp)
      coeff(2)=(2*real(Cggh(2),dp)
     &         +real(Cggh(1)*conjg(Cggh(1)),dp))

      return
      end

