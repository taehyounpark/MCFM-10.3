!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module evolution
      implicit none
      public evolfill_r,evolfill_c

      interface evolfill
      module procedure evolfill_r,evolfill_c
      end interface

      contains

      subroutine evolfill_r(alphasMu,alphasMuh,order,
     & Gamma0,Gamma1,Gamma2,Gamma3,Gamma4,
     & gammai0,gammai1,gammai2,gammai3,
     & beta0,beta1,beta2,beta3,beta4,
     & S,aGamma,agammai)
      implicit none
      include 'types.f'
      include 'constants.f'
      integer,intent(in):: order
      real(dp),intent(in):: alphasMu,alphasMuH
      real(dp),intent(out):: S,aGamma,agammai
      real(dp):: r,lnr,S0,S1,S2,S3,
     &           aGamma0,aGamma1,aGamma2,aGamma3,
     &           agammai0,agammai1,agammai2,agammai3

      real(dp),intent(in):: 
     & Gamma0,Gamma1,Gamma2,Gamma3,Gamma4,
     & gammai0,gammai1,gammai2,gammai3,
     & beta0,beta1,beta2,beta3,beta4
      include 'evolution_include.f'
      end subroutine evolfill_r

      subroutine evolfill_c(alphasMu,alphasMuh,order,
     & Gamma0,Gamma1,Gamma2,Gamma3,Gamma4,
     & gammai0,gammai1,gammai2,gammai3,
     & beta0,beta1,beta2,beta3,beta4,
     & S,aGamma,agammai)
      implicit none
      include 'types.f'
      include 'constants.f'
      integer,intent(in):: order
      complex(dp),intent(in):: alphasMu,alphasMuH
      complex(dp),intent(out):: S,aGamma,agammai
      complex(dp):: r,lnr,S0,S1,S2,S3,
     &              aGamma0,aGamma1,aGamma2,aGamma3,
     &              agammai0,agammai1,agammai2,agammai3

      real(dp),intent(in):: 
     & Gamma0,Gamma1,Gamma2,Gamma3,Gamma4,
     & gammai0,gammai1,gammai2,gammai3,
     & beta0,beta1,beta2,beta3,beta4

      include 'evolution_include.f'
      end subroutine evolfill_c

      end module evolution

