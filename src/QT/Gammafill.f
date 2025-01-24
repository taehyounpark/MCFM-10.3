!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine Gammafill(nfl1)
      implicit none
c     1909.00811v2
      include 'types.f'
      include 'constants.f'
      include 'zeta.f'
      include 'quarkorgluon.f'
      include 'qtconstants.f'
      include 'nfl.f'
      integer:: i,nfl1
      real(dp), parameter:: C(0:1)=[CA,CF]
      nfl=nfl1

c     Eqn.D.2
      beta0=11*CA/3._dp-4*TF*nfl/3._dp
      beta1=34*CA**2/3._dp-2*TF*nfl*(10*CA/3._dp+2*CF)
      beta2=2857*CA**3/54._dp+2*TF*nfl
     & *(CF**2-1415*CA**2/54._dp-205*CF*CA/18._dp)
     & +4*(TF*nfl)**2*(79*CA/54._dp+11*CF/9._dp)


c     Eqn.D.5
      gammaB0(qq)=6*CF
      gammaB1(qq)=2*CF*(CA*(73/9._dp-40*zeta3)
     & +CF*(1.5d0-12*zeta2+24*zeta3)
     & +beta0*(121/18._dp+2*zeta2))

      gammaB2(qq)=2*CF*(CA**2*(
     & 52019/162._dp-1682/27._dp*zeta2-2056/9._dp*zeta3
     & -820/3._dp*zeta4+176/3._dp*zeta2*zeta3+232*zeta5)
     & +CA*CF*(151/4._dp-410/3._dp*zeta2+844/3._dp*zeta3
     & -494/3._dp*zeta4+16*zeta2*zeta3+120*zeta5)
     & +CF**2*(29/2._dp+18*zeta2+68*zeta3+144*zeta4
     & -32*zeta2*zeta3-240*zeta5)
     & +CA*beta0*(-7739/54._dp+650/27._dp*zeta2
     & -1276/9._dp*zeta3+617/3._dp*zeta4)
     & +beta0**2*(-3457/324._dp+10/3._dp*zeta2+16/3._dp*zeta3)
     & +beta1
     & *(1166/27._dp-16/3._dp*zeta2+52/9._dp*zeta3-82/3._dp*zeta4))

c     Eqn.D.6
      gammaB0(gg)=2*beta0
      gammaB1(gg)=2*CA*(CA*(91/9._dp-16*zeta3)
     & +beta0*(47/9._dp-2*zeta2))+2*beta1
      gammaB2(gg)=2*CA*(
     & CA**2*(49373/162._dp-944/27._dp*zeta2-2260/9._dp*zeta3
     & -144*zeta4+128/3._dp*zeta2*zeta3+112*zeta5)
     & +CA*beta0
     & *(-6173/54._dp-376/27._dp*zeta2+140/9._dp*zeta3+117*zeta4)
     & +beta0**2*(-493/81._dp-10/3._dp*zeta2+28/3._dp*zeta3)
     & +beta1*(1765/54._dp-2*zeta2-152/9._dp*zeta3-8*zeta4))+2*beta2

      do i=0,1

c     Eqn. D.7
      gammaS0(i)=0
      gammaS1(i)=2*C(i)
     & *(CA*(-64/9._dp+28*zeta3)+beta0*(-56/9._dp+2*zeta2))
      gammaS2(i)=2*C(i)*(CA**2
     & *(-37871/162._dp+620/27._dp*zeta2+2548/9._dp*zeta3
     & +144*zeta4-176/3._dp*zeta2*zeta3-192*zeta5)
     & +CA*beta0
     & *(4697/54._dp+484/27._dp*zeta2+220/9._dp*zeta3-112*zeta4)
     & +beta0**2*(520/81._dp+10/3._dp*zeta2-28/3._dp*zeta3)
     & +beta1*(-1711/54._dp+2*zeta2+152/9._dp*zeta3+8*zeta4))


c    Eqn. D.9
      tgammaS0(i)=-gammaS0(i)
      tgammaS1(i)=-gammaS1(i)
      tgammaS2(i)=-gammaS2(i)
      tgammaB0(i)=gammaB0(i)+gammaS0(i)
      tgammaB1(i)=gammaB1(i)+gammaS1(i)
      tgammaB2(i)=gammaB2(i)+gammaS2(i)

c     Eqn. D.11
      tildes0(i)=1._dp
      tildes1(i)=-C(i)*2*zeta2
      tildes2(i)=C(i)*(C(i)*5*zeta4+CA*(208/27._dp-4*zeta2+10*zeta4)
     & +beta0*(164/27._dp-5*zeta2-14*zeta3/3._dp))
      tildes3(i)=C(i)*(-C(i)**2*35/6._dp*zeta6
     & +C(i)*CA*(-416/27._dp*zeta2+20*zeta4-35*zeta6)
     & +C(i)*beta0*(-328/27._dp*zeta2+25*zeta4+28/3._dp*zeta2*zeta3)
     & +CA**2*(115895/324._dp-51071/486._dp*zeta2
     & -23396/81._dp*zeta3-58*zeta4+240*zeta2*zeta3-224*zeta5
     & +928/9._dp*zeta3**2-3086/27._dp*zeta6)
     & +CA*beta0*(-363851/2916._dp+2987/486._dp*zeta2-428/81._dp*zeta3
     & +830/9._dp*zeta4-220/3._dp*zeta2*zeta3+1388/9._dp*zeta5)
     & +beta0**2*(-64/729._dp-34/3._dp*zeta2-140/27._dp*zeta3
     & -11/3._dp*zeta4)
     & +beta1*(42727/972._dp-275/18._dp*zeta2-1744/81._dp*zeta3
     & -76/9._dp*zeta4+40/3._dp*zeta2*zeta3-112/9._dp*zeta5))

c     Eqn.D.4
      Gamma0(i)=4*C(i)
      Gamma1(i)=4*C(i)*(CA*(67/9._dp-2*zeta2)-20/9._dp*TF*nfl)
      Gamma2(i)=4*C(i)
     & *(CA**2*(245/6._dp-268*zeta2/9._dp+22*zeta3/3._dp+22*zeta4)
     & +2*TF*nfl*(CA*(-209/27._dp+40*zeta2/9._dp-28*zeta3/3._dp)
     & +CF*(-55/6._dp+8*zeta3))-16/27._dp*(TF*nfl)**2)

c     Eqn.D.10
      tgamman0(i)=0
      tgamman1(i)=2*C(i)*(CA*(28*zeta3-64/9._dp)-56*beta0/9._dp)
      tgamman2(i)=2*C(i)*(
     & CA**2*(-37871/162._dp+620/27._dp*zeta2+2548/9._dp*zeta3
     & +144*zeta4-176/3._dp*zeta2*zeta3-192*zeta5)
     & +CA*beta0*(3865/54._dp+412/27._dp*zeta2+220/9._dp*zeta3-50*zeta4)
     & +beta0**2*(-464/81._dp-8*zeta3)
     & +beta1*(-1711/54._dp+152/9._dp*zeta3+8*zeta4))

c 1602.01829v3, Eq.4.53
      S2(i)=C(i)**2*pi**4/18._dp
     & +C(i)*(CA*(208/27._dp-2*pisq/3._dp+pi**4/9._dp)
     & +beta0*(164/27._dp-5*pisq/6._dp-14*zeta3/3._dp))
      enddo
      end
