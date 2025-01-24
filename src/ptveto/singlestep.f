!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function Cggh(mt,qsq,MuH)
!     single step Wilson coefficient
!     1012.4480, Eqn 2.7
      use LHAPDF, only: getalphas
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp)::mt,qsq,muH,z,ason4piH,alphasMuH
      complex(dp):: Cggh,CggH1,CggH2,ggHf0,ggHf1,ggHf2
      alphasMuH=getalphas(MuH)
      ason4piH=alphasMuH/(4._dp*pi)
      z=qsq/(4._dp*mt**2)
      CggH=alphasMuH*ggHf0(z)*(1._dp
     & +ason4piH*(CggH1(qsq,MuH**2)+ggHf1(z))
     & +ason4piH**2*(CggH2(qsq,MuH**2,z)+ggHf2(z)))
      return
      end

      function CggH1(qsq,muhsq)
!     single step Wilson coefficient
!     1012.4480, Eqn 2.8
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp)::qsq,muhsq
      complex(dp):: CggH1,lnrat
      CggH1=CA*(-lnrat(-qsq,muhsq)**2+pisq/6._dp)
      return
      end

      function CggH2(qsq,muhsq,z)
!     single step Wilson coefficient
!     1012.4480, Eqn 2.9
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'zeta.f'
      include 'nfl.f'
      real(dp)::qsq,muhsq,z,beta0
      complex(dp):: CggH2,ggHf1,ln,lnrat
      ln=lnrat(-qsq,muhsq)
      beta0=(11._dp*CA-4._dp*TF*nfl)/3._dp
      CggH2=0.5_dp*CA**2*ln**4+1/3._dp*CA*beta0*ln**3
     & +CA*((-4/3._dp+pisq/6._dp)*CA-5/3._dp*beta0-ggHf1(z))*ln**2
     & +((59/9._dp-2*zeta3)*CA**2+(19/9._dp-pisq/3)*CA*beta0
     & -ggHf1(z)*beta0)*ln
      return
      end


      function ggHf0(z)
!     1012.4480, Eqn B.1, z=qsq/4/mt**2
      implicit none
      include 'types.f'
      real(dp)::z,realbit
      complex(dp)::arg,ggHf0,fac
      if ((z > 0._dp) .and. (z <= 1._dp)) then
         fac=asin(sqrt(z))**2
      elseif (z > 1) then
         realbit=sqrt(z)+sqrt(1._dp-z)
         arg=cmplx(0._dp,-realbit)
         fac=log(arg)**2
      endif
      ggHf0=1.5_dp/z*(1._dp-abs(1._dp-1._dp/z)*fac)
!      ggHf0=1._dp ! debug for 1-step/2-step comparison
      return
      end

      function ggHf1(z)
!     1012.4480, Eqn B.1, z=qsq/4/mt**2
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp)::z
      complex(dp)::ggHf1
      ggHf1=CA*(5._dp-38._dp/45._dp*z-1289._dp/4725._dp*z**2
     & -155._dp/1134._dp*z**3-5385047._dp/65488500._dp*z**4)
     & +CF*(-3._dp+307._dp/90._dp*z+25813._dp/18900._dp*z**2
     & +3055907._dp/3969000._dp*z**3+659504801._dp/1309770000._dp*z**4)
!      ggHf1=5._dp*CA-3._dp*CF ! debug for 1-step/2-step comparison
      return
      end

      function ggHf2(z)
!     1012.4480, Eqn B.2, z=qsq/4/mt**2
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'zeta.f'
      include 'nfl.f'
      real(dp)::z,beta0
      real(dp),parameter::ln2=0.6931471805599453_dp
      complex(dp)::ggHf2,ln4z,f2(0:2)
      integer j      
      if (z > 0) then 
         ln4z=log(+4._dp*z)-im*pi
      else
         ln4z=log(-4._dp*z)
      endif   
      beta0=(11._dp*CA-4._dp*TF*nfl)/3._dp
      f2(0)=(7*CA**2+11*CA*CF-6*CF*beta0)*ln4z
     & +(-419/27._dp+7*pisq/6._dp+pisq**2/72._dp-44*zeta3)*CA**2
     & +(-217/2._dp-pisq/2._dp+44*zeta3)*CA*CF
     & +(2255/108._dp+5*pisq/12._dp+23*zeta3/3._dp)*CA*beta0
     & -5/6._dp*CA*TF+27/2._dp*CF**2
     & +(41/2._dp-12*zeta3)*CF*beta0-4/3._dp*CF*TF
      f2(1)=CA**2*(-404063/14400._dp+11723/384._dp*zeta3
     & -223/108._dp*ln4z-19/135._dp*pisq)
     & +CF*CA*(-1099453/8100._dp+2297/16._dp*zeta3-242/135._dp*ln4z
     & -953/540._dp*pisq+28/15._dp*pisq*ln2)
     & +CF**2*(-36803/240._dp+13321/96._dp*zeta3+7/3._dp*pisq
     & -56/15._dp*pisq*ln2)
     & +CF*(-4393/405._dp+77/12._dp*zeta3-7337/2700._dp*beta0
     & +39/10._dp*ln4z*beta0+28/45._dp*pisq+7/15._dp*pisq*beta0)
     & +CA*(-64097/129600._dp+77/384._dp*zeta3-269/75._dp*beta0
     & +2/15._dp*ln4z-31/180._dp*ln4z*beta0)
      f2(2)=CA**2*(-3.084463261e9_dp/254016000._dp+110251/9216._dp*zeta3
     & -2869/4536._dp*ln4z-1289/28350._dp*pisq)
     & +CF*CA*(-5.5535378557e10_dp/381024000._dp+2997917/23040._dp*zeta3
     & -18337/28350._dp*ln4z-128447/113400._dp*pisq
     & +1714/1575._dp*pisq*ln2)
     & +CF**2*(-95081911/453600._dp+36173/192._dp*zeta3+857/630._dp*pisq
     & -3428/1575._dp*pisq*ln2)
     & +CA*(+265053121/1524096000._dp-16177/92160._dp*zeta3
     & -45617/47250._dp*beta0+16/315._dp*ln4z-623/5400._dp*ln4z*beta0)
     & +CF*(-8108339/1555200._dp+21973/7680._dp*zeta3
     & -509813/3969000._dp*beta0-8/15._dp*ln4z
     & +29147/18900._dp*ln4z*beta0
     & +1714/4725._dp*pisq+857/3150._dp*pisq*beta0)

      ggHf2=f2(0)
!      return      ! debug for 1-step/2-step comparison
      do j=1,2
         ggHf2=ggHf2+z**j*f2(j)
      enddo

      return
      end
      
