!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqcoeff(Qsq,musq,coeff)
      implicit none
c   Wilson Coefficients for qq->V in units of
c   as/4/pi
      include 'types.f'
      include 'constants.f'
      include 'scet_const.f'
      real(dp),intent(in)::Qsq,musq
      complex(dp)::L,lnrat,HF,HA,HTF
      complex(dp),intent(out)::coeff(2)

      include 'scet_beta.f'

      L=lnrat(-Qsq,musq)

c     0607228, Eq.51
c     Hard functions are simply related to MSbar subtracted
c     on-shell matrix elements
c     see also 0605068, Eq 24
      coeff(1)=CF*(-(L**2-3*L+8._dp-zeta2))

      HF=0.5_dp*L**4-3*L**3+(12.5_dp-pisq/6)*L**2
     & +(-22.5_dp-3*pisq/2+24*zeta3)*L
     & +255/8._dp+3.5_dp*pisq-83*pi**4/360._dp-30*zeta3
      HA=11/9._dp*L**3+(pisq/3._dp-233/18._dp)*L**2
     & +(2545/54._dp+11/9._dp*pisq-26*zeta3)*L
     & -51157/648._dp-337*pisq/108._dp+11*pi**4/45._dp+313/9._dp*zeta3
      HTF=(-4/9._dp*L**3+38/9._dp*L**2
     & -(418/27._dp+4*pisq/9._dp)*L
     & +4085/162._dp+23*pisq/27._dp+4/9._dp*zeta3)
      coeff(2)=CF*(CF*HF+CA*HA+TF*nflav*HTF)
      return
      end

      subroutine qqcoeff3(Qsq,musq,coeff,NfV)
          use types
          implicit complex(dp) (T)
!   Wilson Coefficients for qq->V in units of 
!   as/4/pi
          include 'src/Inc/constants.f'
          include 'src/Inc/nf.f'
          include 'src/Inc/scet_const.f'

          real(dp), intent(in) :: Qsq, musq
          complex(dp)::LL,lnrat
          complex(dp),intent(out) :: coeff(3), NfV
          complex(dp) :: R4,R5
          
          LL=lnrat(-Qsq,musq)

      coeff(1)=CF*(-(LL**2-3*LL+8._dp-zeta2))
      coeff(2)=CF**2*(1/2._dp*(LL**2-3*LL+8._dp-zeta2)**2
     & +(3/2._dp-12*zeta2+24*zeta3)*LL
     & -1/8._dp+29._dp*zeta2-30*zeta3-44/5._dp*zeta2**2)
     & +CF*nf*(-2/9._dp*LL**3+19/9._dp*LL**2-(209/27._dp+4/3._dp*zeta2)*LL
     & +4085/324._dp+23/9._dp*zeta2+2/9._dp*zeta3)
     & +CF*CA*(11/9._dp*LL**3+(2*zeta2-233/18._dp)*LL**2
     & +(2545/54._dp+22/3._dp*zeta2-26*zeta3)*LL
     & -51157/648._dp-zeta2*(337/18._dp-44/5._dp*zeta2)+313/9._dp*zeta3)

      T2=CF*CF
      T3=CF**3
      T4=CA*CA
      T5=1.3333333333333333333D0*CA*CF
      T6=LL*LL
      T7=nf*nf
      T8=LL**3
      T9=LL**4
      T10=LL**5
      T11=((-46.222222222222222222D0*T2-T5)*nf+(-240.0D0*T3+120.0D0*CA*T
     &2+136.0D0*T4*CF)*LL+664.0D0*T3-306.22222222222222222D0*CA*T2-48.22
     &2222222222222222D0*T4*CF)*zeta5+(16.0D0*T3+98.666666666666666667D0
     &*CA*T2-126.22222222222222222D0*T4*CF)*zeta3*zeta3+(((-4.2222222222
     &222222222D0*T2+T5)*nf+(-8.0D0*T3-10.0D0*CA*T2+29.33333333333333333
     &3D0*T4*CF)*LL+250.0D0*T3-416.77777777777777778D0*CA*T2+138.6666666
     &6666666667D0*T4*CF)*zeta2+(0.59259259259259259259D0*CF*LL-1.711934
     &156378600823D0*CF)*T7+((7.7777777777777777778D0*T2-8.0D0*CA*CF)*T6
     &+(-67.777777777777777778D0*T2+80.444444444444444444D0*CA*CF)*LL+16
     &2.7654320987654321D0*T2-158.81481481481481481D0*CA*CF)*nf+(-24.0D0
     &*T3+26.0D0*CA*T2)*T8+(102.0D0*T3-200.77777777777777778D0*CA*T2+88.
     &0D0*T4*CF)*T6+(-214.0D0*T3+813.66666666666666667D0*CA*T2-646.81481
     &481481481481D0*T4*CF)*LL-470.0D0*T3-695.18518518518518519D0*CA*T2+
     &1039.2736625514403292D0*T4*CF)*zeta3+(59.887301587301587301D0*T3-4
     &0.24126984126984127D0*CA*T2-32.550264550264550265D0*T4*CF)*zeta2**
     &3
      T12=T11+(-1.3925925925925925926D0*CF*T7+((5.6D0*T2+2.9333333333333
     &333333D0*CA*CF)*LL-12.259259259259259259D0*T2+0.074074074074074074
     &074D0*CA*CF)*nf+(8.3D0*T3-6.8D0*CA*T2-8.8D0*T4*CF)*T6+(20.7D0*T3+3
     &2.4D0*CA*T2-31.333333333333333333D0*T4*CF)*LL-82.6D0*T3-18.3074074
     &07407407407D0*CA*T2+82.062962962962962963D0*T4*CF)*zeta2*zeta2+((-
     &0.88888888888888888889D0*CF*T6+5.6296296296296296296D0*CF*LL-10.17
     &2839506172839506D0*CF)*T7+((1.1111111111111111111D0*T2+0.888888888
     &88888888889D0*CA*CF)*T8+(-12.444444444444444444D0*T2+5.33333333333
     &33333333D0*CA*CF)*T6+(59.925925925925925926D0*T2-72.39506172839506
     &1729D0*CA*CF)*LL-97.929012345679012346D0*T2+158.51165980795610425D
     &0*CA*CF)*nf+(0.5D0*T3-2.0D0*CA*T2)*T9+(9.0D0*T3-0.1111111111111111
     &1111D0*CA*T2-4.8888888888888888889D0*T4*CF)*T8+(-52.5D0*T3+55.7777
     &77777777777778D0*CA*T2+2.8888888888888888889D0*T4*CF)*T6+(178.5D0*
     &T3-417.03703703703703704D0*CA*T2+214.39506172839506173D0*T4*CF)*LL
     &-268.79166666666666667D0*T3+831.53549382716049383D0*CA*T2-565.5898
     &4910836762688D0*T4*CF)*zeta2
      R4=T12+(-0.074074074074074074074D0*CF*T9+0.93827160493827160494D0
     &*CF*T8-5.012345679012345679D0*CF*T6+13.495198902606310014D0*CF*LL-
     &14.550449626581313824D0*CF)*T7+(0.22222222222222222222D0*T2*T10+(-
     &2.7777777777777777778D0*T2+0.81481481481481481482D0*CA*CF)*T9+(15.
     &185185185185185185D0*T2-12.024691358024691358D0*CA*CF)*T8+(-39.552
     &469135802469136D0*T2+72.543209876543209876D0*CA*CF)*T6+(28.8981481
     &48148148148D0*T2-212.50891632373113855D0*CA*CF)*LL+42.260288065843
     &621399D0*T2+259.13290656912056089D0*CA*CF)*nf-0.166666666666666666
     &67D0*T3*LL**6+(1.5D0*T3-1.2222222222222222222D0*CA*T2)*T10+(-8.5D0
     &*T3+16.611111111111111111D0*CA*T2-2.2407407407407407407D0*T4*CF)*T
     &9+(27.0D0*T3-95.740740740740740741D0*CA*T2+35.419753086419753086D0
     &*T4*CF)*T8+(-63.375D0*T3+318.3904320987654321D0*CA*T2-230.64197530
     &864197531D0*T4*CF)*T6+(98.125D0*T3-575.20833333333333333D0*CA*T2+7
     &17.39026063100137174D0*T4*CF)*LL-211.58333333333333333D0*T3+640.47
     &067901234567901D0*CA*T2-973.22597546105776558D0*T4*CF


      T2=Nc*Nc
      R5=((-400*CF*T2+1600*CF)*zeta5+(70*CF*T2-280*CF)*zeta3+(-6*CF*T2+
     &24*CF)*zeta2*zeta2+(150*CF*T2-600*CF)*zeta2+60*CF*T2-240*CF)/(15*N
     &c)

      coeff(3) = R4
      NfV = R5

      return
      end
