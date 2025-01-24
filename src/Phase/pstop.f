!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine gen2new(r,p,pswt,*)
          use types
          use constants
          implicit none

          include 'mxdim.f'
          include 'mxpart.f'
          include 'x1x2.f'
          include 'phasemin.f'
          include 'energy.f'

          real(dp), intent(in) :: r(mxdim)
          real(dp), intent(out) :: p(mxpart,4)
          real(dp), intent(out) :: pswt

          real(dp) :: lntaum, tau, xjac
          real(dp) :: x1mx2, surd

          real(dp) :: puremass

          real(dp) :: m12, phi, t, sinth, costh
          real(dp) :: p3cm(4)

          real(dp) :: p12(4)

          lntaum=log(taumin)
          tau=exp(lntaum*(1._dp-r(1)))
          xjac=-lntaum*tau

          x1mx2=2._dp*r(2)-1._dp
          surd=sqrt(x1mx2**2+4._dp*tau)

          xx(1)=0.5_dp*(+x1mx2+surd)
          xx(2)=0.5_dp*(-x1mx2+surd)

          xjac=xjac*2._dp/surd

          if   ((xx(1) > 1._dp)
     &    .or. (xx(2) > 1._dp)) return 1

          p(1,4)=-xx(1)*sqrts*0.5_dp
          p(1,1)=0._dp
          p(1,2)=0._dp
          p(1,3)=-xx(1)*sqrts*0.5_dp

          p(2,4)=-xx(2)*sqrts*0.5_dp
          p(2,1)=0._dp
          p(2,2)=0._dp
          p(2,3)=+xx(2)*sqrts*0.5_dp


          p12(:) = -p(1,:)-p(2,:)

          ! we sample directly in t and phi and construct theta from t
          m12 = puremass(p12)
          phi = 2._dp*pi*r(3)
          t = -m12**2 * r(4)
          costh = 1 + 2*t/m12**2
          sinth = sqrt(max(1._dp-costh**2,0._dp))

          p3cm(4) = 1._dp
          p3cm(1) = sinth*sin(phi)
          p3cm(2) = sinth*cos(phi)
          p3cm(3) = costh
          p3cm(:) = p3cm(:) * m12/2

          call boost(m12,p12,p3cm,p(3,:))
          p(4,:) = p12(:) - p(3,:)

          pswt = xjac*(1._dp/8._dp/pi)
          return


          pswt = 0._dp
          return 1

      end subroutine

      subroutine pstop(r,p,pswt,*)
          use types
          use constants
          implicit none

          include 'mxdim.f'
          include 'mxpart.f'
          include 'x1x2.f'
          include 'phasemin.f'
          include 'energy.f'
          include 'masses.f'

          real(dp), intent(in) :: r(mxdim)
          real(dp), intent(out) :: p(mxpart,4)
          real(dp), intent(out) :: pswt

          real(dp) :: wt0, wt1, wt2, wt3
          real(dp) :: lntaum, tau, xjac
          real(dp) :: x1mx2, surd

          real(dp) :: p12(4), p6(4), p126(4), p5(4), p34(4), p3(4), p4(4)

          lntaum=log(taumin)
          tau=exp(lntaum*(1._dp-r(1)))
          xjac=-lntaum*tau

          x1mx2=2._dp*r(2)-1._dp
          surd=sqrt(x1mx2**2+4._dp*tau)

          xx(1)=0.5_dp*(+x1mx2+surd)
          xx(2)=0.5_dp*(-x1mx2+surd)

          xjac=xjac*2._dp/surd

          if   ((xx(1) > 1._dp)
     &    .or. (xx(2) > 1._dp)) return 1


          p(1,4)=-xx(1)*sqrts*0.5_dp
          p(1,1)=0._dp
          p(1,2)=0._dp
          p(1,3)=-xx(1)*sqrts*0.5_dp

          p(2,4)=-xx(2)*sqrts*0.5_dp
          p(2,1)=0._dp
          p(2,2)=0._dp
          p(2,3)=+xx(2)*sqrts*0.5_dp

          wt0 = 1._dp/(2._dp*pi)**2


          p12(:) = -p(1,:)-p(2,:)

          call phi1_2m_bw(0._dp, r(3), r(4), r(5), 0d0, p12, p6, p126, mt, twidth, wt1, *99)
          call phi1_2m_bw(0._dp, r(6), r(7), r(8), 0d0, p126, p5, p34, wmass, wwidth, wt2, *99)
          call phi3m0(r(9), r(10), p34, p3, p4, wt3, *99)

          p(3,:) = p3(:)
          p(4,:) = p4(:)
          p(5,:) = p5(:)
          p(6,:) = p6(:)

          pswt = xjac*wt0*wt1*wt2*wt3
          return


  99      pswt = 0._dp
          return 1

      end subroutine

      subroutine pstopReal(r,p,pswt,*)
          use types
          use constants
          implicit none

          include 'mxdim.f'
          include 'mxpart.f'
          include 'x1x2.f'
          include 'phasemin.f'
          include 'energy.f'
          include 'masses.f'
          include 'ipsgen.f'

          real(dp), intent(in) :: r(mxdim)
          real(dp), intent(out) :: p(mxpart,4)
          real(dp), intent(out) :: pswt

          real(dp) :: wt0, wt1, wt2, wt3, wt4
          real(dp) :: lntaum, tau, xjac
          real(dp) :: x1mx2, surd

          real(dp) :: p12(4), p6(4), p5(4), p34(4), p3(4), p4(4)
          real(dp) :: p345(4), p7(4), p67(4)

          lntaum=log(taumin)
          tau=exp(lntaum*(1._dp-r(1)))
          xjac=-lntaum*tau

          x1mx2=2._dp*r(2)-1._dp
          surd=sqrt(x1mx2**2+4._dp*tau)

          xx(1)=0.5_dp*(+x1mx2+surd)
          xx(2)=0.5_dp*(-x1mx2+surd)

          xjac=xjac*2._dp/surd

          if   ((xx(1) > 1._dp)
     &    .or. (xx(2) > 1._dp)) return 1


          p(1,4)=-xx(1)*sqrts*0.5_dp
          p(1,1)=0._dp
          p(1,2)=0._dp
          p(1,3)=-xx(1)*sqrts*0.5_dp

          p(2,4)=-xx(2)*sqrts*0.5_dp
          p(2,1)=0._dp
          p(2,2)=0._dp
          p(2,3)=+xx(2)*sqrts*0.5_dp

          wt0 = 1._dp/(2._dp*pi)**3


          p12(:) = -p(1,:)-p(2,:)

          select case (ipsgen)
          case (1)
              call phi1_2bw(0.5d0, r(3), r(4), r(5), p12, p345, p67, mt, twidth, wt1, *99)
              call phi3m0(r(6),r(7), p67, p6, p7, wt2, *99)
              call phi1_2m_bw(0d0, 0.5d0,r(9),r(10), 1d-8, p345, p5, p34, wmass, wwidth, wt3, *99)
              call phi3m0(r(11), r(12), p34, p3, p4, wt4, *99)
          case default
              write(6,*) 'Abort in pstop'
              stop
          end select


          p(3,:) = p3(:)
          p(4,:) = p4(:)
          p(5,:) = p5(:)
          p(6,:) = p6(:)
          p(7,:) = p7(:)

          pswt = xjac*wt0*wt1*wt2*wt3*wt4
          return

  99      pswt = 0._dp
          return 1

      end subroutine

      subroutine pstopReal2(r,p,pswt,*)
          use types
          use constants
          implicit none

          include 'mxdim.f'
          include 'mxpart.f'
          include 'x1x2.f'
          include 'phasemin.f'
          include 'energy.f'
          include 'masses.f'
          include 'ipsgen.f'

          real(dp), intent(in) :: r(mxdim)
          real(dp), intent(out) :: p(mxpart,4)
          real(dp), intent(out) :: pswt

          real(dp) :: wt0, wt1, wt2, wt3, wt4, wt5
          real(dp) :: lntaum, tau, xjac
          real(dp) :: x1mx2, surd

          real(dp) :: p12(4), p6(4), p5(4), p34(4), p3(4), p4(4)
          real(dp) :: p345(4), p7(4), p67(4), p678(4), p8(4)

          lntaum=log(taumin)
          tau=exp(lntaum*(1._dp-r(1)))
          xjac=-lntaum*tau

          x1mx2=2._dp*r(2)-1._dp
          surd=sqrt(x1mx2**2+4._dp*tau)

          xx(1)=0.5_dp*(+x1mx2+surd)
          xx(2)=0.5_dp*(-x1mx2+surd)

          xjac=xjac*2._dp/surd

          if   ((xx(1) > 1._dp)
     &    .or. (xx(2) > 1._dp)) return 1


          p(1,4)=-xx(1)*sqrts*0.5_dp
          p(1,1)=0._dp
          p(1,2)=0._dp
          p(1,3)=-xx(1)*sqrts*0.5_dp

          p(2,4)=-xx(2)*sqrts*0.5_dp
          p(2,1)=0._dp
          p(2,2)=0._dp
          p(2,3)=+xx(2)*sqrts*0.5_dp

          wt0 = 1._dp/(2._dp*pi)**4

          p12(:) = -p(1,:)-p(2,:)

          select case (ipsgen)
          case (1)
              call phi1_2bw(0.5d0, r(3), r(4), r(5), p12, p345, p678, mt, twidth, wt1, *99)
              call phi1_2m_nobw(0d0, r(6),r(7),r(8), 1d-8, p678, p8, p67, wt2, *99)
              call phi3m0(r(9),r(10), p67, p6, p7, wt3, *99)
              call phi1_2m_bw(0d0, 0.5d0,r(11),r(12), 1d-8, p345, p5, p34, wmass, wwidth, wt4, *99)
              call phi3m0(r(13), r(14), p34, p3, p4, wt5, *99)
          case default
              write(6,*) 'Abort in pstop'
              stop
          end select


          p(3,:) = p3(:)
          p(4,:) = p4(:)
          p(5,:) = p5(:)
          p(6,:) = p6(:)
          p(7,:) = p7(:)
          p(8,:) = p8(:)

          pswt = xjac*wt0*wt1*wt2*wt3*wt4*wt5
          return

  99      pswt = 0._dp
          return 1

      end subroutine

