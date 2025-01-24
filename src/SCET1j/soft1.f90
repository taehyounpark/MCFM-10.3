!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

    subroutine computeIijm(y12,y23,y31,Iijm)
    implicit none
!--- One-time computation of integrals required for one-jet soft functions
    include 'src/Inc/types.f'
    include 'src/Inc/constants.f'
    real(dp)::y12,y23,y31,y(3,3),Iijm(5),al,be,I0JSTW,I1JSTW,I0,I1
    integer::n
    integer,parameter::i(6)=(/1,2,2,3,3,1/)
    integer,parameter::j(6)=(/2,1,3,2,1,3/)
    integer,parameter::m(6)=(/3,3,1,1,2,2/)

    y(:,:)=zip

    y(1,2)=y12
    y(2,1)=y12
    y(2,3)=y23
    y(3,2)=y23
    y(1,3)=y31
    y(3,1)=y31

    do n=1,5
        al=y(j(n),m(n))/y(i(n),j(n))
        be=y(i(n),m(n))/y(i(n),j(n))
        I0=I0JSTW(al,be)
        I1=I1JSTW(al,be)
        Iijm(n)=I0*log(al)+I1
    enddo

    return
    end subroutine computeIijm
          
          
    function I0JSTW(al,be)
! function I0 from 1102.4344
    implicit none
    include 'src/Inc/types.f'
    include 'src/Inc/constants.f'
    real(dp)::I0JSTW
    real(dp)::al,be,rtal,rtbe,phimax,phicut,dgauss,eps
    common/albe/rtal,rtbe
!$omp threadprivate(/albe/)

    real(dp), external :: I0integrand1,  I0integrand2
    eps=1.e-3_dp
    rtal=sqrt(al)
    rtbe=sqrt(be)
!--setup limits
    if (abs(rtal-rtbe) >= one) then
        phicut=zip
    elseif (rtal+rtbe <= one) then
        phicut=pi
    else
        phicut=acos((al+be-one)/(two*rtal*rtbe))
    endif
!--perform integration
    if (phicut > zip) then
        I0JSTW=dgauss(I0integrand1,zip,phicut,eps)
    else
        I0JSTW=zip
    endif
    if (al-be-one > zip) then
        phimax=asin(one/rtal)
        I0JSTW=I0JSTW &
        +dgauss(I0integrand2,phicut,phimax,eps)
    endif
    I0JSTW=two/pi*I0JSTW
    return
    end function I0JSTW

    function I0integrand1(phi)
    implicit none
    include 'src/Inc/types.f'
    include 'src/Inc/constants.f'
    real(dp)::I0integrand1
    real(dp)::rtal,rtbe,phi,yp,rt
    common/albe/rtal,rtbe
!$omp threadprivate(/albe/)

    rt=one/rtal**2-sin(phi)**2
    if (rt > 0) then
        rt=sqrt(rt)
    else
        rt=zip
    endif
    yp=cos(phi)+rt
    I0integrand1=log(yp/(rtbe/rtal))
    return
    end function I0integrand1

    function I0integrand2(phi)
    implicit none
    include 'src/Inc/types.f'
    include 'src/Inc/constants.f'
    real(dp)::I0integrand2
    real(dp)::rtal,rtbe,phi,yp,ym,rt
    common/albe/rtal,rtbe
!$omp threadprivate(/albe/)
    rt=one/rtal**2-sin(phi)**2
    if (rt > 0) then
        rt=sqrt(rt)
    else
        rt=zip
    endif
    yp=cos(phi)+rt
    ym=cos(phi)-rt
    I0integrand2=log(yp/ym)
    return
    end function I0integrand2


    function I1JSTW(al,be)
! function I1JSTW from 1102.4344
    implicit none
    include 'src/Inc/types.f'
    include 'src/Inc/constants.f'
    real(dp)::I1JSTW
    real(dp)::al,be,rtal,rtbe,phimax,phicut,dgauss,eps
    common/albe/rtal,rtbe
!$omp threadprivate(/albe/)

    real(dp), external :: I1integrand1, I1integrand2

    eps=1.e-3_dp
    rtal=sqrt(al)
    rtbe=sqrt(be)
!--setup limits
    if (abs(rtal-rtbe) >= one) then
        phicut=zip
    elseif (rtal+rtbe <= one) then
        phicut=pi
    else
        phicut=acos((al+be-one)/(two*rtal*rtbe))
    endif
!--perform integration
    if (phicut > zip) then
        I1JSTW=dgauss(I1integrand1,zip,phicut,eps)
    else
        I1JSTW=zip
    endif

    if (al-be-one > zip) then
        phimax=asin(one/rtal)
        I1JSTW=I1JSTW &
        +dgauss(I1integrand2,phicut,phimax,eps)
    endif
    I1JSTW=two/pi*I1JSTW
    return
    end function I1JSTW


    function I1integrand1(phi)
    implicit none
    include 'src/Inc/types.f'
    include 'src/Inc/constants.f'
    real(dp)::I1integrand1
    real(dp)::rtal,rtbe,phi,y,yp,G,rt
    complex(dp)::XSPENZ
    common/albe/rtal,rtbe
!$omp threadprivate(/albe/)
!---statement function
    G(y,phi)=-two*real(XSPENZ(y*exp(im*phi)))
!---end statement function
    rt=one/rtal**2-sin(phi)**2
    if (rt > 0) then
        rt=sqrt(rt)
    else
        rt=zip
    endif
    yp=cos(phi)+rt
    I1integrand1=G(yp,phi)-G(rtbe/rtal,phi)
    return
    end function I1integrand1

    function I1integrand2(phi)
    implicit none
    include 'src/Inc/types.f'
    include 'src/Inc/constants.f'
    real(dp)::I1integrand2
    real(dp)::rtal,rtbe,phi,y,yp,ym,G,rt
    complex(dp)::XSPENZ
    common/albe/rtal,rtbe
!$omp threadprivate(/albe/)
!---statement function
    G(y,phi)=-two*real(XSPENZ(y*exp(im*phi)))
!---end statement function
    rt=one/rtal**2-sin(phi)**2
    if (rt > 0) then
        rt=sqrt(rt)
    else
        rt=zip
    endif
    yp=cos(phi)+rt
    ym=cos(phi)-rt
    I1integrand2=G(yp,phi)-G(ym,phi)
    return
    end function I1integrand2


    subroutine soft_ab_ggg(order,y12,y23,y31, &
    preIijm,j1,j2,j3,j4,j5,j6,soft1,soft2)
    implicit none
!     Abelian piece g(1) g(2) g(3)
!     returns coefficients of [log(tau)^j/tau_+],and for j=-1, delta(tau)
!     in units of as/2/pi
!     y12=n1.n2/2 so that 0<y12<1
!     y23=n2.n3/2 so that 0<y23<1
!     y31=n3.n1/2 so that 0<y31<1
!     order indicates the order in [as/2/pi]
    include 'src/Inc/types.f'
    include 'src/Inc/constants.f'
    include 'src/Inc/nf.f'
    include 'src/Inc/scet_const.f'
    real(dp)::y12,y23,y31,soft1(-1:1),soft2(-1:3),preIijm(5), &
    Iijm(6)
    integer::order,j1,j2,j3,j4,j5,j6
    include 'src/Inc/scet_beta.f'
!      integer,parameter::i(6)=(/1,2,2,3,3,1/)
!      integer,parameter::j(6)=(/2,1,3,2,1,3/)
!      integer,parameter::m(6)=(/3,3,1,1,2,2/)
! d,I12x3=Iijm(1);
! d,I21x3=Iijm(2);
! d,I23x1=Iijm(3);
! d,I32x1=Iijm(4);
! d,I31x2=Iijm(5);
! d,I13x2=Iijm(6);

    soft1(:)=zip
    soft2(:)=zip
    if (order < 1) return

!      y(:,:)=zip

!      y(1,2)=y12
!      y(2,1)=y12
!      y(2,3)=y23
!      y(3,2)=y23
!      y(1,3)=y31
!      y(3,1)=y31

!      do n=1,6
!      al=y(j(n),m(n))/y(i(n),j(n))
!      be=y(i(n),m(n))/y(i(n),j(n))
!      I0=I0JSTW(al,be)
!      I1=I1JSTW(al,be)
!      Iijm(n)=I0*log(al)+I1
!      enddo

    Iijm(1)=preIijm(j1)
    Iijm(2)=preIijm(j2)
    Iijm(3)=preIijm(j3)
    Iijm(4)=preIijm(j4)
    Iijm(5)=preIijm(j5)
    Iijm(6)=preIijm(j6)

    soft1(-1)=CA*(three*half*zeta2 &
    -half*log(y12)**2-half*log(y31)**2-half*log(y23)**2 &
    -Iijm(1)-Iijm(2)-Iijm(3)-Iijm(4)-Iijm(5)-Iijm(6))
    soft1(0)=2*CA*(log(y12)+log(y23)+log(y31))
    soft1(1)=-12*CA

    if (order < 2) return
          
!-------------------
    soft2(-1) = CA**2*( &
    - 243._dp/40._dp*zeta2**2 &
    - 24._dp*log(y12)*zeta3 &
    - 11/4._dp*log(y12)**2*zeta2 &
    + 1._dp/8._dp*log(y12)**4 &
    + quarter*log(y12)**2*log(y31)**2 &
    + quarter*log(y12)**2*log(y23)**2 &
    + half*log(y12)**2*Iijm(1) &
    + half*log(y12)**2*Iijm(2) &
    + half*log(y12)**2*Iijm(3) &
    + half*log(y12)**2*Iijm(4) &
    + half*log(y12)**2*Iijm(5) &
    + half*log(y12)**2*Iijm(6) &
    - 4._dp*log(y12)*log(y31)*zeta2 &
    - 4._dp*log(y12)*log(y23)*zeta2 &
    )
    soft2(-1) = soft2(-1) + CA**2*( &
    - 24._dp*log(y31)*zeta3 &
    - 11/4._dp*log(y31)**2*zeta2 &
    + 1._dp/8._dp*log(y31)**4 &
    + quarter*log(y31)**2*log(y23)**2 &
    + half*log(y31)**2*Iijm(1) &
    + half*log(y31)**2*Iijm(2) &
    + half*log(y31)**2*Iijm(3) &
    + half*log(y31)**2*Iijm(4) &
    + half*log(y31)**2*Iijm(5) &
    + half*log(y31)**2*Iijm(6) &
    - 4._dp*log(y31)*log(y23)*zeta2 &
    - 24._dp*log(y23)*zeta3 &
    - 11/4._dp*log(y23)**2*zeta2 &
    + 1._dp/8._dp*log(y23)**4 &
    )
    soft2(-1) = soft2(-1) + CA**2*( &
    + half*log(y23)**2*Iijm(1) &
    + half*log(y23)**2*Iijm(2) &
    + half*log(y23)**2*Iijm(3) &
    + half*log(y23)**2*Iijm(4) &
    + half*log(y23)**2*Iijm(5) &
    + half*log(y23)**2*Iijm(6) &
    - 3._dp/2._dp*Iijm(1)*zeta2 &
    + half*Iijm(1)**2 &
    + Iijm(1)*Iijm(2) &
    + Iijm(1)*Iijm(3) &
    + Iijm(1)*Iijm(4) &
    + Iijm(1)*Iijm(5) &
    + Iijm(1)*Iijm(6) &
    - 3._dp/2._dp*Iijm(2)*zeta2 &
    )
    soft2(-1) = soft2(-1) + CA**2*( &
    + half*Iijm(2)**2 &
    + Iijm(2)*Iijm(3) &
    + Iijm(2)*Iijm(4) &
    + Iijm(2)*Iijm(5) &
    + Iijm(2)*Iijm(6) &
    - 3._dp/2._dp*Iijm(3)*zeta2 &
    + half*Iijm(3)**2 &
    + Iijm(3)*Iijm(4) &
    + Iijm(3)*Iijm(5) &
    + Iijm(3)*Iijm(6) &
    - 3._dp/2._dp*Iijm(4)*zeta2 &
    + half*Iijm(4)**2 &
    + Iijm(4)*Iijm(5) &
    + Iijm(4)*Iijm(6) &
    )
    soft2(-1) = soft2(-1) + CA**2*( &
    - 3._dp/2._dp*Iijm(5)*zeta2 &
    + half*Iijm(5)**2 &
    + Iijm(5)*Iijm(6) &
    - 3._dp/2._dp*Iijm(6)*zeta2 &
    + half*Iijm(6)**2 &
    )
    soft2(0) = soft2(0) + CA**2*( &
    + 144._dp*zeta3 &
    + 27._dp*log(y12)*zeta2 &
    - log(y12)**3 &
    - log(y12)**2*log(y31) &
    - log(y12)**2*log(y23) &
    - log(y12)*log(y31)**2 &
    - log(y12)*log(y23)**2 &
    - 2._dp*log(y12)*Iijm(1) &
    - 2._dp*log(y12)*Iijm(2) &
    - 2._dp*log(y12)*Iijm(3) &
    - 2._dp*log(y12)*Iijm(4) &
    - 2._dp*log(y12)*Iijm(5) &
    - 2._dp*log(y12)*Iijm(6) &
    + 27._dp*log(y31)*zeta2 &
    )
    soft2(0) = soft2(0) + CA**2*( &
    - log(y31)**3 &
    - log(y31)**2*log(y23) &
    - log(y31)*log(y23)**2 &
    - 2._dp*log(y31)*Iijm(1) &
    - 2._dp*log(y31)*Iijm(2) &
    - 2._dp*log(y31)*Iijm(3) &
    - 2._dp*log(y31)*Iijm(4) &
    - 2._dp*log(y31)*Iijm(5) &
    - 2._dp*log(y31)*Iijm(6) &
    + 27._dp*log(y23)*zeta2 &
    - log(y23)**3 &
    - 2._dp*log(y23)*Iijm(1) &
    - 2._dp*log(y23)*Iijm(2) &
    - 2._dp*log(y23)*Iijm(3) &
    )
    soft2(0) = soft2(0)+CA**2*( &
    - 2._dp*log(y23)*Iijm(4) &
    - 2._dp*log(y23)*Iijm(5) &
    - 2._dp*log(y23)*Iijm(6) &
    )
    soft2(1) = CA**2*( &
    - 162._dp*zeta2 &
    + 10._dp*log(y12)**2 &
    + 8._dp*log(y12)*log(y31) &
    + 8._dp*log(y12)*log(y23) &
    + 10._dp*log(y31)**2 &
    + 8._dp*log(y31)*log(y23) &
    + 10._dp*log(y23)**2 &
    + 12._dp*Iijm(1) &
    + 12._dp*Iijm(2) &
    + 12._dp*Iijm(3) &
    + 12._dp*Iijm(4) &
    + 12._dp*Iijm(5) &
    + 12._dp*Iijm(6) &
    )
    soft2(2) = CA**2*( &
    - 36._dp*log(y12) &
    - 36._dp*log(y31) &
    - 36._dp*log(y23))
    soft2(3) =72._dp*CA**2
    return
    end subroutine soft_ab_ggg
    subroutine soft_ab_qag(order,y12,y23,y31, &
    preIijm,j1,j2,j3,j4,j5,j6,soft1,soft2)
    implicit none
!---   Abelian piece q(1) qbar(2) g(3)
!--- returns coefficients of [log(tau)^j/tau_+],and for j=-1, delta(tau)
!--- in units of as/2/pi
!     y12=n1.n2/2 so that 0<y12<1
!     y23=n2.n3/2 so that 0<y23<1
!     y31=n3.n1/2 so that 0<y31<1
!     order indicates the order in [as/2/pi]
    include 'src/Inc/types.f'
    include 'src/Inc/constants.f'
    include 'src/Inc/nf.f'
    include 'src/Inc/scet_const.f'
    real(dp)::y12,y23,y31,soft1(-1:1),soft2(-1:3),preIijm(5), &
    Iijm(6),Lss,Lsd,Isum12,Isum3456
    integer::order,j1,j2,j3,j4,j5,j6
    include 'src/Inc/scet_beta.f'
!      integer,parameter::i(6)=(/1,2,2,3,3,1/)
!      integer,parameter::j(6)=(/2,1,3,2,1,3/)
!      integer,parameter::m(6)=(/3,3,1,1,2,2/)
! d,I12x3=Iijm(1);
! d,I21x3=Iijm(2);
! d,I23x1=Iijm(3);
! d,I32x1=Iijm(4);
! d,I31x2=Iijm(5);
! d,I13x2=Iijm(6);

    soft1(:)=zip
    soft2(:)=zip
    if (order < 1) return

!      y(:,:)=zip

!      y(1,2)=y12
!      y(2,1)=y12
!      y(2,3)=y23
!      y(3,2)=y23
!      y(1,3)=y31
!      y(3,1)=y31

!      if ((y12 < zip) .or. (y12 > one)) then
!        write(6,*) 'error: y12 outside range, y12 = ',y12
!      endif
!      if ((y23 < zip) .or. (y23 > one)) then
!        write(6,*) 'error: y23 outside range, y23 = ',y23
!      endif
!      if ((y31 < zip) .or. (y31 > one)) then
!        write(6,*) 'error: y31 outside range, y31 = ',y31
!      endif
!      do n=1,6
!      al=y(j(n),m(n))/y(i(n),j(n))
!      be=y(i(n),m(n))/y(i(n),j(n))
!      I0=I0JSTW(al,be)
!      I1=I1JSTW(al,be)
!      Iijm(n)=I0*log(al)+I1
!      enddo

    Iijm(1)=preIijm(j1)
    Iijm(2)=preIijm(j2)
    Iijm(3)=preIijm(j3)
    Iijm(4)=preIijm(j4)
    Iijm(5)=preIijm(j5)
    Iijm(6)=preIijm(j6)

    Lss=log(y31)+log(y23)
    Lsd=log(y31)-log(y23)
    Isum12=Iijm(1)+Iijm(2)
    Isum3456=Iijm(3)+Iijm(4)+Iijm(5)+Iijm(6)

    soft1(1)=-four*CA-eight*CF
    soft1(0)=two*CA*(Lss-log(y12))+four*log(y12)*CF

    soft1(-1)=CA*(-0.25_dp*(Lsd**2+Lss**2)+Isum12-Isum3456 &
    +half*zeta2+half*log(y12)**2) &
    +CF*(-two*Isum12+zeta2-log(y12)**2)


    if (order < 2) return

    soft2(3)=8._dp*(CA+2*CF)**2

    soft2(2)=-12._dp*CA**2*Lss-24._dp*CF*CA*Lss

    soft2(1) = + CA**2*Lsd**2 &
    + 5._dp*CA**2*Lss**2 &
    - 4._dp*CA**2*Isum12 &
    + 4._dp*CA**2*Isum3456 &
    - 18._dp*CA**2*zeta2 &
    + 2._dp*CF*CA*Lsd**2 &
    + 2._dp*CF*CA*Lss**2 &
    + 8._dp*CF*CA*Isum3456 &
    - 72._dp*CF*CA*zeta2 &
    + 16._dp*CF**2*Isum12 &
    - 72._dp*CF**2*zeta2

    soft2(0) = - 1._dp/2._dp*CA**2*Lss*Lsd**2 &
    - 1._dp/2._dp*CA**2*Lss**3 &
    + 2._dp*CA**2*Isum12*Lss &
    - 2._dp*CA**2*Isum3456*Lss &
    + 16._dp*CA**2*zeta3 &
    + 9._dp*CA**2*zeta2*Lss &
    - 4._dp*CF*CA*Isum12*Lss &
    + 64._dp*CF*CA*zeta3 &
    + 18._dp*CF*CA*zeta2*Lss &
    + 64._dp*CF**2*zeta3

    soft2(-1) = + 1._dp/32._dp*CA**2*Lsd**4 &
    + 1._dp/16._dp*CA**2*Lss**2*Lsd**2 &
    + 1._dp/32._dp*CA**2*Lss**4 &
    - 1._dp/4._dp*CA**2*Isum12*Lsd**2 &
    - 1._dp/4._dp*CA**2*Isum12*Lss**2 &
    + 1._dp/2._dp*CA**2*Isum12**2 &
    + 1._dp/4._dp*CA**2*Isum3456*Lsd**2 &
    + 1._dp/4._dp*CA**2*Isum3456*Lss**2 &
    - CA**2*Isum3456*Isum12 &
    + 1._dp/2._dp*CA**2*Isum3456**2 &
    - 8._dp*CA**2*zeta3*Lss &
    - 1._dp/8._dp*CA**2*zeta2*Lsd**2 &
    - 17._dp/8._dp*CA**2*zeta2*Lss**2 &
    + 1._dp/2._dp*CA**2*zeta2*Isum12 &
    - 1._dp/2._dp*CA**2*zeta2*Isum3456
    soft2(-1) = soft2(-1) - 27._dp/40._dp*CA**2*zeta2**2 &
    + 1._dp/2._dp*CF*CA*Isum12*Lsd**2 &
    + 1._dp/2._dp*CF*CA*Isum12*Lss**2 &
    - 2._dp*CF*CA*Isum12**2 &
    + 2._dp*CF*CA*Isum3456*Isum12 &
    - 16._dp*CF*CA*zeta3*Lss &
    - 1._dp/4._dp*CF*CA*zeta2*Lsd**2 &
    - 1._dp/4._dp*CF*CA*zeta2*Lss**2 &
    - CF*CA*zeta2*Isum3456 &
    - 27._dp/10._dp*CF*CA*zeta2**2 &
    + 2._dp*CF**2*Isum12**2 &
    - 2._dp*CF**2*zeta2*Isum12 &
    - 27._dp/10._dp*CF**2*zeta2**2

    return
    end subroutine soft_ab_qag
    subroutine soft_ab_qgq(order,y12,y23,y31, &
    preIijm,j1,j2,j3,j4,j5,j6,soft1,soft2)
    implicit none
!     Abelian piece q(1) g(2) q(3)
!     returns coefficients of [log(tau)^j/tau_+],and for j=-1, delta(tau)
!     in units of as/2/pi
!     y12=n1.n2/2 so that 0<y12<1
!     y23=n2.n3/2 so that 0<y23<1
!     y31=n3.n1/2 so that 0<y31<1
!     order indicates the order in [as/2/pi]
    include 'src/Inc/types.f'
    include 'src/Inc/constants.f'
    include 'src/Inc/nf.f'
    include 'src/Inc/scet_const.f'
    real(dp)::y12,y23,y31,soft1(-1:1),soft2(-1:3),preIijm(5), &
    Iijm(6)
    integer::order,j1,j2,j3,j4,j5,j6
    include 'src/Inc/scet_beta.f'
!      integer,parameter::i(6)=(/1,2,2,3,3,1/)
!      integer,parameter::j(6)=(/2,1,3,2,1,3/)
!      integer,parameter::m(6)=(/3,3,1,1,2,2/)
! d,I12x3=Iijm(1);
! d,I21x3=Iijm(2);
! d,I23x1=Iijm(3);
! d,I32x1=Iijm(4);
! d,I31x2=Iijm(5);
! d,I13x2=Iijm(6);

    soft1(:)=zip
    soft2(:)=zip
    if (order < 1) return

!      y(:,:)=zip
!      y(1,2)=y12
!      y(2,1)=y12
!      y(2,3)=y23
!      y(3,2)=y23
!      y(1,3)=y31
!      y(3,1)=y31

!      if ((y12 < zip) .or. (y12 > one)) then
!        write(6,*) 'error: y12 outside range, y12 = ',y12
!      endif
!      if ((y23 < zip) .or. (y23 > one)) then
!        write(6,*) 'error: y23 outside range, y23 = ',y23
!      endif
!      if ((y31 < zip) .or. (y31 > one)) then
!        write(6,*) 'error: y31 outside range, y31 = ',y31
!      endif

!      do n=1,6
!      al=y(j(n),m(n))/y(i(n),j(n))
!      be=y(i(n),m(n))/y(i(n),j(n))
!      I0=I0JSTW(al,be)
!      I1=I1JSTW(al,be)
!      Iijm(n)=I0*log(al)+I1
!      enddo

    Iijm(1)=preIijm(j1)
    Iijm(2)=preIijm(j2)
    Iijm(3)=preIijm(j3)
    Iijm(4)=preIijm(j4)
    Iijm(5)=preIijm(j5)
    Iijm(6)=preIijm(j6)

    soft1(-1)=CA*(half*zeta2 &
    -half*log(y12)**2+half*log(y31)**2-half*log(y23)**2 &
    -Iijm(1)-Iijm(2)-Iijm(3)-Iijm(4)+Iijm(5)+Iijm(6)) &
    +CF*(zeta2-log(y31)**2-2*Iijm(5)-2*Iijm(6))
    soft1(0)=2*CA*(log(y12)+log(y23)-log(y31))+4*log(y31)*CF
    soft1(1)=-4*CA-8*CF
          
    if (order < 2) return
          
!-------------------
    soft2(3)=8*(CA+2*CF)**2
    soft2(2)=-12*(CA+2*CF)*(CA*(log(y12)+log(y23)-log(y31)) &
    +2*log(y31)*CF)

    soft2(1) =  + CA**2 * ( &
    - 18._dp*zeta2 &
    + 6._dp*log(y12)**2 &
    - 8._dp*log(y12)*log(y31) &
    + 8._dp*log(y12)*log(y23) &
    + 2._dp*log(y31)**2 &
    - 8._dp*log(y31)*log(y23) &
    + 6._dp*log(y23)**2 &
    + 4._dp*Iijm(1) &
    + 4._dp*Iijm(2) &
    + 4._dp*Iijm(3) &
    + 4._dp*Iijm(4) &
    - 4._dp*Iijm(5) &
    - 4._dp*Iijm(6) &
    )
    soft2(1) = soft2(1) + CF*CA * ( &
    - 72._dp*zeta2 &
    + 4._dp*log(y12)**2 &
    + 16._dp*log(y12)*log(y31) &
    - 16._dp*log(y31)**2 &
    + 16._dp*log(y31)*log(y23) &
    + 4._dp*log(y23)**2 &
    + 8._dp*Iijm(1) &
    + 8._dp*Iijm(2) &
    + 8._dp*Iijm(3) &
    + 8._dp*Iijm(4) &
    )
    soft2(1) = soft2(1) + CF**2 * ( &
    - 72._dp*zeta2 &
    + 24._dp*log(y31)**2 &
    + 16._dp*Iijm(5) &
    + 16._dp*Iijm(6) &
    )

    soft2(0) =  + CA**2 * ( &
    + 16._dp*zeta3 &
    + 9._dp*log(y12)*zeta2 &
    - log(y12)**3 &
    + log(y12)**2*log(y31) &
    - log(y12)**2*log(y23) &
    + log(y12)*log(y31)**2 &
    - log(y12)*log(y23)**2 &
    - 2._dp*log(y12)*Iijm(1) &
    - 2._dp*log(y12)*Iijm(2) &
    - 2._dp*log(y12)*Iijm(3) &
    - 2._dp*log(y12)*Iijm(4) &
    + 2._dp*log(y12)*Iijm(5) &
    + 2._dp*log(y12)*Iijm(6) &
    - 9._dp*log(y31)*zeta2 &
    )
    soft2(0) = soft2(0) + CA**2 * ( &
    - log(y31)**3 &
    + log(y31)**2*log(y23) &
    + log(y31)*log(y23)**2 &
    + 2._dp*log(y31)*Iijm(1) &
    + 2._dp*log(y31)*Iijm(2) &
    + 2._dp*log(y31)*Iijm(3) &
    + 2._dp*log(y31)*Iijm(4) &
    - 2._dp*log(y31)*Iijm(5) &
    - 2._dp*log(y31)*Iijm(6) &
    + 9._dp*log(y23)*zeta2 &
    - log(y23)**3 &
    - 2._dp*log(y23)*Iijm(1) &
    - 2._dp*log(y23)*Iijm(2) &
    - 2._dp*log(y23)*Iijm(3) &
    )
    soft2(0) = soft2(0) + CA**2 * ( &
    - 2._dp*log(y23)*Iijm(4) &
    + 2._dp*log(y23)*Iijm(5) &
    + 2._dp*log(y23)*Iijm(6) &
    )
    soft2(0) = soft2(0) + CF*CA * ( &
    + 64._dp*zeta3 &
    + 18._dp*log(y12)*zeta2 &
    - 2._dp*log(y12)**2*log(y31) &
    - 2._dp*log(y12)*log(y31)**2 &
    - 4._dp*log(y12)*Iijm(5) &
    - 4._dp*log(y12)*Iijm(6) &
    + 4._dp*log(y31)**3 &
    - 2._dp*log(y31)**2*log(y23) &
    - 2._dp*log(y31)*log(y23)**2 &
    - 4._dp*log(y31)*Iijm(1) &
    - 4._dp*log(y31)*Iijm(2) &
    - 4._dp*log(y31)*Iijm(3) &
    - 4._dp*log(y31)*Iijm(4) &
    + 8._dp*log(y31)*Iijm(5) &
    )
    soft2(0) = soft2(0) + CF*CA * ( &
    + 8._dp*log(y31)*Iijm(6) &
    + 18._dp*log(y23)*zeta2 &
    - 4._dp*log(y23)*Iijm(5) &
    - 4._dp*log(y23)*Iijm(6) &
    )
    soft2(0) = soft2(0) + CF**2 * ( &
    + 64._dp*zeta3 &
    + 36._dp*log(y31)*zeta2 &
    - 4._dp*log(y31)**3 &
    - 8._dp*log(y31)*Iijm(5) &
    - 8._dp*log(y31)*Iijm(6) &
    )

    soft2(-1) =  + CA**2 * ( &
    - 27._dp/40._dp*zeta2**2 &
    - 8._dp*log(y12)*zeta3 &
    - 9._dp/4._dp*log(y12)**2*zeta2 &
    + 1._dp/8._dp*log(y12)**4 &
    - 1._dp/4._dp*log(y12)**2*log(y31)**2 &
    + 1._dp/4._dp*log(y12)**2*log(y23)**2 &
    + 1._dp/2._dp*log(y12)**2*Iijm(1) &
    + 1._dp/2._dp*log(y12)**2*Iijm(2) &
    + 1._dp/2._dp*log(y12)**2*Iijm(3) &
    + 1._dp/2._dp*log(y12)**2*Iijm(4) &
    - 1._dp/2._dp*log(y12)**2*Iijm(5) &
    - 1._dp/2._dp*log(y12)**2*Iijm(6) &
    + 4._dp*log(y12)*log(y31)*zeta2 &
    - 4._dp*log(y12)*log(y23)*zeta2 &
    )
    soft2(-1) = soft2(-1) + CA**2 * ( &
    + 8._dp*log(y31)*zeta3 &
    - 7._dp/4._dp*log(y31)**2*zeta2 &
    + 1._dp/8._dp*log(y31)**4 &
    - 1._dp/4._dp*log(y31)**2*log(y23)**2 &
    - 1._dp/2._dp*log(y31)**2*Iijm(1) &
    - 1._dp/2._dp*log(y31)**2*Iijm(2) &
    - 1._dp/2._dp*log(y31)**2*Iijm(3) &
    - 1._dp/2._dp*log(y31)**2*Iijm(4) &
    + 1._dp/2._dp*log(y31)**2*Iijm(5) &
    + 1._dp/2._dp*log(y31)**2*Iijm(6) &
    + 4._dp*log(y31)*log(y23)*zeta2 &
    - 8._dp*log(y23)*zeta3 &
    - 9._dp/4._dp*log(y23)**2*zeta2 &
    + 1._dp/8._dp*log(y23)**4 &
    )
    soft2(-1) = soft2(-1) + CA**2 * ( &
    + 1._dp/2._dp*log(y23)**2*Iijm(1) &
    + 1._dp/2._dp*log(y23)**2*Iijm(2) &
    + 1._dp/2._dp*log(y23)**2*Iijm(3) &
    + 1._dp/2._dp*log(y23)**2*Iijm(4) &
    - 1._dp/2._dp*log(y23)**2*Iijm(5) &
    - 1._dp/2._dp*log(y23)**2*Iijm(6) &
    - 1._dp/2._dp*Iijm(1)*zeta2 &
    + 1._dp/2._dp*Iijm(1)**2 &
    + Iijm(1)*Iijm(2) &
    + Iijm(1)*Iijm(3) &
    + Iijm(1)*Iijm(4) &
    - Iijm(1)*Iijm(5) &
    - Iijm(1)*Iijm(6) &
    - 1._dp/2._dp*Iijm(2)*zeta2 &
    )
    soft2(-1) = soft2(-1) + CA**2 * ( &
    + 1._dp/2._dp*Iijm(2)**2 &
    + Iijm(2)*Iijm(3) &
    + Iijm(2)*Iijm(4) &
    - Iijm(2)*Iijm(5) &
    - Iijm(2)*Iijm(6) &
    - 1._dp/2._dp*Iijm(3)*zeta2 &
    + 1._dp/2._dp*Iijm(3)**2 &
    + Iijm(3)*Iijm(4) &
    - Iijm(3)*Iijm(5) &
    - Iijm(3)*Iijm(6) &
    - 1._dp/2._dp*Iijm(4)*zeta2 &
    + 1._dp/2._dp*Iijm(4)**2 &
    - Iijm(4)*Iijm(5) &
    - Iijm(4)*Iijm(6) &
    )
    soft2(-1) = soft2(-1) + CA**2 * ( &
    + 1._dp/2._dp*Iijm(5)*zeta2 &
    + 1._dp/2._dp*Iijm(5)**2 &
    + Iijm(5)*Iijm(6) &
    + 1._dp/2._dp*Iijm(6)*zeta2 &
    + 1._dp/2._dp*Iijm(6)**2 &
    )
    soft2(-1) = soft2(-1) + CF*CA * ( &
    - 27._dp/10._dp*zeta2**2 &
    - 16._dp*log(y12)*zeta3 &
    - 1._dp/2._dp*log(y12)**2*zeta2 &
    + 1._dp/2._dp*log(y12)**2*log(y31)**2 &
    + log(y12)**2*Iijm(5) &
    + log(y12)**2*Iijm(6) &
    - 8._dp*log(y12)*log(y31)*zeta2 &
    + 8._dp*log(y31)**2*zeta2 &
    - 1._dp/2._dp*log(y31)**4 &
    + 1._dp/2._dp*log(y31)**2*log(y23)**2 &
    + log(y31)**2*Iijm(1) &
    + log(y31)**2*Iijm(2) &
    + log(y31)**2*Iijm(3) &
    + log(y31)**2*Iijm(4) &
    )
    soft2(-1) = soft2(-1) + CF*CA * ( &
    - 2._dp*log(y31)**2*Iijm(5) &
    - 2._dp*log(y31)**2*Iijm(6) &
    - 8._dp*log(y31)*log(y23)*zeta2 &
    - 16._dp*log(y23)*zeta3 &
    - 1._dp/2._dp*log(y23)**2*zeta2 &
    + log(y23)**2*Iijm(5) &
    + log(y23)**2*Iijm(6) &
    - Iijm(1)*zeta2 &
    + 2._dp*Iijm(1)*Iijm(5) &
    + 2._dp*Iijm(1)*Iijm(6) &
    - Iijm(2)*zeta2 &
    + 2._dp*Iijm(2)*Iijm(5) &
    + 2._dp*Iijm(2)*Iijm(6) &
    - Iijm(3)*zeta2 &
    )
    soft2(-1) = soft2(-1) + CF*CA * ( &
    + 2._dp*Iijm(3)*Iijm(5) &
    + 2._dp*Iijm(3)*Iijm(6) &
    - Iijm(4)*zeta2 &
    + 2._dp*Iijm(4)*Iijm(5) &
    + 2._dp*Iijm(4)*Iijm(6) &
    - 2._dp*Iijm(5)**2 &
    - 4._dp*Iijm(5)*Iijm(6) &
    - 2._dp*Iijm(6)**2 &
    )
    soft2(-1) = soft2(-1) + CF**2 * ( &
    - 27._dp/10._dp*zeta2**2 &
    - 32._dp*log(y31)*zeta3 &
    - 9._dp*log(y31)**2*zeta2 &
    + 1._dp/2._dp*log(y31)**4 &
    + 2._dp*log(y31)**2*Iijm(5) &
    + 2._dp*log(y31)**2*Iijm(6) &
    - 2._dp*Iijm(5)*zeta2 &
    + 2._dp*Iijm(5)**2 &
    + 4._dp*Iijm(5)*Iijm(6) &
    - 2._dp*Iijm(6)*zeta2 &
    + 2._dp*Iijm(6)**2 &
    )
    return
    end subroutine soft_ab_qgq
    subroutine soft_nab_ggg(order,y12,y23,y31, &
    preIijm,j1,j2,j3,j4,j5,j6,soft2)
    implicit none
!     non-abelian piece
!     returns coefficients of [log(tau)^j/tau_+],and for j=-1, delta(tau)
!     in units of as/2/pi
!     y12=n1.n2/2 so that 0<y12<1
!     y23=n2.n3/2 so that 0<y23<1
!     y31=n3.n1/2 so that 0<y31<1
!     order indicates the order in [as/2/pi]
    include 'src/Inc/types.f'
    include 'src/Inc/constants.f'
    include 'src/Inc/nf.f'
    include 'src/Inc/scet_const.f'
    real(dp)::y12,y23,y31,soft2(-1:3),preIijm(5), &
    Iijm(6),Lss,c(0:3,0:3)
    integer::order,j1,j2,j3,j4,j5,j6,m,n
    logical, parameter:: useBLP= .FALSE. 
    include 'src/Inc/scet_beta.f'
!      integer,parameter::i(6)=(/1,2,2,3,3,1/)
!      integer,parameter::j(6)=(/2,1,3,2,1,3/)
!      integer,parameter::m(6)=(/3,3,1,1,2,2/)
! d,I12x3=Iijm(1);
! d,I21x3=Iijm(2);
! d,I23x1=Iijm(3);
! d,I32x1=Iijm(4);
! d,I31x2=Iijm(5);
! d,I13x2=Iijm(6);
    soft2(:)=zip

    if (order < 2) return

!      y(:,:)=zip

!      y(1,2)=y12
!      y(2,1)=y12
!      y(2,3)=y23
!      y(3,2)=y23
!      y(1,3)=y31
!      y(3,1)=y31

!      do n=1,6
!      al=y(j(n),m(n))/y(i(n),j(n))
!      be=y(i(n),m(n))/y(i(n),j(n))
!      I0=I0JSTW(al,be)
!      I1=I1JSTW(al,be)
!      Iijm(n)=I0*log(al)+I1
!      enddo

    Iijm(1)=preIijm(j1)
    Iijm(2)=preIijm(j2)
    Iijm(3)=preIijm(j3)
    Iijm(4)=preIijm(j4)
    Iijm(5)=preIijm(j5)
    Iijm(6)=preIijm(j6)

    soft2(3)=zip
    soft2(2) =6._dp*CA*be0
    soft2(1) = soft2(1) &
    - 3._dp/2._dp*CA*Ga1 &
    - 2._dp*log(y12)*CA*be0 &
    - 2._dp*log(y31)*CA*be0 &
    - 2._dp*log(y23)*CA*be0

    soft2(0) = &
    - 3._dp/4._dp*CA*gams1 &
    - 3._dp/2._dp*CA*zeta2*be0 &
    + quarter*log(y12)*CA*Ga1 &
    + half*log(y12)**2*CA*be0 &
    + quarter*log(y31)*CA*Ga1 &
    + half*log(y31)**2*CA*be0 &
    + quarter*log(y23)*CA*Ga1 &
    + half*log(y23)**2*CA*be0 &
    + Iijm(1)*CA*be0 &
    + Iijm(2)*CA*be0 &
    + Iijm(3)*CA*be0 &
    + Iijm(4)*CA*be0 &
    + Iijm(5)*CA*be0 &
    + Iijm(6)*CA*be0

    Lss=log(y31)+log(y23)
! BLP version
    soft2(-1)=74.1772_dp+30.2529_dp*Lss-18.7425_dp*Lss**2 &
    -5.61366_dp*Lss**3-0.584101_dp*Lss**4 &
    -0.0428711_dp*Lss**5-0.000301809_dp*Lss**6 &
    +0.00010388_dp*Lss**7+(2.55309e-6_dp)*Lss**8
!      write(6,*) 'BLP  ggg',soft2(-1)

    if (useBLP) return

! CEMW version
    c(:,:)=zip
    c(0,0)=63.1950_dp
    c(1,0)=33.6554_dp
    c(2,0)=-11.0155_dp
    c(3,0)=-2.2668_dp
    c(0,1)=33.5583_dp
    c(0,2)=-11.0910_dp
    c(0,3)=-2.2787_dp
    soft2(-1)=zip
    do m=0,3
        do n=0,3
            if ((m /= 0) .AND. (n /= 0)) cycle
            soft2(-1)=soft2(-1)+c(m,n)*log(y31)**m*log(y23)**n
        enddo
    enddo
!      write(6,*) 'CEMW ggg',soft2(-1)
!      pause

! c_{(0,0)}=63.1950 \pm 0.9619=42.3422 \pm 0.7973=39.1008 \pm 0.6981$\\
! c_{(1,0)}=33.6554 \pm 0.8395=25.2889 \pm 0.6945=13.7263 \pm 0.6150$\\
! c_{(2,0)}=-11.0155 \pm 0.2511=-9.0285 \pm 0.2071=-2.7374 \pm 0.1860$\\
! c_{(3,0)}=-2.2668 \pm 0.0244=-2.1487 \pm 0.0201=0.0164 \pm 0.0183$\\
! c_{(0,1)}=33.5583 \pm 0.8373=25.0011 \pm 0.6949=25.5908 \pm 0.6019$\\
! c_{(0,2)}=-11.0910 \pm 0.2493=-9.1788 \pm 0.2072=-8.7492 \pm 0.1773$\\
! c_{(0,3)}=-2.2787 \pm 0.0241=-2.1676 \pm 0.0201=-2.1255 \pm 0.0170$

    return
    end subroutine soft_nab_ggg
    subroutine soft_nab_qag(order,y12,y23,y31, &
    preIijm,j1,j2,j3,j4,j5,j6,soft2)
    implicit none
!--- Non-abelian piece q(1) qbar(2) g(3)
!--- returns coefficients of [log(tau)^j/tau_+],and for j=-1, delta(tau)
!--- in units of as/2/pi
!     y12=n1.n2/2 so that 0<y12<1
!     y23=n2.n3/2 so that 0<y23<1
!     y31=n3.n1/2 so that 0<y31<1
!     order indicates the order in [as/2/pi]
    include 'src/Inc/types.f'
    include 'src/Inc/constants.f'
    include 'src/Inc/nf.f'
    include 'src/Inc/scet_const.f'
    real(dp)::y12,y23,y31,soft2(-1:3),preIijm(5), &
    Iijm(6),Lss,Lsd,Isum12,Isum3456,c(0:3,0:3)
    integer::order,j1,j2,j3,j4,j5,j6,m,n
    logical, parameter:: useBLP= .FALSE. 
    include 'src/Inc/scet_beta.f'
!      integer,parameter::i(6)=(/1,2,2,3,3,1/)
!      integer,parameter::j(6)=(/2,1,3,2,1,3/)
!      integer,parameter::m(6)=(/3,3,1,1,2,2/)
! d,I12x3=Iijm(1);
! d,I21x3=Iijm(2);
! d,I23x1=Iijm(3);
! d,I32x1=Iijm(4);
! d,I31x2=Iijm(5);
! d,I13x2=Iijm(6);

    soft2(:)=zip
    if (order < 2) return

!      y(:,:)=zip

!      y(1,2)=y12
!      y(2,1)=y12
!      y(2,3)=y23
!      y(3,2)=y23
!      y(1,3)=y31
!      y(3,1)=y31

!      if ((y12 < zip) .or. (y12 > one)) then
!        write(6,*) 'error: y12 outside range, y12 = ',y12
!      endif
!      if ((y23 < zip) .or. (y23 > one)) then
!        write(6,*) 'error: y23 outside range, y23 = ',y23
!      endif
!      if ((y31 < zip) .or. (y31 > one)) then
!        write(6,*) 'error: y31 outside range, y31 = ',y31
!      endif

!      do n=1,6
!      al=y(j(n),m(n))/y(i(n),j(n))
!      be=y(i(n),m(n))/y(i(n),j(n))
!      I0=I0JSTW(al,be)
!      I1=I1JSTW(al,be)
!      Iijm(n)=I0*log(al)+I1
!      enddo

    Iijm(1)=preIijm(j1)
    Iijm(2)=preIijm(j2)
    Iijm(3)=preIijm(j3)
    Iijm(4)=preIijm(j4)
    Iijm(5)=preIijm(j5)
    Iijm(6)=preIijm(j6)

    Lss=log(y31)+log(y23)
    Lsd=log(y31)-log(y23)
    Isum12=Iijm(1)+Iijm(2)
    Isum3456=Iijm(3)+Iijm(4)+Iijm(5)+Iijm(6)


    soft2(3)= 0
    soft2(2)=2._dp*CA*be0 + 4._dp*CF*be0
    soft2(1)= - half*CA*Ga1 - 2._dp*CA*be0*Lss - CF*Ga1
    soft2(0)= - quarter*CA*gams1 + quarter*CA*Ga1*Lss + CA*be0* &
    Isum3456 - CA*be0*Isum12 + quarter*CA*be0*Lsd**2 + quarter* &
    CA*be0*Lss**2 - half*CA*zeta2*be0 - half*CF*gams1 &
    + 2._dp*CF*be0*Isum12 - CF*zeta2*be0

! BLP version
    soft2(-1)=60.1426_dp+41.0237_dp*Lss-2.57381_dp*Lss**2 &
    -0.692668_dp*Lss**3+0.226885_dp*Lss**4 &
    +0.0209989_dp*Lss**5+0.00090579_dp*Lss**6 &
    +0.0000132008_dp*Lss**7
!      write(6,*) 'BLP  qag',soft2(-1)

    if (useBLP) return

! CEMW version
    c(:,:)=zip
    c(0,0)=42.3422_dp
    c(1,0)=25.2889_dp
    c(2,0)=-9.0285_dp
    c(3,0)=-2.1487_dp
    c(0,1)=25.0011_dp
    c(0,2)=-9.1788_dp
    c(0,3)=-2.1676_dp
          
! enforce 1<->2 symmetry by hand
    c(1,0)=(c(1,0)+c(0,1))/two
    c(2,0)=(c(2,0)+c(0,2))/two
    c(3,0)=(c(3,0)+c(0,3))/two
    c(0,1)=c(1,0)
    c(0,2)=c(2,0)
    c(0,3)=c(3,0)

    soft2(-1)=zip
    do m=0,3
        do n=0,3
            if ((m /= 0) .AND. (n /= 0)) cycle
            soft2(-1)=soft2(-1)+c(m,n)*log(y31)**m*log(y23)**n
        enddo
    enddo
!      write(6,*) 'CEMW qag',soft2(-1)
!      pause

! c_{(0,0)}=63.1950 \pm 0.9619=42.3422 \pm 0.7973=39.1008 \pm 0.6981$\\
! c_{(1,0)}=33.6554 \pm 0.8395=25.2889 \pm 0.6945=13.7263 \pm 0.6150$\\
! c_{(2,0)}=-11.0155 \pm 0.2511=-9.0285 \pm 0.2071=-2.7374 \pm 0.1860$\\
! c_{(3,0)}=-2.2668 \pm 0.0244=-2.1487 \pm 0.0201=0.0164 \pm 0.0183$\\
! c_{(0,1)}=33.5583 \pm 0.8373=25.0011 \pm 0.6949=25.5908 \pm 0.6019$\\
! c_{(0,2)}=-11.0910 \pm 0.2493=-9.1788 \pm 0.2072=-8.7492 \pm 0.1773$\\
! c_{(0,3)}=-2.2787 \pm 0.0241=-2.1676 \pm 0.0201=-2.1255 \pm 0.0170$

!      soft2(-1)=soft2(-1)/four ! DEBUG

    return
    end subroutine soft_nab_qag
    subroutine soft_nab_qgq(order,y12,y23,y31, &
    preIijm,j1,j2,j3,j4,j5,j6,soft2)
    implicit none
!     non-abelian piece  q(1) g(2) q(3)
!     returns coefficients of [log(tau)^j/tau_+],and for j=-1, delta(tau)
!     in units of as/2/pi
!     y12=n1.n2/2 so that 0<y12<1
!     y23=n2.n3/2 so that 0<y23<1
!     y31=n3.n1/2 so that 0<y31<1
!     order indicates the order in [as/2/pi]
    include 'src/Inc/types.f'
    include 'src/Inc/constants.f'
    include 'src/Inc/nf.f'
    include 'src/Inc/scet_const.f'
    real(dp)::y12,y23,y31,soft2(-1:3),preIijm(5), &
    Iijm(6),Lss13,Lss23,c(0:3,0:3)
    integer::order,j1,j2,j3,j4,j5,j6,m,n
    logical, parameter:: useBLP= .FALSE. 
    include 'src/Inc/scet_beta.f'
!      integer,parameter::i(6)=(/1,2,2,3,3,1/)
!      integer,parameter::j(6)=(/2,1,3,2,1,3/)
!      integer,parameter::m(6)=(/3,3,1,1,2,2/)
! d,I12x3=Iijm(1);
! d,I21x3=Iijm(2);
! d,I23x1=Iijm(3);
! d,I32x1=Iijm(4);
! d,I31x2=Iijm(5);
! d,I13x2=Iijm(6);

    soft2(:)=zip
    if (order < 2) return

!      y(:,:)=zip

!      y(1,2)=y12
!      y(2,1)=y12
!      y(2,3)=y23
!      y(3,2)=y23
!      y(1,3)=y31
!      y(3,1)=y31

!      do n=1,6
!      al=y(j(n),m(n))/y(i(n),j(n))
!      be=y(i(n),m(n))/y(i(n),j(n))
!      I0=I0JSTW(al,be)
!      I1=I1JSTW(al,be)
!      Iijm(n)=I0*log(al)+I1
!      enddo

    Iijm(1)=preIijm(j1)
    Iijm(2)=preIijm(j2)
    Iijm(3)=preIijm(j3)
    Iijm(4)=preIijm(j4)
    Iijm(5)=preIijm(j5)
    Iijm(6)=preIijm(j6)

    soft2(3)=zip
    soft2(2)=2*be0*(CA+2*CF)
    soft2(1)=+CA*(-half*Ga1 &
    -two*log(y12)*be0+two*log(y31)*be0-two*log(y23)*be0) &
    +CF*(-Ga1-four*log(y31)*be0)
    soft2(0)=CA*( &
    - quarter*gams1 &
    - half*zeta2*be0 &
    + quarter*log(y12)*Ga1 &
    + half*log(y12)**2*be0 &
    - quarter*log(y31)*Ga1 &
    - half*log(y31)**2*be0 &
    + quarter*log(y23)*Ga1 &
    + half*log(y23)**2*be0 &
    + Iijm(1)*be0 &
    + Iijm(2)*be0 &
    + Iijm(3)*be0 &
    + Iijm(4)*be0 &
    - Iijm(5)*be0 &
    - Iijm(6)*be0 &
    )
    soft2(0) = soft2(0) + CF * ( &
    - half*gams1 &
    - zeta2*be0 &
    + half*log(y31)*Ga1 &
    + log(y31)**2*be0 &
    + two*Iijm(5)*be0 &
    + two*Iijm(6)*be0 &
    )

    Lss13=log(y31)
    Lss23=log(y23)
! BLP version
    soft2(-1)=-200.096_dp-223.007_dp*Lss13-210.796_dp*Lss23 &
    -112.644_dp*Lss13**2-118.705_dp*Lss23**2 &
    -30.8755_dp*Lss13**3-32.7046_dp*Lss23**3 &
    -5.57982_dp*Lss13**4-5.43977_dp*Lss23**4 &
    -0.666669_dp*Lss13**5-0.63292_dp*Lss23**5 &
    -0.0533767_dp*Lss13**6-0.0484011_dp*Lss23**6 &
    -0.0028547_dp*Lss13**7-0.00240717_dp*Lss23**7 &
    -0.0000952163_dp*Lss13**8-0.0000728699_dp*Lss23**8 &
    -1.45696e-6_dp*Lss13**9-9.83475e-7_dp*Lss23**9
!      write(6,*) 'BLP  qgq',soft2(-1)

    if (useBLP) return

! CEMW version
    c(:,:)=zip
    c(0,0)=39.1008_dp
    c(1,0)=13.7263_dp
    c(2,0)=-2.7374_dp
    c(3,0)=0.0164_dp
    c(0,1)=25.5908_dp
    c(0,2)=-8.7492_dp
    c(0,3)=-2.1255_dp
    soft2(-1)=zip
    do m=0,3
        do n=0,3
            if ((m /= 0) .AND. (n /= 0)) cycle
            soft2(-1)=soft2(-1)+c(m,n)*log(y31)**m*log(y23)**n
        enddo
    enddo
!      write(6,*) 'CEMW qgq',soft2(-1)
!      pause

! c_{(0,0)}=63.1950 \pm 0.9619=42.3422 \pm 0.7973=39.1008 \pm 0.6981$\\
! c_{(1,0)}=33.6554 \pm 0.8395=25.2889 \pm 0.6945=13.7263 \pm 0.6150$\\
! c_{(2,0)}=-11.0155 \pm 0.2511=-9.0285 \pm 0.2071=-2.7374 \pm 0.1860$\\
! c_{(3,0)}=-2.2668 \pm 0.0244=-2.1487 \pm 0.0201=0.0164 \pm 0.0183$\\
! c_{(0,1)}=33.5583 \pm 0.8373=25.0011 \pm 0.6949=25.5908 \pm 0.6019$\\
! c_{(0,2)}=-11.0910 \pm 0.2493=-9.1788 \pm 0.2072=-8.7492 \pm 0.1773$\\
! c_{(0,3)}=-2.2787 \pm 0.0241=-2.1676 \pm 0.0201=-2.1255 \pm 0.0170$

!      soft2(-1)=soft2(-1)/four ! DEBUG

    return
    end subroutine soft_nab_qgq
