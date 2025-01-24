!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

module SCET_Jet
    use types
    implicit none

    contains

    ! note that disttau = .false. actually returns the distribution in log(tau)
    ! coefficients
    subroutine jetq(order,Qi,J1q,J2q, scale_in, disttau)
!     Formula taken from A10,A11 of 1505.04794v1
    implicit none
    include 'constants.f'
    include 'scale.f'
    include 'nf.f'
    include 'scet_const.f'
    integer:: order
    real(dp):: Qi,logQi,J1q(-1:1),J2q(-1:3)
    real(dp), intent(in), optional :: scale_in
    real(dp) :: scale_used
    logical, intent(in), optional :: disttau
    include 'scet_beta.f'

    if (present(scale_in)) then
        scale_used = scale_in
    else
        scale_used = scale
    endif

!      J0q=1._dp
    J1q(1)=Gaq0
    J1q(0)=-half*gBq0
    J1q(-1)=CF*(7._dp-pisq)

!--- we want result as a distribution in log(taucut/scale)
!--- rather than log(taucut*Qi/musq) => extra +log(Qi/scale)
          
    if ((.not. present(disttau)) .or. disttau .eqv. .true.) then
        logQi=log(Qi/scale_used)
        J1q(-1)=J1q(-1)+J1q(0)*logQi+J1q(1)*logQi**2/2._dp
        J1q(0)=J1q(0)+logQi*J1q(1)
    endif

    J2q(:)=zip
    if (order < 2) return

    J2q(3)=+half*Gaq0**2
    J2q(2)=-half*Gaq0*(3*half*gBq0+be0)
    J2q(1)=Gaq1-Gaq0**2*zeta2+half*gBq0*(half*gBq0+be0) &
    +Gaq0*CF*(7._dp-pisq)
    J2q(0)=Gaq0**2*zeta3+Gaq0*gBq0*pisq/12._dp-half*gBq1 &
    -(half*gBq0+be0)*CF*(7._dp-pisq)
    J2q(-1)=CF*( &
    CF*(205/8._dp-67*pisq/6._dp+14*pisq**2/15._dp-18*zeta3) &
    +CA*(1417/108._dp-7*pisq/9._dp-17*pisq**2/180._dp-18*zeta3) &
    +be0*(4057/216._dp-17*pisq/9._dp-4*zeta3/3._dp))

    if ((.not. present(disttau)) .or. disttau .eqv. .true.) then
!--- add extra logs
        J2q(-1)=J2q(-1) &
        +J2q(0)*logQi &
        +J2q(1)*logQi**2/2._dp &
        +J2q(2)*logQi**3/3._dp &
        +J2q(3)*logQi**4/4._dp
        J2q(0)=J2q(0) &
        +J2q(1)*logQi+J2q(2)*logQi**2+J2q(3)*logQi**3
        J2q(1)=J2q(1) &
        +2._dp*J2q(2)*logQi+3._dp*J2q(3)*logQi**2
        J2q(2)=J2q(2) &
        +3._dp*J2q(3)*logQi
    endif
         
    return
    end subroutine jetq

    subroutine jetg(order,Qi,J1g,J2g, scale_in)
!     Formula taken from A10,A12 of 1505.04794v1
    implicit none
    include 'constants.f'
    include 'scale.f'
    include 'nf.f'
    include 'scet_const.f'
    integer:: order
    real(dp):: Qi,logQi,J1g(-1:1),J2g(-1:3)
    real(dp), intent(in), optional :: scale_in
    real(dp) :: scale_used
    include 'scet_beta.f'
!      J0g=1._dp

    if (present(scale_in)) then
        scale_used = scale_in
    else
        scale_used = scale
    endif

    J1g(1)=Gag0
    J1g(0)=-half*gBg0
    J1g(-1)=CA*(4/3._dp-pisq)+5/3._dp*be0

!--- we want result as a distribution in log(taucut/scale)
!--- rather than log(taucut*Qi/musq) => extra +log(Qi/scale)
    logQi=log(Qi/scale_used)
          
    J1g(-1)=J1g(-1)+J1g(0)*logQi+J1g(1)*logQi**2/2._dp
    J1g(0)=J1g(0)+logQi*J1g(1)

    J2g(:)=zip
    if (order < 2) return

    J2g(3)=+half*Gag0**2
    J2g(2)=-half*Gag0*(3*half*gBg0+be0)
    J2g(1)=Gag1-Gag0**2*zeta2+half*gBg0*(half*gBg0+be0) &
    +Gag0*(CA*(4/3._dp-pisq)+5/3._dp*be0)
    J2g(0)=Gag0**2*zeta3+gag0*gBg0*half*zeta2-half*gBg1 &
    -(half*gBg0+be0)*(CA*(4/3._dp-pisq)+5/3._dp*be0)
    J2g(-1)= &
    +CA**2*(4255/108._dp-26*pisq/9._dp+151*pisq**2/180._dp-72*zeta3) &
    +CA*be0*(-115/108._dp-65*pisq/18._dp+56*zeta3/3._dp) &
    +be0**2*(25/9._dp-pisq/3._dp)+be1*(55/12._dp-4*zeta3)

!--- add extra logs
    J2g(-1)=J2g(-1) &
    +J2g(0)*logQi &
    +J2g(1)*logQi**2/2._dp &
    +J2g(2)*logQi**3/3._dp &
    +J2g(3)*logQi**4/4._dp
    J2g(0)=J2g(0) &
    +J2g(1)*logQi+J2g(2)*logQi**2+J2g(3)*logQi**3
    J2g(1)=J2g(1) &
    +2._dp*J2g(2)*logQi+3._dp*J2g(3)*logQi**2
    J2g(2)=J2g(2) &
    +3._dp*J2g(3)*logQi
         
    return
    end subroutine jetg

end module
