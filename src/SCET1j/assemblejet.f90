!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

    function assemblejet(order,taucut,beama0,beamb0,beama1,beamb1, &
    beama2,beamb2,soft1,soft2,jet1,jet2,hard)
    implicit none
!---- Given beam, soft, jet and hard functions, calculates
!---- O(alphas^order) correction after all convolutions
    include 'src/Inc/types.f'
    include 'src/Inc/constants.f'
    include 'src/Inc/nf.f'
    include 'src/Inc/qcdcouple.f'
    include 'src/Inc/scale.f'
    include 'src/Inc/scet_const.f'
    include 'src/Inc/kpart.f'
    include 'src/Inc/toploops.f'
    integer:: order
    real(dp):: assemblejet,soft1(-1:1),soft2(-1:3),hard(2), &
    beama0,beamb0, &
    beama1(-1:1),beamb1(-1:1), &
    beama2(-1:3),beamb2(-1:3), &
    jet1(-1:1),jet2(-1:3),Lt1cut,full(-1:3),taucut

    Lt1cut = log(taucut/scale)

    if(coeffonly)then
        assemblejet=zip
    else
        assemblejet=beama0*beamb0
    endif
          
    if((order==1) .OR. &
    ((order==2) .AND. (coeffonly.eqv. .FALSE. )))then
        full(1)=beama1(1)*beamb0+beamb1(1)*beama0+beama0*beamb0* &
        soft1(1)+half*beama0*beamb0*jet1(1)

        full(0)=beama1(0)*beamb0+beamb1(0)*beama0+beama0*beamb0* &
        soft1(0)+half*beama0*beamb0*jet1(0)

        full(-1)=beama1(-1)*beamb0+beamb1(-1)*beama0+beama0*beamb0* &
        soft1(-1)+half*beama0*beamb0*jet1(-1)+beama0*beamb0* &
        hard(1)
        assemblejet=assemblejet+ason2pi*( &
        +      full(-1) &
        +full(0)*Lt1cut &
        +full(1)*Lt1cut**2/two)
    endif

    if(order>1)then
        full(3)=beama2(3)*beamb0+beamb2(3)*beama0+beama1(1)*beamb1(1) &
        +beama1(1)*beamb0*soft1(1)+half*beama1(1)*beamb0*jet1(1 &
        )+beamb1(1)*beama0*soft1(1)+half*beamb1(1)*beama0*jet1( &
        &  1)+half*beama0*beamb0*soft1(1)*jet1(1)+beama0*beamb0* &
        soft2(3)+1._dp/4._dp*beama0*beamb0*jet2(3)

        full(2)=beama2(2)*beamb0+beamb2(2)*beama0+3._dp/2._dp*beama1(0) &
        *beamb1(1)+3._dp/2._dp*beama1(0)*beamb0*soft1(1)+3._dp/4._dp* &
        beama1(0)*beamb0*jet1(1)+3._dp/2._dp*beama1(1)*beamb1(0)+3._dp/ &
        &  2._dp*beama1(1)*beamb0*soft1(0)+3._dp/4._dp*beama1(1)*beamb0* &
        jet1(0)+3._dp/2._dp*beamb1(0)*beama0*soft1(1)+3._dp/4._dp* &
        beamb1(0)*beama0*jet1(1)+3._dp/2._dp*beamb1(1)*beama0*soft1(0) &
        +3._dp/4._dp*beamb1(1)*beama0*jet1(0)+3._dp/4._dp*beama0*beamb0* &
        soft1(0)*jet1(1)+3._dp/4._dp*beama0*beamb0*soft1(1)*jet1(0)+ &
        beama0*beamb0*soft2(2)+1._dp/4._dp*beama0*beamb0*jet2(2)

        full(1)=beama2(1)*beamb0+beamb2(1)*beama0+beama1(-1)*beamb1(1 &
        )+beama1(-1)*beamb0*soft1(1)+half*beama1(-1)*beamb0* &
        jet1(1)+2._dp*beama1(0)*beamb1(0)+2._dp*beama1(0)*beamb0* &
        soft1(0)+beama1(0)*beamb0*jet1(0)+beama1(1)*beamb1(-1)-2._dp &
        *beama1(1)*beamb1(1)*zeta2+beama1(1)*beamb0*soft1(-1)-2._dp* &
        beama1(1)*beamb0*soft1(1)*zeta2+half*beama1(1)*beamb0* &
        jet1(-1)-beama1(1)*beamb0*jet1(1)*zeta2+beama1(1)*beamb0* &
        hard(1)+beamb1(-1)*beama0*soft1(1)+half*beamb1(-1)* &
        beama0*jet1(1)+2._dp*beamb1(0)*beama0*soft1(0)+beamb1(0)* &
        beama0*jet1(0)+beamb1(1)*beama0*soft1(-1)-2._dp*beamb1(1)* &
        beama0*soft1(1)*zeta2+half*beamb1(1)*beama0*jet1(-1)- &
        beamb1(1)*beama0*jet1(1)*zeta2+beamb1(1)*beama0*hard(1)+1._dp/ &
        &  2._dp*beama0*beamb0*soft1(-1)*jet1(1)+beama0*beamb0*soft1(0)* &
        jet1(0)+half*beama0*beamb0*soft1(1)*jet1(-1)-beama0* &
        beamb0*soft1(1)*jet1(1)*zeta2+beama0*beamb0*soft1(1)*hard(1) &
        +beama0*beamb0*soft2(1)
        full(1)=full(1)+half*beama0*beamb0*jet1(1)*hard(1)+1._dp &
        /4._dp*beama0*beamb0*jet2(1)

        full(0)=beama2(0)*beamb0+beamb2(0)*beama0+beama1(-1)*beamb1(0 &
        )+beama1(-1)*beamb0*soft1(0)+half*beama1(-1)*beamb0* &
        jet1(0)+beama1(0)*beamb1(-1)-beama1(0)*beamb1(1)*zeta2+ &
        beama1(0)*beamb0*soft1(-1)-beama1(0)*beamb0*soft1(1)*zeta2+ &
        half*beama1(0)*beamb0*jet1(-1)-half*beama1(0)*beamb0 &
        *jet1(1)*zeta2+beama1(0)*beamb0*hard(1)-beama1(1)*beamb1(0)* &
        zeta2+2._dp*beama1(1)*beamb1(1)*zeta3-beama1(1)*beamb0*soft1( &
        &  0)*zeta2+2._dp*beama1(1)*beamb0*soft1(1)*zeta3-half* &
        beama1(1)*beamb0*jet1(0)*zeta2+beama1(1)*beamb0*jet1(1)*zeta3 &
        +beamb1(-1)*beama0*soft1(0)+half*beamb1(-1)*beama0* &
        jet1(0)+beamb1(0)*beama0*soft1(-1)-beamb1(0)*beama0*soft1(1) &
        *zeta2+half*beamb1(0)*beama0*jet1(-1)-half*beamb1( &
        &  0)*beama0*jet1(1)*zeta2+beamb1(0)*beama0*hard(1)-beamb1(1)* &
        beama0*soft1(0)*zeta2+2._dp*beamb1(1)*beama0*soft1(1)*zeta3- &
        half*beamb1(1)*beama0*jet1(0)*zeta2+beamb1(1)*beama0*jet1(1 &
        )*zeta3
        full(0)=full(0)+half*beama0*beamb0*soft1(-1)*jet1(0)+ &
        half*beama0*beamb0*soft1(0)*jet1(-1)-half*beama0* &
        beamb0*soft1(0)*jet1(1)*zeta2+beama0*beamb0*soft1(0)*hard(1) &
        -half*beama0*beamb0*soft1(1)*jet1(0)*zeta2+beama0* &
        beamb0*soft1(1)*jet1(1)*zeta3+beama0*beamb0*soft2(0)+half &
        *beama0*beamb0*jet1(0)*hard(1)+1._dp/4._dp*beama0*beamb0*jet2(0)

        full(-1)=beama2(-1)*beamb0+beamb2(-1)*beama0+beama1(-1)* &
        beamb1(-1)+beama1(-1)*beamb0*soft1(-1)+half*beama1(-1)* &
        beamb0*jet1(-1)+beama1(-1)*beamb0*hard(1)-beama1(0)*beamb1(0 &
        )*zeta2+beama1(0)*beamb1(1)*zeta3-beama1(0)*beamb0*soft1(0)* &
        zeta2+beama1(0)*beamb0*soft1(1)*zeta3-half*beama1(0)* &
        beamb0*jet1(0)*zeta2+half*beama1(0)*beamb0*jet1(1)*zeta3 &
        +beama1(1)*beamb1(0)*zeta3-1._dp/10._dp*beama1(1)*beamb1(1)* &
        zeta2**2+beama1(1)*beamb0*soft1(0)*zeta3-1._dp/10._dp*beama1(1 &
        )*beamb0*soft1(1)*zeta2**2+half*beama1(1)*beamb0*jet1(0)* &
        zeta3-1._dp/20._dp*beama1(1)*beamb0*jet1(1)*zeta2**2+beamb1(-1 &
        )*beama0*soft1(-1)+half*beamb1(-1)*beama0*jet1(-1)+ &
        beamb1(-1)*beama0*hard(1)-beamb1(0)*beama0*soft1(0)*zeta2+ &
        beamb1(0)*beama0*soft1(1)*zeta3-half*beamb1(0)*beama0* &
        jet1(0)*zeta2+half*beamb1(0)*beama0*jet1(1)*zeta3+ &
        beamb1(1)*beama0*soft1(0)*zeta3-1._dp/10._dp*beamb1(1)*beama0* &
        soft1(1)*zeta2**2
        full(-1)=full(-1)+half*beamb1(1)*beama0*jet1(0)*zeta3- &
        &  1._dp/20._dp*beamb1(1)*beama0*jet1(1)*zeta2**2+half*beama0* &
        beamb0*soft1(-1)*jet1(-1)+beama0*beamb0*soft1(-1)*hard(1) &
        -half*beama0*beamb0*soft1(0)*jet1(0)*zeta2+half*beama0* &
        beamb0*soft1(0)*jet1(1)*zeta3+half*beama0*beamb0*soft1(1) &
        *jet1(0)*zeta3-1._dp/20._dp*beama0*beamb0*soft1(1)*jet1(1)* &
        zeta2**2+beama0*beamb0*soft2(-1)+half*beama0*beamb0* &
        jet1(-1)*hard(1)+1._dp/4._dp*beama0*beamb0*jet2(-1)+beama0* &
        beamb0*hard(2)
        
        if (onlyaxial) then
! collect only the terms that contain hard(1)
        full(3)=0._dp

        full(2)=0._dp

        full(1)=+beama1(1)*beamb0*hard(1)+beamb1(1)*beama0*hard(1) &
        +beama0*beamb0*soft1(1)*hard(1)+half*beama0*beamb0*jet1(1)*hard(1)

        full(0)=+beama1(0)*beamb0*hard(1)+beamb1(0)*beama0*hard(1) &
        +beama0*beamb0*soft1(0)*hard(1)+half*beama0*beamb0*jet1(0)*hard(1)

        full(-1)=beama1(-1)*beamb0*hard(1)+beamb1(-1)*beama0*hard(1) &
        +beama0*beamb0*soft1(-1)*hard(1)+half*beama0*beamb0*jet1(-1)*hard(1)
        endif
        
        assemblejet=assemblejet+ason2pi**2*( &
        +full(-1) &
        +full(0)*Lt1cut &
        +full(1)*Lt1cut**2/two &
        +full(2)*Lt1cut**3/three &
        +full(3)*Lt1cut**4/four)
    endif

    return
    end function assemblejet
