!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
       subroutine gg_hgg_mass_tb(p,msq)
        use double_precision
        use sprod_dp
        implicit none
        include 'Inc/hdecaymode.f'
        include 'msq_struc.f'
        include 'masses.f'

        real(dp), intent(in) :: p(mxpart,4)
        real(dp), intent(out) :: msq(-nf:nf,-nf:nf)

        integer, parameter :: iglue1 = 5, iglue2 = 6
        real(dp) :: hdecay,s34
        real(dp) :: dotvec, msqhgamgam

        s34 = dotvec(p(3,:)+p(4,:), p(3,:)+p(4,:))

        if (hdecaymode == 'tlta') then
            call htautaudecay(p,3,4,hdecay)
        elseif (hdecaymode == 'bqba') then
            call hbbdecay(p,3,4,hdecay)
        elseif (hdecaymode == 'gaga') then
            hdecay=msqhgamgam(s34)
        else
        write(6,*) 'Unimplemented process in gg_hgg_v'
        stop
        endif
        hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)

        call gg_hgg_mass_tb_nodecay(p,iglue1,iglue2,msq)

        msq = msq * hdecay
        msq_struc = msq_struc * hdecay

      end subroutine
