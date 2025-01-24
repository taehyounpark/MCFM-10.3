!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      ! only for true constants (no nf or such)
      module constants
          use types
          implicit none

          public

          real(dp), parameter :: pi = 3.14159265358979311599796346854418516_dp
          real(dp), parameter :: pisq = 9.8696044010893586188344909998762_dp
          real(dp), parameter :: fourpi = 12.566370614359172953850573533118_dp
          real(dp), parameter :: zeta2 = 1.6449340668482264364724151666460_dp
          real(dp), parameter :: zeta3 = 1.2020569031595942853997381615114500_dp
          real(dp), parameter :: zeta4 = 1.0823232337111381915160036965411679_dp
          real(dp), parameter :: zeta5 = 1.0369277551433699263313654864570342_dp

          real(dp), parameter :: EulerGamma = 0.57721566490153286060651209008240243_dp

          complex(dp), parameter ::  im = (0._dp, 1._dp)

          real(dp), parameter :: cf = 4._dp/3._dp
          real(dp), parameter :: ca = 3._dp
          real(dp), parameter :: tf = 1._dp/2._dp
          real(dp), parameter :: nc = 3._dp
          real(dp), parameter :: tr = 1._dp/2._dp
      end module
