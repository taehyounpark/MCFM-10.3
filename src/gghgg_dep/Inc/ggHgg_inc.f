!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      program main
      use openloops
      use double_precision
      implicit none

      integer :: id,j,rflav(5)
      real(8) :: psp(0:3,5), m2l2, acc, sqrt_s = 1000,amp2

      call set_parameter("order_ew", 1) ! coupling order

      id = register_process("21 1 -> 25 1 21", 12) !
      rflav=[0,1,25,1,0]


      as = 0.118d0
      mtsq=173d0**2
      vevsq = 246.22056907348590d0**2
      call start() ! initialise
      call set_parameter("verbose", 1)
      call set_parameter("mass(25)", 125.);
      call set_parameter("mass(6)", 173.);
      call set_parameter("alpha_s", 0.118d0)
      call set_parameter("mu", 100)

      do j=1,100
            call phase_space_point(id, sqrt_s, psp) ! obtain a random phase-space point from Rambo

            call evaluate_loop2(id, psp, m2l2, acc) ! calculate tree and loop matrix elements
            !print *, m2l2
            call setreal_mcfm(psp,rflav,amp2)
            !print *, amp2
            write(6,*) j,amp2/m2l2
      if (j == 36) then
            !write(6,*) 'psp',psp
            write(6,*) "p1",psp(:,1)
            write(6,*) "p2",psp(:,2)
            write(6,*) "p3",psp(:,3)
            write(6,*) "p4",psp(:,4)
            write(6,*) "p5",psp(:,5)
      endif
      enddo

      call finish()
      end program main
