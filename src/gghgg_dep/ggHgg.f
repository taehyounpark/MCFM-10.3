!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      program main
      use mod_qcdloop_c
      use openloops
      use setreal_mcfm_generic
      use testreal_generic
      implicit none
      integer, parameter :: dp=selected_real_kind(precision(1.0d0))
      integer, parameter :: qp=selected_real_kind(2*precision(1.0d0))

      integer :: id,j,rflav(5)
      real(dp) :: psp(0:3,5), m2l2, acc, sqrt_s = 1000,amp2,gf,vsq
      real(qp) :: psp_qp(0:3,5),amp2_qp,m2l2_qp,acc_qp,sqrt_s_qp=1000,gf_qp,vsq_qp

      call qlcachesize(200)
      call set_parameter("order_ew", 1) ! coupling order

      id = register_process("2 -2 -> 25 1 -1", 12) !
      rflav=[2, -2,25, 1,-1 ]

      call start() ! initialise
      call set_parameter("verbose", 1)
      call set_parameter("mass(25)", 125._dp);
      call set_parameter("mass(6)", 173._dp);
      call set_parameter("alpha_s", 0.118_dp)
      call set_parameter("mu", 100._dp)

      do j=1,2
            call phase_space_point(id, sqrt_s, psp) ! obtain a random phase-space point from Rambo
            call evaluate_loop2(id, psp, m2l2, acc) ! calculate tree and loop matrix elements
            print *, 'm2l2,acc          ',m2l2,acc
            write(6,*)
            print *, 'Double precision'
            write(6,*) 'rflav',rflav
            call setreal_mcfm(psp,rflav,amp2)
c            write(6,*) 'psp(0,1)   ',psp(0,1)
c            write(6,*) 'psp(0,2)   ',psp(0,2)
c            write(6,*) 'psp(0,3)   ',psp(0,3)
c            write(6,*) 'psp(0,4)   ',psp(0,4)
c            write(6,*) 'psp(0,5)   ',psp(0,5)
            write(6,*) 'j,amp2,amp2/m2l2         ',j,amp2,amp2/m2l2
            write(6,*)
            print *, 'Quad precision'
            psp_qp=real(psp(:,:),kind=qp)
            m2l2_qp=real(m2l2,kind=qp)
            acc_qp=real(acc,kind=qp)
            print *, 'm2l2,acc',m2l2_qp,acc_qp
c            write(6,*) 'psp_qp(0,1)',psp_qp(0,1)
c            write(6,*) 'psp_qp(0,2)',psp_qp(0,2)
c            write(6,*) 'psp_qp(0,3)',psp_qp(0,3)
c            write(6,*) 'psp_qp(0,4)',psp_qp(0,4)
c            write(6,*) 'psp_qp(0,5)',psp_qp(0,5)
            call setreal_mcfm(psp_qp,rflav,amp2_qp)
            write(6,*) 'j,amp2_qp,amp2_qp/m2l2_qp',j,amp2_qp,amp2_qp/m2l2_qp
            write(6,*)
            write(6,*)


c            call testReal(psp)

      enddo

      call finish()
      end program main


