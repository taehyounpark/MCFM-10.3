!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      module evtnum_module
      implicit none
      integer :: evtnum = 0
      end module evtnum_module

      subroutine writecsv(p,msq)
      use evtnum_module
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'nf.f'
      include 'first.f'
      include 'csvfile.f'
      real(dp):: p(mxpart,4),msq(fn:nf,fn:nf)
      integer :: i, j, iostat
      character(len=255)::rundir, filepath
      common/rundir/rundir

      ! csvfile = rundir // '/' // csvfile
      filepath = trim(rundir) // "/" // trim(csvfile)

      if (first) then
      first=.false.
      evtnum = 0
      open(unit=10, file=filepath, status='replace', action='write', iostat=i)
      ! Write headers
      write(10, '(A)', advance='no') 'evtnum,'
      write(10, '(A)', advance='no') 'p1_px,p1_py,p1_pz,p1_E,'
      write(10, '(A)', advance='no') 'p2_px,p2_py,p2_pz,p2_E,'
      write(10, '(A)', advance='no') 'p3_px,p3_py,p3_pz,p3_E,'
      write(10, '(A)', advance='no') 'p4_px,p4_py,p4_pz,p4_E,'
      write(10, '(A)', advance='no') 'p5_px,p5_py,p5_pz,p5_E,'
      write(10, '(A)', advance='no') 'p6_px,p6_py,p6_pz,p6_E,'
      write(10, '(A)', advance='yes') 'msq'
  
      if (iostat /= 0) then
            print*, 'Error opening file: ', csvfile
            stop
      endif
      endif

      ! Write data rows
      ! evetn number
      write(10, '(I10)', advance='no') evtnum
      write(10, '(A)', advance='no') ','
      ! parton momenta
      do i = 1, 6
        do j = 1, 4
          write(10, '(E22.15)', advance='no') p(i, j)
          write(10, '(A)', advance='no') ','
        end do
      end do
      ! matrix element (squared)
      write(10, '(E22.15)', advance='yes') msq(0,0)

      ! incrememt the event number
      evtnum = evtnum+1
  
      ! Close the file
      return
      end
