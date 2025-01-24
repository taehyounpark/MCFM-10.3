!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      module evtnum_module
      implicit none
      end module evtnum_module

      subroutine output_c6(p, msq_sig_sm, msq_bkg_sm, msq_int_sm, msq_sbi_sm, msq_sig_c6, msq_int_c6, msq_sbi_c6, wt)
      use omp_lib 
      use MCFMStorage 
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'nf.f'
      include 'first.f'
      include 'ipsgen.f'
      include 'csvfile.f'
      include 'nplot.f'
      include 'nproc.f'
      include 'higgs_trilinear.f'
      real(dp) :: p(mxpart,4),wt,msq_sig_c6(c6_nval,-nf:nf,-nf:nf),msq_int_c6(c6_nval,-nf:nf,-nf:nf),msq_sbi_c6(c6_nval,-nf:nf,-nf:nf),msq_sig_sm(-nf:nf,-nf:nf),msq_bkg_sm(-nf:nf,-nf:nf),msq_int_sm(-nf:nf,-nf:nf),msq_sbi_sm(-nf:nf,-nf:nf)
      integer :: i, j, k
      character(len=255)::rundir, index_str
      common/rundir/rundir
      integer :: io_status
      integer :: file_unit
      character(len=255), save :: file_path

      write(file_path, '(A,I0)') trim(rundir) // "/" // trim(csvfile)
      file_unit = 10

      if (first) then
      first=.false.
      open(unit=file_unit, file=file_path, status='replace', action='write', iostat=io_status)
      if (io_status /= 0) then
        print*, 'Error opening file: ', file_path
        stop
      endif

      ! header
      write(file_unit, '(A)', advance='no') 'p1_px,p1_py,p1_pz,p1_E,'
      write(file_unit, '(A)', advance='no') 'p2_px,p2_py,p2_pz,p2_E,'
      write(file_unit, '(A)', advance='no') 'p3_px,p3_py,p3_pz,p3_E,'
      write(file_unit, '(A)', advance='no') 'p4_px,p4_py,p4_pz,p4_E,'
      write(file_unit, '(A)', advance='no') 'p5_px,p5_py,p5_pz,p5_E,'
      write(file_unit, '(A)', advance='no') 'p6_px,p6_py,p6_pz,p6_E,'
      write(file_unit, '(A)', advance='no') 'msq_sig_sm,'
      write(file_unit, '(A)', advance='no') 'msq_bkg_sm,'
      write(file_unit, '(A)', advance='no') 'msq_int_sm,'
      write(file_unit, '(A)', advance='no') 'msq_sbi_sm,'
      do k=1,c6_nval
        write(index_str, '(A,I0)') 'msq_sig_c6_',k
        write(file_unit, '(A)', advance='no') trim(adjustl(index_str))
        write(file_unit, '(A)', advance='no') ',' 
        write(index_str, '(A,I0)') 'msq_int_c6_',k
        write(file_unit, '(A)', advance='no') trim(adjustl(index_str))
        write(file_unit, '(A)', advance='no') ','
        write(index_str, '(A,I0)') 'msq_sbi_c6_',k
        write(file_unit, '(A)', advance='no') trim(adjustl(index_str))
        write(file_unit, '(A)', advance='no') ','
      end do
      write(file_unit, '(A)', advance='no') 'wt'
      write(file_unit, '(A)', advance='yes') ''
  
      endif

      ! parton momenta
      do i = 1, 6
        do j = 1, 4
          write(file_unit, '(E15.7)', advance='no') p(i, j)
          write(file_unit, '(A)', advance='no') ','
        end do
      end do

      ! msq
      write(file_unit, '(E15.7)', advance='no') msq_sig_sm(0,0)
      write(file_unit, '(A)', advance='no') ',' 
      write(file_unit, '(E15.7)', advance='no') msq_bkg_sm(0,0)
      write(file_unit, '(A)', advance='no') ','
      write(file_unit, '(E15.7)', advance='no') msq_int_sm(0,0)
      write(file_unit, '(A)', advance='no') ','
      write(file_unit, '(E15.7)', advance='no') msq_sbi_sm(0,0)
      write(file_unit, '(A)', advance='no') ','
      do k=1,c6_nval
        write(file_unit, '(E15.7)', advance='no') msq_sig_c6(k,0,0)
        write(file_unit, '(A)', advance='no') ','
        write(file_unit, '(E15.7)', advance='no') msq_int_c6(k,0,0)
        write(file_unit, '(A)', advance='no') ','
        write(file_unit, '(E15.7)', advance='no') msq_sbi_c6(k,0,0)
        write(file_unit, '(A)', advance='no') ','
      end do 

      ! weight
      write(file_unit, '(E15.7)', advance='no') wt
      !write(file_unit, '(A)', advance='no') ','

      ! newline character
      write(file_unit, '(A)', advance='yes') ''

      ! Close the file
      return
      end