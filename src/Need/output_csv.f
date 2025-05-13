!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      module evtnum_module
      implicit none
      end module evtnum_module

      subroutine output_csv(p, msq_sig_sm, msq_bkg_sm, msq_int_sm, msq_sbi_sm, msq_sig_bsm, msq_int_bsm, msq_sbi_bsm, wt)
      use omp_lib 
      use MCFMStorage 
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'nf.f'
      include 'first.f'
      include 'ipsgen.f'
      include 'nplot.f'
      include 'nproc.f'
      include 'csvfile.f'
      include 'bsm_higgs.f'
      real(dp) :: p(mxpart,4),wt,msq_sig_bsm(c6_nval*ct_nval*cg_nval,-nf:nf,-nf:nf),msq_int_bsm(c6_nval*ct_nval*cg_nval,-nf:nf,-nf:nf),msq_sbi_bsm(c6_nval*ct_nval*cg_nval,-nf:nf,-nf:nf),msq_sig_sm(-nf:nf,-nf:nf),msq_bkg_sm(-nf:nf,-nf:nf),msq_int_sm(-nf:nf,-nf:nf),msq_sbi_sm(-nf:nf,-nf:nf)
      integer :: i, j, k, l, m, n
      character(len=255)::rundir, index_str
      common/rundir/rundir
      integer :: io_status
      integer :: event_file_unit, metadata_file_unit
      character(len=255), save :: event_file_path, metadata_file_path

      write(event_file_path, '(A,I0)') trim(rundir) // "/" // trim(eventfile)
      event_file_unit = 10

      write(metadata_file_path, '(A,I0)') trim(rundir) // "/" // trim(metadatafile)
      metadata_file_unit = 11

      if (first) then

      first=.false.

      open(unit=metadata_file_unit, file=metadata_file_path, status='replace', action='write', iostat=io_status)
      if (io_status /= 0) then
        print*, 'Error opening file: ', metadata_file_path
        stop
      endif
      write(metadata_file_unit, '(A)', advance='yes') 'n,c6,ct,cg'
      n = 1
      do k=1,c6_nval
        do l=1,ct_nval
          do m=1,cg_nval
            write(metadata_file_unit, '(I0)', advance='no') n
            write(metadata_file_unit, '(A)', advance='no') ',' 
            write(metadata_file_unit, '(E15.7)', advance='no') c6_init + (k-1)*c6_step
            write(metadata_file_unit, '(A)', advance='no') ',' 
            write(metadata_file_unit, '(E15.7)', advance='no') ct_init + (l-1)*ct_step
            write(metadata_file_unit, '(A)', advance='no') ',' 
            write(metadata_file_unit, '(E15.7)', advance='no') cg_init + (m-1)*cg_step
            write(metadata_file_unit, '(A)', advance='yes') ''
            n = n+1
          end do
        end do
      end do
      close(unit=metadata_file_unit, status='keep')

      open(unit=event_file_unit, file=event_file_path, status='replace', action='write', iostat=io_status)
      if (io_status /= 0) then
        print*, 'Error opening file: ', event_file_path
        stop
      endif

      ! header
      write(event_file_unit, '(A)', advance='no') 'p1_px,p1_py,p1_pz,p1_E,'
      write(event_file_unit, '(A)', advance='no') 'p2_px,p2_py,p2_pz,p2_E,'
      write(event_file_unit, '(A)', advance='no') 'p3_px,p3_py,p3_pz,p3_E,'
      write(event_file_unit, '(A)', advance='no') 'p4_px,p4_py,p4_pz,p4_E,'
      write(event_file_unit, '(A)', advance='no') 'p5_px,p5_py,p5_pz,p5_E,'
      write(event_file_unit, '(A)', advance='no') 'p6_px,p6_py,p6_pz,p6_E,'
      write(event_file_unit, '(A)', advance='no') 'msq_sig_sm,'
      write(event_file_unit, '(A)', advance='no') 'msq_bkg_sm,'
      write(event_file_unit, '(A)', advance='no') 'msq_int_sm,'
      write(event_file_unit, '(A)', advance='no') 'msq_sbi_sm,'
      n = 1
      do k=1,c6_nval
        do l=1,ct_nval
          do m=1,cg_nval
            write(event_file_unit, '(A)', advance='no') 'msq_sig_bsm_'
            write(event_file_unit, '(I0)', advance='no') n
            write(event_file_unit, '(A)', advance='no') ','
            write(event_file_unit, '(A)', advance='no') 'msq_int_bsm_'
            write(event_file_unit, '(I0)', advance='no') n
            write(event_file_unit, '(A)', advance='no') ','
            write(event_file_unit, '(A)', advance='no') 'msq_sbi_bsm_'
            write(event_file_unit, '(I0)', advance='no') n
            write(event_file_unit, '(A)', advance='no') ','
            n = n+1
          end do
        end do
      end do
      write(event_file_unit, '(A)', advance='yes') 'wt'
  
      endif

      ! parton momenta
      do i = 1, 6
        do j = 1, 4
          write(event_file_unit, '(E15.7)', advance='no') p(i, j)
          write(event_file_unit, '(A)', advance='no') ','
        end do
      end do

      ! msq
      write(event_file_unit, '(E15.7)', advance='no') msq_sig_sm(0,0)
      write(event_file_unit, '(A)', advance='no') ',' 
      write(event_file_unit, '(E15.7)', advance='no') msq_bkg_sm(0,0)
      write(event_file_unit, '(A)', advance='no') ','
      write(event_file_unit, '(E15.7)', advance='no') msq_int_sm(0,0)
      write(event_file_unit, '(A)', advance='no') ','
      write(event_file_unit, '(E15.7)', advance='no') msq_sbi_sm(0,0)
      write(event_file_unit, '(A)', advance='no') ','
      n = 1
      do k=1,c6_nval
        do l=1,ct_nval
          do m=1,cg_nval
            write(event_file_unit, '(E15.7)', advance='no') msq_sig_bsm(n,0,0)
            write(event_file_unit, '(A)', advance='no') ','
            write(event_file_unit, '(E15.7)', advance='no') msq_int_bsm(n,0,0)
            write(event_file_unit, '(A)', advance='no') ','
            write(event_file_unit, '(E15.7)', advance='no') msq_sbi_bsm(n,0,0)
            write(event_file_unit, '(A)', advance='no') ','
            n = n+1
          end do
        end do
      end do

      ! weight
      write(event_file_unit, '(E15.7)', advance='yes') wt
      !write(event_file_unit, '(A)', advance='no') ','

      ! Close the file
      return
      end