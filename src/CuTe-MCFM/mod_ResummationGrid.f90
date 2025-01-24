
!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module qtResummationGrids
      use ieee_arithmetic
      use types
      use LHAPDF
      use PDFerrors
      implicit none

      public :: readAllGrids

      private

    contains

    subroutine readAllGrids()
        use Beamfunctions3L, only: getbeam, getbeam2
        use parseinput
        use iso_fortran_env
#ifdef HAVE_MPI
        use mpi
#endif
        implicit none
#if defined(HAVE_MPI) || defined(HAVE_COARRAY)
        include 'src/Inc/mpicommon.f'
#endif

        integer :: numpdfmembers, pdfentry
        integer :: i,j,k,l,m,n,nimg,bb
        integer :: ierr,iread

        character(len=255) :: dir = "/data/tneumann/local/share/LHAPDF/"
        character(len=255) :: outdir = '/data/tneumann/grids/'
        character(len=512) :: path, imsg, outpath

        character(len=255) :: pdftype, pdfformat, dummy
        character(len=5000) :: nextline

! total number of beam functions generated; different (as, Lperp)
#define NUM_BEAMFUNS 17


#ifdef HAVE_COARRAY
        real(dp), allocatable :: pdfvals(:,:,:,:)[:]
        logical(atomic_logical_kind), allocatable :: work_done(:)[:]
        logical(atomic_logical_kind), save :: prev[*]
        real(dp), save :: xvals(512)[*]
        real(dp), save :: q2vals(512)[*]
        integer, save :: flavors(13)[*]
        integer, save :: curPDF[*]
#else
        real(dp) :: xvals(512)
        real(dp) :: q2vals(512)
        real(dp), allocatable :: pdfvals(:,:,:,:)
        integer :: flavors(13)
        integer :: curPDF
#endif


#ifdef HAVE_COARRAY
        integer, save :: numx[*], numq2[*], numflavors[*]
#else
        integer :: numx, numq2, numflavors
#endif
        integer, save :: fl
!$omp threadprivate(fl)

#if defined(HAVE_COARRAY)

        if (this_image() /= 1) then
            write (*,*) "Worker on image ", this_image(), "launched"

            do while (.true.)
                sync all ! sync numx, numq2, numflavors
                currentPDF = curPDF
                if (allocated(pdfvals)) deallocate(pdfvals)
                allocate(pdfvals(NUM_BEAMFUNS,numflavors,numq2,numx)[*], source=huge(0._dp))
                if (allocated(work_done)) deallocate(work_done)
                allocate(work_done(numx)[*], source=.false.)

                do j=1,numx
                    call atomic_cas(work_done(j)[1], prev, .false., .true.)
                    if (prev .eqv. .false.) then
                        !write (*,*) "Image ", this_image(), " got work x = ", xvals(j)
!$omp parallel do collapse(2)
                        do l=1,numq2
                            do m=1,numflavors
                                fl = flavors(m)
                                if (fl == 21) fl = 0
                                ! we hardcode ih=1 and beam=1
                                ! this means no pdfchannel selection should be done in the input file
                                ! for this grid generation
                                pdfvals(1,m,l,j)[1] =  getbeam(1, fl, 0, 0, xvals(j), q2vals(l))*xvals(j)
                                pdfvals(2,m,l,j)[1] =  getbeam(1, fl, 1, 0, xvals(j), q2vals(l))*xvals(j)
                                pdfvals(3,m,l,j)[1] =  getbeam(1, fl, 1, 1, xvals(j), q2vals(l))*xvals(j)

                                pdfvals(4,m,l,j)[1] =  getbeam(1, fl, 2, 0, xvals(j), q2vals(l))*xvals(j)
                                pdfvals(5,m,l,j)[1] =  getbeam(1, fl, 2, 1, xvals(j), q2vals(l))*xvals(j)
                                pdfvals(6,m,l,j)[1] =  getbeam(1, fl, 2, 2, xvals(j), q2vals(l))*xvals(j)

                                pdfvals(7,m,l,j)[1] =  getbeam(1, fl, 3, 0, xvals(j), q2vals(l))*xvals(j)
                                pdfvals(8,m,l,j)[1] =  getbeam(1, fl, 3, 1, xvals(j), q2vals(l))*xvals(j)
                                pdfvals(9,m,l,j)[1] =  getbeam(1, fl, 3, 2, xvals(j), q2vals(l))*xvals(j)
                                pdfvals(10,m,l,j)[1] =  getbeam(1, fl, 3, 3, xvals(j), q2vals(l))*xvals(j)

                                pdfvals(11,m,l,j)[1] = getbeam2(1, fl, 1, 0, xvals(j), q2vals(l))*xvals(j)

                                ! g^4
                                pdfvals(12,m,l,j)[1] =  getbeam(1, fl, 4, 4, xvals(j), q2vals(l))*xvals(j)

                                ! g^5
                                pdfvals(13,m,l,j)[1] =  getbeam(1, fl, 4, 3, xvals(j), q2vals(l))*xvals(j)
                                pdfvals(14,m,l,j)[1] =  getbeam(1, fl, 5, 5, xvals(j), q2vals(l))*xvals(j)

                                ! g^6
                                pdfvals(15,m,l,j)[1] =  getbeam(1, fl, 4, 2, xvals(j), q2vals(l))*xvals(j)
                                pdfvals(16,m,l,j)[1] =  getbeam(1, fl, 5, 4, xvals(j), q2vals(l))*xvals(j)
                                pdfvals(17,m,l,j)[1] =  getbeam(1, fl, 6, 6, xvals(j), q2vals(l))*xvals(j)

                                do bb=1,NUM_BEAMFUNS
                                    if (ieee_is_nan(pdfvals(bb,m,l,j)[1])) pdfvals(bb,m,l,j)[1] = 0._dp
                                enddo

                                ! otherwise format length is too short for for E+XXX,
                                do bb=1,NUM_BEAMFUNS
                                    if (abs(pdfvals(bb,m,l,j)[1]) < 1d-99) pdfvals(bb,m,l,j)[1] = 0._dp
                                enddo

                            enddo
                        enddo
!$omp end parallel do
                    else
                        continue
                    endif
                enddo
                 write (*,*) "Image ", this_image(), " checked all work"
                sync all
            enddo
        endif

        if (this_image() /= 1) then
            return
        endif

#elif defined(HAVE_MPI)
        logical :: working(world_size-1)
        logical :: taskdone, morework
        integer :: nextx, sendx, recvx, mpistatus(mpi_status_size)
        integer :: wxid
        real(dp) :: wx
        real(dp), allocatable :: pdfvals_slice(:,:,:)

        working(:) = .false.

        if (rank /= 0) then

            write (*,*) "Worker on rank ", rank, " launched"

            do while (.true.)

            call mpi_recv(wx, 1, mpi_double_precision, 0, 1, mpi_comm_world, mpistatus, ierr)
            call mpi_recv(wxid, 1, mpi_integer, 0, 2, mpi_comm_world, mpistatus, ierr)
            call mpi_recv(numq2, 1, mpi_integer, 0, 3, mpi_comm_world, mpistatus, ierr)
            call mpi_recv(numflavors, 1, mpi_integer, 0, 4, mpi_comm_world, mpistatus, ierr)
            call mpi_recv(q2vals, size(q2vals), mpi_double_precision, 0, 5, mpi_comm_world, mpistatus, ierr)
            call mpi_recv(flavors, size(flavors), mpi_integer, 0, 6, mpi_comm_world, mpistatus, ierr)
            call mpi_recv(curPDF, 1, mpi_integer, 0, NUM_BEAMFUNS, mpi_comm_world, mpistatus, ierr)


            if (allocated(pdfvals_slice)) deallocate(pdfvals_slice)
            allocate(pdfvals_slice(NUM_BEAMFUNS,numflavors,numq2), source=huge(0._dp))

!$omp parallel do collapse(2)
            do l=1,numq2
                do m=1,numflavors
                    fl = flavors(m)
                    if (fl == 21) fl = 0
                    ! we hardcode ih=1 and beam=1
                    ! this means no pdfchannel selection should be done in the input file
                    ! for this grid generation
                    pdfvals_slice(1,m,l) = getbeam(1, fl, 0, 0, wx, q2vals(l))*wx
                    pdfvals_slice(2,m,l) = getbeam(1, fl, 1, 0, wx, q2vals(l))*wx
                    pdfvals_slice(3,m,l) = getbeam(1, fl, 1, 1, wx, q2vals(l))*wx
                    pdfvals_slice(4,m,l) = getbeam(1, fl, 2, 0, wx, q2vals(l))*wx
                    pdfvals_slice(5,m,l) = getbeam(1, fl, 2, 1, wx, q2vals(l))*wx
                    pdfvals_slice(6,m,l) = getbeam(1, fl, 2, 2, wx, q2vals(l))*wx

                    pdfvals_slice(7,m,l) = getbeam(1, fl, 3, 0, wx, q2vals(l))*wx
                    pdfvals_slice(8,m,l) = getbeam(1, fl, 3, 1, wx, q2vals(l))*wx
                    pdfvals_slice(9,m,l) = getbeam(1, fl, 3, 2, wx, q2vals(l))*wx
                    pdfvals_slice(10,m,l) = getbeam(1, fl, 3, 3, wx, q2vals(l))*wx

                    pdfvals_slice(11,m,l) = getbeam2(1, fl, 1, 0, wx, q2vals(l))*wx

                    pdfvals_slice(12,m,l) = getbeam(1, fl, 4, 4, wx, q2vals(l))*wx

                    pdfvals_slice(13,m,l) = getbeam(1, fl, 4, 3, wx, q2vals(l))*wx
                    pdfvals_slice(14,m,l) = getbeam(1, fl, 5, 5, wx, q2vals(l))*wx

                    pdfvals_slice(15,m,l) = getbeam(1, fl, 4, 2, wx, q2vals(l))*wx
                    pdfvals_slice(16,m,l) = getbeam(1, fl, 5, 4, wx, q2vals(l))*wx
                    pdfvals_slice(17,m,l) = getbeam(1, fl, 6, 6, wx, q2vals(l))*wx

                    do bb=1,NUM_BEAMFUNS
                        if (ieee_is_nan(pdfvals_slice(bb,m,l))) pdfvals_slice(bb,m,l) = 0._dp
                    enddo

                    do bb=1,NUM_BEAMFUNS
                        if (abs(pdfvals_slice(bb,m,l)) < 1d-99) pdfvals_slice(bb,m,l) = 0._dp
                    enddo

                enddo
            enddo
!$omp end parallel do

            ! signal master work is ready to pick up
            call mpi_send(wxid, 1, mpi_integer, 0, 1, mpi_comm_world, ierr)

            ! this is the signal to send work to master
            call mpi_recv(morework, 1, mpi_logical, 0, 0, mpi_comm_world, mpistatus, ierr)

            call mpi_send(pdfvals_slice, NUM_BEAMFUNS*numflavors*numq2, mpi_double_precision, &
                0, 2, mpi_comm_world, ierr)

            enddo
        endif

#endif

        call cfg_get(cfg, "resummation%gridoutpath", outdir)
        call cfg_get(cfg, "resummation%gridinpath", dir)

        pdfentry = 0
        do i=1,numPDFsets
            if (dopdferrors .eqv. .true.) then
                numpdfmembers = lhapdf_number(trim(PDFnames(i)))
            else
                numpdfmembers = 1
            endif

            write (*,*) "RUNNING ","mkdir -p "//trim(outdir)//trim(PDFnames(i))//'_B00'

            call execute_command_line("mkdir -p "//trim(outdir)//trim(PDFnames(i))//'_B00')
            call execute_command_line("cp "//trim(dir)//trim(PDFnames(i))//'/'//trim(PDFnames(i))//'.info ' &
                //trim(outdir)//trim(PDFnames(i))//'_B00/'//trim(PDFnames(i))//'_B00.info' )

            call execute_command_line("mkdir -p "//trim(outdir)//trim(PDFnames(i))//'_B10')
            call execute_command_line("cp "//trim(dir)//trim(PDFnames(i))//'/'//trim(PDFnames(i))//'.info ' &
                //trim(outdir)//trim(PDFnames(i))//'_B10/'//trim(PDFnames(i))//'_B10.info' )

            call execute_command_line("mkdir -p "//trim(outdir)//trim(PDFnames(i))//'_B11')
            call execute_command_line("cp "//trim(dir)//trim(PDFnames(i))//'/'//trim(PDFnames(i))//'.info ' &
                //trim(outdir)//trim(PDFnames(i))//'_B11/'//trim(PDFnames(i))//'_B11.info' )

            call execute_command_line("mkdir -p "//trim(outdir)//trim(PDFnames(i))//'_B20')
            call execute_command_line("cp "//trim(dir)//trim(PDFnames(i))//'/'//trim(PDFnames(i))//'.info ' &
                //trim(outdir)//trim(PDFnames(i))//'_B20/'//trim(PDFnames(i))//'_B20.info' )

            call execute_command_line("mkdir -p "//trim(outdir)//trim(PDFnames(i))//'_B21')
            call execute_command_line("cp "//trim(dir)//trim(PDFnames(i))//'/'//trim(PDFnames(i))//'.info ' &
                //trim(outdir)//trim(PDFnames(i))//'_B21/'//trim(PDFnames(i))//'_B21.info' )

            call execute_command_line("mkdir -p "//trim(outdir)//trim(PDFnames(i))//'_B22')
            call execute_command_line("cp "//trim(dir)//trim(PDFnames(i))//'/'//trim(PDFnames(i))//'.info ' &
                //trim(outdir)//trim(PDFnames(i))//'_B22/'//trim(PDFnames(i))//'_B22.info' )

            call execute_command_line("mkdir -p "//trim(outdir)//trim(PDFnames(i))//'_B30')
            call execute_command_line("cp "//trim(dir)//trim(PDFnames(i))//'/'//trim(PDFnames(i))//'.info ' &
                //trim(outdir)//trim(PDFnames(i))//'_B30/'//trim(PDFnames(i))//'_B30.info' )

            call execute_command_line("mkdir -p "//trim(outdir)//trim(PDFnames(i))//'_B31')
            call execute_command_line("cp "//trim(dir)//trim(PDFnames(i))//'/'//trim(PDFnames(i))//'.info ' &
                //trim(outdir)//trim(PDFnames(i))//'_B31/'//trim(PDFnames(i))//'_B31.info' )

            call execute_command_line("mkdir -p "//trim(outdir)//trim(PDFnames(i))//'_B32')
            call execute_command_line("cp "//trim(dir)//trim(PDFnames(i))//'/'//trim(PDFnames(i))//'.info ' &
                //trim(outdir)//trim(PDFnames(i))//'_B32/'//trim(PDFnames(i))//'_B32.info' )

            call execute_command_line("mkdir -p "//trim(outdir)//trim(PDFnames(i))//'_B33')
            call execute_command_line("cp "//trim(dir)//trim(PDFnames(i))//'/'//trim(PDFnames(i))//'.info ' &
                //trim(outdir)//trim(PDFnames(i))//'_B33/'//trim(PDFnames(i))//'_B33.info' )

            call execute_command_line("mkdir -p "//trim(outdir)//trim(PDFnames(i))//'_B44')
            call execute_command_line("cp "//trim(dir)//trim(PDFnames(i))//'/'//trim(PDFnames(i))//'.info ' &
                //trim(outdir)//trim(PDFnames(i))//'_B44/'//trim(PDFnames(i))//'_B44.info' )

            call execute_command_line("mkdir -p "//trim(outdir)//trim(PDFnames(i))//'_B43')
            call execute_command_line("cp "//trim(dir)//trim(PDFnames(i))//'/'//trim(PDFnames(i))//'.info ' &
                //trim(outdir)//trim(PDFnames(i))//'_B43/'//trim(PDFnames(i))//'_B43.info' )

            call execute_command_line("mkdir -p "//trim(outdir)//trim(PDFnames(i))//'_B55')
            call execute_command_line("cp "//trim(dir)//trim(PDFnames(i))//'/'//trim(PDFnames(i))//'.info ' &
                //trim(outdir)//trim(PDFnames(i))//'_B55/'//trim(PDFnames(i))//'_B55.info' )

            call execute_command_line("mkdir -p "//trim(outdir)//trim(PDFnames(i))//'_B42')
            call execute_command_line("cp "//trim(dir)//trim(PDFnames(i))//'/'//trim(PDFnames(i))//'.info ' &
                //trim(outdir)//trim(PDFnames(i))//'_B42/'//trim(PDFnames(i))//'_B42.info' )

            call execute_command_line("mkdir -p "//trim(outdir)//trim(PDFnames(i))//'_B54')
            call execute_command_line("cp "//trim(dir)//trim(PDFnames(i))//'/'//trim(PDFnames(i))//'.info ' &
                //trim(outdir)//trim(PDFnames(i))//'_B54/'//trim(PDFnames(i))//'_B54.info' )

            call execute_command_line("mkdir -p "//trim(outdir)//trim(PDFnames(i))//'_B66')
            call execute_command_line("cp "//trim(dir)//trim(PDFnames(i))//'/'//trim(PDFnames(i))//'.info ' &
                //trim(outdir)//trim(PDFnames(i))//'_B66/'//trim(PDFnames(i))//'_B66.info' )

            ! for gg
            call execute_command_line("mkdir -p "//trim(outdir)//trim(PDFnames(i))//'_G10')
            call execute_command_line("cp "//trim(dir)//trim(PDFnames(i))//'/'//trim(PDFnames(i))//'.info ' &
                //trim(outdir)//trim(PDFnames(i))//'_G10/'//trim(PDFnames(i))//'_G10.info' )


            do j=0,numpdfmembers-1
                write (*,*) "PROCESSING ", j+1, " of ", numpdfmembers, " pdf members"
                currentpdf = pdfentry

                write (outpath, '(A,I4.4,A)') trim(outdir)//trim(PDFnames(i))//'_B00/'// &
                    trim(PDFnames(i))//'_B00_', j, '.dat'
                open(unit=15, file=outpath, status='replace', form='formatted')

                write (outpath, '(A,I4.4,A)') trim(outdir)//trim(PDFnames(i))//'_B10/'// &
                    trim(PDFnames(i))//'_B10_', j, '.dat'
                open(unit=16, file=outpath, status='replace', form='formatted')

                write (outpath, '(A,I4.4,A)') trim(outdir)//trim(PDFnames(i))//'_B11/'// &
                    trim(PDFnames(i))//'_B11_', j, '.dat'
                open(unit=17, file=outpath, status='replace', form='formatted')

                write (outpath, '(A,I4.4,A)') trim(outdir)//trim(PDFnames(i))//'_B20/'// &
                    trim(PDFnames(i))//'_B20_', j, '.dat'
                open(unit=18, file=outpath, status='replace', form='formatted')

                write (outpath, '(A,I4.4,A)') trim(outdir)//trim(PDFnames(i))//'_B21/'// &
                    trim(PDFnames(i))//'_B21_', j, '.dat'
                open(unit=19, file=outpath, status='replace', form='formatted')

                write (outpath, '(A,I4.4,A)') trim(outdir)//trim(PDFnames(i))//'_B22/'// &
                    trim(PDFnames(i))//'_B22_', j, '.dat'
                open(unit=20, file=outpath, status='replace', form='formatted')

                write (outpath, '(A,I4.4,A)') trim(outdir)//trim(PDFnames(i))//'_B30/'// &
                    trim(PDFnames(i))//'_B30_', j, '.dat'
                open(unit=21, file=outpath, status='replace', form='formatted')

                write (outpath, '(A,I4.4,A)') trim(outdir)//trim(PDFnames(i))//'_B31/'// &
                    trim(PDFnames(i))//'_B31_', j, '.dat'
                open(unit=22, file=outpath, status='replace', form='formatted')

                write (outpath, '(A,I4.4,A)') trim(outdir)//trim(PDFnames(i))//'_B32/'// &
                    trim(PDFnames(i))//'_B32_', j, '.dat'
                open(unit=23, file=outpath, status='replace', form='formatted')

                write (outpath, '(A,I4.4,A)') trim(outdir)//trim(PDFnames(i))//'_B33/'// &
                    trim(PDFnames(i))//'_B33_', j, '.dat'
                open(unit=24, file=outpath, status='replace', form='formatted')

                write (outpath, '(A,I4.4,A)') trim(outdir)//trim(PDFnames(i))//'_G10/'// &
                    trim(PDFnames(i))//'_G10_', j, '.dat'
                open(unit=25, file=outpath, status='replace', form='formatted')

                write (outpath, '(A,I4.4,A)') trim(outdir)//trim(PDFnames(i))//'_B44/'// &
                    trim(PDFnames(i))//'_B44_', j, '.dat'
                open(unit=26, file=outpath, status='replace', form='formatted')

                write (outpath, '(A,I4.4,A)') trim(outdir)//trim(PDFnames(i))//'_B43/'// &
                    trim(PDFnames(i))//'_B43_', j, '.dat'
                open(unit=27, file=outpath, status='replace', form='formatted')

                write (outpath, '(A,I4.4,A)') trim(outdir)//trim(PDFnames(i))//'_B55/'// &
                    trim(PDFnames(i))//'_B55_', j, '.dat'
                open(unit=28, file=outpath, status='replace', form='formatted')

                write (outpath, '(A,I4.4,A)') trim(outdir)//trim(PDFnames(i))//'_B42/'// &
                    trim(PDFnames(i))//'_B42_', j, '.dat'
                open(unit=29, file=outpath, status='replace', form='formatted')

                write (outpath, '(A,I4.4,A)') trim(outdir)//trim(PDFnames(i))//'_B54/'// &
                    trim(PDFnames(i))//'_B54_', j, '.dat'
                open(unit=30, file=outpath, status='replace', form='formatted')

                write (outpath, '(A,I4.4,A)') trim(outdir)//trim(PDFnames(i))//'_B66/'// &
                    trim(PDFnames(i))//'_B66_', j, '.dat'
                open(unit=31, file=outpath, status='replace', form='formatted')
            

                write (*,'(A,A,A,I3)') "Processing ", trim(PDFnames(i)), " member ", j
                write (path, '(A,I4.4,A)') trim(dir)//trim(PDFnames(i))//'/'//trim(PDFnames(i))//'_', j, '.dat'
                open(unit=11, file=trim(path), status='old', form='formatted', iostat=ierr, iomsg=imsg)

!               read(11,'(A)') pdftype
!               write(unit=15,fmt='(A)') trim(pdftype)
!               write(unit=16,fmt='(A)') trim(pdftype)
!               write(unit=17,fmt='(A)') trim(pdftype)
!               write(unit=18,fmt='(A)') trim(pdftype)
!               write(unit=19,fmt='(A)') trim(pdftype)
!               write(unit=20,fmt='(A)') trim(pdftype)
!               write(unit=21,fmt='(A)') trim(pdftype)

!               read(11,'(A)') pdfformat
!               write(unit=15,fmt='(A)') trim(pdfformat)
!               write(unit=16,fmt='(A)') trim(pdfformat)
!               write(unit=17,fmt='(A)') trim(pdfformat)
!               write(unit=18,fmt='(A)') trim(pdfformat)
!               write(unit=19,fmt='(A)') trim(pdfformat)
!               write(unit=20,fmt='(A)') trim(pdfformat)
!               write(unit=21,fmt='(A)') trim(pdfformat)

!               read(11,'(A)') dummy
!               write(unit=15,fmt='(A)') trim(dummy)
!               write(unit=16,fmt='(A)') trim(dummy)
!               write(unit=17,fmt='(A)') trim(dummy)
!               write(unit=18,fmt='(A)') trim(dummy)
!               write(unit=19,fmt='(A)') trim(dummy)
!               write(unit=20,fmt='(A)') trim(dummy)
!               write(unit=21,fmt='(A)') trim(dummy)

                ierr = 0
                do while (ierr == 0)
                    read(11, '(A)') nextline
                    do bb=15,15+NUM_BEAMFUNS-1
                        write(unit=bb,fmt='(A)') trim(nextline)
                    enddo
                    if (trim(nextline) == "---") then
                        ierr = 1
                    endif
                enddo

                ierr = 0
                read(11,'(A)', iostat=ierr) nextline
                do while (ierr == 0)

                    xvals = huge(0._dp)
                    read(nextline, *, iostat=iread) xvals

                    read(11,'(A)') nextline
                    q2vals = huge(0._dp)
                    read(nextline, *, iostat=iread) q2vals

                    read(11,'(A)') nextline
                    flavors = huge(0)
                    read(nextline, *, iostat=iread) flavors

                    numx = count(xvals /= huge(0._dp))
                    numq2 = count(q2vals /= huge(0._dp))
                    numflavors = count(flavors /= huge(0))

                    write (*,'(A,I3,A,I3,A,I3)') "Processing block with number of x,q2,flavors ", numx, " ", numq2, " ", numflavors

                    do k=1,numx
                        do bb=15,15+NUM_BEAMFUNS-1
                            write(unit=bb, fmt='(ES16.8)', advance="no") xvals(k)
                        enddo
                    enddo
                    
                    do bb=15,15+NUM_BEAMFUNS-1
                        write (unit=bb, fmt='(A)') ''
                    enddo

                    do k=1,numq2
                        do bb=15,15+NUM_BEAMFUNS-1
                            write(unit=bb, fmt='(ES16.8)', advance="no") q2vals(k)
                        enddo
                    enddo

                    do bb=15,15+NUM_BEAMFUNS-1
                        write (unit=bb, fmt='(A)') ''

                    enddo

                    do k=1,numflavors
                        do bb=15,15+NUM_BEAMFUNS-1
                            write(unit=bb, fmt='(I3)', advance="no") flavors(k)
                        enddo
                    enddo

                    ! m: flavors
                    ! l: q2
                    ! k: x
                    ! beam 1:6, m, l, k
#ifndef HAVE_COARRAY
                    if (allocated(pdfvals)) deallocate(pdfvals)
                    allocate(pdfvals(NUM_BEAMFUNS,numflavors,numq2,numx), source=huge(0._dp))
#endif

! with MPI we distribute x values to workers and fill pdfvals
#ifdef HAVE_MPI
                    if (allocated(pdfvals_slice)) deallocate(pdfvals_slice)
                    allocate(pdfvals_slice(NUM_BEAMFUNS,numflavors,numq2), source=huge(0._dp))

                    working(:) = .false.

                    nextx = 1
                    recvx = 0
                    do while (recvx < numx)
                        do n=1,world_size-1
                            if ((.not. working(n)) .and. (nextx <= numx)) then
                                wx = xvals(nextx)
                                wxid = nextx
                                nextx = nextx + 1

                                ! push wx, wxid and other stuff to worker
                                call mpi_send(wx, 1, mpi_double_precision, n, 1, mpi_comm_world, ierr)
                                call mpi_send(wxid, 1, mpi_integer, n, 2, mpi_comm_world, ierr)
                                call mpi_send(numq2, 1, mpi_integer, n, 3, mpi_comm_world, ierr)
                                call mpi_send(numflavors, 1, mpi_integer, n, 4, mpi_comm_world, ierr)
                                call mpi_send(q2vals, size(q2vals), mpi_double_precision, n, 5, mpi_comm_world, ierr)
                                call mpi_send(flavors, size(flavors), mpi_integer, n, 6, mpi_comm_world, ierr)
                                call mpi_send(currentPDF, 1, mpi_integer, n, NUM_BEAMFUNS, mpi_comm_world, ierr)

                                ! skip these in the input file
                                do l=1,numq2
                                    read(11,'(A)', iostat=iread) nextline
                                enddo

                                working(n) = .true.
                            endif
                        enddo

                        ! receive one finished task
                        mpistatus = 0
                        call mpi_recv(wxid, 1, mpi_integer, mpi_any_source, 1, mpi_comm_world, mpistatus, ierr)

                        !if ((j < numpdfmembers-1) .and. (i < numPDFsets) .and. (recvx+1 < numx)) then
                            morework = .true.
                        !else
                            !morework = .false.
                        !endif

                        ! signal worker that it should send stuff
                        call mpi_send(morework, 1, mpi_logical, mpistatus(mpi_source), 0, mpi_comm_world, ierr)

                        write (*,'(A,ES16.8,A,I3)') "MASTER: Receiving x = ", xvals(wxid), " from rank ", mpistatus(mpi_source)
                        call mpi_recv(pdfvals_slice, NUM_BEAMFUNS*numflavors*numq2, &
                            mpi_double_precision, mpistatus(mpi_source), 2, mpi_comm_world, mpistatus, ierr)
                        pdfvals(:,:,:,wxid) = pdfvals_slice(:,:,:)
                        working(mpistatus(mpi_source)) = .false.
                        recvx = recvx + 1
                    enddo

#elif HAVE_COARRAY
                    do nimg=2,num_images()
                        numx[nimg] = numx
                        numq2[nimg] = numq2
                        numflavors[nimg] = numflavors
                        flavors(1:numflavors)[nimg] = flavors(1:numflavors)
                        q2vals(1:numq2)[nimg] = q2vals(1:numq2)
                        xvals(1:numx)[nimg] = xvals(1:numx)
                        curPDF[nimg] = currentPDF
                    enddo
                    sync all
                    if (allocated(pdfvals)) deallocate(pdfvals)
                    allocate(pdfvals(NUM_BEAMFUNS,numflavors,numq2,numx)[*], source=huge(0._dp))
                    if (allocated(work_done)) deallocate(work_done)
                    allocate(work_done(numx)[*], source=.false.)
                    sync all

                    do k=1,numx
                        do l=1,numq2
                            read(11,'(A)', iostat=iread) nextline
                        enddo
                    enddo

! without MPI or Coarray we calculate everything local
#else
                    do k=1,numx
                        do l=1,numq2
                            read(11,'(A)', iostat=iread) nextline
                        enddo
                    enddo

                    do k=1,numx
                        write (*,'(A,ES16.8)') "Processing x = ", xvals(k)
!$omp parallel do collapse(2)
                        do l=1,numq2

                            ! original grid values
                            !read(nextline, *, iostat=iread) pdfvals(:,l,1)

                            ! compute beam for current x and q2
                            do m=1,numflavors
                                fl = flavors(m)
                                if (fl == 21) fl = 0

                                ! we hardcode ih=1 and beam=1
                                ! this means no pdfchannel selection should be done in the input file
                                ! for this grid generation
                                pdfvals(1,m,l,k) = getbeam(1, fl, 0, 0, xvals(k), q2vals(l))*xvals(k)
                                pdfvals(2,m,l,k) = getbeam(1, fl, 1, 0, xvals(k), q2vals(l))*xvals(k)
                                pdfvals(3,m,l,k) = getbeam(1, fl, 1, 1, xvals(k), q2vals(l))*xvals(k)
                                pdfvals(4,m,l,k) = getbeam(1, fl, 2, 0, xvals(k), q2vals(l))*xvals(k)
                                pdfvals(5,m,l,k) = getbeam(1, fl, 2, 1, xvals(k), q2vals(l))*xvals(k)
                                pdfvals(6,m,l,k) = getbeam(1, fl, 2, 2, xvals(k), q2vals(l))*xvals(k)

                                pdfvals(7,m,l,k) = getbeam(1, fl, 3, 0, xvals(k), q2vals(l))*xvals(k)
                                pdfvals(8,m,l,k) = getbeam(1, fl, 3, 1, xvals(k), q2vals(l))*xvals(k)
                                pdfvals(9,m,l,k) = getbeam(1, fl, 3, 2, xvals(k), q2vals(l))*xvals(k)
                                pdfvals(10,m,l,k) = getbeam(1, fl, 3, 3, xvals(k), q2vals(l))*xvals(k)

                                pdfvals(11,m,l,k) = getbeam2(1, fl, 1, 0, xvals(k), q2vals(l))*xvals(k)

                                pdfvals(12,m,l,k) = getbeam(1, fl, 4, 4, xvals(k), q2vals(l))*xvals(k)

                                pdfvals(13,m,l,k) = getbeam(1, fl, 4, 3, xvals(k), q2vals(l))*xvals(k)
                                pdfvals(14,m,l,k) = getbeam(1, fl, 5, 5, xvals(k), q2vals(l))*xvals(k)

                                pdfvals(15,m,l,k) = getbeam(1, fl, 4, 2, xvals(k), q2vals(l))*xvals(k)
                                pdfvals(16,m,l,k) = getbeam(1, fl, 5, 4, xvals(k), q2vals(l))*xvals(k)
                                pdfvals(17,m,l,k) = getbeam(1, fl, 6, 6, xvals(k), q2vals(l))*xvals(k)

                                do bb=1,NUM_BEAMFUNS
                                    if (ieee_is_nan(pdfvals(bb,m,l,k))) pdfvals(bb,m,l,k) = 0._dp

                                enddo

                            enddo
                        enddo
!$omp end parallel do
                    enddo
#endif

                    do k=1,numx
                        do l=1,numq2
                            do bb=15,15+NUM_BEAMFUNS-1
                                write(unit=bb, fmt='(A)') ''
                            enddo

                            do m=1,numflavors
                                do bb=15,15+NUM_BEAMFUNS-1
                                    write (unit=bb, fmt='(ES16.8)', advance="no") pdfvals(bb-14,m,l,k)
                                enddo
                            enddo
                        enddo
                    enddo

                    read(11,'(A)', iostat=ierr) dummy
                    do bb=15,15+NUM_BEAMFUNS-1
                        write (unit=bb, fmt='(A)') ''
                        write (unit=bb, fmt='(A)') '---'
                    enddo

                    read(11,'(A)', iostat=ierr) nextline
                end do

                close(unit=11)

                do bb=15,15+NUM_BEAMFUNS-1
                    close(unit=bb)
                enddo

                pdfentry = pdfentry + 1

            enddo
        enddo

        currentpdf = 0

        write (*,*) "Grid generation is complete."
        write (*,*) "Please copy the generated PDF sets to your $LHAPDF_DATA_PATH"
        write (*,*) "Then set makegrid to .false."
        write (*,*) ""
        write (*,*) "Warning: Please remove any ForcePositive: 1 settings"
        write (*,*) "in the LHAPDF .info files. Most prominently this affects"
        write (*,*) "CT14nnlo."
        stop 0

    end subroutine


end module

