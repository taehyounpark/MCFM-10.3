!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module MCFMPlotting
      use Superhisto
      implicit none

      public :: plot_setup_uniform
      public :: plot_setup_custom
      public :: plots_allocate
      public :: plot_book

      private

      integer, save :: numplots = 0
      type(sh_histogram), save, allocatable :: histobuf(:)

      contains

      ! called once to initialize individual plots
      function plot_setup_uniform(xmin,xmax,dx,title)
          use types
          implicit none

          integer :: plot_setup_uniform
          real(dp), intent(in) :: xmin, xmax, dx
          character(len=*), intent(in) :: title

          include 'mpicommon.f'

          if (rank == 0) then
              write (*,'(A,A,A)') "Uniform histogram initialized for '", title, "'"
          endif

          plot_setup_uniform = plot_setup_internal()
          call histobuf(plot_setup_uniform)%init(title,xmin,xmax,dx)
      end function 

      function plot_setup_custom(xarr,title)
          use types
          implicit none

          integer :: plot_setup_custom
          real(dp), intent(in) :: xarr(:)
          character(len=*), intent(in) :: title

          include 'mpicommon.f'

          if (rank == 0) then
              write (*,'(A,A,A)') "Custom histogram initialized for '", title, "'"
          endif

          plot_setup_custom = plot_setup_internal()
          call histobuf(plot_setup_custom)%init_custom(title,xarr)
      end function 


      function plot_setup_internal()
          use types
          implicit none

          ! to make threadStorageOp, iteration accumulation etc. work
          ! eventually will be adjusted
          integer :: nplotmax
          common/nplotmax/nplotmax

          integer :: plot_setup_internal
          type(sh_histogram), allocatable :: newhistos(:)

          if (allocated(histobuf) .and. (numplots == 0)) then
              write (*,*) __FILE__//" ERROR: Histograms already allocated but first plot"
              error stop
          endif

          numplots = numplots + 1
          plot_setup_internal = numplots
          nplotmax = numplots

          if (.not. allocated(histobuf)) then
              allocate(histobuf(numplots))
          else
              allocate(newhistos(numplots))
              newhistos(1:numplots-1) = histobuf
              call move_alloc(newhistos, histobuf)
          endif

      end function

      subroutine plots_allocate()
          use MCFMStorage
          use Scalevar
          use PDFerrors
          use SCET
          use Integration
          implicit none

          integer :: j,k,m

          ! for current currentPart, and currentIPS
          ! perform histogram allocation, if not already done so

          j = currentPart
          k = currentIps

!$omp master
          do j=1,maxParts; do k=1,maxIPS
              allocate(masterStorage(j,k)%histCentral%histos, source=histobuf)
              allocate(iterationStorage(j,k)%histCentral%histos, source=histobuf)

              if (doScalevar) then
                  allocate(iterationStorage(j,k)%histScalevar(maxscalevar+extrascalevar))
                  do m=1,maxscalevar+extrascalevar
                      allocate(iterationStorage(j,k)%histScalevar(m)%histos, source=histobuf)
                  enddo

                  allocate(masterStorage(j,k)%histScalevar(maxscalevar+extrascalevar))
                  do m=1,maxscalevar+extrascalevar
                      allocate(masterStorage(j,k)%histScalevar(m)%histos, source=histobuf)
                  enddo
              endif

              if (maxPDFsets > 0) then
                  allocate(iterationStorage(j,k)%histPDFerrors(maxPDFsets))
                  do m=1,maxPDFsets
                      allocate(iterationStorage(j,k)%histPDFerrors(m)%histos, source=histobuf)
                  enddo

                  allocate(masterStorage(j,k)%histPDFerrors(maxPDFsets))
                  do m=1,maxPDFsets
                      allocate(masterStorage(j,k)%histPDFerrors(m)%histos, source=histobuf)
                  enddo
              endif

              allocate(iterationStorage(j,k)%histTaucut(size(tcutarray)))
              do m=1,size(tcutarray)
                  allocate(iterationStorage(j,k)%histTaucut(m)%histos, source=histobuf)
              enddo
              allocate(masterStorage(j,k)%histTaucut(size(tcutarray)))
              do m=1,size(tcutarray)
                  allocate(masterStorage(j,k)%histTaucut(m)%histos, source=histobuf)
              enddo
          enddo; enddo
!$omp end master

!$omp parallel
          do j=1,maxParts; do k=1,maxIPS

              allocate(threadStorage(j,k)%histCentral%histos, source=histobuf)

              if (doScalevar) then
                  allocate(threadStorage(j,k)%histScalevar(maxscalevar+extrascalevar))
                  do m=1,maxscalevar+extrascalevar
                      allocate(threadStorage(j,k)%histScalevar(m)%histos, source=histobuf)
                  enddo
              endif

              if (maxPDFsets > 0) then
                  allocate(threadStorage(j,k)%histPDFerrors(maxPDFsets))
                  do m=1,maxPDFsets
                      allocate(threadStorage(j,k)%histPDFerrors(m)%histos, source=histobuf)
                  enddo
              endif

              allocate(threadStorage(j,k)%histTaucut(size(tcutarray)))
              do m=1,size(tcutarray)
                  allocate(threadStorage(j,k)%histTaucut(m)%histos, source=histobuf)
              enddo
          enddo; enddo
!$omp end parallel

      end subroutine

      ! replaces bookplot
      subroutine plot_book(ids, vals, wt0, wts)
          use types
          use MCFMStorage
          use Scalevar
          use PDFerrors
          use SCET
          use Integration
          implicit none
          include 'kpart.f'
          include 'nproc.f'

          integer, intent(in) :: ids(:)
          real(dp), intent(in) :: vals(:)
          real(dp), intent(in) :: wt0
          real(dp), intent(in) :: wts(:)

          real(dp) :: nomTau
          integer :: i,n,nplots

          nplots = size(ids)

          ! for 161 and 165 we always call tmpbook, since these
          ! modified routines always call threadStorageTmpCommit to properly book everything
          ! once all partial contributions have been added

          if (.not. ((maxPDFsets > 0) .and. kpart == kreal .and. currentPDF > 0)) then
              ! whether or not to actually include the weight for the nominal taucut
              if (includeTaucutgrid(currentNd)) then
                  if (kpart == kreal .or. nproc == 1610 .or. nproc == 1650) then
                      do n=1,nplots
                          call threadStorage(currentPart,currentIps)%histCentral%histos(n)%tmpbook(vals(n),wts(n))
                      enddo
                  else
                      do n=1,nplots
                          call threadStorage(currentPart,currentIps)%histCentral%histos(n)%book(vals(n),wts(n))
                      enddo
                  endif

                  if (doScalevar) then
                      if (kpart == kreal .or. nproc == 1610 .or. nproc == 1650) then
                          do i=1,maxscalevar+extrascalevar
                              do n=1,nplots
                                  call threadStorage(currentPart,currentIps)%histScalevar(i)% &
                                      histos(n)%tmpbook(vals(n),wts(n) * (scalereweight(i) - 1._dp))
                              enddo
                          enddo
                      else
                          do i=1,maxscalevar+extrascalevar
                              do n=1,nplots
                                  call threadStorage(currentPart,currentIps)%histScalevar(i)% &
                                      histos(n)%book(vals(n),wts(n) * (scalereweight(i) - 1._dp))
                              enddo
                          enddo
                      endif
                  endif
              endif

              if (size(tcutarray) > 0) then
                  if (includeTaucutgrid(currentNd)) then
                      nomTau = 1._dp
                  else
                      nomTau = 0._dp
                  endif

                   if (kpart == kreal .or. nproc == 1610 .or. nproc == 1650) then
                       do i=1,size(tcutarray)
                           do n=1,nplots
                               call threadStorage(currentPart,currentIps)%histTaucut(i)% &
                                   histos(n)%tmpbook(vals(n), wts(n) * (nomTau - scetreweight(i)))
                           enddo
                       enddo
                   else
                       do i=1,size(tcutarray)
                           do n=1,nplots
                               call threadStorage(currentPart,currentIps)%histTaucut(i)% &
                                   histos(n)%book(vals(n), wts(n) * (nomTau - scetreweight(i)))
                           enddo
                       enddo
                   endif

              endif

          endif ! exclusion of pdferrors with real and currentPDF > 0

          if (nproc == 1610 .or. nproc == 1650) then
              ! new style pdf uncertainties (singletop realint)
              if ((maxPDFsets > 0) .and. kpart == kreal .and. includeTaucutgrid(currentNd)) then
                  do i=1,maxPDFsets
                      do n=1,nplots
                          if (wt0 /= 0._dp) then
                          call threadStorage(currentPart,currentIps)%histPDFerrors(i)% &
                              histos(n)%tmpbook(vals(n), wts(n)/wt0 * pdfreweight(i))
                          else
                              call threadStorage(currentPart,currentIps)%histPDFerrors(i)% &
                                  histos(n)%tmpbook(vals(n), pdfreweight(i))
                          endif
                      enddo
                  enddo
              endif
          else
              ! old style pdf uncertainties (general realint)
              if ((maxPDFsets > 0) .and. kpart == kreal .and. currentPDF > 0 .and. includeTaucutgrid(currentNd)) then
                  do n=1,nplots
                      if (wt0 /= 0._dp) then
                      call threadStorage(currentPart,currentIps)%histPDFerrors(currentPDF)% &
                          histos(n)%tmpbook(vals(n), wts(n)/wt0 * pdfreweight(currentPDF))
                      else
                          call threadStorage(currentPart,currentIps)%histPDFerrors(currentPDF)% &
                              histos(n)%tmpbook(vals(n), pdfreweight(currentPDF))
                      endif
                  enddo
              endif
          endif


          if ((maxPDFsets > 0) .and. kpart /= kreal .and. includeTaucutgrid(0)) then
              if (nproc == 1610 .or. nproc == 1650) then
                  do i=1,maxPDFsets
                      do n=1,nplots
                          if (wt0 /= 0._dp) then
                          call threadStorage(currentPart,currentIps)%histPDFerrors(i)% &
                              histos(n)%tmpbook(vals(n), wts(n)/wt0 * pdfreweight(i))
                          else
                              call threadStorage(currentPart,currentIps)%histPDFerrors(i)% &
                                  histos(n)%tmpbook(vals(n), pdfreweight(i))
                          endif
                      enddo
                  enddo
              else
                  do i=1,maxPDFsets
                      do n=1,nplots
                          if (wt0 /= 0._dp) then
                          call threadStorage(currentPart,currentIps)%histPDFerrors(i)% &
                              histos(n)%book(vals(n), wts(n)/wt0 * pdfreweight(i))
                          else
                              call threadStorage(currentPart,currentIps)%histPDFerrors(i)% &
                                  histos(n)%book(vals(n), pdfreweight(i))
                          endif
                      enddo
                  enddo
              endif
          endif

      end subroutine


end module
