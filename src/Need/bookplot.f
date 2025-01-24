!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine bookplot(n,tag,titlex,var,wt,wt2,xmin,xmax,dx,llplot)
          use types
          use Superhisto
          use PDFerrors
          use Scalevar
          use SCET
          use MCFMStorage
          use omp_lib
          implicit none
      include 'nplot.f'
      include 'kpart.f'
      include 'vegas_common.f'
      include 'scalevar.f'
      include 'nproc.f'

      integer, intent(in) :: n,tag
      character(len=3), intent(in) ::  llplot
      character*(*), intent(in) ::  titlex
      real(dp), intent(in) :: var, wt, wt2, xmin, xmax, dx
      real(dp) :: nomTau

      integer, parameter:: tagbook=1, tagplot=2
      integer :: i

      character(len=*), parameter, dimension(8) :: scalevarlabels =
     &      [": scale uu", ": scale dd",
     &       ": scale u-", ": scale d-",
     &       ": scale -u", ": scale -d",
     &       ": scale ud", ": scale du"]

      if (storageAllocated) then
        ! initialize thread local histograms
        associate (chist => threadStorage(currentPart,currentIps)%histCentral%histos(n))
            ! use central histogram as marker that central
            ! as well as scalevariation and pdf error histograms
            ! have been initialized
            if (shinitialized(chist) .eqv. .false.) then
                !write (*,*) "thread storage initialization ", trim(titlex), n, "on thread", omp_get_thread_num()
                call chist%init(trim(titlex), xmin, xmax, dx)
                threadStorage(currentPart,currentIps)%used = .true.

                if (doScalevar) then
                    do i=1,maxscalevar+extrascalevar
                        associate (shist => threadStorage(currentPart,currentIps)%histScalevar(i)%histos(n))
                            if (shinitialized(shist) .eqv. .false.) then
                                call shist%init(trim(titlex)//scalevarlabels(i), xmin, xmax, dx)
                                threadStorage(currentPart,currentIps)%used = .true.
c                               write (*,*) "thread storage scalevar", i, "initialized ", n,
c    &                                      "on thread", omp_get_thread_num()
                            endif
                        end associate
                    enddo
                endif

                if (maxPDFsets > 0) then
                    do i=1,maxPDFsets
                        associate (phist => threadStorage(currentPart,currentIps)%histPDFerrors(i)%histos(n))
                            if (shinitialized(phist) .eqv. .false.) then
                                call phist%init(trim(titlex), xmin, xmax, dx)
                                threadStorage(currentPart,currentIps)%used = .true.
c                               write (*,*) "thread storage PDFerrors", i, "initialized ", n,
c    &                                      "on thread", omp_get_thread_num()
                            endif
                        end associate
                    enddo
                endif

                do i=1,size(tcutarray)
                    associate (thist => threadStorage(currentPart,currentIps)%histTaucut(i)%histos(n))
                        if (shinitialized(thist) .eqv. .false.) then
                            call thist%init(trim(titlex), xmin, xmax, dx)
                            threadStorage(currentPart,currentIps)%used = .true.
                        endif
                    end associate
                enddo


            endif
        end associate

        ! initialize accumulator master histogram
        ! and iteration accumulator histograms

!$omp master
        ! use histCentral as marker that everything has been properly allocated
        if (shinitialized(masterStorage(currentPart,currentIps)%histCentral%histos(n)) .eqv. .false.) then

          masterStorage(currentPart,currentIps)%histCentral%histos(n) =
     &          threadStorage(currentPart,currentIps)%histCentral%histos(n)
          masterStorage(currentPart,currentIps)%used = .true.

          if (doScalevar) then
              do i=1,maxscalevar+extrascalevar
                  masterStorage(currentPart,currentIps)%histScalevar(i)%histos(n) =
     &              threadStorage(currentPart,currentIps)%histScalevar(i)%histos(n)
              enddo
          endif

          if (maxPDFsets > 0) then
              do i=1,maxPDFsets
                  masterStorage(currentPart,currentIps)%histPDFerrors(i)%histos(n) =
     &              threadStorage(currentPart,currentIps)%histPDFerrors(i)%histos(n)
              enddo
          endif

          do i=1,size(tcutarray)
              masterStorage(currentPart,currentIps)%histTaucut(i)%histos(n) =
     &          threadStorage(currentPart,currentIps)%histTaucut(i)%histos(n)
          enddo

        endif

        if (shinitialized(iterationStorage(currentPart,currentIps)%histCentral%histos(n)) .eqv. .false.) then
          iterationStorage(currentPart,currentIps)%histCentral%histos(n) =
     &          threadStorage(currentPart,currentIps)%histCentral%histos(n)
          iterationStorage(currentPart,currentIps)%used = .true.

          if (doScalevar) then
              do i=1,maxscalevar+extrascalevar
                  iterationStorage(currentPart,currentIps)%histScalevar(i)%histos(n) =
     &              threadStorage(currentPart,currentIps)%histScalevar(i)%histos(n)
              enddo
          endif

          if (maxPDFsets > 0) then
              do i=1,maxPDFsets
                  iterationStorage(currentPart,currentIps)%histPDFerrors(i)%histos(n) =
     &              threadStorage(currentPart,currentIps)%histPDFerrors(i)%histos(n)
              enddo
          endif

          do i=1,size(tcutarray)
              iterationStorage(currentPart,currentIps)%histTaucut(i)%histos(n) =
     &          threadStorage(currentPart,currentIps)%histTaucut(i)%histos(n)
          enddo
        endif

!$omp end master
      endif

c do the booking

      if (tag == tagplot) then
          ! For the real emission bookplot is called for each dipole contribution.
          ! Since we loop over the PDFs, we only book once for the central value
          ! and for the scale varied results, but for currentPDF > 0 below.

          if (.not. ((maxPDFsets > 0) .and. kpart == kreal .and. currentPDF > 0)) then

              ! whether or not to actually include the weight for the nominal taucut
              if (includeTaucutgrid(currentNd)) then
                  associate (chist => threadStorage(currentPart,currentIps)%histCentral%histos(n) )
                      if (kpart == kreal .or. nproc == 1610 .or. nproc == 1650) then
                          ! realint will call threadStorageTmpCommit to properly book it
                          ! once all dipole contributions have been added
                          call chist%tmpbook(var,wt)
                      else
                          call chist%book(var,wt)
                      endif
                  end associate

                  if (doScalevar) then
                      do i=1,maxscalevar+extrascalevar
                          associate (shist => threadStorage(currentPart,currentIps)%histScalevar(i)%histos(n) )
                              if (kpart == kreal .or. nproc == 1610 .or. nproc == 1650) then
                                  call shist%tmpbook(var,wt * (scalereweight(i) - 1._dp))
                              else
                                  call shist%book(var,wt * (scalereweight(i) - 1._dp))
                              endif
                          end associate
                      enddo
                  endif
              endif

              if (size(tcutarray) > 0) then
                  if (includeTaucutgrid(currentNd)) then
                      nomTau = 1._dp
                  else
                      nomTau = 0._dp
                  endif
                  do i=1,size(tcutarray)
                      associate(thist => threadStorage(currentPart,currentIps)%histTaucut(i)%histos(n))
                          if (kpart == kreal .or. nproc == 1610 .or. nproc == 1650) then
                              call thist%tmpbook(var, wt * (nomTau - scetreweight(i)))
                          else
                              call thist%book(var, wt * (nomTau - scetreweight(i)))
                          endif
                      end associate
                  enddo
              endif

          endif ! exclusion of pdferrors with real and currentPDF > 0

          if (nproc == 1610 .or. nproc == 1650) then
              ! new style pdf uncertainties (singletop realint)
              if ((maxPDFsets > 0) .and. kpart == kreal .and. includeTaucutgrid(currentNd)) then
                  do i=1,maxPDFsets
                      associate (phist => threadStorage(currentPart,currentIps)%histPDFerrors(i)%histos(n) )
                          call phist%tmpbook(var, pdfreweight(i))
                      end associate
                  enddo
              endif
          else
              ! old style pdf uncertainties (general realint)
              if ((maxPDFsets > 0) .and. kpart == kreal .and. currentPDF > 0 .and. includeTaucutgrid(currentNd)) then
                  associate (phist => threadStorage(currentPart,currentIps)%histPDFerrors(currentPDF)%histos(n) )
                      call phist%tmpbook(var, pdfreweight(currentPDF))
                  end associate
              endif
          endif


          if ((maxPDFsets > 0) .and. kpart /= kreal .and. includeTaucutgrid(0)) then
              do i=1,maxPDFsets
                  associate (phist => threadStorage(currentPart,currentIps)%histPDFerrors(i)%histos(n) )
                      if (nproc == 1610 .or. nproc == 1650) then
                          call phist%tmpbook(var, pdfreweight(i))
                      else
                          call phist%book(var, pdfreweight(i))
                      endif
                  end associate
              enddo
          endif

      endif

      end subroutine
