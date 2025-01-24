!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine mcfm_exit(xinteg,xinteg_err)
          use omp_lib
          use types
          use MCFMStorage
          use MCFMPrint
          use SCET, only : doMultitaucut
          use PDFerrors, only : doPDFerrors
          use Scalevar, only : doScalevar
          use SingletopPrint, only: finalizeStorageFixedOrder,
     &        printcross_singletop, writeAllHistogramsSingletop
      implicit none
      include 'mpicommon.f'
      include 'nproc.f'

      real(dp), intent(in) :: xinteg, xinteg_err

      if (rank == 0) then
          call finalizeStorage

          if (nprocbelow == 1610) then
              call finalizeStorageFixedOrder
          endif

          if (nprocbelow == 1610) then
              call printcross_singletop()
          endif

          call printallcross()

          if (doPDFerrors) then
            call printPDFuncertainties()
          endif

          if (doPDFerrors) then
            call printPDFuncertainties()
          endif

          if (doScalevar) then
              call printScaleuncertainties()
          endif

          if (doMultitaucut) then
              call printTaucuts()
          endif

          if (nprocbelow == 1610) then
              call writeAllHistogramsSingletop()
              !call writeAllHistograms()
          else
          call writeAllHistograms()
          endif

      endif

#ifdef HAVE_COARRAY
      sync all
#endif

      end
