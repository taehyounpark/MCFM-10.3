!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
module LHAPDF
      use iso_c_binding
      use types
      implicit none

      type pdf
          private
          character(:), allocatable :: setname
          type(c_ptr) :: ptr
          contains
          procedure :: xfxq2 => pdf_xfxq2
          procedure :: nmem => lhapdf_number_pdf
          procedure :: alphas => pdf_alphas
          procedure :: numflavors => pdf_numflavors
          procedure :: quarkThreshold => pdf_quarkThreshold
          procedure :: quarkMass => pdf_quarkMass
          procedure :: orderqcd => pdf_orderqcd
      end type

      type, bind(c) :: pdfuncertainty
          real(c_double) :: central
          real(c_double) :: errplus
          real(c_double) :: errminus
          real(c_double) :: errsymm
          real(c_double) :: scale
          real(c_double) :: errplus_pdf
          real(c_double) :: errminus_pdf
          real(c_double) :: errsymm_pdf
          real(c_double) :: err_par
      end type

      public :: pdf
      public :: pdf_xfxq2
      public :: lhapdf_loadMember
      public :: lhapdf_info
      public :: lhapdf_number
      public :: pdfuncertainty
      public :: pdf_readgrid

      ! these routines work on the global pdfs
      public :: fdist
      public :: fdist_one
      public :: initcentral, initall
      public :: getalphas, getorderqcd
      public :: getnumflavors
      public :: getquarkThreshold, getquarkMass
      public :: computeUncertainty

      public :: lhapdf_pathsPrepend

      public :: initcentralResummation
      public :: initallResummation
      public :: fdist_one_beam
      public :: fdist_one_beam2

      interface lhapdf_number
          module procedure lhapdf_number_pdf
          module procedure lhapdf_number_name
      end interface

      private

      type(pdf), allocatable, save :: pdfs(:)

      type(pdf), allocatable, save :: respdfs(:,:)

      interface
          subroutine lhapdf_info() bind(C, name="lhapdf_info")
          end subroutine

          subroutine c_lhapdf_pathsPrepend(path) bind(C, name="lhapdf_pathsPrepend")
              use iso_c_binding
              implicit none
              character(kind=c_char), dimension(*), intent(in) :: path
          end subroutine

          function c_lhapdf_loadMember(setname, member) bind(C, name="lhapdf_loadMember")
              use iso_c_binding
              implicit none
              type(c_ptr) :: c_lhapdf_loadMember 
              character(kind=c_char), dimension(*), intent(in) :: setname
              integer(c_int), value, intent(in) :: member
          end function

          function c_lhapdf_evolve(pdfptr, x, q2, flav) bind(C, name="lhapdf_evolve")
              use iso_c_binding
              implicit none
              real(c_double) :: c_lhapdf_evolve
              type(c_ptr), value, intent(in) :: pdfptr
              real(c_double), value, intent(in) :: x, q2
              integer(c_int), value, intent(in) :: flav
          end function

          function c_lhapdf_number(setname) bind(C, name="lhapdf_number")
              use iso_c_binding
              implicit none
              integer(c_int) :: c_lhapdf_number
              character(kind=c_char), dimension(*), intent(in) :: setname
          end function

          function c_lhapdf_alphas(pdfptr, q2) bind(C, name="lhapdf_alphas")
              use iso_c_binding
              implicit none
              real(c_double) :: c_lhapdf_alphas
              type(c_ptr), value, intent(in) :: pdfptr
              real(c_double), value, intent(in) :: q2
          end function

          function c_lhapdf_numflavors(pdfptr, q) bind(C, name="lhapdf_numFlavors")
              use iso_c_binding
              implicit none
              integer(c_int) :: c_lhapdf_numflavors
              type(c_ptr), value, intent(in) :: pdfptr
              real(c_double), value, intent(in) :: q
          end function

          function c_lhapdf_quarkThreshold(pdfptr, i) bind(C, name="lhapdf_quarkThreshold")
              use iso_c_binding
              implicit none
              real(c_double) :: c_lhapdf_quarkThreshold
              type(c_ptr), value, intent(in) :: pdfptr
              integer(c_int), value, intent(in) :: i
          end function

          function c_lhapdf_quarkMass(pdfptr, i) bind(C, name="lhapdf_quarkMass")
              use iso_c_binding
              implicit none
              real(c_double) :: c_lhapdf_quarkMass
              type(c_ptr), value, intent(in) :: pdfptr
              integer(c_int), value, intent(in) :: i
          end function

          function c_lhapdf_orderqcd(pdfptr) bind(C, name="lhapdf_orderqcd")
              use iso_c_binding
              implicit none
              integer(c_int) :: c_lhapdf_orderqcd
              type(c_ptr), value, intent(in) :: pdfptr
          end function

          function c_lhapdf_computeUncertainty(pdfptr, array, num) bind(C, name="lhapdf_computeUncertainty")
              use iso_c_binding
              import :: pdfuncertainty
              implicit none
              type(pdfuncertainty) :: c_lhapdf_computeUncertainty
              type(c_ptr), value, intent(in) :: pdfptr
              type(c_ptr), value, intent(in) :: array
              integer, value, intent(in) :: num
          end function

          subroutine c_lhapdf_getconfig(pdfptr, key, outvalue, olen) bind(C, name="lhapdf_getconfig")
              use iso_c_binding
              implicit none
              type(c_ptr), value, intent(in) :: pdfptr
              character(kind=c_char), dimension(*), intent(in) :: key
              character(kind=c_char), dimension(*), intent(inout) :: outvalue
              integer(c_int), value, intent(in) :: olen
          end subroutine

        function c_lhapdf_numxKnots(pdfptr) bind(C, name="lhapdf_numxKnots")
            use iso_c_binding
            implicit none
            type(c_ptr), value, intent(in) :: pdfptr
            integer(c_int) :: c_lhapdf_numxKnots
        end function

        function c_lhapdf_numq2Knots(pdfptr) bind(C, name="lhapdf_numq2Knots")
            use iso_c_binding
            implicit none
            type(c_ptr), value, intent(in) :: pdfptr
            integer(c_int) :: c_lhapdf_numq2Knots
        end function

        subroutine c_lhapdf_getxKnots(pdfptr, xknots) bind(C, name="lhapdf_getxKnots")
            use iso_c_binding
            implicit none
            type(c_ptr), value, intent(in) :: pdfptr
            type(c_ptr), value, intent(in) :: xknots
        end subroutine

        subroutine c_lhapdf_getq2Knots(pdfptr, q2knots) bind(C, name="lhapdf_getq2Knots")
            use iso_c_binding
            implicit none
            type(c_ptr), value, intent(in) :: pdfptr
            type(c_ptr), value, intent(in) :: q2knots
        end subroutine

      end interface

      contains

      subroutine lhapdf_getconfig(pdf_in, key, outvalue)
          implicit none
          class(pdf), intent(in) :: pdf_in
          character(len=*), intent(in) :: key
          character(len=*), intent(inout) :: outvalue

          call c_lhapdf_getconfig(pdf_in%ptr, key//c_null_char, outvalue, len(outvalue))
      end subroutine

      function computeUncertainty(array, pdfnum)
          implicit none
          type(pdfuncertainty) :: computeUncertainty
          integer, intent(in) :: pdfnum
          real(dp), target, intent(in) :: array(:)

          computeUncertainty = c_lhapdf_computeUncertainty(pdfs(pdfnum)%ptr, c_loc(array), size(array))
      end function

      function pdf_orderqcd(pdf_in)
          implicit none
          integer :: pdf_orderqcd
          class(pdf), intent(in) :: pdf_in

          pdf_orderqcd = c_lhapdf_orderqcd(pdf_in%ptr)
      end function

    function getpdf(numpdf)
        implicit none
        integer, intent(in) :: numpdf
        type(pdf) :: getpdf

        getpdf = pdfs(numpdf)
    end function
      
    subroutine pdf_readgrid(pdfnum, xvals, q2vals)
        implicit none
        integer, intent(in) :: pdfnum

        integer :: numx, numq2
        real(dp), allocatable, target, intent(out) :: xvals(:)
        real(dp), allocatable, target, intent(out) :: q2vals(:)

        numx = c_lhapdf_numxKnots(pdfs(pdfnum)%ptr)
        numq2 = c_lhapdf_numq2Knots(pdfs(pdfnum)%ptr)

        allocate(xvals(numx))
        allocate(q2vals(numq2))

        call c_lhapdf_getxKnots(pdfs(pdfnum)%ptr, c_loc(xvals(1)))
        call c_lhapdf_getq2Knots(pdfs(pdfnum)%ptr, c_loc(q2vals(1)))

    end subroutine

      function lhapdf_number_pdf(pdf_in)
          implicit none
          integer :: lhapdf_number_pdf
          class(pdf), intent(in) :: pdf_in

          lhapdf_number_pdf = c_lhapdf_number(trim(pdf_in%setname)//c_null_char)
      end function

      function lhapdf_number_name(setname)
          implicit none
          integer :: lhapdf_number_name
          character(*), intent(in) :: setname

          lhapdf_number_name = c_lhapdf_number(trim(setname)//c_null_char)
      end function

      subroutine lhapdf_pathsPrepend(path)
          implicit none
          character(*), intent(in) :: path

          call c_lhapdf_pathsPrepend(trim(path)//c_null_char)

      end subroutine

      function lhapdf_loadMember(setname,member)
          implicit none
          type(pdf) :: lhapdf_loadMember

          character(*), intent(in) :: setname
          integer, intent(in) :: member

          lhapdf_loadMember%setname = trim(setname)
          lhapdf_loadMember%ptr = c_lhapdf_loadMember(trim(setname)//c_null_char, member)

      end function

      function pdf_xfxq2(pdf_in, x, q2, flav)
          use ieee_arithmetic
          implicit none
          real(dp) :: pdf_xfxq2
          class(pdf), intent(in) :: pdf_in
          real(dp), intent(in) :: x,q2
          integer, intent(in) :: flav

          real(dp) :: q2_use

          if (ieee_is_nan(q2) .or. (.not. ieee_is_finite(q2))) then
              write (*,*) "xfxq2 called with q2 = ", q2
              q2_use = 1d0
          else
              q2_use = q2
          endif

          pdf_xfxq2 = c_lhapdf_evolve(pdf_in%ptr, x, q2_use, flav)
      end function

      function pdf_alphas(pdf_in, q)
          use ieee_arithmetic
          implicit none
          real(dp) :: pdf_alphas
          class(pdf), intent(in) :: pdf_in
          real(dp), intent(in) :: q

          if (ieee_is_nan(q) .or. (.not. ieee_is_finite(q)) .or. &
                q < 1._dp) then
              pdf_alphas = c_lhapdf_alphas(pdf_in%ptr, 1._dp)
          else
              pdf_alphas = c_lhapdf_alphas(pdf_in%ptr, q**2)
          endif

      end function

      function pdf_numflavors(pdf_in, q)
          implicit none
          integer :: pdf_numflavors
          class(pdf), intent(in) :: pdf_in
          real(dp), intent(in) :: q

          pdf_numflavors = c_lhapdf_numflavors(pdf_in%ptr, q)
      end function

      function pdf_quarkThreshold(pdf_in, i)
          implicit none
          real(dp) :: pdf_quarkThreshold
          class(pdf), intent(in) :: pdf_in
          integer, intent(in) :: i

          pdf_quarkThreshold = c_lhapdf_quarkThreshold(pdf_in%ptr, i)
      end function

      function pdf_quarkMass(pdf_in, i)
          implicit none
          real(dp) :: pdf_quarkMass
          class(pdf), intent(in) :: pdf_in
          integer, intent(in) :: i

          pdf_quarkMass = c_lhapdf_quarkMass(pdf_in%ptr, i)
      end function

!!!! routines affecting global pdfs state

      function getorderqcd()
          use PDFerrors, only: currentPDF
          implicit none
          integer :: getorderqcd

          getorderqcd = pdfs(currentPDF)%orderqcd()
      end function

      ! get alphas from central pdf
      function getalphas(q)
          use PDFerrors, only: currentPDF
          use qtResummation_params
          implicit none
          real(dp) :: getalphas
          real(dp), intent(in) :: q

          real(dp) :: T
          include 'masses.f'

          if (fix_alphas_nf5) then
              T=2._dp*log(Q/ZMASS)
!$OMP CRITICAL
              call NEWTON1(T, 0.118d0, getalphas, 3, 5)
!$OMP END CRITICAL
          else
              getalphas = pdfs(currentPDF)%alphas(q)
          endif
      end function

      function getnumflavors(q)
          use PDFerrors, only: currentPDF
          use qtResummation_params
          implicit none
          integer :: getnumflavors
          real(dp), intent(in) :: q

          if (fix_alphas_nf5) then
              getnumflavors = 5
          else
              getnumflavors = pdfs(currentPDF)%numflavors(q)
          endif
      end function

      function getquarkThreshold(i)
          use PDFerrors, only: currentPDF
          real(dp) :: getquarkThreshold
          integer, intent(in) :: i

          getquarkThreshold = pdfs(currentPDF)%quarkThreshold(i)
      end function

      function getquarkMass(i)
          use PDFerrors, only: currentPDF
          real(dp) :: getquarkMass
          integer, intent(in) :: i

          getquarkMass = pdfs(currentPDF)%quarkMass(i)
      end function

      subroutine initcentral(setnames, setmembers)
          use omp_lib
          implicit none
          character(len=*), intent(in) :: setnames(:)
          integer, intent(in) :: setmembers(:)
          integer :: j

          if (.not. allocated(pdfs)) then
              allocate(pdfs(0:size(setnames)-1))
              do j=1,size(setnames)
                  pdfs(j-1) = lhapdf_loadMember(trim(setnames(j))//c_null_char,setmembers(j))
              enddo
          endif
      end subroutine

      subroutine initcentralResummation(setnames, setmembers)
          use omp_lib
          implicit none
          character(len=*), intent(in) :: setnames(:)
          integer, intent(in) :: setmembers(:)
          integer :: j

          if (.not. allocated(respdfs)) then
              allocate(respdfs(16,0:size(setnames)-1))
              do j=1,size(setnames)
                  respdfs(1,j-1) = lhapdf_loadMember(trim(setnames(j))//'_B10'//c_null_char,setmembers(j))
                  respdfs(2,j-1) = lhapdf_loadMember(trim(setnames(j))//'_B11'//c_null_char,setmembers(j))
                  respdfs(3,j-1) = lhapdf_loadMember(trim(setnames(j))//'_B20'//c_null_char,setmembers(j))
                  respdfs(4,j-1) = lhapdf_loadMember(trim(setnames(j))//'_B21'//c_null_char,setmembers(j))
                  respdfs(5,j-1) = lhapdf_loadMember(trim(setnames(j))//'_B22'//c_null_char,setmembers(j))
                  respdfs(6,j-1) = lhapdf_loadMember(trim(setnames(j))//'_B30'//c_null_char,setmembers(j))
                  respdfs(7,j-1) = lhapdf_loadMember(trim(setnames(j))//'_B31'//c_null_char,setmembers(j))
                  respdfs(8,j-1) = lhapdf_loadMember(trim(setnames(j))//'_B32'//c_null_char,setmembers(j))
                  respdfs(9,j-1) = lhapdf_loadMember(trim(setnames(j))//'_B33'//c_null_char,setmembers(j))
                  respdfs(10,j-1) = lhapdf_loadMember(trim(setnames(j))//'_G10'//c_null_char,setmembers(j))
                  respdfs(11,j-1) = lhapdf_loadMember(trim(setnames(j))//'_B44'//c_null_char,setmembers(j))
                  respdfs(12,j-1) = lhapdf_loadMember(trim(setnames(j))//'_B43'//c_null_char,setmembers(j))
                  respdfs(13,j-1) = lhapdf_loadMember(trim(setnames(j))//'_B55'//c_null_char,setmembers(j))

                  respdfs(14,j-1) = lhapdf_loadMember(trim(setnames(j))//'_B42'//c_null_char,setmembers(j))
                  respdfs(15,j-1) = lhapdf_loadMember(trim(setnames(j))//'_B54'//c_null_char,setmembers(j))
                  respdfs(16,j-1) = lhapdf_loadMember(trim(setnames(j))//'_B66'//c_null_char,setmembers(j))
              enddo
          endif
      end subroutine

      subroutine initallResummation(setnames, setmembers)
          use omp_lib
          implicit none
          character(len=*), intent(in) :: setnames(:)
          integer, intent(in) :: setmembers(:)

          integer :: nmax
          integer :: j,k
          integer :: accum
          integer, allocatable :: memberCounts(:)

          allocate(memberCounts(size(setnames)))
          do j=1,size(setnames)
              memberCounts(j) = lhapdf_number(trim(setnames(j)))
          enddo
          nmax = sum(memberCounts(:))

          if (.not. allocated(respdfs)) then
              allocate(respdfs(16,0:nmax-1))

              accum = 0
              do j=1,size(setnames)
                  do k=0, memberCounts(j) - 1
                      respdfs(1,accum) = lhapdf_loadMember(trim(setnames(j))//'_B10'//c_null_char, k)
                      respdfs(2,accum) = lhapdf_loadMember(trim(setnames(j))//'_B11'//c_null_char, k)
                      respdfs(3,accum) = lhapdf_loadMember(trim(setnames(j))//'_B20'//c_null_char, k)
                      respdfs(4,accum) = lhapdf_loadMember(trim(setnames(j))//'_B21'//c_null_char, k)
                      respdfs(5,accum) = lhapdf_loadMember(trim(setnames(j))//'_B22'//c_null_char, k)
                      respdfs(6,accum) = lhapdf_loadMember(trim(setnames(j))//'_B30'//c_null_char, k)
                      respdfs(7,accum) = lhapdf_loadMember(trim(setnames(j))//'_B31'//c_null_char, k)
                      respdfs(8,accum) = lhapdf_loadMember(trim(setnames(j))//'_B32'//c_null_char, k)
                      respdfs(9,accum) = lhapdf_loadMember(trim(setnames(j))//'_B33'//c_null_char, k)
                      respdfs(10,accum) = lhapdf_loadMember(trim(setnames(j))//'_G10'//c_null_char, k)
                      respdfs(11,accum) = lhapdf_loadMember(trim(setnames(j))//'_B44'//c_null_char, k)
                      respdfs(12,accum) = lhapdf_loadMember(trim(setnames(j))//'_B43'//c_null_char, k)
                      respdfs(13,accum) = lhapdf_loadMember(trim(setnames(j))//'_B55'//c_null_char, k)

                      respdfs(14,accum) = lhapdf_loadMember(trim(setnames(j))//'_B42'//c_null_char, k)
                      respdfs(15,accum) = lhapdf_loadMember(trim(setnames(j))//'_B54'//c_null_char, k)
                      respdfs(16,accum) = lhapdf_loadMember(trim(setnames(j))//'_B66'//c_null_char, k)
                      accum = accum + 1
                  enddo
              enddo
          endif
      end subroutine initallResummation

      subroutine initall(setnames, setmembers)
          use omp_lib
          implicit none
          character(len=*), intent(in) :: setnames(:)
          integer, intent(in) :: setmembers(:)

          integer :: nmax
          integer :: j,k
          integer :: accum
          integer, allocatable :: memberCounts(:)

          allocate(memberCounts(size(setnames)))
          do j=1,size(setnames)
              memberCounts(j) = lhapdf_number(trim(setnames(j)))
          enddo
          nmax = sum(memberCounts(:))

          if (.not. allocated(pdfs)) then
              allocate(pdfs(0:nmax-1))

              accum = 0
              do j=1,size(setnames)
                  do k=0, memberCounts(j) - 1
                      pdfs(accum) = lhapdf_loadMember(trim(setnames(j))//c_null_char, k)
                      accum = accum + 1
                  enddo
              enddo
          endif
      end subroutine

      function fdist_one_beam2(ih, powAs, powLperp, x, xmu, parton, ibeam_in)
          use MCFMStorage, only: selectpdfs
          use PDFerrors, only: currentPDF
          implicit none
          real(dp) :: fdist_one_beam2
          integer, intent(in) :: powAs, powLperp
          integer, intent(in) :: ih
          ! ih = +1 for proton
          ! ih = -1 for anti-proton
          real(dp), intent(in) :: x, xmu
          integer, intent(in) :: parton
          integer, intent(in), optional :: ibeam_in

          integer :: ibeam_use
          integer :: iresuse

          if (present(ibeam_in)) then
              ibeam_use = ibeam_in
          else
              ibeam_use = 1
          endif

          fdist_one_beam2 = 0._dp

          if (powAs == 1 .and. powLperp == 0) then
              if (ih == 1) then
                  if (selectpdfs(ibeam_use,parton)) then
                      fdist_one_beam2 = respdfs(10,currentPDF)%xfxq2(x,xmu**2,parton)/x
                  endif
              else
                  if (selectpdfs(ibeam_use,parton)) then
                      fdist_one_beam2 = respdfs(10,currentPDF)%xfxq2(x,xmu**2,-parton)/x
                  endif
              endif
          else
              fdist_one_beam2 = 0._dp
          endif

      end function

      function fdist_one_beam(ih, powAs, powLperp, x, xmu, parton, ibeam_in)
          use MCFMStorage, only: selectpdfs
          use PDFerrors, only: currentPDF
          implicit none
          real(dp) :: fdist_one_beam
          integer, intent(in) :: powAs, powLperp
          integer, intent(in) :: ih
          ! ih = +1 for proton
          ! ih = -1 for anti-proton
          real(dp), intent(in) :: x, xmu
          integer, intent(in) :: parton
          integer, intent(in), optional :: ibeam_in

          integer :: ibeam_use
          integer :: iresuse

          if (present(ibeam_in)) then
              ibeam_use = ibeam_in
          else
              ibeam_use = 1
          endif

          fdist_one_beam = 0._dp

          if (powAs == 0 .and. powLperp == 0) then
              if (ih == 1) then
                  if (selectpdfs(ibeam_use,parton)) then
                      fdist_one_beam = pdfs(currentPDF)%xfxq2(x,xmu**2,parton)/x
                  endif
              else
                  if (selectpdfs(ibeam_use,parton)) then
                      fdist_one_beam = pdfs(currentPDF)%xfxq2(x,xmu**2,-parton)/x
                  endif
              endif
          else
              if (powAs == 1 .and. powLperp == 0) then
                  iresuse = 1
              elseif (powAs == 1 .and. powLperp == 1) then
                  iresuse = 2
              elseif (powAs == 2 .and. powLperp == 0) then
                  iresuse = 3
              elseif (powAs == 2 .and. powLperp == 1) then
                  iresuse = 4
              elseif (powAs == 2 .and. powLperp == 2) then
                  iresuse = 5
              elseif (powAs == 3) then
                  if (powLperp == 0) then
                      iresuse = 6
                  elseif (powLperp == 1) then
                      iresuse = 7
                  elseif (powLperp == 2) then
                      iresuse = 8
                  elseif (powLperp == 3) then
                      iresuse = 9
                  endif
              elseif (powAs == 4) then
                  if (powLperp == 4) then
                      iresuse = 11
                  elseif (powLperp == 3) then
                      iresuse = 12
                  elseif (powLperp == 2) then
                      iresuse = 14
                  endif
              elseif (powAs == 5) then
                  if (powLperp == 5) then
                      iresuse = 13
                  elseif (powLperp == 4) then
                      iresuse = 15
                  endif
              elseif (powAs == 6) then
                  if (powLperp == 6) then
                      iresuse = 16
                  endif
              endif

              if (ih == 1) then
                  if (selectpdfs(ibeam_use,parton)) then
                      fdist_one_beam = respdfs(iresuse,currentPDF)%xfxq2(x,xmu**2,parton)/x
                  endif
              else
                  if (selectpdfs(ibeam_use,parton)) then
                      fdist_one_beam = respdfs(iresuse,currentPDF)%xfxq2(x,xmu**2,-parton)/x
                  endif
              endif

          endif

      end function

      function fdist_one(ih, x, xmu, parton, ibeam_in)
          use MCFMStorage, only: selectpdfs
          use PDFerrors, only: currentPDF
          implicit none
          real(dp) :: fdist_one
          integer, intent(in) :: ih
          ! ih = +1 for proton
          ! ih = -1 for anti-proton
          real(dp), intent(in) :: x, xmu
          integer, intent(in) :: parton
          integer, intent(in), optional :: ibeam_in

          integer :: ibeam_use

          logical :: dummypdf
          common/dummypdf/dummypdf

          if (present(ibeam_in)) then
              ibeam_use = ibeam_in
          else
              ibeam_use = 1
          endif

          fdist_one = 0._dp

          if (x > 1._dp) then
              fdist_one = 0._dp
              return
          endif

          if (dummypdf) then
              if (selectpdfs(ibeam_use,parton)) then
                  !fdist_one = x**1d0 * (1d0-x)**1d0
                  fdist_one = x*(1d0-x)**(2d0+sin(real(parton,dp))) / x
              endif
          else
              if (ih == 1) then
                  if (selectpdfs(ibeam_use,parton)) then
                      fdist_one = pdfs(currentPDF)%xfxq2(x,xmu**2,parton)/x
                  endif
              else
                  if (selectpdfs(ibeam_use,parton)) then
                      fdist_one = pdfs(currentPDF)%xfxq2(x,xmu**2,-parton)/x
                  endif
              endif
          endif

      end function

      subroutine fdist_lhapdf_photonforgluon(ih, x, xmu, fx)
          use PDFerrors, only: currentPDF
          implicit none
          integer, intent(in) :: ih
          ! ih = +1 for proton
          ! ih = -1 for anti-proton
          real(dp), intent(in) :: x, xmu
          real(dp), intent(out) :: fx(-5:5)

          integer :: i

          if (x > 1d0) then
              fx(:) = 0._dp
              return
          endif

          if (ih ==1) then
              do i=-5,5
                  fx(i) = pdfs(currentPDF)%xfxq2(x, xmu**2, i) / x
              enddo
          else
              do i=-5,5
                  fx(i) = pdfs(currentPDF)%xfxq2(x, xmu**2, -i) / x
              enddo
          endif

          fx(0) = pdfs(currentPDF)%xfxq2(x, xmu**2, 22) / x

      end subroutine

      subroutine fdist_lhapdf_photonseparately(ih, x, xmu, fx, fxa)
          use PDFerrors, only: currentPDF
          implicit none
          integer, intent(in) :: ih
          ! ih = +1 for proton
          ! ih = -1 for anti-proton
          real(dp), intent(in) :: x, xmu
          real(dp), intent(out) :: fx(-5:5),fxa

          integer :: i

          if (x > 1d0) then
              fx(:) = 0._dp
              return
          endif

          if (ih ==1) then
              do i=-5,5
                  fx(i) = pdfs(currentPDF)%xfxq2(x, xmu**2, i) / x
              enddo
          else
              do i=-5,5
                  fx(i) = pdfs(currentPDF)%xfxq2(x, xmu**2, -i) / x
              enddo
          endif

          fxa = pdfs(currentPDF)%xfxq2(x, xmu**2, 22) / x

      end subroutine

      subroutine fdist(ih, x, xmu, fx, ibeam_in, mask)
          use MCFMStorage, only: selectpdfs
          use PDFerrors, only: currentPDF
          use singletop2_nnlo_vars, only: usemask, maskb1, maskb2
          implicit none
          include 'Cabibbo.f'
          include 'kprocess.f'
          integer, intent(in) :: ih
          ! ih = +1 for proton
          ! ih = -1 for anti-proton
          real(dp), intent(in) :: x, xmu
          real(dp), intent(out) :: fx(-5:5)
          integer, intent(in), optional :: ibeam_in
          logical, intent(in), optional :: mask(-5:5)
          real(dp) :: fxtmp(-5:5)

          logical :: dummypdf
          common/dummypdf/dummypdf

          logical :: mask_use(-5:5)

          integer :: ibeam_use

          integer :: i

          real(dp) :: fxa
          common/photonpdf/fxa
!$omp threadprivate(/photonpdf/)
      
          if (present(ibeam_in)) then
              ibeam_use = ibeam_in
          else
              ibeam_use = 1
          endif

          if (.not. present(mask)) then
              if (usemask) then
                  if (ibeam_in == 1) then
                      mask_use = maskb1
                  else
                      mask_use = maskb2
                  endif
              else
                  mask_use = .true.
              endif
          else
              mask_use = mask
          endif


          fx(:) = 0._dp

          if (dummypdf) then
              do i=-5,5
                  if (selectpdfs(ibeam_use,i)) then
                      fx(i) = x*(1d0-x)**(2d0+sin(real(i,dp))) / x
                  endif
              enddo
          else
              if (ih == 1) then
                  do i=-5,5
                      if (selectpdfs(ibeam_use,i) .and. mask_use(i)) then
                          !fx(i) = x**0.1d0*(1-x)
                          fx(i) = pdfs(currentPDF)%xfxq2(x, xmu**2, i) / x
                      endif
                  enddo
              else
                  do i=-5,5
                      if (selectpdfs(ibeam_use,i) .and. mask_use(i)) then
                          !fx(i) = x**0.1d0*(1-x)
                          fx(i) = pdfs(currentPDF)%xfxq2(x, xmu**2, -i) / x
                      endif
                  enddo
              endif
          endif

! special catches for process with a photon in the initial state
          if ((kcase == kWgaj_a) .or. (kcase == kWln_aq)) then
            call fdist_lhapdf_photonforgluon(ih,x,xmu,fx)
          elseif (kcase == kWgajja) then
            call fdist_lhapdf_photonseparately(ih,x,xmu,fx,fxa)
          endif
      
          if (CKMrotate) then
            fxtmp(:)=fx(:)
            fx(+1)=Vud_rotate**2*fxtmp(+1)+Vus_rotate**2*fxtmp(+3)+Vub_rotate**2*fxtmp(+5)
            fx(+3)=Vcd_rotate**2*fxtmp(+1)+Vcs_rotate**2*fxtmp(+3)+Vcb_rotate**2*fxtmp(+5)
            fx(-1)=Vud_rotate**2*fxtmp(-1)+Vus_rotate**2*fxtmp(-3)+Vub_rotate**2*fxtmp(-5)
            fx(-3)=Vcd_rotate**2*fxtmp(-1)+Vcs_rotate**2*fxtmp(-3)+Vcb_rotate**2*fxtmp(-5)
          endif

      end subroutine

end module
