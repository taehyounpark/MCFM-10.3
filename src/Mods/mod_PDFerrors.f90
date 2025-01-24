!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module PDFerrors
    use types
    implicit none

    public

    logical, public, save :: doPDFerrors

    ! this is set to .true. when there are different as(mZ) values
    ! among the PDF sets and members
    logical, public, save :: doPDFAlphas

    ! support for multiple pdf sets
    character(len=256), allocatable :: PDFnames(:)
    integer, allocatable :: PDFmembers(:)
    ! this is the size of PDFnames and PDFmembers
    integer :: numPDFsets

    integer, public, save :: currentPDF
!$omp threadprivate(currentPDF)


    ! this is just the size of pdfreweight
    integer, save :: maxPDFsets

    real(dp), public, save, allocatable :: pdfreweight(:)
!$omp threadprivate(pdfreweight)

end module
