!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module MCFMTaufit
    use Minpack
    use types
    implicit none

    public :: fittau, autofit

    private

    real(dp), allocatable :: xdat(:), ydat(:)
!$omp threadprivate(xdat, ydat)
    real(dp), allocatable :: weights(:)
!$omp threadprivate(weights)

    contains

    subroutine fcn(m, n, x, fvec, fjac, ldfjac, iflag)
        integer, intent(in) :: m, n
        integer, intent(in) :: ldfjac
        real(dp) :: fvec(m), fjac(ldfjac,n)
        real(dp), intent(inout) :: x(n)
        integer, intent(out) :: iflag

        include 'kpart.f'
        include 'taucut.f'
        integer :: i

        if (iflag == 0) then
            write (*,'(a)') ''
            do i=1,n
                write (*,'(g14.6)') x(i)
            enddo
        elseif (iflag == 1) then
            ! calculate the function at x and return this vector in fvec
            if (origKpart == ksnlo) then
                ! did not include linear term, since there is generally not enough precision
                ! to limit it
                fvec(:) = (x(1) + x(2)*xdat*log(xdat) + x(3)*xdat - ydat) * weights
            elseif (origKpart == knnlo) then
!                fvec(:) = (x(1) + x(2)*xdat*log(xdat)**3 + x(3)*xdat*log(xdat)**2 + x(4)*xdat - ydat) * weights
                fvec(:) = (x(1) + x(2)*xdat*log(xdat)**3 + x(3)*xdat*log(xdat)**2 - ydat) * weights
            endif
        elseif (iflag == 2) then
            ! calculate jacobian at x(i) and return this matrix in fjac
            if (origKpart == ksnlo) then
                fjac(:, 1) = (1._dp) * weights
                fjac(:, 2) = (xdat * log(xdat)) * weights
                fjac(:, 3) = (xdat) * weights
            elseif (origKpart == knnlo) then
                fjac(:, 1) = (1._dp) * weights
                fjac(:, 2) = (xdat * log(xdat)**3) * weights
                fjac(:, 3) = (xdat * log(xdat)**2) * weights
!                fjac(:, 4) = (xdat) * weights
           endif
        else
            write (*,*) "Called with unexpected value of iflag = ", iflag
            stop
        endif

    end subroutine

    subroutine fcn2(m, n, x, fvec, fjac, ldfjac, iflag)
        integer, intent(in) :: m, n
        integer, intent(in) :: ldfjac
        real(dp) :: fvec(m), fjac(ldfjac,n)
        real(dp), intent(inout) :: x(n)
        integer, intent(out) :: iflag

        include 'kpart.f'
        include 'taucut.f'
        integer :: i

        if (iflag == 0) then
            write (*,'(a)') ''
            do i=1,n
                write (*,'(g14.6)') x(i)
            enddo
        elseif (iflag == 1) then
            ! calculate the function at x and return this vector in fvec
            if (origKpart == ksnlo) then
                ! did not include linear term, since there is generally not enough precision
                ! to limit it
                fvec(:) = (x(1) + x(2)*xdat*log(xdat) + x(3)*xdat - ydat) * weights
            elseif (origKpart == knnlo) then
                fvec(:) = (x(1) + x(2)*xdat*log(xdat)**3 + x(3)*xdat*log(xdat)**2 + x(4)*xdat - ydat) * weights
            endif
        elseif (iflag == 2) then
            ! calculate jacobian at x(i) and return this matrix in fjac
            if (origKpart == ksnlo) then
                fjac(:, 1) = (1._dp) * weights
                fjac(:, 2) = (xdat * log(xdat)) * weights
                fjac(:, 3) = (xdat) * weights
            elseif (origKpart == knnlo) then
                fjac(:, 1) = (1._dp) * weights
                fjac(:, 2) = (xdat * log(xdat)**3) * weights
                fjac(:, 3) = (xdat * log(xdat)**2) * weights
                fjac(:, 4) = (xdat) * weights
           endif
        else
            write (*,*) "Called with unexpected value of iflag = ", iflag
            stop
        endif

    end subroutine

    subroutine fittau(xi, yi, weights_in, params, errs, resvar, info)
        real(dp), intent(in) :: xi(:), yi(:), weights_in(:)
        real(dp), intent(inout) :: params(:)
        real(dp), intent(inout) :: errs(:)
        real(dp), intent(out) :: resvar
        integer, intent(out) :: info
        
        real(dp), allocatable :: fvec(:), fjac(:,:)
        real(dp), allocatable :: diag(:), qtf(:)
        integer, allocatable :: ipvt(:)
        real(dp), parameter :: tol = 1e-4_dp

        real(dp), allocatable :: perm(:,:), r(:,:), hess(:,:)
        real(dp), allocatable :: hess_compressed(:), hess_inv_compressed(:)
        real(dp), allocatable :: workspace(:)
        integer :: nullty, ifault

        integer :: numEqns, n, ldfjac
        integer :: nfev, njev
        integer :: i,j,k,l

        integer :: rdf
        real(dp) :: deviance

        logical, parameter :: weightedFit = .true.

        n = size(params)

        numEqns = size(xi)
        ldfjac = numEqns
        if (allocated(xdat)) then
            deallocate(xdat)
        endif
        if (allocated(ydat)) then
            deallocate(ydat)
        endif
        if (allocated(weights)) then
            deallocate(weights)
        endif

        allocate(diag(numEqns))
        allocate(qtf(n))
        allocate(ipvt(n))

        allocate(perm(n,n))
        allocate(r(n,n))
        allocate(hess(n,n))

        allocate(xdat(numEqns))
        allocate(ydat(numEqns))
        allocate(weights(numEqns))
        xdat(:) = xi(:)
        ydat(:) = yi(:)
        weights(:) = weights_in(:)

        !for unweighted fit:
        if (.not. weightedFit) then
            !weights(:) = 1._dp
        endif

        allocate(fvec(numEqns))
        allocate(fjac(numEqns,n))

        if (n == 3) then
            call lmder(fcn, numEqns, n, params, fvec, fjac, ldfjac, tol, tol, 0._dp, &
                        100*(n+1), diag, 1, 100._dp, 0, info, nfev, njev, ipvt, qtf)
        else
            call lmder(fcn2, numEqns, n, params, fvec, fjac, ldfjac, tol, tol, 0._dp, &
                        100*(n+1), diag, 1, 100._dp, 0, info, nfev, njev, ipvt, qtf)
        endif
        
        ! column j of perm is column ipvt(j) of the identity matrix
        perm(:,:) = 0._dp
        do j=1,n
            perm(ipvt(j), j) = 1._dp
        enddo

        do i=1,n
            do j=1,n
                if (i <= j) then
                    r(i,j) = fjac(i,j)
                else
                    r(i,j) = 0._dp
                endif
            enddo
        enddo

        hess = matmul(perm, matmul(matmul(transpose(r),r), transpose(perm)))

        allocate(hess_compressed(n*(n+1)/2))
        k = 1
        do i=1,n
            do j=1,n
                if (i >= j) then
                    hess_compressed(k) = hess(i,j)
                    k = k + 1
                endif
            enddo
        enddo

        allocate(hess_inv_compressed(n*(n+1)/2))
        allocate(workspace(n*(n+1)/2))
        call syminv(hess_compressed, n, hess_inv_compressed, workspace, nullty, ifault)

        rdf = numEqns - n
        deviance = sum(fvec**2)
        resvar = deviance / real(rdf,dp)

        !allocate(errs(n))
        k = 1
        l = 1
        do i=1,n
            do j=1,n
                if (i >= j) then
                    if (i == j) then
                        if (weightedFit) then
                            errs(l) = sqrt(hess_inv_compressed(k))
                        else
                            errs(l) = sqrt(hess_inv_compressed(k) * resvar)
                        endif
                        l = l + 1
                    endif
                    k = k + 1
                endif
            enddo
        enddo

        ! for unweighted fit
        !write (*,*) "Residual standard error: ", sqrt(resvar)

    end subroutine

    subroutine autofit(ydat, weights, fitcross, fiterr, resvar, info)
        use SCET, only : tcutarray
        implicit none
        real(dp), intent(in) :: ydat(:), weights(:)
        real(dp), intent(out) :: fitcross, fiterr
        real(dp), intent(out) :: resvar
        integer, intent(out) :: info

        real(dp) :: resvar2
        real(dp), allocatable :: params(:), errs(:), params2(:), errs2(:)
        integer :: info1,info2

        include 'kpart.f'

        if (origKpart == ksnlo .and. size(tcutarray) < 3) then
            error stop "Need at least 3 taucut values for a NLO fit"
        elseif (origKpart == knnlo .and. size(tcutarray) < 4) then
            error stop "Need at least 4 taucut values for a NNLO fit"
        endif

        if (origKpart == ksnlo) then
            allocate(params(3))
            allocate(errs(3))
            params = [0d0, 0d0, 0d0]
        elseif (origKpart == knnlo) then
            allocate(params(3))
            allocate(errs(3))
            params = [0d0, 0d0, 0d0]
        endif

        call fittau(tcutarray, ydat, weights, params, errs, resvar, info1)
        fitcross = params(1)
        fiterr = errs(1)
        info = info1
        
        if (origKpart == knnlo) then
            allocate(params2(4))
            allocate(errs2(4))
            params2 = [0d0, 0d0, 0d0, 0d0]
            ! At NNLO: use fit including subleading term if it results in an improvement
            call fittau(tcutarray, ydat, weights, params2, errs2, resvar2, info2)
            if (resvar > 1._dp .and. (resvar2 < resvar)) then
                fitcross = params2(1)
                fiterr = errs2(1)
                resvar = resvar2
                info = info2
            endif
        endif

    end subroutine

end module
