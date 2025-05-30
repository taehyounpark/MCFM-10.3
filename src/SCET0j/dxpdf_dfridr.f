!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine dxpdf_dfridr(ih,ibeam,x,h,err,dxfx)
          use LHAPDF
c Routine for computing derivative wrt x of pdf routines,
c based on original NR dfridr routine;
c Calls to pdfs are cached for optimization.

c Returns the derivative of a function func at a point x by Ridders'
c method of polynomial extrapolation. The value h is input as an
c estimated initial stepsize; it need not be small, but rather should be
c an increment in x over which func changes substantially. An estimate
c of the error in the derivative is returned as err.
c Parameters: Stepsize is decreased by CON at each iteration. Max size of
c tableau is set by NTAB. Return when error is SAFE worse than the best so far.

c Adapted from Numerical Recipes in Fortran 77, Section 5.7
c  (C) Copr. 1986-92 Numerical Recipes Software .
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'facscale.f'
      include 'nf.f'
      real(dp):: dxfx(-5:5),err,h,x
      real(dp), parameter:: CON=1.4_dp, CON2=CON*CON,
     &  BIG=1.e30_dp, SAFE=2._dp
      integer, parameter:: NTAB=10
      integer i,j,k,ih,ibeam
      real(dp):: errt,fac,hh,a(NTAB,NTAB)
      real(dp):: pdfup(-5:5),pdfdn(-5:5),
     & pdfup_save(NTAB,-5:5),pdfdn_save(NTAB,-5:5)
      logical:: calledpdf(2:NTAB)

c Catch exceptional point, at which pdfs vanish anyway
      if (x == 1._dp) then
        dxfx=0._dp
        return
      endif

c Catch improper usage
      if(h == 0._dp) then
        write(6,*) 'h must be nonzero in dxpdf_dfridr'
        stop
      endif

      calledpdf(:)=.false.

      call fdist(ih,x+h,facscale,pdfup)
      call fdist(ih,x-h,facscale,pdfdn)
      pdfup_save(1,:)=pdfup(:)
      pdfdn_save(1,:)=pdfdn(:)

      do k=-nf,nf

      hh=h
      a(1,1)=(pdfup_save(1,k)-pdfdn_save(1,k))/(two*hh)
c Restoring these two lines turns this routine into simple midpoint one
c      dxfx(k)=a(1,1)
c      cycle
      err=BIG
c Loop over successively smaller stepsizes and higher orders of extrapolation
      do i=2,NTAB
        hh=hh/CON
        if (calledpdf(i) .eqv. .false.) then
          call fdist(ih,x+hh,facscale,pdfup)
          call fdist(ih,x-hh,facscale,pdfdn)
          pdfup_save(i,:)=pdfup(:)
          pdfdn_save(i,:)=pdfdn(:)
          calledpdf(i)=.true.
        endif
        a(1,i)=(pdfup_save(i,k)-pdfdn_save(i,k))/(two*hh)
        fac=CON2
        do j=2,i
c Compute extrapolations of various orders with no new function evaluations
          a(j,i)=(a(j-1,i)*fac-a(j-1,i-1))/(fac-1._dp)
          fac=CON2*fac
c Compare new extrapolation to one order lower, both at present stepsize and previous one
          errt=max(abs(a(j,i)-a(j-1,i)),abs(a(j,i)-a(j-1,i-1)))
c If error is decreased, save the improved answer
          if (errt <= err) then
            err=errt
            dxfx(k)=a(j,i)
          endif
        enddo
c If higher order is significantly worse (by a factor SAFE), quit early
        if (abs(a(i,i)-a(i-1,i-1)) >= SAFE*err) exit
      enddo

      enddo

      return
      end

