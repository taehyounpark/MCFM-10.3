C=============================================================================
C---  HPLs  of Rank 4
C=============================================================================
c--- main forking function
      double complex function HPL4(n1, n2, n3, n4, x)
      implicit none
      integer n1,n2,n3,n4
      double complex x,ris
      double complex HPL4at0,HPL4at1,HPL4atm1
      double complex HPL4ar1,HPL4arm1,HPL4ar0,HPL4else
      double precision rad1,radm1,rad0

      rad1 = 0.01d0
      radm1 = 0.025d0
      rad0 = 0.025d0

      if ((abs(n1).gt.1).or.(abs(n2).gt.1).or.(abs(n3).gt.1)
     &     .or.(abs(n4).gt.1)) then
         print*, ""
         print*, "****************"
         print*, "Error in HPL4:"
         print*, "Indices",n1,n2,n3,n4," out of range !"
         print*, "Aborting..."
         print*,"****************"
         stop
      endif

      ris = dcmplx(0d0,0d0)
      
      if (x.eq.dcmplx(0d0,0d0)) then
c         print*, "I'm in 1"
         ris = HPL4at0(n1,n2,n3,n4)
      elseif (x.eq.dcmplx(1d0,0d0)) then
c         print*, "I'm in 2"
         ris = HPL4at1(n1,n2,n3,n4)
      elseif (x.eq.dcmplx(-1d0,0d0)) then
c         print*, "I'm in 3"
         ris = HPL4atm1(n1,n2,n3,n4)
      elseif (abs(x-dcmplx(1d0,0d0)).lt.rad1) then
c         print*, "I'm in 4"
         ris = HPL4ar1(n1,n2,n3,n4,x)
      elseif (abs(x+dcmplx(1d0,0d0)).lt.radm1) then
c         print*, "I'm in 5"
         ris = HPL4arm1(n1,n2,n3,n4,x)
      elseif (abs(x-dcmplx(0d0,0d0)).lt.rad0) then
c         print*, "I'm in 6"
         ris = HPL4ar0(n1,n2,n3,n4,x)
      else
c         print*, "I'm in 7"
         ris = HPL4else(n1,n2,n3,n4,x)
      endif
      HPL4=ris

c      write(*,'(I1,I1,I1,I1,D16.16,F16.16,F16.16)') 
c     &     n1,n2,n3,n4,x,abs(x),abs(x-dcmplx(1d0,0d0))

c      print*, n1,n2,n3,n4,x,abs(x),ris

      return
      end function
c ------------------------------------------------
      double complex function HPL4at0(n1, n2, n3, n4)
      implicit none
      integer n1,n2,n3,n4,j
      double complex ris
    
      j=1+(n4+1)*1+(n3+1)*3+(n2+1)*9+(n1+1)*27
      ris = dcmplx(0d0,0d0)

      if (j.eq.41) then             
         print*, ""
         print*, "****************"
         print*, "ERROR in HPL4: "
         print*, "HPL4(",n1,",",n2,",",n3,",",n4
     &        ,",0) is divergent!"
         print*, "Aborting..."
         print*,"****************"
         stop
      endif
      HPL4at0=ris
      return
      end function
c ------------------------------------------------
c --- Real part of HPL4     
      double precision function HPL4real(n1,n2,n3,n4,xr,xi)
      implicit none
      double precision xr,xi
      integer n1,n2,n3,n4
      double complex x,HPL4
      x=dcmplx(xr,xi)
      HPL4real = dreal(HPL4(n1,n2,n3,n4,x))
      return
      end

c --- Imaginary part of HPL4 
      double precision function HPL4im(n1,n2,n3,n4,xr,xi)
      implicit none
      double precision xr,xi
      integer n1,n2,n3,n4
      double complex x,HPL4
      x=dcmplx(xr,xi)
      HPL4im = dimag(HPL4(n1,n2,n3,n4,x))
      return
      end
