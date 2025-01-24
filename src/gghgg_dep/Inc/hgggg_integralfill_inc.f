!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      use mod_qcdloop_c
c      use iso_c_binding
c      use iso_fortran_env
      implicit none
c calculation of the integrals needed order results but into arrays (e.g. Dint)
      include 'Inc/hgggglabels.f'
      include 'Inc/IntResults.f'
      integer p1,p2,p3,p4
      real(dp)::mtsq
      real(dp)::zip,s12,s13,s14,s23,s24,s34,s123,s124,s134,s234,Mhsq,mu2

c Dummy scale: there is no dependence on this in the result (integrals are finite)
      zip=0._dp
      mu2 = mtsq
      s12=s(p1,p2)
      s13=s(p1,p3)
      s14=s(p1,p4)
      s23=s(p2,p3)
      s24=s(p2,p4)
      s34=s(p3,p4)
      s123=s12+s13+s23
      s124=s12+s14+s24
      s234=s23+s24+s34
      s134=s13+s14+s34
      Mhsq=s12+s13+s14+s23+s24+s34

c Boxes
c      write(6,*) "zip = ",zip
c      write(6,*) "s34 = ",s34
c     write(6,*) "Mhsq = ",Mhsq
c      write(6,*) "s12 = ",s12
c      write(6,*) "s234 = ",s234
c      write(6,*) "mtsq = ",mtsq
c      write(6,*) "mu2 = ",mu2
      Dint(1)= qlI4(zip,zip,s34,Mhsq,s12,s234,mtsq,mtsq,mtsq,mtsq,mu2,0)
      Dint(2)= qlI4(zip,zip,s23,Mhsq,s14,s234,mtsq,mtsq,mtsq,mtsq,mu2,0)
      Dint(3)= qlI4(zip,zip,s34,Mhsq,s12,s134,mtsq,mtsq,mtsq,mtsq,mu2,0)
      Dint(4)= qlI4(zip,zip,s14,Mhsq,s23,s134,mtsq,mtsq,mtsq,mtsq,mu2,0)

      Dint(5)= qlI4(zip,zip,s14,Mhsq,s23,s124,mtsq,mtsq,mtsq,mtsq,mu2,0)
      Dint(6)= qlI4(Mhsq,zip,zip,s12,s124,s34,mtsq,mtsq,mtsq,mtsq,mu2,0)
      Dint(7)= qlI4(zip,zip,s23,Mhsq,s14,s123,mtsq,mtsq,mtsq,mtsq,mu2,0)
      Dint(8)= qlI4(Mhsq,zip,zip,s12,s123,s34,mtsq,mtsq,mtsq,mtsq,mu2,0)

      Dint(9)=qlI4(zip,s34,zip,Mhsq,s134,s234,mtsq,mtsq,mtsq,mtsq,mu2,0)
      Dint(10)=qlI4(zip,s23,zip,Mhsq,s123,s234,mtsq,mtsq,mtsq,mtsq,mu2,0)
      Dint(11)=qlI4(zip,s14,zip,Mhsq,s124,s134,mtsq,mtsq,mtsq,mtsq,mu2,0)
      Dint(12)=qlI4(zip,Mhsq,zip,s12,s123,s124,mtsq,mtsq,mtsq,mtsq,mu2,0)

      Dint(13)=qlI4(zip,zip,zip,s124,s14,s12,mtsq,mtsq,mtsq,mtsq,mu2,0)
      Dint(14)=qlI4(zip,zip,zip,s123,s12,s23,mtsq,mtsq,mtsq,mtsq,mu2,0)
      Dint(15)=qlI4(zip,zip,zip,s134,s34,s14,mtsq,mtsq,mtsq,mtsq,mu2,0)
      Dint(16)=qlI4(zip,zip,zip,s234,s23,s34,mtsq,mtsq,mtsq,mtsq,mu2,0)

      !write(6,*) "Dint",Dint(1)

c Triangles
      Cint(1)=qlI3(s234,zip,Mhsq,mtsq,mtsq,mtsq,mu2,0)
      Cint(2)=qlI3(s134,zip,Mhsq,mtsq,mtsq,mtsq,mu2,0)
      Cint(3)=qlI3(s124,zip,Mhsq,mtsq,mtsq,mtsq,mu2,0)
      Cint(4)=qlI3(s123,zip,Mhsq,mtsq,mtsq,mtsq,mu2,0)

      Cint(5)=qlI3(s34,Mhsq,s12,mtsq,mtsq,mtsq,mu2,0)
      Cint(6)=qlI3(s23,Mhsq,s14,mtsq,mtsq,mtsq,mu2,0)

      Cint(7)=qlI3(s23,s123,zip,mtsq,mtsq,mtsq,mu2,0)
      Cint(8)=qlI3(s34,s134,zip,mtsq,mtsq,mtsq,mu2,0)
      Cint(9)=qlI3(s34,s234,zip,mtsq,mtsq,mtsq,mu2,0)
      Cint(10)=qlI3(s14,s124,zip,mtsq,mtsq,mtsq,mu2,0)

      Cint(11)=qlI3(s12,s123,zip,mtsq,mtsq,mtsq,mu2,0)
      Cint(12)=qlI3(s14,s134,zip,mtsq,mtsq,mtsq,mu2,0)
      Cint(13)=qlI3(s12,s124,zip,mtsq,mtsq,mtsq,mu2,0)
      Cint(14)=qlI3(s23,s234,zip,mtsq,mtsq,mtsq,mu2,0)

      Cint(15)=qlI3(s12,zip,zip,mtsq,mtsq,mtsq,mu2,0)
      Cint(16)=qlI3(s23,zip,zip,mtsq,mtsq,mtsq,mu2,0)
      Cint(17)=qlI3(s34,zip,zip,mtsq,mtsq,mtsq,mu2,0)
      Cint(18)=qlI3(s14,zip,zip,mtsq,mtsq,mtsq,mu2,0)

c Bubbles
      Bint(1)=qlI2(s123,mtsq,mtsq,mu2,0)
      Bint(2)=qlI2(s234,mtsq,mtsq,mu2,0)
      Bint(3)=qlI2(s134,mtsq,mtsq,mu2,0)
      Bint(4)=qlI2(s124,mtsq,mtsq,mu2,0)

      Bint(5)=qlI2(s12,mtsq,mtsq,mu2,0)
      Bint(6)=qlI2(s23,mtsq,mtsq,mu2,0)
      Bint(7)=qlI2(s34,mtsq,mtsq,mu2,0)
      Bint(8)=qlI2(s14,mtsq,mtsq,mu2,0)

      Bint(9)=qlI2(Mhsq,mtsq,mtsq,mu2,0)

      return

