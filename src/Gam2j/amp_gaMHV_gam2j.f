!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!


c==========gamma MHV LEADING COLOR ===========


      function amp_qqbggga_lc_gaMHV(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: amp_qqbggga_lc_gaMHV
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f'
      integer i1,i2,i3,i4,i5
      complex(dp) :: ALO,Vpole,rat
      complex(dp) :: L0,L1,Bubs,Boxes,Lsm1,lnrat
      complex(dp) :: Boxint(4),boxcoeff(4)
      integer i

      ALO=za(i1,i5)**2/(za(i2,i3)*za(i3,i4)*za(i4,i1))
c      write(6,*) 'Tree = ',ALO


c-----pole expansion
      Vpole =   3*epinv*epinv2 + epinv*
     &   (1.5 + lnrat(musq,-s(i1,i4)) +
     &     lnrat(musq,-s(i2,i3)) + lnrat(musq,-s(i3,i4))) +
     &  (6 + 3*lnrat(musq,-s(i1,i4)) +
     &     lnrat(musq,-s(i1,i4))**2 +
     &     lnrat(musq,-s(i2,i3))**2 + lnrat(musq,-s(i3,i4))**2
     &     )/2.


      Vpole=Vpole*ALO

      Boxint(1)=Lsm1(-s(i2,i5),-s(i3,i4),-s(i1,i5),-s(i3,i4))
      Boxint(2)=Lsm1(-s(i3,i4),-s(i1,i5),-s(i2,i3),-s(i1,i5))
      Boxint(3)=Lsm1(-s(i1,i4),-s(i2,i5),-s(i3,i4),-s(i2,i5))
c     -----s34, s25, s15;
      boxcoeff(1) = ALO
c-----s15, s34, s23;
      boxcoeff(2)= ALO
c-----s25, s34, s14;
      boxcoeff(3)=ALO

c------Boxes
      Boxes=czip
      do i=1,3
         Boxes=Boxes+boxcoeff(i)*Boxint(i)
      enddo


c======triangles (Get from overall V pole pieces)

c======bubbles
c===== s25
      Bubs=(L0(-s(i2,i5),-s(i1,i4))*za(i1,i3)*za(i1,i5)*za(i2,i5)*
     &     zb(i3,i2))/(s(i1,i4)*za(i1,i4)*za(i2,i3)*za(i3,i4)) +
     &  (L1(-s(i2,i5),-s(i1,i4))*za(i1,i3)**2*za(i2,i5)**2*
     &     zb(i3,i2)**2)/
     &   (2.*s(i1,i4)**2*za(i1,i4)*za(i2,i3)*za(i3,i4))


      rat=  za(i1,i5)**2/(2.*za(i1,i4)*za(i2,i3)*za(i3,i4)) +
     &  (za(i1,i3)**2*za(i2,i5)*zb(i3,i2)**2)/
     &   (2.*za(i1,i4)**2*za(i2,i3)*za(i3,i4)*zb(i4,i1)*
     &     zb(i5,i2)) -
     &  (zb(i3,i2)*zb(i4,i2))/
     &   (2.*za(i3,i4)*zb(i5,i1)*zb(i5,i2)) +
     &  (za(i2,i4)*zb(i4,i2)**2)/
     &     (2.*za(i2,i3)*za(i3,i4)*zb(i5,i1)*zb(i5,i2))
c---- Additional rational pieces (D+2 bubs)
      rat=rat+(2*zb(i3,i2)*zb(i4,i2))/(3.*za(i3,i4)*zb(i5,i1)*zb(i5,i2))

      amp_qqbggga_lc_gaMHV=Bubs+Boxes+Vpole+rat


c      write(6,*) 'CC  = ',amp_qqbggga_lc_gaMHV*im
c      write(6,*) 'Boxes =',Boxes
c      write(6,*) 'Bubs = ',Bubs
c      write(6,*) 'Vpole = ',Vpole
c      write(6,*) 'RAT = ',(rat+Bubs)*im

      return
      end


c=========gamma MHV subleading color

      function amp_qqbggga_slc_gaMHV(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: amp_qqbggga_slc_gaMHV
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f'
      integer i1,i2,i3,i4,i5
      complex(dp) :: rat
      complex(dp) :: ALO,Vpole
      complex(dp) :: L0,L1,Bubs,Boxes,Lsm1,lnrat
      complex(dp) :: Boxint(5),boxcoeff(5)
      integer i

      ALO=za(i1,i5)**2/(za(i2,i3)*za(i3,i4)*za(i4,i1))
c      write(6,*) 'Tree = ',ALO

c-----pole expansion
      Vpole =  (-3._dp - epinv*epinv2
     &     + epinv*(-1.5 - lnrat(musq,-s(i1,i2))) -
     &     lnrat(musq,-s(i1,i2))**2/2.
     &     - (3*lnrat(musq,-s(i2,i5)))/2._dp)
      Vpole=Vpole*ALO

      Boxint(1)=Lsm1(-s(i3,i5),-s(i1,i4),-s(i2,i5),-s(i1,i4))
      Boxint(2)=Lsm1(-s(i4,i5),-s(i2,i3),-s(i1,i5),-s(i2,i3))
      Boxint(3)=Lsm1(-s(i4,i5),-s(i1,i2),-s(i3,i5),-s(i1,i2))
      Boxint(4)=Lsm1(-s(i1,i2),-s(i3,i5),-s(i1,i4),-s(i3,i5))
      Boxint(5)=Lsm1(-s(i2,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))

      boxcoeff(1)= (za(i1,i2)**2*za(i3,i5)**2)/
     &  (za(i1,i4)*za(i2,i3)**3*za(i3,i4))
      boxcoeff(2)=(za(i1,i5)**2/(za(i1,i4)*za(i2,i3)*za(i3,i4)))
      boxcoeff(3)= (za(i1,i3)**2*za(i4,i5)**2)/
     &  (za(i1,i4)*za(i2,i3)*za(i3,i4)**3)
      boxcoeff(4)=(za(i1,i2)**2*za(i4,i5)**2)/
     &  (za(i1,i4)*za(i2,i3)*za(i2,i4)**2*za(i3,i4))
      boxcoeff(5)=(za(i1,i5)**2/(za(i1,i4)*za(i2,i3)*za(i3,i4)))

c------Boxes
      Boxes=czip
      do i=1,5
         Boxes=Boxes+boxcoeff(i)*Boxint(i)

      enddo



c===== completed bubbles
      bubs= -((lnrat(-s(i3,i5),-s(i2,i5))*za(i1,i2)*za(i1,i5)*za(i2,i5))/
     &     (za(i1,i4)*za(i2,i3)**2*za(i2,i4))) -
     &  (lnrat(-s(i4,i5),-s(i2,i5))*za(i1,i5)*za(i3,i5))/
     &   (za(i2,i3)*za(i3,i4)**2) +
     &  (lnrat(-s(i3,i5),-s(i2,i5))*za(i1,i5)*za(i4,i5))/
     &   (za(i2,i4)*za(i3,i4)**2) +
     &  (L0(-s(i1,i2),-s(i4,i5))*za(i1,i2)*za(i3,i5)**2*zb(i3,i2))/
     &   (s(i4,i5)*za(i2,i3)*za(i3,i4)**2) -
     &  (L0(-s(i1,i4),-s(i3,i5))*za(i1,i2)*za(i2,i5)**2*zb(i4,i2))/
     &   (s(i3,i5)*za(i2,i3)**2*za(i2,i4)) -
     &  (L0(-s(i1,i2),-s(i3,i5))*za(i1,i2)*za(i4,i5)**2*zb(i4,i2))/
     &   (s(i3,i5)*za(i2,i4)*za(i3,i4)**2) +
     &  (L0(-s(i1,i4),-s(i2,i5))*za(i1,i2)*za(i3,i5)**2*zb(i4,i3))/
     &   (s(i2,i5)*za(i2,i3)**2*za(i3,i4)) -
     &  (L1(-s(i1,i4),-s(i2,i5))*za(i1,i4)*za(i3,i5)**2*zb(i4,i3)**2)/
     &   (2.*s(i2,i5)**2*za(i2,i3)*za(i3,i4))

      rat=  (zb(i4,i3)*(-(za(i2,i5)*za(i3,i4)*zb(i4,i2)) +
     &      za(i1,i5)*za(i3,i5)*zb(i5,i1)))/
     &  (2.*za(i2,i3)*za(i2,i5)*za(i3,i4)*zb(i5,i1)*zb(i5,i2))


      amp_qqbggga_slc_gaMHV=Bubs+Boxes+Vpole+rat
c      write(6,*) 'CC  = ',amp_qqbggga_slc_gaMHV*im
c      write(6,*) 'Boxes =',Boxes*im
c      write(6,*) 'Bubs = ',Bubs*im
c      write(6,*) 'Vpole = ',Vpole*im
c      write(6,*) 'RAT = ',(rat)*im

c      stop
      return
      end
