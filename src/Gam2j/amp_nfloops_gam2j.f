!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

c==== These are the A_{61} N_F loops which go with the
c==== tree-level coupling structure (A_s in BDK four parton notation)
c===== only non-zero for ga-MHV and purely rational
      function amp_qqbggga_nf_A51_gaMHV(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: amp_qqbggga_nf_A51_gaMHV
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer i1,i2,i3,i4,i5


      amp_qqbggga_nf_A51_gaMHV=
     &     (zb(i3,i2)*zb(i4,i2))/(3._dp*za(i3,i4)*zb(i5,i1)*zb(i5,i2))

      return
      end


c===== This are the routines for the closed fermion amplitudes
c==== q(i1)^- qb(i2)^+ g(i3), g(i4) , gamma(i5)
c=====Note that the photon only couples to the loop fermion, and that this
c==== is for a single charge.
c=====C.W March 2016



c===== g MHV
      function amp_qqbggga_nf_gMHV(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: amp_qqbggga_nf_gMHV
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer i1,i2,i3,i4,i5
      complex(dp) :: rat
      complex(dp) :: L0,L1,Bubs,Boxes,Lsm1
      complex(dp) :: Boxint(4),boxcoeff(4)


c---- boxes
      Boxint(1)=Lsm1(-s(i3,i4),-s(i1,i2),-s(i3,i5),-s(i1,i2))
      boxcoeff(1)=-((za(i1,i5)**2*za(i3,i4)**2 +
     &      za(i1,i4)**2*za(i3,i5)**2)/
     &    (za(i1,i2)*za(i4,i5)**4))

      Boxes=boxint(1)*boxcoeff(1)

c==== bubs
      Bubs=  (-2*L0(-s(i3,i4),-s(i1,i2))*za(i1,i4)*za(i1,i5)*
     &     za(i3,i4)*za(i3,i5)*zb(i5,i4))/
     &   (s(i1,i2)*za(i1,i2)*za(i4,i5)**3) -
     &  (2*L0(-s(i3,i5),-s(i1,i2))*za(i1,i4)*za(i1,i5)*
     &     za(i3,i4)*za(i3,i5)*zb(i5,i4))/
     &   (s(i1,i2)*za(i1,i2)*za(i4,i5)**3) -
     &  (L1(-s(i3,i4),-s(i1,i2))*za(i1,i5)**2*za(i3,i4)**2*
     &     zb(i5,i4)**2)/
     &   (s(i1,i2)**2*za(i1,i2)*za(i4,i5)**2) -
     &  (L1(-s(i3,i5),-s(i1,i2))*za(i1,i4)**2*za(i3,i5)**2*
     &     zb(i5,i4)**2)/
     &   (s(i1,i2)**2*za(i1,i2)*za(i4,i5)**2)

c==== rat (after completion)

      Rat=(za(i1,i3)**2*zb(i5,i4))/
     &  (za(i1,i2)**2*za(i4,i5)*zb(i2,i1))

      amp_qqbggga_nf_gMHV=Boxes+Bubs+Rat
      return
      end

c===== Adjacent g-MHV
      function amp_qqbggga_nf_gMHVadj(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: amp_qqbggga_nf_gMHVadj
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer i1,i2,i3,i4,i5
      complex(dp) :: rat
      complex(dp) :: L0,L1,Bubs,Boxes,Lsm1
      complex(dp) :: Boxint(4),boxcoeff(4)



c---- boxes
      Boxint(1)=Lsm1(-s(i3,i4),-s(i1,i2),-s(i4,i5),-s(i1,i2))
      boxcoeff(1)= -((za(i1,i5)**2*za(i3,i4)**2 +
     &      za(i1,i3)**2*za(i4,i5)**2)/
     &    (za(i1,i2)*za(i3,i5)**4))


      Boxes=boxint(1)*boxcoeff(1)

c==== bubs
      Bubs= (2*L0(-s(i3,i4),-s(i1,i2))*za(i1,i3)*za(i1,i5)*
     &     za(i3,i4)*za(i4,i5)*zb(i5,i3))/
     &   (s(i1,i2)*za(i1,i2)*za(i3,i5)**3) +
     &  (2*L0(-s(i4,i5),-s(i1,i2))*za(i1,i3)*za(i1,i5)*
     &     za(i3,i4)*za(i4,i5)*zb(i5,i3))/
     &   (s(i1,i2)*za(i1,i2)*za(i3,i5)**3) -
     &  (L1(-s(i3,i4),-s(i1,i2))*za(i1,i5)**2*za(i3,i4)**2*
     &     zb(i5,i3)**2)/
     &   (s(i1,i2)**2*za(i1,i2)*za(i3,i5)**2) -
     &  (L1(-s(i4,i5),-s(i1,i2))*za(i1,i3)**2*za(i4,i5)**2*
     &     zb(i5,i3)**2)/
     &   (s(i1,i2)**2*za(i1,i2)*za(i3,i5)**2)


c==== rat (after completion)

      Rat= (za(i1,i4)**2*zb(i5,i3))/(za(i1,i2)**2*za(i3,i5)*zb(i2,i1))

      amp_qqbggga_nf_gMHVadj=Boxes+Bubs+Rat
c      write(6,*) amp_qqbggga_nf_gMHVadj*im
      return
      end


c===== GAMMA MHV
      function amp_qqbggga_nf_gaMHV(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: amp_qqbggga_nf_gaMHV
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer i1,i2,i3,i4,i5
      complex(dp) :: rat
      complex(dp) :: ALO,Vpole
      complex(dp) :: L0,L1,Bubs,Boxes,Lsm1
      complex(dp) :: Boxint(4),boxcoeff(4)

      ALO=za(i1,i5)**2/(za(i2,i3)*za(i3,i4)*za(i4,i1))
c      write(6,*) 'Tree = ',ALO

c-----pole expansion
      Vpole =czip
      Vpole=Vpole*ALO

c---- boxes
      Boxint(1)=Lsm1(-s(i3,i5),-s(i1,i2),-s(i4,i5),-s(i1,i2))
      boxcoeff(1)=-((za(i1,i4)**2*za(i3,i5)**2 +
     &      za(i1,i3)**2*za(i4,i5)**2)/
     &     (za(i1,i2)*za(i3,i4)**4))


      Boxes=boxint(1)*boxcoeff(1)
c      write(6,*) boxcoeff(1)/2._dp*s(i3,i5)*s(i4,i5)*im
c==== bubs
      Bubs=(-2*L0(-s(i3,i5),-s(i1,i2))*za(i1,i3)*za(i1,i4)*
     &     za(i3,i5)*za(i4,i5)*zb(i4,i3))/
     &   (s(i1,i2)*za(i1,i2)*za(i3,i4)**3) -
     &  (2*L0(-s(i4,i5),-s(i1,i2))*za(i1,i3)*za(i1,i4)*
     &     za(i3,i5)*za(i4,i5)*zb(i4,i3))/
     &   (s(i1,i2)*za(i1,i2)*za(i3,i4)**3) -
     &  (L1(-s(i3,i5),-s(i1,i2))*za(i1,i4)**2*
     &     za(i3,i5)**2*zb(i4,i3)**2)/
     &   (s(i1,i2)**2*za(i1,i2)*za(i3,i4)**2) -
     &  (L1(-s(i4,i5),-s(i1,i2))*za(i1,i3)**2*
     &     za(i4,i5)**2*zb(i4,i3)**2)/
     &   (s(i1,i2)**2*za(i1,i2)*za(i3,i4)**2)


c==== rat (after completion)

      Rat= (za(i1,i5)**2*zb(i4,i3))/
     &     (za(i1,i2)**2*za(i3,i4)*zb(i2,i1))

      amp_qqbggga_nf_gaMHV=Boxes+Bubs+Rat

c      write(6,*) Boxes*im
c      write(6,*) bubs*im
c      write(6,*) rat*im
c      write(6,*)  amp_qqbggga_nf_gaMHV*im
      return
      end
