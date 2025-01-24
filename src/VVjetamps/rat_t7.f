!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine rat_t7(j1,j2,j3,j4,j5,j6,j7,za,zb,res)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4,j5,j6,j7
      complex(dp):: zab,zaba,res,zab2345,zab2346,zab2347,
     & zab5342,zab6342,zab1342,zab7342,zab2156,zab3156
      real(dp):: s23p24,s156,s234

c--- statement function
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
c--- end statement function
      s23p24=s(j2,j3)+s(j2,j4)
      s156=s(j1,j5)+s(j1,j6)+s(j5,j6)
      s234=s(j2,j3)+s(j2,j4)+s(j3,j4)
      zab2345=(zab(j2,j4,j5)+zab(j2,j3,j5))
      zab2346=(zab(j2,j4,j6)+zab(j2,j3,j6))
      zab2156=(zab(j2,j1,j6)+zab(j2,j5,j6))
      zab3156=(zab(j3,j1,j6)+zab(j3,j5,j6))
      zab2347=(zab(j2,j4,j7)+zab(j2,j3,j7))
      zab1342=(zab(j1,j4,j2)+zab(j1,j3,j2))
      zab5342=(zab(j5,j4,j2)+zab(j5,j3,j2))
      zab6342=(zab(j6,j4,j2)+zab(j6,j3,j2))
      zab7342=(zab(j7,j4,j2)+zab(j7,j3,j2))
      zaba=-za(j2,j6)*zab2346-za(j2,j5)*zab2345

      res=-za(j2,j1)**2*za(j2,j3)*zab(j2,j3,j4)*za(j2,j5)*zb(j7,j2)
     & *zab2346*zab7342/(za(j1,j7)*za(j2,j7)**2*s23p24*zab2347*zaba)
      res=res
     & +(za(j2,j1)*za(j2,j3)*(za(j2,j1)*zb(j4,j2)*zab2346*zab5342
     & +zab(j2,j3,j4)*za(j2,j5)*zb(j6,j2)*zab1342))
     & /(za(j1,j7)*za(j2,j7)*s23p24*zaba)

      res=res
     & +za(j2,j1)*za(j2,j3)*zab(j2,j3,j4)*za(j2,j5)*zb(j7,j2)
     & *zab1342*zab2346/(za(j1,j7)*za(j2,j7)*s23p24*zab2347*zaba)

      res=res
     & +(za(j2,j1)*(za(j2,j3)*zb(j4,j2)*za(j5,j1)*zab2346
     & +za(j2,j1)*za(j3,j5)*zb(j4,j2)*zab2346
     & +za(j2,j3)*za(j5,j1)*zb(j6,j4)*s23p24))
     & /(za(j1,j7)*za(j2,j7)*zaba)

      res=res
     & +(za(j2,j1)**2*za(j2,j3)*zab(j2,j3,j4)*za(j2,j5)*zab2346
     & *(zb(j6,j2)*zab6342+zb(j5,j2)*zab5342))
     & /(za(j1,j7)*za(j2,j7)*s23p24*zaba**2)

      res=res
     & +(zab2347*zab1342**2*za(j1,j5)*za(j2,j3)*za(j2,j7)**2
     & *zb(j1,j6)*zb(j2,j4)
     & +zab2347*zab1342*s23p24*za(j1,j5)*za(j2,j7)**2*zb(j2,j4)*zab3156
     & -zab2347*zab1342*za(j1,j5)*zb(j2,j4)*za(j2,j7)*za(j2,j3)
     & *(za(j1,j2)*zb(j1,j6)*zab7342-za(j2,j7)*zb(j5,j6)*zab5342)
     & -zab2347*za(j1,j2)*za(j1,j5)*za(j2,j7)*zb(j2,j4)*zab7342
     & *(s23p24*zab3156+za(j2,j3)*zb(j5,j6)*zab5342)
     & +zab1342*zab7342*zb(j2,j7)*za(j1,j5)*zab2156
     & *za(j2,j3)*za(j2,j7)*zab(j2,j3,j4)
     & +za(j1,j2)*za(j1,j5)*za(j2,j3)*zab(j2,j3,j4)*zb(j2,j7)
     & *zab2156*zab7342**2)
     & /(za(j1,j7)*za(j2,j7)**2*s156*s23p24*zab2347*zab7342)

      res=res
     & -za(j2,j3)*zb(j4,j2)*zb(j6,j2)*zab1342**2*zab5342
     & /(za(j1,j7)*s23p24*(zb(j6,j2)*zab6342+zb(j5,j2)*zab5342)*zab7342)

      res=res
     & -zb(j4,j2)*zb(j6,j2)*zab1342
     & *(za(j3,j5)*zab1342+zab(j3,j4,j2)*za(j5,j1))
     & /(za(j1,j7)*(zb(j6,j2)*zab6342+zb(j5,j2)*zab5342)*zab7342)

      res=res
     & -2*zab(j1,j5,j6)*zb(j4,j2)*za(j5,j1)
     & *(-zab(j3,j4,j7)*s23p24+za(j2,j3)*zab(j4,j3,j2)*zb(j7,j4)
     & +za(j2,j3)*zab(j3,j4,j2)*zb(j7,j3))
     & /(za(j1,j7)*s156*s234*s23p24)

      res=res/(2d0*s(j3,j4)*s(j5,j6))

      return
      end

