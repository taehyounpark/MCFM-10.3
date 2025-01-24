!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine rat_t13(j1,j2,j3,j4,j5,j6,j7,za,zb,res)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4,j5,j6,j7
      complex(dp):: zab,res,zab2176,zab1567,zab6174,zab3567,zab7562,
     & zab2567,zab6172,zab3176,zab5176,zab7564
      real(dp):: y13, nnn,s57p67 ,s16p67 ,s234,sum,a,b,c,lamp,lamm
c--- statement function
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
      sum=0.5d0*(s(j6,j1)+s(j6,j7)+s(j5,j1)+s(j5,j7))

      a=1;b=2*sum;c=s(j6,j5)*s(j1,j7)
      call solvequadratic(a,b,c,lamp,lamm)
      y13=-lamm
      nnn=y13**2-s(j6,j5)*s(j1,j7)

      zab2176=(zab(j2,j1,j6)+zab(j2,j7,j6))
      zab1567=(zab(j1,j6,j7)+zab(j1,j5,j7))
      zab6174=(zab(j6,j1,j4)+zab(j6,j7,j4))
      zab3567=(zab(j3,j6,j7)+zab(j3,j5,j7))
      zab7562=(zab(j7,j6,j2)+zab(j7,j5,j2))
      zab2567=(zab(j2,j6,j7)+zab(j2,j5,j7))
      zab6172=(zab(j6,j1,j2)+zab(j6,j7,j2))
      zab3176=(zab(j3,j1,j6)+zab(j3,j7,j6))
      zab5176=(zab(j5,j1,j6)+zab(j5,j7,j6))
      zab7564=(zab(j7,j6,j4)+zab(j7,j5,j4))
      s57p67 =s(j5,j7)+s(j6,j7)
      s16p67 =s(j1,j6)+s(j6,j7)
      s234   =s(j2,j3)+s(j2,j4)+s(j3,j4)


      res=(zb(j2,j4)*za(j5,j6)*zb(j5,j6)*y13*zab1567**2*zab5176
     & *(zb(j1,j7)*za(j3,j4)*zab7564-zb(j1,j7)*za(j2,j3)*zab7562
     & -zb(j1,j4)*za(j3,j4)*y13
     & +zb(j1,j2)*za(j2,j3)*y13))
     & /(nnn**2*s234*(y13-s57p67)**2)
      res=res
     & +(zb(j2,j4)*za(j5,j6)*zb(j5,j6)*y13
     & *(za(j1,j5)*y13-za(j1,j7)*za(j5,j6)*zb(j6,j7))*zab1567*zab5176
     & *(zb(j1,j7)*za(j3,j4)*zab7564-zb(j1,j7)*za(j2,j3)*zab7562
     & -zb(j1,j4)*za(j3,j4)*y13
     & +zb(j1,j2)*za(j2,j3)*y13))
     & /(nnn**2*s234*(y13-s57p67)
     & *(za(j5,j7)*y13+zb(j1,j6)*za(j1,j7)*za(j5,j6)))
      res=res
     & +(zb(j2,j4)*za(j5,j6)*y13
     & *(za(j1,j5)*y13-za(j1,j7)*za(j5,j6)*zb(j6,j7))
     & *(zb(j1,j6)*y13+zb(j1,j7)*zb(j5,j6)*za(j5,j7))
     & *(zb(j5,j7)*y13+za(j1,j6)*zb(j1,j7)*zb(j5,j6))*zab1567*zab5176
     & *(zb(j1,j7)*za(j3,j4)*zab7564-zb(j1,j7)*za(j2,j3)*zab7562
     & -zb(j1,j4)*za(j3,j4)*y13
     & +zb(j1,j2)*za(j2,j3)*y13))
     & /(zb(j1,j7)*nnn**2*s234*(y13-s16p67)*(y13-s57p67)**2
     & *(za(j5,j7)*y13+zb(j1,j6)*za(j1,j7)*za(j5,j6)))
       res=res
     & -(zb(j2,j4)*za(j5,j6)*y13
     & *(za(j1,j5)*y13-za(j1,j7)*za(j5,j6)*zb(j6,j7))**2
     & *(zb(j1,j6)*y13+zb(j1,j7)*zb(j5,j6)*za(j5,j7))
     & *(zb(j5,j7)*y13+za(j1,j6)*zb(j1,j7)*zb(j5,j6))*zab5176
     & *(zb(j1,j7)*za(j3,j4)*zab7564-zb(j1,j7)*za(j2,j3)*zab7562
     & -zb(j1,j4)*za(j3,j4)*y13
     & +zb(j1,j2)*za(j2,j3)*y13))
     & /(zb(j1,j7)*nnn**2*s234*(y13-s16p67)*(y13-s57p67)
     & *(za(j5,j7)*y13+zb(j1,j6)*za(j1,j7)*za(j5,j6))**2)
      res=res
     & +(za(j5,j6)*y13*(za(j1,j5)*y13-za(j1,j7)*za(j5,j6)*zb(j6,j7))**2
     & *(zb(j1,j6)*y13+zb(j1,j7)*zb(j5,j6)*za(j5,j7))
     & *(za(j1,j7)*zab3567-za(j1,j3)*y13)*zab5176
     & *(zb(j5,j6)*zab6172+zb(j2,j5)*y13)
     & *(zb(j1,j7)*zab7564-zb(j1,j4)*y13))
     & /(za(j1,j7)*nnn**2*(y13-s16p67)*(y13-s57p67)
     & *(za(j5,j7)*y13+zb(j1,j6)*za(j1,j7)*za(j5,j6))
     & *(za(j5,j6)*zab2176+za(j2,j5)*y13)
     & *(zb(j1,j7)*zab7562-zb(j1,j2)*y13))
     & -(za(j5,j6)*y13*(za(j1,j5)*y13-za(j1,j7)*za(j5,j6)*zb(j6,j7))**2
     & *(zb(j1,j6)*y13+zb(j1,j7)*zb(j5,j6)*za(j5,j7))
     & *(za(j1,j7)*zab2567-za(j1,j2)*y13)
     & *(za(j5,j6)*zab3176+za(j3,j5)*y13)*zab5176
     & *(zb(j5,j6)*zab6172+zb(j2,j5)*y13)
     & *(zb(j1,j7)*zab7564-zb(j1,j4)*y13))
     & /(za(j1,j7)*nnn**2*(y13-s16p67)*(y13-s57p67)
     & *(za(j5,j7)*y13+zb(j1,j6)*za(j1,j7)*za(j5,j6))
     & *(za(j5,j6)*zab2176+za(j2,j5)*y13)**2
     & *(zb(j1,j7)*zab7562-zb(j1,j2)*y13))
      res=res
     & +(za(j5,j6)*zb(j5,j6)*y13*zab1567**2*(za(j1,j7)*zab3567-za(j1,j3)*y13)*zab5176
     & *(zb(j1,j7)*zab7564-zb(j1,j4)*y13))
     & /(nnn**2*(y13-s57p67)**2*(za(j1,j7)*zab2567-za(j1,j2)*y13))
      res=res
     & +(za(j5,j6)*zb(j5,j6)*y13*(za(j1,j5)*y13-za(j1,j7)*za(j5,j6)*zb(j6,j7))
     & *zab1567*(za(j5,j6)*zab3176+za(j3,j5)*y13)*zab5176
     & *(zb(j1,j7)*zab7564-zb(j1,j4)*y13))
     & /(nnn**2*(y13-s57p67)*(za(j5,j7)*y13+zb(j1,j6)*za(j1,j7)*za(j5,j6))
     & *(za(j5,j6)*zab2176+za(j2,j5)*y13))
     & +(za(j5,j6)*y13*(za(j1,j5)*y13-za(j1,j7)*za(j5,j6)*zb(j6,j7))
     & *(zb(j1,j6)*y13+zb(j1,j7)*zb(j5,j6)*za(j5,j7))
     & *(zb(j5,j7)*y13+za(j1,j6)*zb(j1,j7)*zb(j5,j6))*zab1567
     & *(za(j5,j6)*zab3176+za(j3,j5)*y13)*zab5176
     & *(zb(j1,j7)*zab7564-zb(j1,j4)*y13))
     & /(zb(j1,j7)*nnn**2*(y13-s16p67)*(y13-s57p67)**2
     & *(za(j5,j7)*y13+zb(j1,j6)*za(j1,j7)*za(j5,j6))
     & *(za(j5,j6)*zab2176+za(j2,j5)*y13))
      res=res
     & -(za(j5,j6)*y13*(za(j1,j5)*y13-za(j1,j7)*za(j5,j6)*zb(j6,j7))**2
     & *(zb(j1,j6)*y13+zb(j1,j7)*zb(j5,j6)*za(j5,j7))
     & *(zb(j5,j7)*y13+za(j1,j6)*zb(j1,j7)*zb(j5,j6))
     & *(za(j5,j6)*zab3176+za(j3,j5)*y13)*zab5176
     & *(zb(j1,j7)*zab7564-zb(j1,j4)*y13))
     & /(zb(j1,j7)*nnn**2*(y13-s16p67)*(y13-s57p67)
     & *(za(j5,j7)*y13+zb(j1,j6)*za(j1,j7)*za(j5,j6))**2
     & *(za(j5,j6)*zab2176+za(j2,j5)*y13))
     & +(zb(j2,j4)*za(j5,j6)*(za(j1,j5)*y13-za(j1,j7)*za(j5,j6)*zb(j6,j7))
     & *(zb(j1,j6)*y13+zb(j1,j7)*zb(j5,j6)*za(j5,j7))
     & *(za(j1,j5)*zab3176+za(j1,j7)*za(j3,j5)*zb(j6,j7)))
     & /(za(j1,j7)*(za(j5,j7)*y13+zb(j1,j6)*za(j1,j7)*za(j5,j6))
     & *(za(j5,j6)*zab2176+za(j2,j5)*y13)
     & *(zb(j1,j7)*zab7562-zb(j1,j2)*y13))
     & +(zb(j5,j6)*y13*(za(j1,j5)*y13-za(j1,j7)*za(j5,j6)*zb(j6,j7))*zab1567**2
     & *(za(j5,j6)*zab2176+za(j2,j5)*y13)
     & *(za(j1,j7)*zab3567-za(j1,j3)*y13)
     & *(zb(j5,j6)*zab6174+zb(j4,j5)*y13)
     & *(zb(j1,j7)*zab7562-zb(j1,j2)*y13))
     & /(nnn**2*(y13-s57p67)**2*(za(j1,j7)*zab2567-za(j1,j2)*y13)**2
     & *(zb(j5,j6)*zab6172+zb(j2,j5)*y13))
      res=res
     & -(zb(j5,j6)*y13*(za(j1,j5)*y13-za(j1,j7)*za(j5,j6)*zb(j6,j7))*zab1567**2
     & *(za(j5,j6)*zab3176+za(j3,j5)*y13)
     & *(zb(j5,j6)*zab6174+zb(j4,j5)*y13)
     & *(zb(j1,j7)*zab7562-zb(j1,j2)*y13))
     & /(nnn**2*(y13-s57p67)**2*(za(j1,j7)*zab2567-za(j1,j2)*y13)
     & *(zb(j5,j6)*zab6172+zb(j2,j5)*y13))
     & +(zb(j1,j7)*zb(j2,j4)*zb(j5,j6)*y13
     & *(za(j1,j5)*y13-za(j1,j7)*za(j5,j6)*zb(j6,j7))
     & *(za(j5,j7)*y13+zb(j1,j6)*za(j1,j7)*za(j5,j6))*zab1567**2
     & *(za(j3,j4)*zb(j5,j6)*zab6174-za(j2,j3)*zb(j5,j6)*zab6172
     & +za(j3,j4)*zb(j4,j5)*y13
     & -za(j2,j3)*zb(j2,j5)*y13))
     & /(za(j1,j7)*nnn**2*s234*(y13-s57p67)**2
     & *(zb(j5,j7)*y13+za(j1,j6)*zb(j1,j7)*zb(j5,j6)))
     & -(zb(j1,j7)*zb(j2,j4)*zb(j5,j6)*y13
     & *(za(j1,j5)*y13-za(j1,j7)*za(j5,j6)*zb(j6,j7))**2*zab1567
     & *(za(j3,j4)*zb(j5,j6)*zab6174-za(j2,j3)*zb(j5,j6)*zab6172
     & +za(j3,j4)*zb(j4,j5)*y13
     & -za(j2,j3)*zb(j2,j5)*y13))
     & /(za(j1,j7)*nnn**2*s234*(y13-s57p67)
     & *(zb(j5,j7)*y13+za(j1,j6)*zb(j1,j7)*zb(j5,j6)))
     & -(zb(j2,j4)*y13*(za(j1,j5)*y13-za(j1,j7)*za(j5,j6)*zb(j6,j7))**2
     & *(zb(j1,j6)*y13+zb(j1,j7)*zb(j5,j6)*za(j5,j7))*zab1567
     & *(za(j3,j4)*zb(j5,j6)*zab6174-za(j2,j3)*zb(j5,j6)*zab6172
     & +za(j3,j4)*zb(j4,j5)*y13
     & -za(j2,j3)*zb(j2,j5)*y13))
     & /(za(j1,j7)*nnn**2*s234*(y13-s16p67)*(y13-s57p67)**2)
     & -(zb(j2,j4)*y13*(za(j1,j5)*y13-za(j1,j7)*za(j5,j6)*zb(j6,j7))**3
     & *(zb(j1,j6)*y13+zb(j1,j7)*zb(j5,j6)*za(j5,j7))
     & *(za(j3,j4)*zb(j5,j6)*zab6174-za(j2,j3)*zb(j5,j6)*zab6172
     & +za(j3,j4)*zb(j4,j5)*y13
     & -za(j2,j3)*zb(j2,j5)*y13))
     & /(za(j1,j7)*nnn**2*s234*(y13-s16p67)*(y13-s57p67)
     & *(za(j5,j7)*y13+zb(j1,j6)*za(j1,j7)*za(j5,j6)))
      res=res
     & +(zb(j1,j7)*zb(j5,j6)*y13*(za(j1,j5)*y13-za(j1,j7)*za(j5,j6)*zb(j6,j7))
     & *(za(j5,j7)*y13+zb(j1,j6)*za(j1,j7)*za(j5,j6))*zab1567**2
     & *(za(j1,j7)*zab3567-za(j1,j3)*y13)
     & *(zb(j5,j6)*zab6174+zb(j4,j5)*y13))
     & /(za(j1,j7)*nnn**2*(y13-s57p67)**2
     & *(zb(j5,j7)*y13+za(j1,j6)*zb(j1,j7)*zb(j5,j6))
     & *(za(j1,j7)*zab2567-za(j1,j2)*y13))
     & -(zb(j1,j7)*zb(j5,j6)*y13*(za(j1,j5)*y13-za(j1,j7)*za(j5,j6)*zb(j6,j7))**2
     & *zab1567*(za(j1,j7)*zab3567-za(j1,j3)*y13)
     & *(zb(j5,j6)*zab6174+zb(j4,j5)*y13))
     & /(za(j1,j7)*nnn**2*(y13-s57p67)*(zb(j5,j7)*y13+za(j1,j6)*zb(j1,j7)*zb(j5,j6))
     & *(za(j1,j7)*zab2567-za(j1,j2)*y13))
      res=res
     & -(y13*(za(j1,j5)*y13-za(j1,j7)*za(j5,j6)*zb(j6,j7))**2
     & *(zb(j1,j6)*y13+zb(j1,j7)*zb(j5,j6)*za(j5,j7))*zab1567
     & *(za(j1,j7)*zab3567-za(j1,j3)*y13)*(zb(j5,j6)*zab6174+zb(j4,j5)*y13))
     & /(za(j1,j7)*nnn**2*(y13-s16p67)*(y13-s57p67)**2
     & *(za(j1,j7)*zab2567-za(j1,j2)*y13))
     & -(y13*(za(j1,j5)*y13-za(j1,j7)*za(j5,j6)*zb(j6,j7))**3
     & *(zb(j1,j6)*y13+zb(j1,j7)*zb(j5,j6)*za(j5,j7))
     & *(za(j5,j6)*zab3176+za(j3,j5)*y13)*(zb(j5,j6)*zab6174+zb(j4,j5)*y13))
     & /(za(j1,j7)*nnn**2*(y13-s16p67)*(y13-s57p67)
     & *(za(j5,j7)*y13+zb(j1,j6)*za(j1,j7)*za(j5,j6))
     & *(za(j5,j6)*zab2176+za(j2,j5)*y13))
     & +(zb(j2,j4)*zb(j5,j6)*(y13-s16p67)*zab1567
     & *(za(j1,j5)*za(j1,j7)*zab3567-za(j1,j7)*za(j3,j5)*zab1567
     & -za(j1,j3)*za(j1,j5)*y13))
     & /(za(j1,j7)*(y13-s57p67)*(za(j1,j7)*zab2567-za(j1,j2)*y13)
     & *(zb(j5,j6)*zab6172+zb(j2,j5)*y13))
      res=res
     & -(za(j1,j5)*zb(j4,j6)*zab1567*(za(j1,j7)*zab3567-za(j1,j3)*y13))
     & /(za(j1,j7)*(y13-s57p67)*(za(j1,j7)*zab2567-za(j1,j2)*y13))
     & -(za(j1,j5)*zb(j4,j6)*(za(j1,j5)*y13-za(j1,j7)*za(j5,j6)*zb(j6,j7))
     & *(za(j5,j6)*zab3176+za(j3,j5)*y13))
     & /(za(j1,j7)*(za(j5,j7)*y13+zb(j1,j6)*za(j1,j7)*za(j5,j6))
     & *(za(j5,j6)*zab2176+za(j2,j5)*y13))
     & -(za(j1,j5)*zb(j2,j4)*(za(j3,j4)*zb(j4,j6)-za(j2,j3)*zb(j2,j6))*zab1567)
     & /(za(j1,j7)*s234*(y13-s57p67))
     & -(za(j1,j5)*zb(j2,j4)*(za(j3,j4)*zb(j4,j6)-za(j2,j3)*zb(j2,j6))
     & *(za(j1,j5)*y13-za(j1,j7)*za(j5,j6)*zb(j6,j7)))
     & /(za(j1,j7)*s234*(za(j5,j7)*y13+zb(j1,j6)*za(j1,j7)*za(j5,j6)))

      res=-res/(2d0*s(j4,j3)*s(j6,j5))
      return
      end
