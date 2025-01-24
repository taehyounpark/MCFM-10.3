!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine rat_t14(j1,j2,j3,j4,j5,j6,j7,za,zb,res)
      implicit none
      include 'types.f'
      integer j1,j2,j3,j4,j5,j6,j7
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      complex(dp):: zab,res,zab2346,zab5273,zab2347,zab3276,zab3274,
     & zab2341,zab5342,zab3271,zab7342,zab1273,zab1342,zab3156,
     & den1,den2,den3,den4,den5,den6
      real(dp):: sum,soln,s156,s23p24,s23p37,a,b,c,lamp,lamm
c--- statement function
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
c--- end statement function

      sum=0.5d0*(s(j2,j4)+s(j2,j3)+s(j4,j7)+s(j3,j7))
      a=1;b=2*sum;c=s(j2,j7)*s(j4,j3)
      call solvequadratic(a,b,c,lamp,lamm)
      soln=-lamm

      zab2346=(zab(j2,j4,j6)+zab(j2,j3,j6))
      zab5273=(zab(j5,j2,j3)+zab(j5,j7,j3))
      zab2347=(zab(j2,j4,j7)+zab(j2,j3,j7))
      zab3276=(zab(j3,j2,j6)+zab(j3,j7,j6))
      zab3274=(zab(j3,j2,j4)+zab(j3,j7,j4))
      zab2341=(zab(j2,j4,j1)+zab(j2,j3,j1))
      zab5342=(zab(j5,j4,j2)+zab(j5,j3,j2))
      zab3271=(zab(j3,j2,j1)+zab(j3,j7,j1))
      zab7342=(zab(j7,j4,j2)+zab(j7,j3,j2))
      zab1273=(zab(j1,j2,j3)+zab(j1,j7,j3))
      zab1342=(zab(j1,j4,j2)+zab(j1,j3,j2))
      zab3156=(zab(j3,j1,j6)+zab(j3,j5,j6))
      s23p24=s(j2,j3)+s(j2,j4)
      s23p37=s(j2,j3)+s(j3,j7)
      s156=s(j1,j5)+s(j1,j6)+s(j5,j6)
      den1=(zab(j3,j4,j7)*soln-zab(j3,j2,j7)*s(j3,j4))
      den2=(s(j3,j4)*zab3271-zab(j3,j4,j1)*soln)
      den3=(zab(j7,j4,j3)*soln-s(j3,j4)*zab(j7,j2,j3))
      den4=(s(j2,j7)*zab2341-zab(j2,j7,j1)*soln)
      den5=(s(j3,j4)*zab1273-zab(j1,j4,j3)*soln)
      den6=(s(j2,j7)*zab1342-zab(j1,j7,j2)*soln)

      res=(-(s(j3,j4)*za(j5,j1)*zb(j7,j2)*zab3274
     & *(za(j2,j3)*s(j2,j7)*zb(j5,j6)*zab5342
     & +zb(j1,j6)*za(j2,j3)*s(j2,j7)*zab1342
     & -za(j2,j3)*zb(j5,j6)*zab(j5,j7,j2)*soln
     & -zab3156*s(j2,j7)*soln
     & -zb(j1,j6)*zab(j1,j7,j2)*za(j2,j3)*soln
     & +s(j2,j7)*zab(j3,j5,j6)*s23p24
     & +s(j2,j7)*zab(j3,j1,j6)*s23p24))
     & /(den1*s156*zab7342))

     & -s(j3,j4)*zb(j7,j2)*zab3274
     & *(s(j3,j4)*zab3276-zab(j3,j4,j6)*soln)
     & *(za(j2,j3)*s(j2,j7)*zab5342-za(j2,j3)*zab(j5,j7,j2)*soln
     & -s(j2,j7)*za(j3,j5)*soln
     & +s(j2,j7)*za(j3,j5)*s23p24)
     & /(den1*den2*zab7342)

     & -(s(j3,j4)*za(j5,j1)*zb(j7,j2)*(soln-s23p24)
     & *(zab(j3,j7,j2)*soln-s(j2,j7)*zab(j3,j4,j2))*zab3274
     & *(s(j3,j4)*zab3276-zab(j3,j4,j6)*soln))*s(j2,j7)
     & /(den1*den2*den6*zab7342)

     & +(za(j5,j1)*zb(j7,j2)*(zab(j2,j7,j4)*soln-zab(j2,j3,j4)*s(j2,j7))
     & *(za(j2,j3)*s(j3,j4)*zb(j5,j6)*zab5273
     & +zb(j1,j6)*za(j2,j3)*s(j3,j4)*zab1273
     & -za(j2,j3)*zab(j5,j4,j3)*zb(j5,j6)*soln
     & -zab(j2,j4,j3)*zab3156*soln
     & -zab(j1,j4,j3)*zb(j1,j6)*za(j2,j3)*soln
     & +zab(j2,j7,j3)*s(j3,j4)*zab(j3,j5,j6)
     & +zab(j2,j7,j3)*zab(j3,j1,j6)*s(j3,j4)))
     & /(den3*s156*zab2347)

     & +(zb(j7,j2)*(zab(j2,j7,j4)*soln-zab(j2,j3,j4)*s(j2,j7))
     & *(s(j2,j7)*zab2346-zab(j2,j7,j6)*soln)
     & *(za(j2,j3)*s(j3,j4)*zab5273-za(j2,j3)*zab(j5,j4,j3)*soln
     & -zab(j2,j4,j3)*za(j3,j5)*soln
     & +zab(j2,j7,j3)*s(j3,j4)*za(j3,j5)))
     & /(den3*den4*zab2347)

     & +(za(j5,j1)*zb(j7,j2)*(zab(j2,j4,j3)*soln-zab(j2,j7,j3)*s(j3,j4))
     & *(zab(j2,j7,j4)*soln-zab(j2,j3,j4)*s(j2,j7))
     & *(zab(j3,j4,j3)*soln-s(j3,j4)*s23p37)
     & *(s(j2,j7)*zab2346-zab(j2,j7,j6)*soln))
     & /(den3*den4*den5*zab2347)

      res=-res/(2d0*s(j3,j4)*s(j5,j6)*za(j2,j7)*s(j2,j7))

      return
      end

