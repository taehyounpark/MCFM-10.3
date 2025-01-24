!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine rat_t15(j1,j2,j3,j4,j5,j6,j7,za,zb,res)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer::j1,j2,j3,j4,j5,j6,j7
      complex(dp):: zab,res,zab1564,zab1346,zab2564,zab2346,
     & zab3564,zab4567,zab4562,zab4561,zab5346,zab6342,zab6341,zab6347,
     & zab7564,zab7346
      real(dp):: y13, nn,sum,s45p46,s36p46,a,b,c,lamp,lamm
c--- statement function
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
c--- statement function
      sum=0.5d0*s(j4,j6)+0.5d0*s(j4,j5)+0.5d0*s(j3,j6)+0.5d0*s(j3,j5)
      a=1;b=2*sum;c=s(j4,j3)*s(j6,j5)
      call solvequadratic(a,b,c,lamp,lamm)
      y13=-lamm
      nn=y13**2-s(j4,j3)*s(j6,j5)

      zab1564=(zab(j1,j6,j4)+zab(j1,j5,j4))
      zab1346=(zab(j1,j4,j6)+zab(j1,j3,j6))
      zab2564=(zab(j2,j6,j4)+zab(j2,j5,j4))
      zab2346=(zab(j2,j4,j6)+zab(j2,j3,j6))

      zab3564=(zab(j3,j6,j4)+zab(j3,j5,j4))
      zab4567=(zab(j4,j6,j7)+zab(j4,j5,j7))
      zab4562=(zab(j4,j6,j2)+zab(j4,j5,j2))

      zab4561=(zab(j4,j6,j1)+zab(j4,j5,j1))
      zab5346=(zab(j5,j4,j6)+zab(j5,j3,j6))
      zab6342=(zab(j6,j4,j2)+zab(j6,j3,j2))
      zab6341=(zab(j6,j4,j1)+zab(j6,j3,j1))
      zab6347=(zab(j6,j4,j7)+zab(j6,j3,j7))

      zab7564=(zab(j7,j6,j4)+zab(j7,j5,j4))
      zab7346=(zab(j7,j4,j6)+zab(j7,j3,j6))
      s45p46=s(j4,j5)+s(j4,j6)
      s36p46=s(j3,j6)+s(j4,j6)

      res=(za(j2,j7)*zb(j5,j6)*y13
     & *(za(j3,j5)*y13+za(j3,j4)*zb(j4,j6)*za(j5,j6))**2
     & *(zb(j4,j5)*y13-zb(j3,j4)*za(j3,j6)*zb(j5,j6))
     & *(za(j3,j4)*zab1564+za(j1,j3)*y13)
     & *(zb(j3,j4)*zab4567-zb(j3,j7)*y13)
     & *(zb(j5,j6)*zab6342+zb(j2,j5)*y13))
     & /(nn**2*(y13-s45p46)*(zb(j5,j6)*zab6341+zb(j1,j5)*y13)
     & *(zb(j5,j6)*zab6347-zb(j5,j7)*y13)
     & *(za(j3,j4)*zab7564-za(j3,j7)*y13))
     & +(za(j2,j7)*zb(j4,j7)*zb(j5,j6)*(y13-s36p46)
     & *(za(j1,j5)*za(j3,j4)*zab3564-za(j3,j4)*za(j3,j5)*zab1564
     & -za(j1,j3)*za(j3,j5)*y13)
     & *(zb(j5,j6)*zab6342+zb(j2,j5)*y13))
     & /((zb(j5,j6)*zab6341+zb(j1,j5)*y13)
     & *(zb(j5,j6)*zab6347-zb(j5,j7)*y13)
     & *(za(j3,j4)*zab7564-za(j3,j7)*y13))
     & +(za(j2,j7)*za(j3,j4)*zb(j5,j6)*y13
     & *(za(j3,j5)*y13+za(j3,j4)*zb(j4,j6)*za(j5,j6))
     & *(zb(j4,j5)*y13-zb(j3,j4)*za(j3,j6)*zb(j5,j6))
     & *(za(j5,j6)*zab1346+za(j1,j5)*y13)
     & *(za(j3,j4)*zab1564+za(j1,j3)*y13)*zab3564
     & *(zb(j3,j4)*zab4567-zb(j3,j7)*y13))
     & /(nn**2*(y13-s45p46)*(za(j3,j4)*zab2564+za(j2,j3)*y13)
     & *(zb(j5,j6)*zab6347-zb(j5,j7)*y13)
     & *(za(j3,j4)*zab7564-za(j3,j7)*y13))
     & -(za(j2,j7)*za(j3,j4)*zb(j5,j6)*y13
     & *(za(j3,j5)*y13+za(j3,j4)*zb(j4,j6)*za(j5,j6))
     & *(zb(j4,j5)*y13-zb(j3,j4)*za(j3,j6)*zb(j5,j6))
     & *(za(j3,j4)*zab1564+za(j1,j3)*y13)**2
     & *(za(j5,j6)*zab2346+za(j2,j5)*y13)*zab3564
     & *(zb(j3,j4)*zab4562+zb(j2,j3)*y13))
     & /(nn**2*(y13-s45p46)*(za(j3,j4)*zab2564+za(j2,j3)*y13)**2
     & *(zb(j5,j6)*zab6342+zb(j2,j5)*y13)
     & *(za(j3,j4)*zab7564-za(j3,j7)*y13))
      res=res
     & +(za(j2,j7)*zb(j5,j6)*y13
     & *(za(j3,j5)*y13+za(j3,j4)*zb(j4,j6)*za(j5,j6))**2
     & *(zb(j4,j5)*y13-zb(j3,j4)*za(j3,j6)*zb(j5,j6))
     & *(za(j3,j4)*zab1564+za(j1,j3)*y13)**2
     & *(zb(j3,j4)*zab4562+zb(j2,j3)*y13))
     & /(nn**2*(y13-s45p46)*(za(j3,j4)*zab2564+za(j2,j3)*y13)
     & *(zb(j5,j6)*zab6342+zb(j2,j5)*y13)
     & *(za(j3,j4)*zab7564-za(j3,j7)*y13))
     & -(zb(j2,j4)*za(j2,j7)*zb(j5,j6)*(y13-s36p46)
     & *(za(j3,j4)*zab1564+za(j1,j3)*y13)
     & *(za(j1,j5)*za(j3,j4)*zab3564-za(j3,j4)*za(j3,j5)*zab1564
     & -za(j1,j3)*za(j3,j5)*y13))
     & /((za(j3,j4)*zab2564+za(j2,j3)*y13)
     & *(zb(j5,j6)*zab6342+zb(j2,j5)*y13)
     & *(za(j3,j4)*zab7564-za(j3,j7)*y13))
      res=res
     & +(za(j2,j7)*za(j3,j4)*zb(j3,j4)*za(j5,j6)*zb(j5,j6)*y13
     & *(za(j3,j4)*zab1564+za(j1,j3)*y13)*zab3564*zab5346
     & *(zb(j5,j6)*zab6342+zb(j2,j5)*y13))
     & /(nn**2*(zb(j5,j6)*zab6341+zb(j1,j5)*y13)
     & *(za(j3,j4)*zab7564-za(j3,j7)*y13))
      res=res
     & +(za(j2,j7)*za(j3,j4)*y13
     & *(za(j3,j5)*y13+za(j3,j4)*zb(j4,j6)*za(j5,j6))
     & *(zb(j3,j6)*y13-zb(j3,j4)*za(j4,j5)*zb(j5,j6))
     & *(zb(j4,j5)*y13-zb(j3,j4)*za(j3,j6)*zb(j5,j6))
     & *(za(j5,j6)*zab1346+za(j1,j5)*y13)*zab3564
     & *(zb(j5,j6)*zab6342+zb(j2,j5)*y13))
     & /(nn**2*(y13-s36p46)*(y13-s45p46)
     & *(zb(j5,j6)*zab6341+zb(j1,j5)*y13)
     & *(za(j3,j4)*zab7564-za(j3,j7)*y13))
     & +(za(j1,j5)*za(j2,j7)*za(j3,j4)*zb(j4,j6)*zab3564
     & *(zb(j5,j6)*zab6342+zb(j2,j5)*y13))
     & /((zb(j5,j6)*zab6341+zb(j1,j5)*y13)
     & *(za(j3,j4)*zab7564-za(j3,j7)*y13))
      res=res
     & -(za(j2,j7)*za(j3,j4)*zb(j5,j6)*y13
     & *(za(j3,j5)*y13+za(j3,j4)*zb(j4,j6)*za(j5,j6))
     & *(zb(j4,j5)*y13-zb(j3,j4)*za(j3,j6)*zb(j5,j6))
     & *(za(j5,j6)*zab1346+za(j1,j5)*y13)*zab3564
     & *(zb(j3,j4)*zab4561+zb(j1,j3)*y13)
     & *(zb(j5,j6)*zab6342+zb(j2,j5)*y13))
     & /(nn**2*(y13-s45p46)*(zb(j5,j6)*zab6341+zb(j1,j5)*y13)**2
     & *(za(j3,j4)*zab7564-za(j3,j7)*y13))
      res=res
     & +(za(j2,j7)*za(j3,j4)*zb(j5,j6)*y13
     & *(za(j3,j5)*y13+za(j3,j4)*zb(j4,j6)*za(j5,j6))
     & *(zb(j4,j5)*y13-zb(j3,j4)*za(j3,j6)*zb(j5,j6))
     & *(za(j5,j6)*zab1346+za(j1,j5)*y13)*zab3564
     & *(zb(j3,j4)*zab4562+zb(j2,j3)*y13))
     & /(nn**2*(y13-s45p46)*(zb(j5,j6)*zab6341+zb(j1,j5)*y13)
     & *(za(j3,j4)*zab7564-za(j3,j7)*y13))
     & +(za(j2,j7)*za(j3,j4)*zb(j3,j4)*za(j5,j6)*zb(j5,j6)*y13
     & *(za(j3,j4)*zab1564+za(j1,j3)*y13)**2*zab3564*zab5346)
     & /(nn**2*(za(j3,j4)*zab2564+za(j2,j3)*y13)
     & *(za(j3,j4)*zab7564-za(j3,j7)*y13))
     & +(za(j2,j7)*za(j3,j4)*y13
     & *(za(j3,j5)*y13+za(j3,j4)*zb(j4,j6)*za(j5,j6))
     & *(zb(j3,j6)*y13-zb(j3,j4)*za(j4,j5)*zb(j5,j6))
     & *(zb(j4,j5)*y13-zb(j3,j4)*za(j3,j6)*zb(j5,j6))
     & *(za(j5,j6)*zab1346+za(j1,j5)*y13)
     & *(za(j3,j4)*zab1564+za(j1,j3)*y13)*zab3564)
     & /(nn**2*(y13-s36p46)*(y13-s45p46)
     & *(za(j3,j4)*zab2564+za(j2,j3)*y13)
     & *(za(j3,j4)*zab7564-za(j3,j7)*y13))
      res=res
     & +(za(j1,j5)*za(j2,j7)*za(j3,j4)*zb(j4,j6)
     & *(za(j3,j4)*zab1564+za(j1,j3)*y13)*zab3564)
     & /((za(j3,j4)*zab2564+za(j2,j3)*y13)
     & *(za(j3,j4)*zab7564-za(j3,j7)*y13))
     & -(za(j2,j7)*za(j3,j4)*zb(j5,j6)*y13
     & *(za(j3,j5)*y13+za(j3,j4)*zb(j4,j6)*za(j5,j6))
     & *(zb(j4,j5)*y13-zb(j3,j4)*za(j3,j6)*zb(j5,j6))
     & *(za(j3,j4)*zab1564+za(j1,j3)*y13)*zab3564
     & *(zb(j3,j4)*zab4567-zb(j3,j7)*y13)
     & *(zb(j5,j6)*zab6342+zb(j2,j5)*y13)
     & *(za(j5,j6)*zab7346-za(j5,j7)*y13))
     & /(nn**2*(y13-s45p46)*(zb(j5,j6)*zab6341+zb(j1,j5)*y13)
     & *(zb(j5,j6)*zab6347-zb(j5,j7)*y13)
     & *(za(j3,j4)*zab7564-za(j3,j7)*y13)**2)
      res=res
     & -(za(j2,j7)*za(j3,j4)*zb(j5,j6)*y13
     & *(za(j3,j5)*y13+za(j3,j4)*zb(j4,j6)*za(j5,j6))
     & *(zb(j4,j5)*y13-zb(j3,j4)*za(j3,j6)*zb(j5,j6))
     & *(za(j3,j4)*zab1564+za(j1,j3)*y13)**2*zab3564
     & *(zb(j3,j4)*zab4567-zb(j3,j7)*y13)
     & *(za(j5,j6)*zab7346-za(j5,j7)*y13))
     & /(nn**2*(y13-s45p46)*(za(j3,j4)*zab2564+za(j2,j3)*y13)
     & *(zb(j5,j6)*zab6347-zb(j5,j7)*y13)
     & *(za(j3,j4)*zab7564-za(j3,j7)*y13)**2)
      res=res
     & -(za(j2,j7)*zb(j3,j4)*za(j5,j6)*y13
     & *(za(j3,j5)*y13+za(j3,j4)*zb(j4,j6)*za(j5,j6))
     & *(zb(j3,j6)*y13-zb(j3,j4)*za(j4,j5)*zb(j5,j6))
     & *(za(j5,j6)*zab1346+za(j1,j5)*y13)
     & *(zb(j3,j4)*zab4562+zb(j2,j3)*y13)*zab5346
     & *(zb(j5,j6)*zab6347-zb(j5,j7)*y13)
     & *(za(j3,j4)*zab7564-za(j3,j7)*y13))
     & /(nn**2*(y13-s36p46)*(zb(j3,j4)*zab4561+zb(j1,j3)*y13)
     & *(zb(j3,j4)*zab4567-zb(j3,j7)*y13)
     & *(za(j5,j6)*zab7346-za(j5,j7)*y13)**2)
      res=res
     & -(za(j2,j7)*zb(j3,j4)*za(j5,j6)*y13
     & *(za(j3,j5)*y13+za(j3,j4)*zb(j4,j6)*za(j5,j6))
     & *(zb(j3,j6)*y13-zb(j3,j4)*za(j4,j5)*zb(j5,j6))
     & *(za(j5,j6)*zab1346+za(j1,j5)*y13)**2*zab5346
     & *(zb(j5,j6)*zab6347-zb(j5,j7)*y13)
     & *(za(j3,j4)*zab7564-za(j3,j7)*y13))
     & /(nn**2*(y13-s36p46)*(za(j5,j6)*zab2346+za(j2,j5)*y13)
     & *(zb(j3,j4)*zab4567-zb(j3,j7)*y13)
     & *(za(j5,j6)*zab7346-za(j5,j7)*y13)**2)
     & +(za(j2,j7)*za(j3,j4)*zb(j3,j4)*za(j5,j6)*y13
     & *(zb(j3,j6)*y13-zb(j3,j4)*za(j4,j5)*zb(j5,j6))
     & *(za(j5,j6)*zab1346+za(j1,j5)*y13)*zab3564
     & *(zb(j3,j4)*zab4562+zb(j2,j3)*y13)*zab5346
     & *(zb(j5,j6)*zab6347-zb(j5,j7)*y13))
     & /(nn**2*(y13-s36p46)*(zb(j3,j4)*zab4561+zb(j1,j3)*y13)
     & *(zb(j3,j4)*zab4567-zb(j3,j7)*y13)
     & *(za(j5,j6)*zab7346-za(j5,j7)*y13))
      res=res
     & +(za(j2,j7)*zb(j3,j4)*za(j5,j6)*y13
     & *(za(j3,j5)*y13+za(j3,j4)*zb(j4,j6)*za(j5,j6))
     & *(zb(j3,j6)*y13-zb(j3,j4)*za(j4,j5)*zb(j5,j6))
     & *(za(j5,j6)*zab1346+za(j1,j5)*y13)
     & *(za(j3,j4)*zab1564+za(j1,j3)*y13)*zab5346
     & *(zb(j5,j6)*zab6347-zb(j5,j7)*y13))
     & /(nn**2*(y13-s36p46)*(za(j5,j6)*zab2346+za(j2,j5)*y13)
     & *(zb(j3,j4)*zab4567-zb(j3,j7)*y13)
     & *(za(j5,j6)*zab7346-za(j5,j7)*y13))
      res=res
     & +(za(j2,j7)*za(j3,j4)*zb(j3,j4)*za(j5,j6)*y13
     & *(zb(j3,j6)*y13-zb(j3,j4)*za(j4,j5)*zb(j5,j6))
     & *(za(j5,j6)*zab1346+za(j1,j5)*y13)**2*zab3564*zab5346
     & *(zb(j5,j6)*zab6342+zb(j2,j5)*y13))
     & /(nn**2*(y13-s36p46)*(za(j5,j6)*zab2346+za(j2,j5)*y13)
     & *(zb(j3,j4)*zab4562+zb(j2,j3)*y13)
     & *(za(j5,j6)*zab7346-za(j5,j7)*y13))
     & -(za(j2,j7)*zb(j3,j4)*za(j5,j6)*y13
     & *(za(j3,j5)*y13+za(j3,j4)*zb(j4,j6)*za(j5,j6))
     & *(zb(j3,j6)*y13-zb(j3,j4)*za(j4,j5)*zb(j5,j6))
     & *(za(j5,j6)*zab1346+za(j1,j5)*y13)**2
     & *(za(j3,j4)*zab2564+za(j2,j3)*y13)*zab5346
     & *(zb(j5,j6)*zab6342+zb(j2,j5)*y13))
     & /(nn**2*(y13-s36p46)*(za(j5,j6)*zab2346+za(j2,j5)*y13)**2
     & *(zb(j3,j4)*zab4562+zb(j2,j3)*y13)
     & *(za(j5,j6)*zab7346-za(j5,j7)*y13))
      res=res
     & +(za(j2,j7)*zb(j3,j4)*za(j5,j6)*y13
     & *(za(j3,j5)*y13+za(j3,j4)*zb(j4,j6)*za(j5,j6))
     & *(zb(j3,j6)*y13-zb(j3,j4)*za(j4,j5)*zb(j5,j6))
     & *(za(j3,j4)*zab1564+za(j1,j3)*y13)*zab5346
     & *(zb(j5,j6)*zab6342+zb(j2,j5)*y13))
     & /(nn**2*(y13-s36p46)*(zb(j3,j4)*zab4561+zb(j1,j3)*y13)
     & *(za(j5,j6)*zab7346-za(j5,j7)*y13))
      res=res
     & -(za(j2,j7)*zb(j3,j4)*za(j5,j6)*y13
     & *(za(j3,j5)*y13+za(j3,j4)*zb(j4,j6)*za(j5,j6))
     & *(zb(j3,j6)*y13-zb(j3,j4)*za(j4,j5)*zb(j5,j6))
     & *(za(j3,j4)*zab1564+za(j1,j3)*y13)
     & *(zb(j3,j4)*zab4562+zb(j2,j3)*y13)*zab5346
     & *(zb(j5,j6)*zab6341+zb(j1,j5)*y13))
     & /(nn**2*(y13-s36p46)*(zb(j3,j4)*zab4561+zb(j1,j3)*y13)**2
     & *(za(j5,j6)*zab7346-za(j5,j7)*y13))
      res=res
     & +(za(j2,j7)*zb(j3,j4)*za(j5,j6)*zb(j5,j6)*y13
     & *(za(j3,j5)*y13+za(j3,j4)*zb(j4,j6)*za(j5,j6))
     & *(za(j3,j4)*zab1564+za(j1,j3)*y13)
     & *(zb(j3,j4)*zab4562+zb(j2,j3)*y13)*zab5346)
     & /(nn**2*(zb(j3,j4)*zab4561+zb(j1,j3)*y13)
     & *(za(j5,j6)*zab7346-za(j5,j7)*y13))
      res=res
     & +(za(j2,j7)*zb(j3,j4)*za(j5,j6)*zb(j5,j6)*y13
     & *(za(j3,j5)*y13+za(j3,j4)*zb(j4,j6)*za(j5,j6))
     & *(za(j5,j6)*zab1346+za(j1,j5)*y13)
     & *(za(j3,j4)*zab1564+za(j1,j3)*y13)*zab5346)
     & /(nn**2*(za(j5,j6)*zab2346+za(j2,j5)*y13)
     & *(za(j5,j6)*zab7346-za(j5,j7)*y13))
      res=res
     & -(za(j2,j7)*zb(j4,j7)*za(j5,j6)
     & *(zb(j3,j6)*y13-zb(j3,j4)*za(j4,j5)*zb(j5,j6))
     & *(za(j3,j5)*zab1346-za(j1,j5)*za(j3,j4)*zb(j4,j6))
     & *(zb(j3,j4)*zab4562+zb(j2,j3)*y13))
     & /((zb(j3,j4)*zab4561+zb(j1,j3)*y13)
     & *(zb(j3,j4)*zab4567-zb(j3,j7)*y13)
     & *(za(j5,j6)*zab7346-za(j5,j7)*y13))
      res=res
     & +(zb(j2,j4)*za(j2,j7)*za(j5,j6)
     & *(zb(j3,j6)*y13-zb(j3,j4)*za(j4,j5)*zb(j5,j6))
     & *(za(j3,j5)*zab1346-za(j1,j5)*za(j3,j4)*zb(j4,j6))
     & *(za(j5,j6)*zab1346+za(j1,j5)*y13))
     & /((za(j5,j6)*zab2346+za(j2,j5)*y13)
     & *(zb(j3,j4)*zab4562+zb(j2,j3)*y13)
     & *(za(j5,j6)*zab7346-za(j5,j7)*y13))
      res=res
     & +(za(j2,j7)*y13*(za(j3,j5)*y13+za(j3,j4)*zb(j4,j6)*za(j5,j6))**2
     & *(zb(j3,j6)*y13-zb(j3,j4)*za(j4,j5)*zb(j5,j6))
     & *(zb(j4,j5)*y13-zb(j3,j4)*za(j3,j6)*zb(j5,j6))
     & *(za(j5,j6)*zab1346+za(j1,j5)*y13)
     & *(zb(j3,j4)*zab4562+zb(j2,j3)*y13))
     & /(nn**2*(y13-s36p46)*(y13-s45p46)
     & *(zb(j3,j4)*zab4561+zb(j1,j3)*y13)
     & *(za(j5,j6)*zab7346-za(j5,j7)*y13))
      res=res
     & +(za(j1,j5)*za(j2,j7)*zb(j4,j6)
     & *(za(j3,j5)*y13+za(j3,j4)*zb(j4,j6)*za(j5,j6))
     & *(zb(j3,j4)*zab4562+zb(j2,j3)*y13))
     & /((zb(j3,j4)*zab4561+zb(j1,j3)*y13)
     & *(za(j5,j6)*zab7346-za(j5,j7)*y13))
      res=res
     & +(za(j2,j7)*y13*(za(j3,j5)*y13+za(j3,j4)*zb(j4,j6)*za(j5,j6))**2
     & *(zb(j3,j6)*y13-zb(j3,j4)*za(j4,j5)*zb(j5,j6))
     & *(zb(j4,j5)*y13-zb(j3,j4)*za(j3,j6)*zb(j5,j6))
     & *(za(j5,j6)*zab1346+za(j1,j5)*y13)**2)
     & /(nn**2*(y13-s36p46)*(y13-s45p46)
     & *(za(j5,j6)*zab2346+za(j2,j5)*y13)
     & *(za(j5,j6)*zab7346-za(j5,j7)*y13))
      res=res
     & +(za(j1,j5)*za(j2,j7)*zb(j4,j6)
     & *(za(j3,j5)*y13+za(j3,j4)*zb(j4,j6)*za(j5,j6))
     & *(za(j5,j6)*zab1346+za(j1,j5)*y13))
     & /((za(j5,j6)*zab2346+za(j2,j5)*y13)
     & *(za(j5,j6)*zab7346-za(j5,j7)*y13))
     & +(za(j1,j5)*za(j3,j4)*zab3564
     & *(za(j1,j7)*zb(j4,j6)*zb(j5,j6)*zab6347
     & +za(j1,j2)*zb(j4,j6)*zb(j5,j6)*zab6342
     & -za(j1,j7)*zb(j4,j6)*zb(j5,j7)*y13
     & +za(j1,j7)*zb(j4,j7)*zb(j5,j6)*y13
     & -za(j1,j2)*zb(j2,j4)*zb(j5,j6)*y13
     & +za(j1,j2)*zb(j2,j5)*zb(j4,j6)*y13
     & -za(j1,j7)*zb(j4,j7)*zb(j5,j6)*s36p46
     & +za(j1,j2)*zb(j2,j4)*zb(j5,j6)*s36p46))
     & /((za(j3,j4)*zab1564+za(j1,j3)*y13)
     & *(zb(j5,j6)*zab6341+zb(j1,j5)*y13))
     & +(za(j3,j4)*zb(j3,j4)*za(j5,j6)*zb(j5,j6)*y13*zab3564*zab5346
     & *(za(j1,j7)*zb(j5,j6)*zab6347+za(j1,j2)*zb(j5,j6)*zab6342
     & -za(j1,j7)*zb(j5,j7)*y13
     & +za(j1,j2)*zb(j2,j5)*y13))
     & /(nn**2*(zb(j5,j6)*zab6341+zb(j1,j5)*y13))
      res=res
     & +(za(j3,j4)*y13*(za(j3,j5)*y13+za(j3,j4)*zb(j4,j6)*za(j5,j6))
     & *(zb(j3,j6)*y13-zb(j3,j4)*za(j4,j5)*zb(j5,j6))
     & *(zb(j4,j5)*y13-zb(j3,j4)*za(j3,j6)*zb(j5,j6))
     & *(za(j5,j6)*zab1346+za(j1,j5)*y13)*zab3564
     & *(za(j1,j7)*zb(j5,j6)*zab6347+za(j1,j2)*zb(j5,j6)*zab6342
     & -za(j1,j7)*zb(j5,j7)*y13
     & +za(j1,j2)*zb(j2,j5)*y13))
     & /(nn**2*(y13-s36p46)*(y13-s45p46)
     & *(za(j3,j4)*zab1564+za(j1,j3)*y13)
     & *(zb(j5,j6)*zab6341+zb(j1,j5)*y13))
      res=res
     & -(za(j3,j4)*zb(j5,j6)*y13
     & *(za(j3,j5)*y13+za(j3,j4)*zb(j4,j6)*za(j5,j6))
     & *(zb(j4,j5)*y13-zb(j3,j4)*za(j3,j6)*zb(j5,j6))
     & *(za(j5,j6)*zab1346+za(j1,j5)*y13)*zab3564
     & *(zb(j3,j4)*zab4561+zb(j1,j3)*y13)
     & *(za(j1,j7)*zb(j5,j6)*zab6347+za(j1,j2)*zb(j5,j6)*zab6342
     & -za(j1,j7)*zb(j5,j7)*y13
     & +za(j1,j2)*zb(j2,j5)*y13))
     & /(nn**2*(y13-s45p46)*(za(j3,j4)*zab1564+za(j1,j3)*y13)
     & *(zb(j5,j6)*zab6341+zb(j1,j5)*y13)**2)
     & +(za(j3,j4)*zb(j3,j4)*za(j5,j6)*y13
     & *(zb(j3,j6)*y13-zb(j3,j4)*za(j4,j5)*zb(j5,j6))*zab3564*zab5346
     & *(za(j1,j7)*zb(j5,j6)*zab6347+za(j1,j2)*zb(j5,j6)*zab6342
     & -za(j1,j7)*zb(j5,j7)*y13
     & +za(j1,j2)*zb(j2,j5)*y13))
     & /(nn**2*(y13-s36p46)*(zb(j3,j4)*zab4561+zb(j1,j3)*y13))
      res=res
     & +(zb(j5,j6)*y13*(za(j3,j5)*y13+za(j3,j4)*zb(j4,j6)*za(j5,j6))**2
     & *(zb(j4,j5)*y13-zb(j3,j4)*za(j3,j6)*zb(j5,j6))
     & *(za(j1,j7)*zb(j3,j4)*zab4567+za(j1,j2)*zb(j3,j4)*zab4562
     & -za(j1,j7)*zb(j3,j7)*y13
     & +za(j1,j2)*zb(j2,j3)*y13))
     & /(nn**2*(y13-s45p46)*(zb(j5,j6)*zab6341+zb(j1,j5)*y13))
     & -(za(j3,j5)*(za(j1,j7)*zb(j4,j7)-za(j1,j2)*zb(j2,j4))
     & *zb(j5,j6)*(y13-s36p46))/(zb(j5,j6)*zab6341+zb(j1,j5)*y13)
     & -(zb(j3,j4)*za(j5,j6)*y13
     & *(za(j3,j5)*y13+za(j3,j4)*zb(j4,j6)*za(j5,j6))
     & *(zb(j3,j6)*y13-zb(j3,j4)*za(j4,j5)*zb(j5,j6))
     & *(za(j3,j4)*zab1564+za(j1,j3)*y13)
     & *(za(j1,j7)*zb(j3,j4)*zab4567+za(j1,j2)*zb(j3,j4)*zab4562
     & -za(j1,j7)*zb(j3,j7)*y13
     & +za(j1,j2)*zb(j2,j3)*y13)*zab5346
     & *(zb(j5,j6)*zab6341+zb(j1,j5)*y13))
     & /(nn**2*(y13-s36p46)*(za(j5,j6)*zab1346+za(j1,j5)*y13)
     & *(zb(j3,j4)*zab4561+zb(j1,j3)*y13)**2)
      res=res
     & +(zb(j3,j4)*za(j5,j6)*zb(j5,j6)*y13
     & *(za(j3,j5)*y13+za(j3,j4)*zb(j4,j6)*za(j5,j6))
     & *(za(j3,j4)*zab1564+za(j1,j3)*y13)
     & *(za(j1,j7)*zb(j3,j4)*zab4567+za(j1,j2)*zb(j3,j4)*zab4562
     & -za(j1,j7)*zb(j3,j7)*y13
     & +za(j1,j2)*zb(j2,j3)*y13)*zab5346)
     & /(nn**2*(za(j5,j6)*zab1346+za(j1,j5)*y13)
     & *(zb(j3,j4)*zab4561+zb(j1,j3)*y13))
      res=res
     & +(za(j1,j5)*(za(j3,j5)*y13+za(j3,j4)*zb(j4,j6)*za(j5,j6))
     & *(za(j1,j7)*zb(j3,j4)*zb(j4,j6)*zab4567
     & +za(j1,j2)*zb(j3,j4)*zb(j4,j6)*zab4562
     & +za(j1,j7)*zb(j3,j6)*zb(j4,j7)*y13
     & -za(j1,j7)*zb(j3,j7)*zb(j4,j6)*y13
     & +za(j1,j2)*zb(j2,j3)*zb(j4,j6)*y13
     & -za(j1,j2)*zb(j2,j4)*zb(j3,j6)*y13
     & -za(j1,j7)*zb(j3,j4)*za(j4,j5)*zb(j4,j7)*zb(j5,j6)
     & +za(j1,j2)*zb(j2,j4)*zb(j3,j4)*za(j4,j5)*zb(j5,j6)))
     & /((za(j5,j6)*zab1346+za(j1,j5)*y13)
     & *(zb(j3,j4)*zab4561+zb(j1,j3)*y13))
     & +(y13*(za(j3,j5)*y13+za(j3,j4)*zb(j4,j6)*za(j5,j6))**2
     & *(zb(j3,j6)*y13-zb(j3,j4)*za(j4,j5)*zb(j5,j6))
     & *(zb(j4,j5)*y13-zb(j3,j4)*za(j3,j6)*zb(j5,j6))
     & *(za(j1,j7)*zb(j3,j4)*zab4567+za(j1,j2)*zb(j3,j4)*zab4562
     & -za(j1,j7)*zb(j3,j7)*y13
     & +za(j1,j2)*zb(j2,j3)*y13))
     & /(nn**2*(y13-s36p46)*(y13-s45p46)
     & *(zb(j3,j4)*zab4561+zb(j1,j3)*y13))
     & -(za(j3,j5)*(za(j1,j7)*zb(j4,j7)-za(j1,j2)*zb(j2,j4))
     & *(zb(j3,j6)*y13-zb(j3,j4)*za(j4,j5)*zb(j5,j6)))
     & /(zb(j3,j4)*zab4561+zb(j1,j3)*y13)
      res=res/(2d0*s(j4,j3)*s(j6,j5)*za(j1,j7)*za(j2,j7))
      return
      end
