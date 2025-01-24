!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine AP1(z,p1)
      implicit none
c     See, for example, 1405.1044v2, Eq.(A.7,A.8)
      include 'types.f'
      include 'constants.f'
      include 'zeta.f'
      include 'qtconstants.f'
      include 'nfl.f'
      include 'singletlabels.f'
c     gg=0,qqV=1,qbV=2,qqS=3,qg=4,gq=5,qqDS=6

      include 'distributions.f'
c     dmin=-2,dmax=1,rglr=-2,delt=-1,plus=0,lpls=1
c     general definition of functions
c     -2 distribution-free piece
c     -1 delta function
c      0 coefficient of L0(z)
c      1 coefficient of L1(z)
      real(dp), intent(in) :: z
      real(dp) :: omz,opz,opzsq,p1(0:6,dmin:dmax),
     & fpgg,fpggmz,fpgq,fpgqmz,fpqg,fpqgmz,s2z,splits2,s2zesw,ddilog

c     Initialize all to zero
      p1(:,:)=0

      s2z = splits2(z)
      omz=1._dp-z
      opz=1._dp+z
      opzsq=1._dp+z**2
c gg
c This code now deprecated in favor of results from ESW below
c      pggA(delt) = CA*(-1._dp+3*zeta3)+beta0
c      pggA(plus) = Gamma1(0)/CA/4.0_dp            !Coeff of 1/(1-z)

c      if (z == 1._dp) then
c         pggA(rglr)=0
c      else
c         fpgg = 2*(omz+z**2)**2/z/omz
c         fpggmz = -2*(opz+z**2)**2/z/opz
c!         s2z = splits2(z)

c         pggA(rglr) = CA*(fpgg*(-2*log(omz)+log(z)/2.0_dp)*log(z)
c     &            +fpggmz*(s2z+log(z)**2/2.0_dp)
c     &            +4*opz*log(z)**2
c     &            -4/3._dp*(9._dp+11*z**2)*log(z)
c     &            -277/18.0_dp/z + 19*omz + 277/18.0_dp*z**2 )
c     &            +beta0*(13/6.0_dp/z-1.5_dp*omz
c     &            -13/6.0_dp*z**2+opz*log(z))
c         pggA(rglr) = pggA(rglr)+Gamma1(0)/CA/4.0_dp*omz/z
c      endif

c      pggF(delt) = -CF
c      pggF(plus) = 0.0_dp
c      pggF(rglr) = CF*(4/3._dp/z-16._dp+8*z+20/3._dp*z**2
c     &              -2*opz*log(z)**2 - 2*(3._dp+5*z)*log(z) )

c      p1(gg,:) = 4*(CA*pggA(:) + TR*nfl*pggF(:))

c--- Taken from ESW Page 111
c      fpgg = 2*(omz+z**2)**2/z/omz
c      fpggmz = -2*(opz+z**2)**2/z/opz

      fpgg=1._dp/omz+1._dp/z-2._dp+z*omz
      fpggmz=1._dp/opz-1._dp/z-2._dp-z*opz
      S2zesw=-2*ddilog(-z)+log(z)**2/2._dp-2*log(z)*log(1._dp+z)-pi**2/6

      p1(gg,delt)=+CA*TR*nfl*(-4._dp/3._dp)+CA**2*(+8._dp/3._dp+3*zeta3)
     & +CF*TR*nfl*(-1)
      p1(gg,lpls)=zip
      p1(gg,plus)=-20._dp/9._dp*CA*TR*nfl+CA**2*(67._dp/9._dp-1._dp/3._dp*pi**2)
      p1(gg,rglr)=CF*TR*nfl*(-16._dp+8*z+20._dp/3._dp*z**2+4._dp/3._dp/z
     & -(6._dp+10*z)*log(z)-(2._dp+2*z)*log(z)**2)
     & +CA*TR*nfl*(2._dp-2*z+26._dp/9._dp*(z**2-1._dp/z)
     & -4._dp/3._dp*(1._dp+z)*log(z)-20._dp/9._dp*fpgg)
     & +CA**2*(27._dp/2*(1._dp-z)+67._dp/9._dp*(z**2-1._dp/z)
     & -(25._dp/3._dp-11._dp/3._dp*z+44._dp/3._dp*z**2)*log(z)
     & +4*(1._dp+z)*log(z)**2+2*fpggmz*S2zesw
     & +(67._dp/9._dp-4*log(z)*log(1._dp-z)+log(z)**2-pi**2/3._dp)*fpgg)
     & -p1(gg,plus)/(1._dp-z)
c---

      fpgq=(1._dp+omz**2)/z
      fpgqmz=(1._dp+opz**2)/(-z)
      fpqg=omz**2+z**2
      fpqgmz=opz**2+z**2

cgq
      p1(gq,rglr) =
     & CA*(fpgq*(log(omz)**2-2*log(omz)*log(z)-101/18.0_dp-pisq/6.0_dp)
     &  +fpgqmz*s2z + 2*z*log(omz)+(2._dp+z)*log(z)**2
     &  -(36._dp+15*z+8*z**2)/3._dp*log(z)
     &  +(56._dp-z+88*z**2)/18.0_dp )
     &  -CF*(fpgq*log(omz)**2+(3*fpgq+2*z)*log(omz)
     &   +(2._dp-z)/2.0_dp*log(z)**2
     &   -(4._dp+7*z)/2.0_dp*log(z) + (5._dp+7*z)/2.0_dp )
     &  +beta0*(fpgq*(log(omz)+5/3._dp)+z)
       p1(gq,rglr)=p1(gq,rglr)*CF

cqg
      p1(qg,rglr)=
     & CF*(fpqg*(log((omz)/z)**2-2*log(omz/z)-pisq/3._dp+5._dp)
     &             + 2*log(omz) - (1._dp-2*z)/2.0_dp*log(z)**2
     &             - (1._dp-4*z)/2.0_dp*log(z) + 2._dp - 9/2.0_dp*z )
     &        +CA*(fpqg*(-log(omz)**2+2*log(omz)+22/3._dp*log(z)
     &             -109/9._dp+pisq/6.0_dp) +fpqgmz*s2z
     &             -2*log(omz)-(1._dp+2*z)*log(z)**2
     &             +(68*z-19._dp)/3._dp*log(z)
     &             +20/9._dp/z+91/9._dp+7/9._dp*z )
      p1(qg,rglr)=p1(qg,rglr)*TR

cqqV
c delta function
      p1(qqV,delt) = CF*(3/8.0_dp-pisq/2.0_dp+6*zeta3)
     &           +CA*(1/4.0_dp-3*zeta3)
     &           +beta0*(1/8.0_dp+pisq/6.0_dp)
c L0 piece
      p1(qqV,plus) = Gamma1(1)/CF/4.0_dp  
      if (z == 1._dp) then
         p1(qqV,rglr)=0
      else
         p1(qqV,rglr) = -CF*((opzsq)/omz*(2*log(omz)+3/2.0_dp)*log(z)
     &              +opz/2.0_dp*log(z)**2
     &              +(3._dp+7*z)/2.0_dp*log(z)+5*omz )
     &              +CA*(1/2.0_dp*(opzsq)/omz*log(z)**2
     &              +opz*log(z) + 3*omz )
     &              +beta0*(1/2.0_dp*(opzsq)/omz*log(z)+omz)
         p1(qqV,rglr) = p1(qqV,rglr)-Gamma1(1)/CF/8._dp*(1._dp+z)
      endif


cqbV
      p1(qbV,rglr)=(2*CF-CA)*(opzsq/opz*(s2z+log(z)**2/2.0_dp)
     &                    +opz*log(z)+2*(omz)  )

cqqS
      p1(qqS,rglr)=TR*(-opz*log(z)**2+(1._dp+5*z+8/3._dp*z**2)*log(z)
     &            +20/9._dp/z-2._dp+6*z-56/9._dp*z**2 )

      p1(qqV,:)=p1(qqV,:)*CF
      p1(qbV,:)=p1(qbV,:)*CF
      p1(qqS,:)=p1(qqS,:)*CF

c     correct overall normalization for (as/4/pi)
      p1(:,:)=4*p1(:,:)
      return
      end
