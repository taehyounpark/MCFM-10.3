!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine I2slash(z,I2slsh)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'zeta.f'
      include 'singletlabels.f'
c     gg=0,qqV=1,qbV=2,qqS=3,qg=4,gq=5,qqDS=6
      include 'distributions.f'
c     dmin=-2,dmax=1,rglr=-2,delt=-1,plus=0,lpls=1
c     -2=finite,delta=-1,0=L0(z),1=L1(z)
      include 'nfl.f'
      real(dp):: tp0(0:6,dmin:dmax),mp0(0:6,dmin:dmax)
      real(dp):: z,mz,I2slsh(0:6,dmin:dmax),Li2,Li3,omz,opz,
     & Li2z,Li2mz,Li2omz,Li3z,Li3mz,Li3omz,Li3zonopz,
     & factorplus,factorplus1,factorrglr,
     & ggfactorplus,ggfactorplus1,ggfactorrglr
      integer,parameter::n1=-1,n2=1
      real(dp):: Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2)

c Explicit expressions for hplog (in terms of Li2, Li3) are now included
c below, but the following may be uncommented for comparison
c     Using Hplog program, hep-ph/0107173
c      call hplog(z,3,Hc1,Hc2,Hc3,Hc4,Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)

      omz=one-z
      opz=one+z
      Li2z=Li2(z); Li2mz=Li2(-z); Li2omz=Li2(omz)
      Li3z=Li3(z); Li3mz=Li3(-z); Li3omz=Li3(omz)
      Li3zonopz=Li3(z/(opz))

      Hr1(0)=(log(z))
      Hr1(-1)=(log(opz))
      Hr1(1)=(-log(omz))

      Hr2(0,0)=(1d0/2d0*log(z)**2)
      Hr2(0,1)=(Li2z)
      Hr2(0,-1)=(-Li2mz)
      Hr2(1,0)=(-log(z)*log(omz)-Li2z)
      Hr2(1,1)=(0.5_dp*log(omz)**2)
      Hr2(-1,0)=(log(z)*log(opz)+Li2mz)

      Hr3(0,0,0)=(1/6d0*log(z)**3)
      Hr3(1,0,0)=(-0.5_dp*log(omz)*log(z)**2-log(z)*Li2z+Li3z)
      Hr3(-1,0,0)=(0.5_dp*log(opz)*log(z)**2+log(z)*Li2mz-Li3mz)
      Hr3(-1,-1,0)=(Li3zonopz-log(z+one)**3/6+(log(z)*log(z+one)**2)/2d0+Li3mz)
      Hr3(0,1,0)=(log(z)*Li2z-2*Li3z)
      Hr3(0,-1,0)=(-log(z)*Li2mz+2*Li3mz)
      Hr3(1,0,1)=(-log(omz)*Li2z
     & -2*(0.5_dp*log(omz)**2*log(z)+log(omz)*Li2omz-Li3omz+zeta3))
      Hr3(1,1,0)=(log(omz)**2*log(z)+log(omz)*Li2omz
     & +log(omz)*Li2z-Li3omz+zeta3)
      Hr3(0,1,1)=(0.5_dp*log(omz)**2*log(z)+log(omz)*Li2omz-Li3omz+zeta3)
      Hr3(1,1,1)=(-1./6d0*log(omz)**3)

c     initialize to zero
      I2slsh(:,:)=0
c     1403.6451v2, Eqns on pages 16-18
      call tildep0(+z,tp0)
      mz=-z
      call tildep0(mz,mp0)
      I2slsh(gg,delt)=CA**2*(25*zeta4/4._dp-77*zeta3/9._dp
     & -67*zeta2/6._dp+1214/81._dp)
     & +CA*TF*nfl*(28*zeta3/9._dp+10*zeta2/3._dp-328/81._dp)
      I2slsh(gg,rglr)=CA**2*(-16*(opz)*Hr3(0,0,0)
     & +(8*(omz)*(11._dp-z+11*z**2))/3._dp/z*(Hr2(1,0)+zeta2)
     & +2*(25._dp-11*z+44*z**2)/3._dp*Hr2(0,0)-2*z/3._dp*Hr1(1)
     & -(701._dp+149*z+536*z**2)/9._dp*Hr1(0)
     & +4*(-196._dp+174*z-186*z**2+211*z**3)/9._dp/z)
     & +CA*TF*nfl*((8*(opz))/3._dp*Hr2(0,0)+4*z/3._dp*Hr1(1)
     & +(4*(13._dp+10*z))/9._dp*Hr1(0)
     & -(4*(-65._dp+54*z-54*z**2+83*z**3))/(27._dp*z))
     &  +CF*TF*nfl*(8*(opz)*Hr3(0,0,0)+4*(3._dp+z)*Hr2(0,0)
     & +24*(opz)*Hr1(0)-(8*(omz)*(one-23*z+z**2))/(3._dp*z))

c Old code
c     tp0(gg) has only plus and rglr components
c      factor=CA**2*(-4*Hr3(0,0,0)+8*Hr3(0,1,0)+8*Hr3(0,1,1)-8*Hr3(1,0,0)
c     & +8*Hr3(1,0,1)+8*Hr3(1,1,0)+52*zeta3-808/27._dp)
c     & +CA*TF*nfl*224/27._dp
c      I2slsh(gg,plus)=I2slsh(gg,plus)+factor*tp0(gg,plus)
c      I2slsh(gg,rglr)=I2slsh(gg,rglr)+factor*tp0(gg,rglr)
c     mp0 has only rglr component
c      I2slsh(gg,rglr)=I2slsh(gg,rglr)
c     & +mp0(gg,rglr)*CA**2*(-16*Hr3(-1,-1,0)+8*Hr3(-1,0,0)+16*Hr3(0,-1,0)
c     & -4*Hr3(0,0,0)-8*Hr3(0,1,0)-8*Hr1(-1)*zeta2+4*zeta3)
c New code -- need to make logs explicit and remove z-dependence

      ggfactorplus =  + CA**2 * (
     &     - 808.D0/27.D0
     &     + 52*zeta3
     &     + 4*log(z)*log(one - z)**2
     &     + 4*log(z)**2*log(one - z)
     &     - 2.D0/3.D0*log(z)**3
     &     + 16*log(z)*Li2(z)
     &     - 24*Li3(z)
     &     )
      ggfactorplus = ggfactorplus + TF*nfl*CA * (
     &     + 224.D0/27.D0
     &     )

      ggfactorrglr =  + CA**2 * (
     &     + 1616.D0/27.D0
     &     - 4/z/(opz)*zeta3
     &     - 808.D0/27.D0/z
     &     + 52/z*zeta3
     &     - 8/(opz)*zeta3
     &     - 104*zeta3
     &     - 12*z/(opz)*zeta3
     &     - 808.D0/27.D0*z
     &     + 52*z*zeta3
     &     - 8*z**2/(opz)*zeta3
     &     + 808.D0/27.D0*z**2
     &     - 52*z**2*zeta3
     &     - 4*z**3/(opz)*zeta3
     &     + 4.D0/3.D0*log(one + z)/z*pi**2/(opz)
     &     )
      ggfactorrglr = ggfactorrglr + CA**2 * (
     &     + 8.D0/3.D0*log(one + z)*pi**2/(opz)
     &     + 4*log(one + z)*z*pi**2/(opz)
     &     + 8.D0/3.D0*log(one + z)*z**2*pi**2/(opz)
     &     + 4.D0/3.D0*log(one + z)*z**3*pi**2/(opz)
     &     - 8.D0/3.D0*log(one + z)**3/z/(opz)
     &     - 16.D0/3.D0*log(one + z)**3/(opz)
     &     - 8*log(one + z)**3*z/(opz)
     &     - 16.D0/3.D0*log(one + z)**3*z**2/(opz)
     &     - 8.D0/3.D0*log(one + z)**3*z**3/(opz)
     &     - 8*log(z)*log(one - z)**2
     &     + 4*log(z)*log(one - z)**2/z
     &     + 4*log(z)*log(one - z)**2*z
     &     - 4*log(z)*log(one - z)**2*z**2
     &     + 8*log(z)*log(one + z)**2/z/(opz)
     &     )
      ggfactorrglr = ggfactorrglr + CA**2 * (
     &     + 16*log(z)*log(one + z)**2/(opz)
     &     + 24*log(z)*log(one + z)**2*z/(opz)
     &     + 16*log(z)*log(one + z)**2*z**2/(opz)
     &     + 8*log(z)*log(one + z)**2*z**3/(opz)
     &     - 8*log(z)**2*log(one - z)
     &     + 4*log(z)**2*log(one - z)/z
     &     + 4*log(z)**2*log(one - z)*z
     &     - 4*log(z)**2*log(one - z)*z**2
     &     - 4*log(z)**2*log(one + z)/z/(opz)
     &     - 8*log(z)**2*log(one + z)/(opz)
     &     - 12*log(z)**2*log(one + z)*z/(opz)
     &     - 8*log(z)**2*log(one + z)*z**2/(opz)
     &     - 4*log(z)**2*log(one + z)*z**3/(opz)
     &     + 4.D0/3.D0*log(z)**3
     &     )
      ggfactorrglr = ggfactorrglr + CA**2 * (
     &     + 2.D0/3.D0*log(z)**3/z/(opz)
     &     - 2.D0/3.D0*log(z)**3/z
     &     + 4.D0/3.D0*log(z)**3/(opz)
     &     + 2*log(z)**3*z/(opz)
     &     - 2.D0/3.D0*log(z)**3*z
     &     + 4.D0/3.D0*log(z)**3*z**2/(opz)
     &     + 2.D0/3.D0*log(z)**3*z**2
     &     + 2.D0/3.D0*log(z)**3*z**3/(opz)
     &     + 8*log(z)*Li2( - z)/z/(opz)
     &     + 16*log(z)*Li2( - z)/(opz)
     &     + 24*log(z)*Li2( - z)*z/(opz)
     &     + 16*log(z)*Li2( - z)*z**2/(opz)
     &     + 8*log(z)*Li2( - z)*z**3/(opz)
     &     - 32*log(z)*Li2(z)
     &     )
      ggfactorrglr = ggfactorrglr + CA**2 * (
     &     + 8*log(z)*Li2(z)/z/(opz)
     &     + 16*log(z)*Li2(z)/z
     &     + 16*log(z)*Li2(z)/(opz)
     &     + 24*log(z)*Li2(z)*z/(opz)
     &     + 16*log(z)*Li2(z)*z
     &     + 16*log(z)*Li2(z)*z**2/(opz)
     &     - 16*log(z)*Li2(z)*z**2
     &     + 8*log(z)*Li2(z)*z**3/(opz)
     &     - 8*Li3( - z)/z/(opz)
     &     - 16*Li3( - z)/(opz)
     &     - 24*Li3( - z)*z/(opz)
     &     - 16*Li3( - z)*z**2/(opz)
     &     - 8*Li3( - z)*z**3/(opz)
     &     + 16*Li3(1/(one + z)*z)/z/(opz)
     &     )
      ggfactorrglr = ggfactorrglr + CA**2 * (
     &     + 32*Li3(1/(one + z)*z)/(opz)
     &     + 48*Li3(1/(one + z)*z)*z/(opz)
     &     + 32*Li3(1/(one + z)*z)*z**2/(opz)
     &     + 16*Li3(1/(one + z)*z)*z**3/(opz)
     &     + 48*Li3(z)
     &     - 16*Li3(z)/z/(opz)
     &     - 24*Li3(z)/z
     &     - 32*Li3(z)/(opz)
     &     - 48*Li3(z)*z/(opz)
     &     - 24*Li3(z)*z
     &     - 32*Li3(z)*z**2/(opz)
     &     + 24*Li3(z)*z**2
     &     - 16*Li3(z)*z**3/(opz)
     &     )
      ggfactorrglr = ggfactorrglr + TF*nfl*CA * (
     &     - 448.D0/27.D0
     &     + 224.D0/27.D0/z
     &     + 224.D0/27.D0*z
     &     - 224.D0/27.D0*z**2
     &     )

      ggfactorplus1 =  + CA**2 * (
     &     - 808.D0/27.D0
     &     + 28*zeta3
     &     )
      ggfactorplus1 = ggfactorplus1 + TF*nfl*CA * (
     &     + 224.D0/27.D0
     &     )

      I2slsh(gg,plus)=ggfactorplus1

      I2slsh(gg,rglr)=I2slsh(gg,rglr)
     & +ggfactorrglr
     & +(ggfactorplus-ggfactorplus1)/(omz)


      I2slsh(gq,rglr)=CF*CA*(-4*(2._dp+z)*Hr3(0,0,0)+16*Hr3(0,1,0)
     & +4*z*Hr2(-1,0)+4*z*Hr2(0,1)+4*z*Hr2(1,1)
     & -(8*(one+z+2*z**2))/3._dp*Hr2(1,0)-22*z/3._dp*Hr1(1)
     & +(2*(36._dp+9*z+8*z**2))/3._dp*Hr2(0,0)
     & -(2*(249._dp-6*z+88*z**2))/9._dp*Hr1(0)
     & -8*zeta3-(2*(4._dp+13*z+8*z**2))/3._dp*zeta2
     & +(4*(one+127*z+152*z**2))/27._dp)
     & +CF**2*(2*(2._dp-z)*Hr3(0,0,0)-(4._dp+3*z)*Hr2(0,0)-4*z*Hr2(1,1)
     & +6*z*Hr1(1)-5*(3._dp-z)*Hr1(0)+(10._dp-z))
     & +CF*TF*nfl*(8*z/3._dp*Hr1(1)-40*z/9._dp)

c     tp0(gq) and mp0(gq) have only the rglr component
      I2slsh(gq,rglr)=I2slsh(gq,rglr)
     & +tp0(gq,rglr)*(CF*CA*(4*Hr3(1,1,1)+4*Hr3(0,1,1)+4*Hr3(1,0,1)
     & +4*Hr3(1,1,0)+8*Hr3(0,1,0)-4*Hr3(1,0,0)+44/3._dp*(Hr2(1,0)+zeta2)
     & -22/3._dp*Hr2(1,1)+152/9._dp*Hr1(1)+24*zeta3-1580/27._dp)
     & +CF**2*(-4*Hr3(1,1,1)+6*Hr2(1,1)-16*Hr1(1))
     & +CF*TF*nfl*(8/3._dp*Hr2(1,1)-40/9._dp*Hr1(1)+224/27._dp))
     & +mp0(gq,rglr)*CF*CA
     & *(-8*Hr3(-1,-1,0)+4*Hr3(-1,0,0)+8*Hr3(0,-1,0)-4*Hr1(-1)*zeta2)

      I2slsh(qqS,rglr)=CF*TF*(4*(opz)*Hr3(0,0,0)
     & -(2*(3._dp+3*z+8*z**2))/3._dp*Hr2(0,0)
     & -(8*(omz)*(2._dp-z+2*z**2))/(3._dp*z)*(Hr2(1,0)+zeta2)
     & +(4*(21._dp-30*z+32*z**2))/9._dp*Hr1(0)
     & +(2*(omz)*(172._dp-143*z+136*z**2))/(27._dp*z))

      I2slsh(qqV,delt)=CF*CA*(5*zeta4-77*zeta3/9._dp-67*zeta2/6._dp
     & +1214/81._dp)+CF**2*5/4._dp*zeta4
     & +CF*TF*nfl*(28*zeta3/9._dp+10*zeta2/3._dp-328/81._dp)

      I2slsh(qqV,rglr)=CF*CA*(-4*(omz)*Hr2(1,0)-4*z*Hr2(0,0)-2*z*Hr1(1)
     & +2*(1._dp+5*z)*Hr1(0)-6*(omz)*zeta2+44/3._dp*(omz))
     & +CF**2*(2*(opz)*Hr3(0,0,0)+(3._dp+7*z)*Hr2(0,0)+4*(omz)*Hr2(0,1)
     & +12*(omz)*Hr2(1,0)+2*z*Hr1(1)+2*(1._dp-12*z)*Hr1(0)
     & +6*(omz)*zeta2-22*(omz))
     & -CF*TF*nfl*4/3._dp*(omz)
c     tp0(qqV) has only plus and rglr components
c Old code
c      factor=CF*CA*(-2*Hr3(0,0,0)-4*Hr3(0,1,0)-4*Hr3(1,0,1)-4*Hr3(1,1,0)
c     & -11/3._dp*Hr2(0,0)-76/9._dp*Hr1(0)+2*zeta3-404/27._dp)
c     & +CF**2*(8*Hr3(0,1,0)+4*Hr3(0,1,1)-4*Hr3(1,0,0)+8*Hr3(1,0,1)
c     & +8*Hr3(1,1,0)+3*Hr2(0,0)+8*Hr1(0)+24*zeta3)
c     & +CF*TF*nfl*(4/3._dp*Hr2(0,0)+20/9._dp*Hr1(0)+112/27._dp)
c      I2slsh(qqV,plus)=I2slsh(qqV,plus)+factor*tp0(qqV,plus)
c      I2slsh(qqV,rglr)=I2slsh(qqV,rglr)+factor*tp0(qqV,rglr)
c New code -- need to make logs explicit and remove z-dependence
      factorplus =  + CF*CA * (
     &     - 808.D0/27.D0
     &     + 12*zeta3
     &     + 8*log(omz)*Li2omz
     &     - 152.D0/9.D0*log(z)
     &     - 11.D0/3.D0*log(z)**2
     &     - 2.D0/3.D0*log(z)**3
     &     - 8*log(z)*Li2z
     &     - 8*Li3omz
     &     + 16*Li3z
     &     )
      factorplus = factorplus + CF*TF*nfl * (
     &     + 224.D0/27.D0
     &     + 40.D0/9.D0*log(z)
     &     + 4.D0/3.D0*log(z)**2
     &     )
      factorplus = factorplus + CF**2 * (
     &     + 40*zeta3
     &     - 8*log(omz)*Li2omz
     &     + 16*log(z)
     &     + 4*log(z)*log(omz)**2
     &     + 3*log(z)**2
     &     + 4*log(z)**2*log(omz)
     &     + 24*log(z)*Li2z
     &     + 8*Li3omz
     &     - 40*Li3z
     &     )

      factorrglr =  + CF*CA * (
     &     + 404.D0/27.D0
     &     - 6*zeta3
     &     + 404.D0/27.D0*z
     &     - 6*z*zeta3
     &     - 4*log(omz)*Li2omz
     &     - 4*log(omz)*Li2omz*z
     &     + 76.D0/9.D0*log(z)
     &     + 76.D0/9.D0*log(z)*z
     &     + 11.D0/6.D0*log(z)**2
     &     + 11.D0/6.D0*log(z)**2*z
     &     + 1.D0/3.D0*log(z)**3
     &     + 1.D0/3.D0*log(z)**3*z
     &     + 4*log(z)*Li2z
     &     + 4*log(z)*Li2z*z
     &     )
      factorrglr = factorrglr + CF*CA * (
     &     + 4*Li3omz
     &     + 4*Li3omz*z
     &     - 8*Li3z
     &     - 8*Li3z*z
     &     )
      factorrglr = factorrglr + CF*TF*nfl * (
     &     - 112.D0/27.D0
     &     - 112.D0/27.D0*z
     &     - 20.D0/9.D0*log(z)
     &     - 20.D0/9.D0*log(z)*z
     &     - 2.D0/3.D0*log(z)**2
     &     - 2.D0/3.D0*log(z)**2*z
     &     )
      factorrglr = factorrglr + CF**2 * (
     &     - 20*zeta3
     &     - 20*z*zeta3
     &     + 4*log(omz)*Li2omz
     &     + 4*log(omz)*Li2omz*z
     &     - 8*log(z)
     &     - 8*log(z)*z
     &     - 2*log(z)*log(omz)**2
     &     - 2*log(z)*log(omz)**2*z
     &     - 3.D0/2.D0*log(z)**2
     &     - 3.D0/2.D0*log(z)**2*z
     &     - 2*log(z)**2*log(omz)
     &     - 2*log(z)**2*log(omz)*z
     &     - 12*log(z)*Li2z
     &     - 12*log(z)*Li2z*z
     &     )
      factorrglr = factorrglr + CF**2 * (
     &     - 4*Li3omz
     &     - 4*Li3omz*z
     &     + 20*Li3z
     &     + 20*Li3z*z
     &     )

      factorplus1 =  + CF*CA * (
     &     - 808.D0/27.D0
     &     + 28*zeta3
     &     )
      factorplus1 = factorplus1 + CF*TF*nfl * (
     &     + 224.D0/27.D0
     &     )

      I2slsh(qqV,plus)=factorplus1

      I2slsh(qqV,rglr)=I2slsh(qqV,rglr)
     & +factorrglr
     & +(factorplus-factorplus1)/(omz)

c Note that there is an error in the expression for I(2){q/q} on pp.17-18
c of 1403.6451, from which the above expression has been taken.
c This term, +I(2){q'/q}, has been omitted, but it is
c present in both the accompanying FORM file (TPDFsAtNNLO.frm) and the
c original paper, Ref. [36] ibid (Eq. (12) of 1209.0682)
c However, in the end it is not needed here, c.f. comment below
c      I2slsh(qqV,rglr)=I2slsh(qqV,rglr)+I2slsh(qqS,rglr)

      I2slsh(qg,rglr)=CA*TF*(4*(one+2*z)*Hr3(0,0,0)-16*z*Hr3(0,1,0)
     & +2*(19._dp-32*z)/3._dp*Hr2(0,0)-4*Hr2(-1,0)-4*Hr2(1,1)
     & -4*(4._dp+5*z+2*z**2)/3/z*(Hr2(1,0)+zeta2)
     & +2*(-2._dp+z)*Hr1(1)-4*(13._dp-38*z)/9._dp*Hr1(0)+8*z*(zeta3+zeta2)
     & +2*(172._dp-166*z+89*z**2)/(27*z))
     & +CF*TF*(2*(1._dp-2*z)*Hr3(0,0,0)+(5._dp+4*z)*Hr2(0,0)+4*Hr2(0,1)
     & +4*Hr2(1,0)+4*Hr2(1,1)+2*(2._dp-z)*Hr1(1)+(12._dp+7*z)*Hr1(0)-6*zeta2
     & +(23._dp+3*z))
c     tp0(qg) and mp0(qg) have only the 2 component
      I2slsh(qg,rglr)=I2slsh(qg,rglr)+tp0(qg,rglr)
     & *(CA*TF*(4*Hr3(1,0,1)+4*Hr3(1,1,0)-4*Hr3(1,1,1)+4*Hr2(1,1)
     & -44/3._dp*Hr2(0,0)+44/3._dp*(Hr2(1,0)+zeta2)+136/9._dp*Hr1(0)
     & +4*Hr1(1)-298/27._dp)
     & +CF*TF*(4*Hr3(1,1,1)-4*Hr3(1,0,0)+4*Hr3(0,1,1)-4*Hr3(0,0,0)
     & -4*Hr2(1,1)-4*Hr2(1,0)-4*Hr2(0,1)-4*Hr2(0,0)-4*Hr1(1)-4*Hr1(0)
     & +28*zeta3+6*zeta2-36._dp))
     & +mp0(qg,rglr)
     & *CA*TF*(-8*Hr3(-1,-1,0)+4*Hr3(-1,0,0)+8*Hr3(0,-1,0)+4*Hr2(-1,0)
     & -4*Hr1(-1)*zeta2)

      I2slsh(qbV,rglr)=CF*(CA-2*CF)
     & *(4*(omz)*Hr2(1,0)+4*(opz)*Hr2(-1,0)-(3._dp+11*z)*Hr1(0)
     & +2*(3._dp-z)*zeta2-15*(omz))
      I2slsh(qbV,rglr)=I2slsh(qbV,rglr)+CF*(CA-2*CF)*mp0(qqV,rglr)
     & *(8*Hr3(-1,-1,0)-4*Hr3(-1,0,0)-8*Hr3(0,-1,0)+4*Hr3(0,1,0)
     & +2*Hr3(0,0,0)+4*Hr1(-1)*zeta2-2*zeta3)
c This term is not needed, c.f. comment below
c     & +I2slsh(qqS,rglr)

c as written above, the qqV and qbV entries are actually for
c quarks of the same flavor, i.e. they correspond to qqV+qqS and qbV+qqS
c      I2slsh(qqV,:)=I2slsh(qqV,:)-I2slsh(qqS,:)
c      I2slsh(qbV,:)=I2slsh(qbV,:)-I2slsh(qqS,:)
c Rather than perform this subtraction, these terms have just been removed above

c      return

c      write(6,*)'Hr1(0)',Hr1(0)-(log(z))
c      write(6,*)'Hr1(-1)',Hr1(-1)-(log(1+z))
c      write(6,*)'Hr1(+1)',Hr1(1)-(-log(1-z))

c      write(6,*)'Hr2(0,0)',Hr2(0,0)-(0.5_dp*log(z)**2)
c      write(6,*)'Hr2(0,+1)',Hr2(0,1)-(Li2(z))
c      write(6,*)'Hr2(0,-1)',Hr2(0,-1)-(-Li2(-z))
c      write(6,*)'Hr2(+1,0)',Hr2(1,0)-(-log(z)*log(1-z)-Li2(z))
c      write(6,*)'Hr2(+1,+1)',Hr2(1,1)-(0.5_dp*log(1-z)**2)
c      write(6,*)'Hr2(-1,0)',Hr2(-1,0)-(log(z)*log(1+z)+Li2(-z))

c      write(6,*)'Hr3(0,0,0)',
c     & Hr3(0,0,0)-(1/6d0*log(z)**3)
c      write(6,*)'Hr3(1,0,0)',
c     & Hr3(1,0,0)-(-0.5_dp*log(1-z)*log(z)**2-log(z)*Li2(z)+Li3(z))
c      write(6,*)'Hr3(-1,0,0)',
c     & Hr3(-1,0,0)-(0.5_dp*log(1+z)*log(z)**2+log(z)*Li2(-z)-Li3(-z))
c      write(6,*)'Hr3(-1,-1,0)',
c     & Hr3(-1,-1,0)
c     & -(Li3(z/(z+1))-log(z+1)**3/6+(log(z)*log(z+1)**2)/2d0+Li3(-z))
c      write(6,*)'Hr3(0,1,0)',
c     & Hr3(0,1,0)-(log(z)*Li2(z)-2*Li3(z))
c      write(6,*)'Hr3(0,-1,0)',
c     & Hr3(0,-1,0)-(-log(z)*Li2(-z)+2*Li3(-z))
c      write(6,*)'Hr3(1,0,1)',
c     & Hr3(1,0,1)-(-log(1-z)*Li2(z)
c     & -2*(0.5_dp*log(1-z)**2*log(z)+log(1-z)*Li2(1-z)-Li3(1-z)+zeta3))
c      write(6,*)'Hr3(1,1,0)',
c     & Hr3(1,1,0)-(log(1-z)**2*log(z)+log(1-z)*Li2(1-z)
c     & +log(1-z)*Li2(z)-Li3(1-z)+zeta3)
c      write(6,*)'Hr3(0,1,1)',Hr3(0,1,1)
c     & -(0.5_dp*log(1-z)**2*log(z)+log(1-z)*Li2(1-z)-Li3(1-z)+zeta3)
c      write(6,*)'Hr3(1,1,1)',
c     & Hr3(1,1,1)-(-1/6d0*log(1-z)**3)

      return
      end

