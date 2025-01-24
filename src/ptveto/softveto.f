!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine softveto(j,order,pt,R,S)
!     returns renormalized soft function from 2204.02987v1 in exponential regularization
!     S = 1+ason4pi*S(1)+ason4pi^2*S(2)
!     j=0 -- gluon      
!     j=1 -- quark
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'zeta.f'
      include 'Lnu.f'
      include 'scale.f'
      include 'taucut.f'
      real(dp)::pt,R,mu,S(0:2),CR,softperp,softab,softnab,softferm,
     & ln2,ln4,ln8,ln64
      real(dp),parameter::Li4_half=0.5174790616738991_dp
      integer::j,order

      if (j == 0) then
         CR=CA
      elseif (j==1) then
         CR=CF
      else
         write(6,*) 'softveto:j must be 0 or 1'
         stop
      endif
      ln2=log(2._dp)
      ln4=2*ln2
      ln8=3*ln2
      ln64=6*ln2

! Set mu from renormalization scale
      mu=scale

      S(0)=1
      S(1)=-CR*(8*log(pt/mu)*(2*Lnu+log(pt/mu))+pisq/3._dp)

      if (order == 1) return

      softperp=(-2*CR*nf*TR*(328-45*pisq+576*log(pt/mu)**3+
     & 72*log(pt/mu)*(pisq-20*Lnu)+672*Lnu+
     & 144*log(pt/mu)**2*(-5+6*Lnu)+36*zeta3))/81+
     & (CA*CR*(4856-603*pisq+18*pi**4+6336*log(pt/mu)**3+
     & 144*log(pt/mu)**2*(-67+3*pisq+66*Lnu)+
     & 72*log(pt/mu)*(11*pisq+4*(-67+3*pisq)*Lnu)+
     & Lnu*(9696-9072*zeta3)+396*zeta3))/162+
     & (CR**2*(pi**4+576*log(pt/mu)**4+2304*log(pt/mu)**3*Lnu+
     & 48*log(pt/mu)**2*(pisq+48*Lnu**2)+
     & 1152*Lnu*zeta3+96*log(pt/mu)*
     & (pisq*Lnu+12*zeta3)))/18

      softab=((-(CR**2*R**2*(295680*R**4+2744*R**6+149*R**8+7741440*R**2*(-1+ln2)-
     & 15482880*zeta3))/60480-
     & (128*CR**2*(-log(pt/mu)-Lnu)*(4*pisq*R**2-3*(R**4+16*zeta3)))/3))/32

      softnab=(((103744*CA*CR)/9-(4384*CA*CR*pisq)/27-(1504*CA*CR*pi**4)/15+
     & (2948264*CA*CR*R**2)/10125-(1108793639*CA*CR*R**4)/83349000+
     & (100116883541d0*CA*CR*R**6)/60011280000d0-(153495372888217d0*CA*CR*R**8)/
     & 2556000437760000d0-(206848*CA*CR*ln2)/27-
     & (1408*CA*CR*pisq*ln2)/9-(187048*CA*CR*R**2*ln2)/675+
     & (1706191*CA*CR*R**4*ln2)/198450-(25763593*CA*CR*R**6*ln2)/47628000+
     & (48471314251d0*CA*CR*R**8*ln2)/1106493696000d0-(18688*CA*CR*ln2**2)/9-
     & (1024*CA*CR*pisq*ln2**2)/3-(3328*CA*CR*R**2*ln2**2)/45-
     & (4696*CA*CR*R**4*ln2**2)/945+(184*CA*CR*R**6*ln2**2)/525-
     & (28529*CA*CR*R**8*ln2**2)/1871100-(5632*CA*CR*ln2**3)/9+
     & (1024*CA*CR*ln2**4)/3-(34304*CA*CR*ln2*ln8)/27-
     & (5632*CA*CR*ln2**2*ln8)/9+8192*CA*CR*Li4_half-
     & 64*CA*CR*R**2*zeta3+(8*CA*CR*R**4*zeta3)/9+
     & (CA*CR*R**8*zeta3)/8100+7168*CA*CR*ln2*zeta3+
     & (128*CA*CR*log(R)*(-811+822*ln2+396*ln2**2+
     & (786-72*pisq-792*ln2)*(-log(pt/mu)-Lnu)+108*zeta3))/27-
     & (CA*CR*(-log(pt/mu)-Lnu)*(91069440*pisq*(-3801600-518400*R**2+7200*R**4+
     & R**8)+R**8*(50937246539d0-33741818880d0*ln2)-
     & 13113999360d0*R**2*(1429+12480*ln2)-
     & 5575680*R**4*(-9383279+1972320*ln2)+
     & 23232*R**6*(-74801417+33384960*ln2)-5245599744000d0*
     & (-1610-411*ln4+6*ln2*(277+66*ln4)+216*zeta3)))/
     & 1106493696000d0))/32

      softferm=(TR*((-20800*CR*nf)/9+(928*CR*nf*pisq)/27-(104264*CR*nf*R**2)/10125+
     & (13532213*CR*nf*R**4)/83349000+(651413779*CR*nf*R**6)/60011280000d0-
     & (30571386962947d0*CR*nf*R**8)/5112000875520000d0+(40960*CR*nf*ln2)/27+
     & (256*CR*nf*pisq*ln2)/9+(7048*CR*nf*R**2*ln2)/675-
     & (57311*CR*nf*R**4*ln2)/99225-(1272767*CR*nf*R**6*ln2)/47628000+
     & (2692634479d0*CR*nf*R**8*ln2)/553246848000d0+(4864*CR*nf*ln2**2)/9+
     & (448*CR*nf*R**2*ln2**2)/45-(848*CR*nf*R**4*ln2**2)/945+
     & (88*CR*nf*R**6*ln2**2)/1575-(4001*CR*nf*R**8*ln2**2)/935550+
     & (1024*CR*nf*ln2**3)/9+(6656*CR*nf*ln2*ln8)/27+
     & (512*CR*nf*ln2*ln8*ln64)/27-
     & (CR*nf*(31473598464000d0*pisq+R**8*(5664846191d0-4732062720d0*ln2)+
     & 6556999680d0*R**2*(-3071+1680*ln2)-5575680*R**4*
     & (-168401+178080*ln2)+11616*R**6*(-7001023+5322240*ln2)+
     & 2622799872000d0*(-314-87*ln4+6*ln2*(61+12*ln4)))*
     & (-log(pt/mu)-Lnu))/553246848000d0-
     & (128*CR*nf*log(R)*(-163+174*ln2+72*ln2**2-
     & 6*(-23+24*ln2)*(-log(pt/mu)-Lnu)))/27))/16

!      write(6,*) 'softperp',softperp
!      write(6,*) 'softab  ',softab
!      write(6,*) 'softnab ',softnab
!      write(6,*) 'softferm',softferm

      S(2)=softperp+softab+softnab+softferm

      return
      end
