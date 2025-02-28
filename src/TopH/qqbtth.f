!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqbtth(denr,denb,wtqqb)
      implicit none
      include 'types.f'

c***********************************************************************
c     Author: R.K. Ellis                                               *
c     December, 1999.                                                  *
c     calculate the element squared and subtraction terms              *
c     for the process                                                  *
c     My notation                                                      *
c     q(-p1) +qbar(-p2)=bbar(p6)+e-(p7)+nubar(p8)+b(p5)+nu(p3)+e+(p4)  *
c     +b(p9)+bbar(p10)                                                 *
c                                                                      *
c q=3+4+5
c r=3+4+5+h
c a=-6-7-8
c b=-6-7-8-h
c since momenta 3,5,6,8,9,10 are not needed they are re-used to represent
c the the de-massified vectors. Thus q4 = q de-massified wrt p4 etc: this
c hase been stored in wrapper routine (qqb_tth) in p(3,mu)
c q4(mu)=q(mu)-qsq/2/p4Dq*p4(mu)-->p(3,mu)
c a7(mu)=a(mu)-asq/2/p7Da*p7(mu)-->p(5,mu)
c r1(mu)=r(mu)-rsq/2/p1Dr*p1(mu)-->p(6,mu)
c r2(mu)=r(mu)-rsq/2/p2Dr*p2(mu)-->p(8,mu)
c b1(mu)=b(mu)-bsq/2/p1Db*p1(mu)-->p(9,mu)
c b2(mu)=b(mu)-bsq/2/p2Db*p2(mu)-->p(10,mu)
c***********************************************************************
      include 'mxpart.f'
      include 'masses.f'
      include 'zprods_com.f'
      integer:: q4,a7,r1,r2,b1,b2
      parameter(q4=3,a7=5,r1=6,r2=8,b1=9,b2=10)
      real(dp):: wtqqb,denr,denb
      complex(dp):: ampp,ampm
c--diagram 1-Polarity plus
      ampp=mt/denr*((zb(4,q4)*za(q4,r2)*zb(r2,2)+mt**2*zb(4,2))*za(1,7)
     & +(zb(4,q4)*za(q4,1)+zb(4,r1)*za(r1,1))*zb(2,a7)*za(a7,7))
c--diagram 2-Polarity plus
      ampp=ampp
     & +mt/denb*(zb(4,2)*(za(1,b1)*zb(b1,a7)*za(a7,7)+mt**2*za(1,7))
     & +zb(4,q4)*za(q4,1)*(zb(2,b2)*za(b2,7)+zb(2,a7)*za(a7,7)))

c--diagram 1-Polarity minus
      ampm=mt/denr*((zb(4,q4)*za(q4,r1)*zb(r1,1)+mt**2*zb(4,1))*za(2,7)
     & +(zb(4,q4)*za(q4,2)+zb(4,r2)*za(r2,2))*zb(1,a7)*za(a7,7))
c--diagram 2-Polarity minus
      ampm=ampm
     & +mt/denb*(zb(4,1)*(za(2,b2)*zb(b2,a7)*za(a7,7)+mt**2*za(2,7))
     & +zb(4,q4)*za(q4,2)*(zb(1,b1)*za(b1,7)+zb(1,a7)*za(a7,7)))

      wtqqb=(abs(ampp)**2+abs(ampm)**2)/abs(za(1,2)*zb(2,1))**2

      return
      end
