!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine tdecayg(p,pqq,pqb,pb,pg,m)
      implicit none
      include 'types.f'

c***********************************************************************
c     Author: R.K. Ellis, January 2012                                 *
c     top decay  t --> q(pqq)+qb(pb)+b(pb)+g(pg)                       *
c     with bottom and top masses (and radiation form t-b line)         *
c     in massless spinor notation                                      *
c     pqq,pqb,pb,pg are integer::s that point to                         *
c     the appropriate four-momenta in p                                *
c     pqq=quark                                                        *
c     pqb=antiquark                                                    *
c     pb=bottom quark                                                  *
c     pg=gluon                                                         *
c     q(t) is rendered massless wrt to pqb                             *
c     q(b) is rendered massless wrt to pqq                             *
c     returned m(bpol,gpol,tpol)                                       *
c***********************************************************************
      include 'constants.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      real(dp):: p(mxpart,4),q(mxpart,4),dot,sw,twopbDg,twoptDg,
     & altg,alb
      complex(dp):: m(2,2,2),cprop
      integer:: qq,qb,b,g,tg,nu,pb,pqq,pqb,pg
      parameter(tg=1,g=2,qq=3,qb=4,b=5)
      do nu=1,4
      q(tg,nu)=p(pqq,nu)+p(pqb,nu)+p(pb,nu)+p(pg,nu)
      q(g,nu)=p(pg,nu)
      q(qq,nu)=p(pqq,nu)
      q(qb,nu)=p(pqb,nu)
      q(b,nu)=p(pb,nu)
      enddo
      twoptDg=2._dp*dot(q,tg,g)
      twopbDg=2._dp*dot(q,b,g)
      altg=mt**2/(2._dp*dot(q,tg,qb))
      alb=mb**2/(2._dp*dot(q,b,qq))
      do nu=1,4
      q(tg,nu)=q(tg,nu)-altg*q(qb,nu)
      q(b,nu)=q(b,nu)-alb*q(qq,nu)
      enddo
      call spinoru(5,q,za,zb)
      sw=s(qq,qb)
      cprop=cplx2(sw-wmass**2,wmass*wwidth)
c---order of polarizations is the m(bpol,gpol,tpol)
      m(1,1,1)= + cprop**(-1)*twoptDg**(-1) * ( za(qq,b)*za(g,tg)*zb(qb
     &    ,tg)*zb(b,tg)/zb(b,g) )
      m(1,1,1) = m(1,1,1) + cprop**(-1)*twopbDg**(-1) * ( za(qq,g)*za(b
     &    ,g)*zb(qb,tg) )
      m(1,1,1) = m(1,1,1) + cprop**(-1)*twopbDg**(-1)*mb**2 * ( za(qq,g
     &    )*zb(qb,tg)/zb(b,g) )
      m(1,1,1) = m(1,1,1) + cprop**(-1)*mt**2*twoptDg**(-1) * (  - za(
     &    qq,b)*za(qb,g)*zb(qb,b)/za(qb,tg)/zb(b,g) )

      m(1,2,1)= + cprop**(-1)*twoptDg**(-1) * (  - za(qq,b)*za(b,tg)*
     &    zb(qb,tg)*zb(g,tg)/za(b,g) + za(qq,b)*zb(qb,g)*zb(g,tg) )
      m(1,2,1) = m(1,2,1) + cprop**(-1)*twopbDg**(-1)*mb**2 * (  - za(
     &    qq,b)*zb(qq,g)*zb(qb,tg)/za(b,g)/zb(qq,b) )
      m(1,2,1) = m(1,2,1) + cprop**(-1)*mt**2*twoptDg**(-1) * ( za(qq,b
     &    )*za(qb,b)*zb(qb,g)/za(qb,tg)/za(b,g) )

      m(2,1,1)= + cprop**(-1)*twopbDg**(-1)*mb * ( za(qq,g)**2*zb(qb,tg
     &    )/za(qq,b) )

      m(2,2,1)= czip

      m(1,1,2)= czip

      m(1,2,2)= + cprop**(-1)*mt*twoptDg**(-1) * ( za(qq,b)*zb(qb,g)**2
     &    /zb(qb,tg) )

      m(2,1,2)= czip

      m(2,2,2)= czip

      return
      end
