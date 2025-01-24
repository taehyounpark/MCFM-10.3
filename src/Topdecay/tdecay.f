!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine tdecay(p,pqq,pqb,pb,m)
      implicit none
      include 'types.f'

c***********************************************************************
c     Author: R.K. Ellis, January 2012                                 *
c     top decay  t --> q(pqq)+qb(pb)+b(pb)                             *
c     with bottom and top masses (and no radiation)                    *
c     in massless spinor notation                                      *
c     pe,pnb,pc are integer::s that point to                             *
c     the appropriate four-momenta in p                                *
c     pqq=quark                                                        *
c     pqb=antiquark                                                    *
c     pb=bottom quark                                                  *
c     q(t) is rendered massless wrt to pqb                             *
c     q(b) is rendered massless wrt to pqq                             *
c     returned m(bpol,tpol)                                            *
c***********************************************************************
      include 'constants.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      real(dp):: p(mxpart,4),q(mxpart,4),dot,sw,alt,alb
      complex(dp):: m(2,2),cprop
      integer:: qq,qb,b,t,si,pb,pqq,pqb
      parameter(t=1,qq=3,qb=4,b=2)
      do si=1,4
      q(t,si)=p(pqq,si)+p(pqb,si)+p(pb,si)
      q(qq,si)=p(pqq,si)
      q(qb,si)=p(pqb,si)
      q(b,si)=p(pb,si)
      enddo
      alt=mt**2/(2._dp*dot(q,t,qb))
      alb=mb**2/(2._dp*dot(q,b,qq))
      do si=1,4
      q(t,si)=q(t,si)-alt*q(qb,si)
      q(b,si)=q(b,si)-alb*q(qq,si)
      enddo
      call spinoru(4,q,za,zb)
      sw=s(qq,qb)
      cprop=cplx2(sw-wmass**2,wmass*wwidth)
c---order of polarizations is the m(bpol,tpol)
      m(1,1)= - za(qq,b)*zb(qb,t)*cprop**(-1)

      m(1,2)=czip

      m(2,1)=czip

      m(2,2)=czip

      return
      end
