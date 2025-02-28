!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

c--- File written by FORM program adecayW_v.frm on Thu Mar  1 12:30:14 CST 2012
      subroutine adecayW_v(p,pqq,pqb,pc,m)
      implicit none
      include 'types.f'

c***********************************************************************
c     Author: R.K. Ellis, January 2012                                 *
c     Virtual corrections to                                           *
c     antitop decay  a --> qq(pqq)+qb(pqb)+bbar(pc)                    *
c     with bottom and top masses (and no radiation)                    *
c     in massless spinor notation                                      *
c     pqq,pqb,pc are integer::s that point to                            *
c     the appropriate four-momenta in p                                *
c     pqq=electron                                                     *
c     pqb=antineutrino                                                 *
c     pc=anti-bottom quark                                             *
c     q(c) is rendered massless wrt to pqb                             *
c     q(a) is rendered massless wrt to pqq                             *
c     returned m(apol,cpol)                                            *
c***********************************************************************
      include 'constants.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      real(dp):: p(mxpart,4),q(mxpart,4),dot,sw,ala,alc,c0,ct
      complex(dp):: m(2,2),cprop,c0L
      integer:: qb,a,c,qq,si,pc,pqq,pqb
      parameter(a=1,qq=3,qb=4,c=2)
      do si=1,4
      q(a,si)=p(pqq,si)+p(pqb,si)+p(pc,si)
      q(qq,si)=p(pqq,si)
      q(qb,si)=p(pqb,si)
      q(c,si)=p(pc,si)
      enddo
      ala=mt**2/(2._dp*dot(q,a,qq))
      alc=mb**2/(2._dp*dot(q,c,qb))
      do si=1,4
      q(a,si)=q(a,si)-ala*q(qq,si)
      q(c,si)=q(c,si)-alc*q(qb,si)
      enddo
      call spinoru(4,q,za,zb)
      sw=s(qq,qb)
      call coefswdk(sw,ct,c0)
      c0L=cplx1(c0+ct)
      cprop=cplx2(sw-wmass**2,wmass*wwidth)
c---order of polarizations is m(apol,cpol)
      m(1,1)= czip

      m(1,2)= czip

      m(2,1)= czip

      m(2,2)= + cprop**(-1) * ( za(qq,a)*zb(qb,c)*c0L )

      return
      end
