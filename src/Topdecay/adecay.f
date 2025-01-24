!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine adecay(p,pe,pnb,pc,m)
      implicit none
      include 'types.f'

c***********************************************************************
c     Author: R.K. Ellis, January 2012                                 *
c     antitop decay  a --> e(pe)+nb(pnb)+bbar(pc)                      *
c     with bottom and top masses (and no radiation)                    *
c     in massless spinor notation                                      *
c     pe,pnb,pc are integer::s that point to                             *
c     the appropriate four-momenta in p                                *
c     pe=electron                                                      *
c     pnb=antineutrino                                                 *
c     pc=anti-bottom quark                                             *
c     q(c) is rendered massless wrt to pnb                             *
c     q(a) is rendered massless wrt to pe                              *
c     returned m(apol,cpol)                                            *
c***********************************************************************
      include 'constants.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      real(dp):: p(mxpart,4),q(mxpart,4),dot,sw,ala,alc
      complex(dp):: m(2,2),cprop
      integer:: nb,a,c,e,si,pc,pe,pnb
      parameter(a=1,e=3,nb=4,c=2)
      do si=1,4
      q(a,si)=p(pe,si)+p(pnb,si)+p(pc,si)
      q(e,si)=p(pe,si)
      q(nb,si)=p(pnb,si)
      q(c,si)=p(pc,si)
      enddo
      ala=mt**2/(2._dp*dot(q,a,e))
      alc=mb**2/(2._dp*dot(q,c,nb))
      do si=1,4
      q(a,si)=q(a,si)-ala*q(e,si)
      q(c,si)=q(c,si)-alc*q(nb,si)
      enddo

      call spinoru(4,q,za,zb)
      sw=s(e,nb)
      cprop=cplx2(sw-wmass**2,wmass*wwidth)
c---order of polarizations is the m(apol,cpol)
      m(1,1)=czip
      m(1,2)=czip
      m(2,1)=czip
      m(2,2)=-za(e,a)*zb(c,nb)/cprop
      return
      end
