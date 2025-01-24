!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine hqqgg(p1,p2,p3,p4,ampsq,ampsq_ab,ampsq_ba,ampsq_sym)
      implicit none
      include 'types.f'

c     Taken from Kauffman,Desai,Risal
c     PRD 55 1997 (4009)
c     and checked with hep-ph/9903330
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      integer:: p1,p2,p3,p4,j1,j2,j3
      complex(dp):: ab(2,2,2),ba(2,2,2),abppp,abppm,bappm
c      complex(dp):: abmmm,abmmp,bammp
      real(dp):: ampsq,ampsq_ab,ampsq_ba,ampsq_sym
c====statement functions
c  Eq 25
      abppp(p1,p2,p3,p4)=
     & +(za(p2,p1)*zb(p1,p3)+za(p2,p4)*zb(p4,p3))**2*zb(p1,p4)
     & /((s(p1,p2)+s(p1,p4)+s(p2,p4))*za(p2,p4))
     & *(1._dp/s(p1,p2)+1._dp/s(p1,p4))
     & -(za(p2,p1)*zb(p1,p4)+za(p2,p3)*zb(p3,p4))**2
     & *zb(p1,p3)/((s(p1,p2)+s(p1,p3)+s(p2,p3))*s(p1,p2)*za(p2,p3))
     & +(+za(p2,p3)*zb(p3,p1)+za(p2,p4)*zb(p4,p1))**2
     & /(zb(p1,p2)*za(p2,p3)*za(p2,p4)*za(p3,p4))

c      abmmm(p1,p2,p3,p4)=
c     & +(zb(p2,p1)*za(p1,p3)+zb(p2,p4)*za(p4,p3))**2
c     & *za(p1,p4)/((s(p1,p2)+s(p1,p4)+s(p2,p4))
c     & *zb(p2,p4))*(1._dp/s(p1,p2)+1._dp/s(p1,p4))
c     & -(+zb(p2,p1)*za(p1,p4)+zb(p2,p3)*za(p3,p4))**2
c     & *za(p1,p3)/((s(p1,p2)+s(p1,p3)+s(p2,p3))*s(p1,p2)*zb(p2,p3))
c     & +(+zb(p2,p3)*za(p3,p1)+zb(p2,p4)*za(p4,p1))**2
c     & /(za(p1,p2)*zb(p2,p3)*zb(p2,p4)*zb(p3,p4))
c--Eq 26
      abppm(p1,p2,p3,p4)=
     & -za(p2,p4)**3/(za(p1,p2)*za(p2,p3)*za(p3,p4))
     & +zb(p1,p3)**3/(zb(p1,p2)*zb(p1,p4)*zb(p3,p4))
c      abmmp(p1,p2,p3,p4)=
c     & -zb(p2,p4)**3/(zb(p1,p2)*zb(p2,p3)*zb(p3,p4))
c     & +za(p1,p3)**3/(za(p1,p2)*za(p1,p4)*za(p3,p4))

c--Eq 27
      bappm(p1,p2,p3,p4)=
     & -zb(p1,p3)**2*zb(p2,p3)/(zb(p1,p2)*zb(p2,p4)*zb(p3,p4))
     & +za(p1,p4)*za(p2,p4)**2/(za(p1,p2)*za(p1,p3)*za(p3,p4))

c      bammp(p1,p2,p3,p4)=
c     & -za(p1,p3)**2*za(p2,p3)/(za(p1,p2)*za(p2,p4)*za(p3,p4))
c     & +zb(p1,p4)*zb(p2,p4)**2/(zb(p1,p2)*zb(p1,p3)*zb(p3,p4))
c====end statement functions

c====It has been checked that taking the complex conjugate
c====only gets the answer different by an overall (and hence irrelevant)
c====phase (in some crossings).


cORIG      ab(2,2,2)=abppp(p1,p2,p3,p4)
cORIGc      ab(1,1,1)=abmmm(p1,p2,p3,p4)
cORIG      ab(1,1,1)=conjg(ab(2,2,2))
cORIG
cORIG
cORIG      ba(2,2,2)=abppp(p1,p2,p4,p3)
cORIGc      ba(1,1,1)=abmmm(p1,p2,p4,p3)
cORIG      ba(1,1,1)=conjg(ba(2,2,2))
cORIG
cORIG      ab(1,2,2)=abppp(p2,p1,p3,p4)
cORIGc      ab(2,1,1)=abmmm(p2,p1,p3,p4)
cORIG      ab(2,1,1)=conjg(ab(1,2,2))
cORIG
cORIG      ba(1,2,2)=abppp(p2,p1,p4,p3)
cORIGc      ba(2,1,1)=abmmm(p2,p1,p4,p3)
cORIG      ba(2,1,1)=conjg(ba(1,2,2))
cORIG
cORIG      ab(2,2,1)=abppm(p1,p2,p3,p4)
cORIGc      ab(1,1,2)=abmmp(p1,p2,p3,p4)
cORIG      ab(1,1,2)=conjg(ab(2,2,1))
cORIG
cORIG      ba(2,1,2)=abppm(p1,p2,p4,p3)
cORIGc      ba(1,2,1)=abmmp(p1,p2,p4,p3)
cORIG      ba(1,2,1)=conjg(ba(2,1,2))
cORIG
cORIG      ba(2,2,1)=bappm(p1,p2,p3,p4)
cORIGc      ba(1,1,2)=bammp(p1,p2,p3,p4)
cORIG      ba(1,1,2)=conjg(ba(2,2,1))
cORIG
cORIG      ab(2,1,2)=bappm(p1,p2,p4,p3)
cORIGc      ab(1,2,1)=bammp(p1,p2,p4,p3)
cORIG      ab(1,2,1)=conjg(ab(2,1,2))
cORIG
cORIG



      ba(2,2,2)=abppp(p1,p2,p3,p4)
      ba(1,1,1)=conjg(ba(2,2,2))

      ab(2,2,2)=abppp(p1,p2,p4,p3)
      ab(1,1,1)=conjg(ab(2,2,2))

      ab(1,2,2)=abppp(p2,p1,p3,p4)
      ab(2,1,1)=conjg(ab(1,2,2))

      ba(1,2,2)=abppp(p2,p1,p4,p3)
      ba(2,1,1)=conjg(ba(1,2,2))

      ba(2,2,1)=abppm(p1,p2,p3,p4)
      ba(1,1,2)=conjg(ba(2,2,1))

      ab(2,1,2)=abppm(p1,p2,p4,p3)
      ab(1,2,1)=conjg(ab(2,1,2))

      ab(2,2,1)=bappm(p1,p2,p3,p4)
      ab(1,1,2)=conjg(ab(2,2,1))

      ba(2,1,2)=bappm(p1,p2,p4,p3)
      ba(1,2,1)=conjg(ba(2,1,2))


c--- calculate the matrix element as the sum of two colour orderings,
c--- plus a colour-suppressed QE.e-_dplike piece which is symmetric
c--- in the ordering of the two gluons
      ampsq_ab=0._dp
      ampsq_ba=0._dp
      ampsq_sym=0._dp
      do j1=1,2
      do j2=1,2
      do j3=1,2
c         write(*,*) 'j1,j2,j3',j1,j2,j3
c         write(*,*) cf*xn**2/2._dp*abs(ab(j1,j2,j3))**2,
c     &      cf*xn**2/2._dp*abs(ba(j1,j2,j3))**2,
c     &        -cf/2._dp*abs(ab(j1,j2,j3)+ba(j1,j2,j3))**2
      ampsq_ab=ampsq_ab
     &  +cf*xn**2/2._dp*abs(ab(j1,j2,j3))**2
      ampsq_ba=ampsq_ba
     &  +cf*xn**2/2._dp*abs(ba(j1,j2,j3))**2
      ampsq_sym=ampsq_sym
     &  -cf/2._dp*abs(ab(j1,j2,j3)+ba(j1,j2,j3))**2
      enddo
      enddo
      enddo

      ampsq=ampsq_ab+ampsq_ba+ampsq_sym

      return
      end
