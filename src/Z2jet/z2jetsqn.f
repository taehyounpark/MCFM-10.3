!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine z2jetsqn(i1,i2,i5,i6,p,n,za,zb,zab,zba,msq)
      implicit none
      include 'types.f'
c-----Apart from overall factors returns the matrix element squared
c-----msq dependent on the helicities pq and pl of the quark and
c-----lepton lines for
c-----q(-p1)+qbar(-p2)-->l(p3)+al(p4)+g(p5)+g(p6) where
c-----where gluon 6 has been contracted with the vector n
cargument 1-4 represent (i1) incoming quark line
c                       (i2) incoming quark line
c                       (i5) outgoing gluon line
c                       (i6) outgoing gluon line contracted with n
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'mmsqv_cs.f'
      complex(dp):: qcdabn(2,2,2),qcdban(2,2,2),qedn(2,2,2)
      complex(dp):: zab(mxpart,mxpart),zba(mxpart,mxpart)
      real(dp):: msq(2,2),n(4),p(mxpart,4),nDp5
      integer:: i1,i2,i3,i4,i5,i6,pg,pq,pl,icol

      i3=3
      i4=4

      nDp5=n(4)*p(i5,4)-n(3)*p(i5,3)-n(2)*p(i5,2)-n(1)*p(i5,1)

      call checkndotp(p,n,i6)

      call subqcdn(i1,i2,i3,i4,i5,i6,nDp5,za,zb,zab,zba,qcdabn,qcdban)

c--first argument is gluon line
c--second argument is polarization of i5 line pq
c--third argument is polarization of lepton line pl
c  1=L,2=R
      do pq=1,2
      do pl=1,2
      msq(pq,pl)=zip
      do icol=0,2
      mmsqv_cs(icol,pq,pl)=zip
      enddo
      enddo
      enddo

      do pq=1,2
      do pl=1,2
      do pg=1,2
      qedn(pg,pq,pl)=qcdabn(pg,pq,pl)+qcdban(pg,pq,pl)
      mmsqv_cs(0,pq,pl)=mmsqv_cs(0,pq,pl)-ninth*abs(qedn(pg,pq,pl))**2
      mmsqv_cs(1,pq,pl)=mmsqv_cs(1,pq,pl)+abs(qcdabn(pg,pq,pl))**2
      mmsqv_cs(2,pq,pl)=mmsqv_cs(2,pq,pl)+abs(qcdban(pg,pq,pl))**2
      enddo
      msq(pq,pl)=mmsqv_cs(0,pq,pl)+mmsqv_cs(1,pq,pl)+mmsqv_cs(2,pq,pl)
      enddo
      enddo

      return
      end

