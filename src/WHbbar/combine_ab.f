!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine combine_ab(pjet,i,j)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'jetlabel.f'
      integer:: i,j
      real(dp):: pjet(mxpart,4)

c--Run II prescription
      pjet(i,1)=pjet(i,1)+pjet(j,1)
      pjet(i,2)=pjet(i,2)+pjet(j,2)
      pjet(i,3)=pjet(i,3)+pjet(j,3)
      pjet(i,4)=pjet(i,4)+pjet(j,4)

c--- special combination tag for W+heavy quarks
c      if ((kcase==kWbbmas) .or. (kcase==kW_bjet)
c     & .or.(kcase==kWHbbar) .or. (kcase==kZHbbar)
c     &) then
c      if (((jetlabel(i) == 'bq') .and. (jetlabel(j) == 'ba'))
c     &.or.((jetlabel(j) == 'bq') .and. (jetlabel(i) == 'ba'))) then
c        jetlabel(i)='bb'
c        return
c      endif
c      if (((jetlabel(i) == 'bb') .and. (jetlabel(j) == 'pp'))
c     &.or.((jetlabel(j) == 'bb') .and. (jetlabel(i) == 'pp'))) then
c        jetlabel(i)='bb'
c        return
c      endif
c      endif

      if (((jetlabel(i) == 'qb') .and. (jetlabel(j) == 'pp'))
     &.or.((jetlabel(j) == 'qb') .and. (jetlabel(i) == 'pp'))) then
        jetlabel(i)='qb'
        return
      endif
      if (((jetlabel(i) == 'ab') .and. (jetlabel(j) == 'pp'))
     &.or.((jetlabel(j) == 'ab') .and. (jetlabel(i) == 'pp'))) then
        jetlabel(i)='ab'
        return
      endif
      if (((jetlabel(i) == 'qb') .and. (jetlabel(j) == 'ab'))
     &.or.((jetlabel(j) == 'qb') .and. (jetlabel(i) == 'ab'))) then
        jetlabel(i)='pp'
        return
      endif
      if (((jetlabel(i) == 'qb') .and. (jetlabel(j) == 'qj'))
     &.or.((jetlabel(j) == 'qb') .and. (jetlabel(i) == 'qj'))) then
        jetlabel(i)='qb'
        return
      endif
      if (((jetlabel(i) == 'ab') .and. (jetlabel(j) == 'qj'))
     &.or.((jetlabel(j) == 'ab') .and. (jetlabel(i) == 'qj'))) then
        jetlabel(i)='ab'
        return
      endif
      if (((jetlabel(i) == 'qj') .and. (jetlabel(j) == 'pp'))
     &.or.((jetlabel(j) == 'pp') .and. (jetlabel(i) == 'qj'))) then
        jetlabel(i)='qj'
        return
      endif


      if (((jetlabel(i) == 'bq') .and. (jetlabel(j) == 'pp'))
     &.or.((jetlabel(j) == 'bq') .and. (jetlabel(i) == 'pp'))) then
        jetlabel(i)='bq'
        return
      endif
      if (((jetlabel(i) == 'ba') .and. (jetlabel(j) == 'pp'))
     &.or.((jetlabel(j) == 'ba') .and. (jetlabel(i) == 'pp'))) then
        jetlabel(i)='ba'
        return
      endif
      if (((jetlabel(i) == 'bq') .and. (jetlabel(j) == 'ba'))
     &.or.((jetlabel(j) == 'bq') .and. (jetlabel(i) == 'ba'))) then
        jetlabel(i)='pp'
        return
      endif
      if (((jetlabel(i) == 'bq') .and. (jetlabel(j) == 'qj'))
     &.or.((jetlabel(j) == 'bq') .and. (jetlabel(i) == 'qj'))) then
        jetlabel(i)='bq'
        return
      endif
      if (((jetlabel(i) == 'ba') .and. (jetlabel(j) == 'qj'))
     &.or.((jetlabel(j) == 'ba') .and. (jetlabel(i) == 'qj'))) then
        jetlabel(i)='ba'
        return
      endif
      if (((jetlabel(i) == 'qj') .and. (jetlabel(j) == 'pp'))
     &.or.((jetlabel(j) == 'pp') .and. (jetlabel(i) == 'qj'))) then
        jetlabel(i)='qj'
        return
      endif



      return
      end

