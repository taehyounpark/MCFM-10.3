!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function filterWbbmas()
       implicit none
      include 'types.f'
      logical:: filterWbbmas
c--- this routine is specific to the "Wbbmas" processes 401-408;
c--F also for processes 520-529
c--- it inspects the jets to check whether an event should be included
c--- routine returns FALSE if event does not pass the process-specific cuts

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'clustering.f'
      include 'jetlabel.f'
      include 'nproc.f'
      logical:: veto3jets

      if (jets<1) then
         filterWbbmas=.false.
         return
      endif


      filterWbbmas=.true.

c--F  3-jet veto to match with ATLAS analysis. In general we do not want such a veto
      veto3jets = .false.

c---  20 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +b(p5)+b~(p6) [massive]'
      if ((nproc == 20) .or. (nproc == 25)) then
        if (jets < 2) then
          filterWbbmas=.false.
          return
c--- note: jets should already contain bq and ba because of bbproc
        endif
      endif

c--  New final state slicing: Wb, W(bb~), Wbb

c--  401 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+b(p5) [1, 2 or 3 jets, 4FNS]'
      if ((nproc == 401) .or. (nproc == 406)) then
         if ((inclusive .eqv. .false.) .and. (jets  /=  1)) then
            filterWbbmas=.false.
            return
         endif
         if (jets == 1) then
            if((jetlabel(1) /= 'bq').and.(jetlabel(1) /= 'ba')) then
               filterWbbmas=.false.
               return
            endif
         endif
         if (jets == 2) then
            if ((jetlabel(1)=='bb').or.(jetlabel(2)=='bb')) then
               filterWbbmas=.false.
               return
            endif
         endif
         if ((jets == 3) .and. (veto3jets .eqv. .true.)) then
            filterWbbmas=.false.
            return
         endif
      endif

c--  402 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+(b+b~)(p5) [1 or 2 jets, 4FNS]'

      if ((nproc == 402) .or. (nproc == 407)) then
         if ((inclusive .eqv. .false.) .and. (jets  /=  1)) then
            filterWbbmas=.false.
            return
         endif
         if ((jets == 1) .and. (jetlabel(1)  /=  'bb')) then
            filterWbbmas=.false.
            return
         endif
         if (jets == 2) then
            if ((jetlabel(1) /= 'bb').and.(jetlabel(2) /= 'bb')) then
               filterWbbmas=.false.
               return
            endif
         endif
         if (jets == 3) then
            filterWbbmas=.false.
            return
         endif
      endif

c--  403 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+b(p5)+b~(p6) [2 or 3 jets, 4FNS]'
      if ((nproc == 403) .or. (nproc == 408)) then
         if ((jets == 3) .and.
     &       ((inclusive.eqv..false.).or.(veto3jets.eqv..true.))) then
            filterWbbmas=.false.
            return
         endif
         if (jets == 1) then
            filterWbbmas=.false.
            return
         endif
         if (jets == 2) then
            if ((jetlabel(1)=='pp').or.(jetlabel(2)=='pp')) then
               filterWbbmas=.false.
               return
            endif
         endif
      endif


      return
      end


