!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c***********************************************************************
c This routine checks the singularity cancellation, given the two      *
c appropriate routines - real_g, sub_gs - and the momentum array p.    *
c A given sub-process (j,k) is determined to be singular if            *
c  (npart=4) msq(j,k) > 1d3 * (smallest value of msq(j,k))             *
c  (npart=5) msq(j,k) > 1d4 * (smallest value of msq(j,k))             *
c This seems to be a reasonable criterion for this set of singular     *
c points. The singularity is considered uncancelled if                 *
c  msq(j,k)/(sum of dipoles) is more than 2% different from 1          *
c***********************************************************************
      subroutine singcheck(real_g,sub_gs,p)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'nflav.f'
      include 'npart.f'
      include 'ptilde.f'
      include 'kprocess.f'
      include 'jetlabel.f'
      include 'nqcdjets.f'
      include 'incldip.f'
      include 'notag.f'
      include 'Rcut.f'
      integer:: j,k,jj,kk,jmax,kmax,nd,icount
      real(dp):: p(mxpart,4),q(mxpart,4),pjet(mxpart,4),
     & msq(-nf:nf,-nf:nf),msqs(-nf:nf,-nf:nf),msqc(maxd,-nf:nf,-nf:nf),
     & s(mxpart,mxpart),debugsmall,debuglarge,xtoler,smin,dot
      character(len=32):: debugmsg
      logical:: includedipole
      external real_g,sub_gs
      data icount/0/
      save icount

      if     (npart == 3) then
        jmax=5
c        if (kcase==kdirgam) then
c          kmax=5
c        elseif (kcase==ktwojet) then
c          kmax=9
        if ((kcase==ktt_tot) .or. (kcase==kcc_tot)
     &   .or. (kcase==kbb_tot)) then
          kmax=4
        else
        jmax=1
        kmax=1
        xtoler=1d3
c        smin=min(-dot(p,1,5),-dot(p,2,5))
c        if (smin. gt. 1.e-3_dp) return
        endif
        xtoler=1d3
      elseif (npart == 4) then
        jmax=5
        kmax=5
        xtoler=1d3
        if ( (kcase==kqg_tbq) .or.(kcase==kqq_tbg)
     &   .or.(kcase==kepem3j) ) then
          jmax=1
          kmax=1
        endif
        jmax=1
        kmax=1
      elseif (npart == 5) then
        if ((kcase==kt_bbar) .or. (kcase==kbq_tpq)) then
        jmax=5
        kmax=6
        xtoler=1d3
        elseif (kcase==kW_bjet) then
          jmax=1
          kmax=1
          xtoler=1d3
        else
        jmax=1
        kmax=1
        xtoler=1d3
        smin=min(-dot(p,1,5),-dot(p,2,5))
        smin=min(smin,-dot(p,1,6))
        smin=min(smin,-dot(p,2,6))
        smin=min(smin,-dot(p,1,7))
        smin=min(smin,-dot(p,2,7))
        smin=min(smin,+dot(p,5,6))
        smin=min(smin,+dot(p,5,7))
        smin=min(smin,+dot(p,6,7))
c        if (smin. gt. 1.e-1_dp) return
        endif
      elseif (npart == 6) then
        continue
      elseif (npart == 7) then
        jmax=5
        kmax=5
        xtoler=1d3
        if ((kcase==ktt_bbl) .or.(kcase==ktt_bbh)) then
          jmax=1
          kmax=1
        endif
      else
        write(6,*) 'No singularity check implemented yet'
        stop
      endif

      jmax=1
      kmax=1

      do j=1,jmax
      do k=1,kmax
         if     (npart == 3) then
           if ((kcase==ktt_tot) .or. (kcase==kcc_tot)
     &   .or. (kcase==kbb_tot)) then
c             call coll3m(p,k,j)
           elseif (kcase==kdirgam) then
c             call coll3gam(p,k,j)
           elseif (kcase==ktwojet) then
c             call coll3jet(p,k,j)
           else
c             call coll3(p,k,j)
           endif
         elseif (npart == 4) then
           if     (kcase==kW_cjet) then
c             call coll4c(p,k,j)
           elseif (kcase==kW_tndk) then
             if (k == 3) goto 68
c             call coll4t(p,k,j)
           elseif ( (kcase==kqg_tbq) .or. (kcase==kqq_tbg)
     &         .or. (kcase==kepem3j) ) then
             write(6,*) 'p(5,4)=',p(5,4)
             write(6,*) 'p(6,4)=',p(6,4)
           else
c             call coll4a(p,k,j)
           endif
         elseif (npart == 5) then
           if ((kcase==kt_bbar) .or. (kcase==kbq_tpq)) then
c             call coll_stop(p,k,j)
           elseif (kcase==kW_bjet) then
             write(6,*) 'p(5,4)=',p(5,4)
             write(6,*) 'p(6,4)=',p(6,4)
             write(6,*) 'p(7,4)=',p(7,4)
           else
c             call coll5(p,k,j)
           endif
         elseif (npart == 7) then
           if ((kcase==ktt_bbl) .or. (kcase==ktt_bbh)) then
             write(6,*) 'p(9,4)=',p(9,4)
           else
c             call coll5(p,k,j)
           endif
         endif
         msq = 0._dp
         if (includedipole(0,p)) then
             call real_g(p,msq)
         endif

         call sub_gs(p,msqc)

         msqs = 0._dp
         call dotem(9,p,s)
         call genclust2(p,Rcut,pjet,0)
c         write(6,*) 'rl,pt34',0,pttwo(3,4,p)
c         write(6,*) 'nd,jets, msq(5,0)',0,nd,jets,msq(5,0)
         if (jets == -1) then
           write(6,*) 'This point does not have a final state b-jet'
           msq = 0._dp
           goto 68
         endif
         if (jets < nqcdjets-notag) then
c            write(6,*) 'This point does not have enough jets'
           msq = 0._dp
           goto 68
         endif

         write(*,*) 'Point ',j,':'
         write(6,*) '2._dp*p1.p4',s(1,4)
         write(6,*) '2._dp*p2.p4',s(2,4)
         if (npart > 2) then
         write(6,*) '2._dp*p1.p5',s(1,5)
         write(6,*) '2._dp*p2.p5',s(2,5)
         write(6,*) '2._dp*p3.p4',s(3,4)
         write(6,*) '2._dp*p3.p5',s(3,5)
         write(6,*) '2._dp*p4.p5',s(4,5)
         endif
         if (npart > 3) then
         write(6,*) '2._dp*p1.p6',s(1,6)
         write(6,*) '2._dp*p2.p6',s(2,6)
         write(6,*) '2._dp*p3.p6',s(3,6)
         write(6,*) '2._dp*p4.p6',s(4,6)
         write(6,*) '2._dp*p5.p6',s(5,6)
         endif
         if (npart > 4) then
         write(6,*) '2._dp*p1.p7',s(1,7)
         write(6,*) '2._dp*p2.p7',s(2,7)
         write(6,*) '2._dp*p5.p7',s(5,7)
         write(6,*) '2._dp*p6.p7',s(6,7)
         endif
         if (npart > 4) then
         write(6,*) '2._dp*p1.p8',s(1,8)
         write(6,*) '2._dp*p2.p8',s(2,8)
         write(6,*) '2._dp*p6.p8',s(6,8)
         write(6,*) '2._dp*p7.p8',s(7,8)
         endif

         write (6,*) '2._dp*p5.p8',s(5,8)

          if (npart == 3) write(6,*) 'p(5,4)',p(5,4)
          if (npart == 4) write(6,*) 'p(5,4),p(6,4)',p(5,4),p(6,4)
          if (npart == 5) write(6,*) 'p(5,4),p(6,4),p(7,4)',p(5,4),p(6,4),p(7,4)
          if (npart == 6) write(6,*) 'p(7,4),p(8,4)',p(7,4),p(8,4)

c         write(6,*) '0',jets,msq(2,0),(pjet(5,4)+pjet(6,4))**2
c     &                      -(pjet(5,1)+pjet(6,1))**2
c     &                      -(pjet(5,2)+pjet(6,2))**2
c     &                      -(pjet(5,3)+pjet(6,3))**2
         do nd=1,ndmax
            q(:,:) = ptilde(nd,:,:)
            if (incldip(nd) .eqv. .false.) then
              write (*,*) "not including ", nd
              goto 72
            endif
            call dotem(9,q,s)

            if (includedipole(nd,q)) then
               msqs(:,:) = msqs(:,:) + msqc(nd,:,:)
               write (*,*) "ND contrib", nd, includedipole(nd,q), msqc(nd,2,5)
            endif
   72    continue
         enddo

c--- find smallest value of msq
         debugsmall = minval(msq, mask=msq>0._dp)
         debuglarge = maxval(msq, mask=msq>0._dp)

         if ((debuglarge/debugsmall < xtoler)
     &  .or. (debuglarge < 1.e-10_dp)) then
c           write(6,*) ' OK - no singular configurations'
c           goto 68
         endif

         if (debugsmall > 1.e-3_dp) debugsmall=debugsmall/xtoler/2._dp
c         if (debuglarge/debugsmall < xtoler*1d2)
c     &       debugsmall=debuglarge/xtoler

         write (*,*) "       msq-msqs   msq/msqs   msq           msqs"

         do jj=-nflav,nflav
         do kk=-nflav,nflav

c         if ((gqonly) .and. (jj /= 0) .and. (kk /= 0)) goto 69
c         if ((jj > 0) .or. (kk > 0)) goto 69 ! DEBUG: QBG & GQB ONLY
c         if ((jj <= 0) .or. (kk <= 0)) goto 69 ! DEBUG: QQ only

         if (msq(jj,kk) == 0._dp) then
            if (msqs(jj,kk) == 0._dp) then
               goto 69
               debugmsg='   OK   zero msq and subtraction'
            else
               debugmsg=' FAILED subtraction with msq=0'
            endif
            write(*,22) jj,kk,'    n/a    ', '    n/a   ',msq(jj,kk),msqs(jj,kk),
     &      debugmsg
            goto 69
            endif

         if ((msq(jj,kk)/debuglarge) < 1._dp/xtoler) then
            if (abs(abs(msq(jj,kk)/msqs(jj,kk))-1._dp) < 0.02_dp)then
              debugmsg='   OK   cancelled'
            else
              debugmsg='   OK   not singular'
            endif
            write(*,21) jj,kk,msq(jj,kk)-msqs(jj,kk),msq(jj,kk)/msqs(jj,kk),
     &                    msq(jj,kk),msqs(jj,kk),debugmsg
         else
             if (msqs(jj,kk) == 0._dp) then
                debugmsg=' FAILED singular, no subtraction'
                write(*,22) jj,kk,'    n/a   ','    n/a   ',msq(jj,kk),msqs(jj,kk),
     & debugmsg
                goto 69
             endif
             if (abs(abs(msq(jj,kk)/msqs(jj,kk))-1._dp) > 0.02_dp)then
                debugmsg=' FAILED singularity uncancelled'
             else
                debugmsg='   OK   cancelled'
             endif
             write(*,21) jj,kk,msq(jj,kk)-msqs(jj,kk),msq(jj,kk)/msqs(jj,kk),
     &                  msq(jj,kk),msqs(jj,kk),debugmsg
         endif

   69    continue
         enddo
         enddo

c         pause

   68 continue
      enddo
      enddo

c---c--- This block of code is useful for checking the symmetry of
c---c--- the real matrix elements and subtractions
c---      icount=icount+1
c---      if (icount == 1) then
c---c--- save contents of arrays for reals and subtractions
c---        do j=-nf,nf
c---       do k=-nf,nf
c---       premsq(j,k)=msq(j,k)
c---       premsqs(j,k)=msqs(j,k)
c---       enddo
c---       enddo
c---      else
c---c--- compare reals and subtractions with previous iteration
c---        do j=-nf,nf
c---        do k=-nf,nf
c---c        if (j*k > 0) goto 67 ! NO QQ OR QBQB (for Z+2jet 1<>2 symm)
c---c        if (j*k <= 0) goto 67 ! ONLY QQ/QBQB (for Z+2jet 1<>2,5<>6 symm)
c---c       if (j*k < 0) goto 67 ! DEBUG: NO QQB OR QBQ
c---c       if (j*k >= 0) goto 67 ! DEBUG: QQB & QBQ only
c---c      if (j*k <= 0) goto 67 ! DEBUG: QQ & QBQB only
c---c      if ((j <= 0) .or. (k <= 0)) goto 67 ! DEBUG: QQ only
c---c       if ((j*k > 0).or.(abs(j)==5).or.(abs(k)==5)) goto 67 ! NO QQ OR QBQB, OR B QUARK
c---        if ((abs(msq(j,k)) < 1.e-15_dp) .and. (abs(msqs(j,k))<1.e-15_dp))
c---     &                    goto 67 ! skip zero contributions
c---        write(6,23) j,k,premsq(j,k),msq(k,j),premsqs(j,k),msqs(k,j),
c---     &                  premsq(j,k)/msq(k,j),premsqs(j,k)/msqs(k,j)
c---   67        continue
c---        enddo
c---        enddo
c---        icount=0
c---        pause
c---      endif

c      write(6,*) 'Singularity check completed'
c      write(6,*) 'Do ''mcfm | grep -v OK'' to look for failures'
c      stop

      return

   21 format(1x,2i3,f12.6,f10.6,2e14.6,1x,a32)
   22 format(1x,2i3,a10,a10,2e14.6,1x,a32)
c   23 format(2i4,4e12.5,2f10.6)

      end
