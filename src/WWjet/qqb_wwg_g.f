!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_wwg_g(p,msq)
c---  Author: R.K. Ellis
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  + n(p3)+ebar(p4)+e(p5)+nubar(p6)+ g(p7)+ g(p8)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'zcouple_cms.f'
      include 'ewcharge.f'
      include 'pchoice.f'

      integer:: jk,tjk,polg1,polg2,polq,minus,mplus,jtype
      real(dp):: P(mxpart,4),msq(-nf:nf,-nf:nf),
     & ave,t3456,fac,fac1,facqq,msq4q
      real(dp):: junk,
     & uc_ucsq,ubc_ubcsq,ucb_ucbsq,ubcb_ubcbsq,
     & ds_dssq,dbs_dbssq,dsb_dsbsq,dbsb_dbsbsq,
     & us_ussq,ubs_ubssq,usb_usbsq,ubsb_ubsbsq,
     & su_susq,sub_subsq,sbu_sbusq,sbub_sbubsq,
     & uu_uusq,dd_ddsq,ubub_ububsq,dbdb_dbdbsq,
     & ud_udsq,du_dusq,ubdb_ubdbsq,dbub_dbubsq,
     & udb_udbsq,dbu_dbusq,ubd_ubdsq,dub_dubsq,
     & udb_csbsq,dbu_sbcsq,ubd_cbssq,dub_scbsq,
     & uub_ddbsq,uub_ssbsq,uub_ccbsq,uub_uubsq,
     & ddb_ddbsq,ddb_ssbsq,ddb_ccbsq,ddb_uubsq,
     & ubu_dbdsq,ubu_sbssq,ubu_cbcsq,ubu_ubusq,
     & dbd_dbdsq,dbd_sbssq,dbd_cbcsq,dbd_ubusq
      complex(dp)::
     & cd_su,cub_sdb,su_cd,sdb_cub,
     & us_dc,ucb_dsb,dc_us,dsb_ucb,
     & dbs_ubc,dbcb_ubsb,ubc_dbs,ubsb_dbcb,
     & sbd_cbu,sbub_cbdb,cbu_sbd,cbdb_sbub
      complex(dp):: wpwmamp,C(2,2),
     & u_ubAB(3,2,2,2),d_dbAB(3,2,2,2),ub_uAB(3,2,2,2),db_dAB(3,2,2,2),
     & u_gAB(3,2,2,2), d_gAB(3,2,2,2), g_ubAB(3,2,2,2),g_dbAB(3,2,2,2),
     & ub_gAB(3,2,2,2),db_gAB(3,2,2,2),g_uAB(3,2,2,2),g_dAB(3,2,2,2),
     & u_ubBA(3,2,2,2),d_dbBA(3,2,2,2),ub_uBA(3,2,2,2),db_dBA(3,2,2,2),
     & u_gBA(3,2,2,2), d_gBA(3,2,2,2), g_ubBA(3,2,2,2),g_dbBA(3,2,2,2),
     & ub_gBA(3,2,2,2),db_gBA(3,2,2,2),g_uBA(3,2,2,2),g_dBA(3,2,2,2),
     & g_guAB(3,2,2,2),g_guBA(3,2,2,2),g_gdAB(3,2,2,2),g_gdBA(3,2,2,2),
     & ampAB(5),ampBA(5),propwp,propwm,prop12,cprop,
     & AB(2,2,2),BA(2,2,2),ampABx(5),ampBAx(5),ABx(2,2,2),BAx(2,2,2)
      parameter(minus=1,mplus=2)
      real(dp), parameter :: mp(nf)=(/-1d0,+1d0,-1d0,+1d0,-1d0/)

      msq(:,:)=0d0

c      fac=gw**4/4d0
      fac=esq**2/abs(zxw)**2/4d0
      fac1=32d0*gsq**2*cf*xn**2
      facqq=aveqq*(8d0*gsq*esq**2/abs(zxw)**2/4d0)**2*V/4d0

      call spinoru(8,p,za,zb)
c--   s returned from sprodx (common block) is 2*dot product

c--   calculate propagators
      t3456=s(3,4)+s(3,5)+s(3,6)+s(4,5)+s(4,6)+s(5,6)
      propwp=cmplx(s(3,4),kind=dp)
     &      /cmplx(s(3,4)-wmass**2,wmass*wwidth,kind=dp)
      propwm=cmplx(s(5,6))
     &      /cmplx(s(5,6)-wmass**2,wmass*wwidth,kind=dp)

c-----we can put in the real Z propagator and overall factor is one
      prop12=t3456/(t3456-zmass**2+im*zmass*zwidth)
      cprop=propwp*propwm

c-- couplings according to hep-ph/9803250 Eqs. 3.4 and 3.6
c-- First argument is left or right, second is d(1) or u(2)
      do j=1,2
      C(minus,j)=mp(j)*2d0*Q(j)*zxw+prop12*(1d0-2d0*mp(j)*Q(j)*zxw)
      C(mplus,j)=mp(j)*2d0*Q(j)*zxw*(1d0-prop12)
      enddo

c--- DEBUG: 4-quark only for now
c      goto 33

c---remember u-ub is the basic process.
c---case u-ubar
      call wwggamps(1,2,3,4,5,6,7,8,u_ubAB,u_ubBA,2)
c---case u-g
      call wwggamps(1,7,3,4,5,6,2,8,u_gAB,u_gBA,2)
c---case g-u
      call wwggamps(2,7,3,4,5,6,1,8,g_uAB,g_uBA,2)

c---case ubar-u
      call wwggamps(2,1,3,4,5,6,7,8,ub_uAB,ub_uBA,2)
c---case ubar-g
      call wwggamps(7,1,3,4,5,6,2,8,ub_gAB,ub_gBA,2)
c---case g-ubar
      call wwggamps(7,2,3,4,5,6,1,8,g_ubAB,g_ubBA,2)

c---case d-dbar
      call wwggamps(1,2,5,6,3,4,7,8,d_dbAB,d_dbBA,1)
c---case d-g
      call wwggamps(1,7,5,6,3,4,2,8,d_gAB,d_gBA,1)
c---case g-d
      call wwggamps(2,7,5,6,3,4,1,8,g_dAB,g_dBA,1)

c---case dbar-d
      call wwggamps(2,1,5,6,3,4,7,8,db_dAB,db_dBA,1)
c---case dbar-g
      call wwggamps(7,1,5,6,3,4,2,8,db_gAB,db_gBA,1)
c---case g-dbar
      call wwggamps(7,2,5,6,3,4,1,8,g_dbAB,g_dbBA,1)

c---case ggu
      call wwggamps(8,7,3,4,5,6,1,2,g_guAB,g_guBA,2)
c---case ggd
      call wwggamps(8,7,5,6,3,4,1,2,g_gdAB,g_gdBA,1)

c--- first populate the basic set of matrix elements
c--- (for the first generation of quarks)
      do j=-2,2
      do k=-2,2

c-- skip non-diagonal quark flavors except gluon
      if((j  /=  0 .and. k  /=  0) .and. (j  /=  -k)) goto 19

      if     ((j == 0) .and. (k == 0)) then
c--- gg
      jk=0
        ave=avegg
      elseif ((j == 0) .or. (k == 0)) then
c--- gq
        jk=j+k
        ave=aveqg
      else
c--- qq - remember final-state identical particle factor of 1/2
        jk=max(j,k)
        ave=aveqq*half
      endif

      do polg1=1,2
      do polg2=1,2
      do polq=1,2

c---sum is over diagram types (a), (b) and (c)
      do jtype=1,3
        if    (j  <  0 .and. tau(jk)  ==  -1d0 .and. k  /=  0) then
          ampAB(jtype)=db_dAB(jtype,polg1,polg2,polq)
          ampBA(jtype)=db_dBA(jtype,polg1,polg2,polq)
        elseif(j  <  0 .and. tau(jk)  ==   1d0 .and. k  /=  0) then
          ampAB(jtype)=ub_uAB(jtype,polg1,polg2,polq)
          ampBA(jtype)=ub_uBA(jtype,polg1,polg2,polq)
        elseif(j  >  0 .and. tau(jk)  ==  -1d0 .and. k  /=  0) then
          ampAB(jtype)=d_dbAB(jtype,polg1,polg2,polq)
          ampBA(jtype)=d_dbBA(jtype,polg1,polg2,polq)
        elseif(j  >  0 .and. tau(jk)  ==   1d0 .and. k  /=  0) then
          ampAB(jtype)=u_ubAB(jtype,polg1,polg2,polq)
          ampBA(jtype)=u_ubBA(jtype,polg1,polg2,polq)
        elseif(j  ==  0 .and. tau(jk)  ==   1d0 .and. jk  >  0) then
          ampAB(jtype)=g_uAB(jtype,polg1,polg2,polq)
          ampBA(jtype)=g_uBA(jtype,polg1,polg2,polq)
        elseif(j  ==  0 .and. tau(jk)  ==  -1d0 .and. jk  >  0) then
          ampAB(jtype)=g_dAB(jtype,polg1,polg2,polq)
          ampBA(jtype)=g_dBA(jtype,polg1,polg2,polq)
        elseif(j  ==  0 .and. tau(jk)  ==  -1d0 .and. jk  <  0) then
          ampAB(jtype)=g_ubAB(jtype,polg1,polg2,polq)
          ampBA(jtype)=g_ubBA(jtype,polg1,polg2,polq)
        elseif(j  ==  0 .and. tau(jk)  ==   1d0 .and. jk  <  0) then
          ampAB(jtype)=g_dbAB(jtype,polg1,polg2,polq)
          ampBA(jtype)=g_dbBA(jtype,polg1,polg2,polq)
        elseif(k  ==  0 .and. tau(jk)  ==   1d0 .and. jk  >  0) then
          ampAB(jtype)=u_gAB(jtype,polg1,polg2,polq)
          ampBA(jtype)=u_gBA(jtype,polg1,polg2,polq)
        elseif(k  ==  0 .and. tau(jk)  ==  -1d0 .and. jk  >  0) then
          ampAB(jtype)=d_gAB(jtype,polg1,polg2,polq)
          ampBA(jtype)=d_gBA(jtype,polg1,polg2,polq)
        elseif(k  ==  0 .and. tau(jk)  ==  -1d0 .and. jk  <  0) then
          ampAB(jtype)=ub_gAB(jtype,polg1,polg2,polq)
          ampBA(jtype)=ub_gBA(jtype,polg1,polg2,polq)
        elseif(k  ==  0 .and. tau(jk)  ==   1d0 .and. jk  <  0) then
          ampAB(jtype)=db_gAB(jtype,polg1,polg2,polq)
          ampBA(jtype)=db_gBA(jtype,polg1,polg2,polq)
        elseif ((j  ==  0).and.(k  ==  0)) then
          ampAB(jtype)=g_guAB(jtype,polg1,polg2,polq)
          ampBA(jtype)=g_guBA(jtype,polg1,polg2,polq)
          ampABx(jtype)=g_gdAB(jtype,polg1,polg2,polq)
          ampBAx(jtype)=g_gdBA(jtype,polg1,polg2,polq)
        endif
      enddo

c---tjk is equal to 2 (u,c) or 1 (d,s,b)
      tjk=2-mod(abs(jk),2)

c-- include coupling in l1 to account for non-leptonic W decays
      AB(polg1,polg2,polq)=l1*cmplx(fac,kind=dp)*(
     & cprop*(ampAB(1)+C(polq,tjk)*ampAB(2))+ampAB(3)*zxw)

      BA(polg1,polg2,polq)=l1*cmplx(fac,kind=dp)*(
     & cprop*(ampBA(1)+C(polq,tjk)*ampBA(2))+ampBA(3)*zxw)

c--- extra amplitudes for gg case (for gg->ddb)
      if ((j == 0) .and. (k == 0)) then
      ABx(polg1,polg2,polq)=l1*cmplx(fac,kind=dp)*(
     & cprop*(ampABx(1)+C(polq,1)*ampABx(2))+ampABx(3)*zxw)

      BAx(polg1,polg2,polq)=l1*cmplx(fac,kind=dp)*(
     & cprop*(ampBAx(1)+C(polq,1)*ampBAx(2))+ampBAx(3)*zxw)
      endif

      enddo
      enddo
      enddo

      msq(j,k)=fac1*ave*(
     & +cdabs(AB(1,1,1))**2+cdabs(AB(1,1,2))**2
     & +cdabs(AB(1,2,1))**2+cdabs(AB(1,2,2))**2
     & +cdabs(AB(2,2,2))**2+cdabs(AB(2,2,1))**2
     & +cdabs(AB(2,1,2))**2+cdabs(AB(2,1,1))**2

     & +cdabs(BA(1,1,1))**2+cdabs(BA(1,1,2))**2
     & +cdabs(BA(1,2,1))**2+cdabs(BA(1,2,2))**2
     & +cdabs(BA(2,2,2))**2+cdabs(BA(2,2,1))**2
     & +cdabs(BA(2,1,2))**2+cdabs(BA(2,1,1))**2
     & -(
     & +cdabs(AB(1,1,1)+BA(1,1,1))**2
     & +cdabs(AB(1,1,2)+BA(1,1,2))**2
     & +cdabs(AB(1,2,1)+BA(1,2,1))**2
     & +cdabs(AB(1,2,2)+BA(1,2,2))**2
     & +cdabs(AB(2,2,2)+BA(2,2,2))**2
     & +cdabs(AB(2,2,1)+BA(2,2,1))**2
     & +cdabs(AB(2,1,2)+BA(2,1,2))**2
     & +cdabs(AB(2,1,1)+BA(2,1,1))**2)/xn**2)

c--- extra contribution for gg->ddb
      if ((j == 0) .and. (k == 0)) then
      msq(j,k)=msq(j,k)+fac1*ave*(
     & +cdabs(ABx(1,1,1))**2+cdabs(ABx(1,1,2))**2
     & +cdabs(ABx(1,2,1))**2+cdabs(ABx(1,2,2))**2
     & +cdabs(ABx(2,2,2))**2+cdabs(ABx(2,2,1))**2
     & +cdabs(ABx(2,1,2))**2+cdabs(ABx(2,1,1))**2

     & +cdabs(BAx(1,1,1))**2+cdabs(BAx(1,1,2))**2
     & +cdabs(BAx(1,2,1))**2+cdabs(BAx(1,2,2))**2
     & +cdabs(BAx(2,2,2))**2+cdabs(BAx(2,2,1))**2
     & +cdabs(BAx(2,1,2))**2+cdabs(BAx(2,1,1))**2
     & -(
     & +cdabs(ABx(1,1,1)+BAx(1,1,1))**2
     & +cdabs(ABx(1,1,2)+BAx(1,1,2))**2
     & +cdabs(ABx(1,2,1)+BAx(1,2,1))**2
     & +cdabs(ABx(1,2,2)+BAx(1,2,2))**2
     & +cdabs(ABx(2,2,2)+BAx(2,2,2))**2
     & +cdabs(ABx(2,2,1)+BAx(2,2,1))**2
     & +cdabs(ABx(2,1,2)+BAx(2,1,2))**2
     & +cdabs(ABx(2,1,1)+BAx(2,1,1))**2)/xn**2)
c--- multiply the whole answer by 2, to account for generations
      msq(j,k)=msq(j,k)*2d0
      endif

c      if (j == 2 .and. k == -2) then
c      write(6,*) 'uub 1 1 1 :',fac1*ave*(
c     & +cdabs(AB(1,1,1))**2+cdabs(BA(1,1,1))**2
c     & -cdabs(AB(1,1,1)+BA(1,1,1))**2/xn**2)
c      write(6,*) fac1*ave*cdabs(AB(1,1,1))**2,
c     &           fac1*ave*cdabs(BA(1,1,1))**2
c      write(6,*) 'uub 1 2 1 :',fac1*ave*(
c     & +cdabs(AB(1,2,1))**2+cdabs(BA(1,2,1))**2
c     & -cdabs(AB(1,2,1)+BA(1,2,1))**2/xn**2)
c      write(6,*) fac1*ave*cdabs(AB(1,2,1))**2,
c     &           fac1*ave*cdabs(BA(1,2,1))**2
c      write(6,*) 'uub 2 1 1 :',fac1*ave*(
c     & +cdabs(AB(2,1,1))**2+cdabs(BA(2,1,1))**2
c     & -cdabs(AB(2,1,1)+BA(2,1,1))**2/xn**2)
c      write(6,*) fac1*ave*cdabs(AB(2,1,1))**2,
c     &           fac1*ave*cdabs(BA(2,1,1))**2
c      write(6,*) 'uub 2 2 1 :',fac1*ave*(
c     & +cdabs(AB(2,2,1))**2+cdabs(BA(2,2,1))**2
c     & -cdabs(AB(2,2,1)+BA(2,2,1))**2/xn**2)
c      write(6,*) fac1*ave*cdabs(AB(2,2,1))**2,
c     &           fac1*ave*cdabs(BA(2,2,1))**2
c      write(6,*) 'uub 1 1 2 :',fac1*ave*(
c     & +cdabs(AB(1,1,2))**2+cdabs(BA(1,1,2))**2
c     & -cdabs(AB(1,1,2)+BA(1,1,2))**2/xn**2)
c      write(6,*) fac1*ave*cdabs(AB(1,1,2))**2,
c     &           fac1*ave*cdabs(BA(1,1,2))**2
c      write(6,*) 'uub 1 2 2 :',fac1*ave*(
c     & +cdabs(AB(1,2,2))**2+cdabs(BA(1,2,2))**2
c     & -cdabs(AB(1,2,2)+BA(1,2,2))**2/xn**2)
c      write(6,*) fac1*ave*cdabs(AB(1,2,2))**2,
c     &           fac1*ave*cdabs(BA(1,2,2))**2
c      write(6,*) 'uub 2 1 2 :',fac1*ave*(
c     & +cdabs(AB(2,1,2))**2+cdabs(BA(2,1,2))**2
c     & -cdabs(AB(2,1,2)+BA(2,1,2))**2/xn**2)
c      write(6,*) fac1*ave*cdabs(AB(2,1,2))**2,
c     &           fac1*ave*cdabs(BA(2,1,2))**2
c      write(6,*) 'uub 2 2 2 :',fac1*ave*(
c     & +cdabs(AB(2,2,2))**2+cdabs(BA(2,2,2))**2
c     & -cdabs(AB(2,2,2)+BA(2,2,2))**2/xn**2)
c      write(6,*) fac1*ave*cdabs(AB(2,2,2))**2,
c     &           fac1*ave*cdabs(BA(2,2,2))**2
c      endif

   19 continue
      enddo
      enddo

c--- fill out the 2nd-generation contributions
      do j=3,4
      msq(j,-j)=msq(j-2,-(j-2))
      msq(-j,j)=msq(-(j-2),j-2)
      msq(j,0)=msq(j-2,0)
      msq(0,j)=msq(0,j-2)
      msq(-j,0)=msq(-(j-2),0)
      msq(0,-j)=msq(0,-(j-2))
      enddo

c      write(6,*) 'DEBUG: only 2-quark contribution for now'
c--- DEBUG: 2-quark contribution only for now
c      return

c 33   continue

c---*****************4 quark diagrams********************
c---Charge in initial state must sum to that in the final state
c-- Same Weak Isospin

c     1) uc --> uc  (upper line, lower line) (14 diagrams)
c     2) uu --> uu  (upper line lower line + interference)

c---Different weak isospin

c     3) ud --> ud  Upper line, lower line + one W on each line
c     4) us -->us   Upper line, lower line (14 diagrams)
c     5) us -->cd   One W on each line (4 diagrams)

c--  start with four quark processes with both WW's attached to
c--- the same quark line

c--- Processes 1 and 4
      call a8fourqsq(1,2,3,4,5,6,7,8,C,facqq,
     & ds_dssq,su_susq,us_ussq,uc_ucsq)
      call a8fourqsq(1,8,3,4,5,6,7,2,C,facqq,
     & dsb_dsbsq,sub_subsq,usb_usbsq,ucb_ucbsq)
      call a8fourqsq(7,8,3,4,5,6,1,2,C,facqq,
     & dbsb_dbsbsq,sbub_sbubsq,ubsb_ubsbsq,ubcb_ubcbsq)
      call a8fourqsq(7,2,3,4,5,6,1,8,C,facqq,
     & dbs_dbssq,sbu_sbusq,ubs_ubssq,ubc_ubcsq)

c---Processs 2
      call a8qid(1,2,3,4,5,6,7,8,C,2,2,facqq,uu_uusq)
      call a8qid(1,2,5,6,3,4,7,8,C,1,1,facqq,dd_ddsq)
      call a8qid(7,8,3,4,5,6,1,2,C,2,2,facqq,ubub_ububsq)
      call a8qid(7,8,5,6,3,4,1,2,C,1,1,facqq,dbdb_dbdbsq)

c*************************************************************
c  We now have to add the diagrams with four quarks
c  with the W's on different lines
c --
c --
c --                 (3)---<-------(4)
c --                       $  W+
c --                        $
c --        (7)----------<---$----<----u(1)
c --                   o
c --                   o
c --                   o
c--         (8)----------<--------<----d(2)
c--                            $
c--                           $ W-
c--                          $
c --                 (5)---<-------(6)
c--
c----Amplitude for
c     u(-p1)+c(-p2)-->W^+ + W^- d(p7)+s(p8)
c                     |     |
c                     |     --> e^-(p5)+nu~(p6)
c                     |
c                     |---> num(p3)+mu^+(p4)

c    with a facqq of 8*gsq*gw**4/4 and color facqq V/4 removed

c--- Process 5
      us_dc  =  wpwmamp(1,2,3,4,5,6,7,8)
      ucb_dsb=  wpwmamp(1,8,3,4,5,6,7,2)

      dc_us  =  wpwmamp(1,2,5,6,3,4,7,8)
      dsb_ucb=  wpwmamp(1,8,5,6,3,4,7,2)

      ubc_dbs=wpwmamp(7,2,5,6,3,4,1,8)
      ubsb_dbcb=wpwmamp(7,8,5,6,3,4,1,2)

      dbs_ubc=  wpwmamp(7,2,3,4,5,6,1,8)
      dbcb_ubsb=wpwmamp(7,8,3,4,5,6,1,2)

      cd_su=us_dc
      cub_sdb=ucb_dsb

      su_cd=dc_us
      sdb_cub=dsb_ucb

      sbd_cbu=dbs_ubc
      sbub_cbdb=dbcb_ubsb

      cbu_sbd=ubc_dbs
      cbdb_sbub=ubsb_dbcb

c--- contributions that have diagrams from both sets (process 3)
      call a8ududsq(1,2,3,4,5,6,7,8,C,2,1,facqq,ud_udsq,junk)
      call a8ududsq(7,8,3,4,5,6,1,2,C,2,1,facqq,ubdb_ubdbsq,junk)
      call a8ududsq(2,1,3,4,5,6,8,7,C,2,1,facqq,du_dusq,junk)
      call a8ududsq(8,7,3,4,5,6,2,1,C,2,1,facqq,dbub_dbubsq,junk)

      call a8ududsq(1,8,3,4,5,6,7,2,C,2,1,facqq,udb_udbsq,udb_csbsq)
      call a8ududsq(2,7,3,4,5,6,8,1,C,2,1,facqq,dbu_dbusq,dbu_sbcsq)
      call a8ududsq(8,1,3,4,5,6,2,7,C,2,1,facqq,dub_dubsq,dub_scbsq)
      call a8ududsq(7,2,3,4,5,6,1,8,C,2,1,facqq,ubd_ubdsq,ubd_cbssq)

c---- quark-antiquark initial states (combination of processes 1,2,3,4)
c--- annihilation diagrams
      call a8fourqsq(1,8,3,4,5,6,2,7,C,facqq,
     & ddb_ssbsq,ddb_ccbsq,uub_ssbsq,uub_ccbsq)
      call a8fourqsq(2,7,3,4,5,6,1,8,C,facqq,
     & dbd_sbssq,dbd_cbcsq,ubu_sbssq,ubu_cbcsq)

c--- annihilation + one W off each line
      call a8ududsq(1,8,3,4,5,6,2,7,C,2,1,facqq,uub_ddbsq,junk)
      call a8ududsq(8,1,3,4,5,6,7,2,C,2,1,facqq,ddb_uubsq,junk)
      call a8ududsq(2,7,3,4,5,6,1,8,C,2,1,facqq,ubu_dbdsq,junk)
      call a8ududsq(7,2,3,4,5,6,8,1,C,2,1,facqq,dbd_ubusq,junk)

c--- annihilation + two W's off each line
      call a8qid(1,8,3,4,5,6,7,2,C,2,2,facqq,uub_uubsq)
      call a8qid(1,8,5,6,3,4,7,2,C,1,1,facqq,ddb_ddbsq)
      call a8qid(7,2,3,4,5,6,1,8,C,2,2,facqq,ubu_ubusq)
      call a8qid(7,2,5,6,3,4,1,8,C,1,1,facqq,dbd_dbdsq)

c---- QUARKS IN DIFFERENT GENERATIONS
c--- Quarks with charges of opposite sign have 2 contributions:
c---  (1) a W emitted from each line
c---  (2) both W's emitted on the same line

c--- Quarks with the same-sign charge only have contribution (2)

c---- QUARKS IN THE SAME GENERATION
c--- Identical quarks only have contribution (2)

c--- Opposite-sign, different flavour (ud) have one contribution,
c---- the sum of diagrams (1) and (2) above

c--- Same-sign, different flavours (udb) have a similar contribution
c---  (udb->udb) plus a contribution from diagrams (1) alone (udb->csb)

c--- Opposite-sign, same flavours (uub) have a contribution with
c---  W's emitted from the same line and the gluon in s- and t-
c---  channels (uub->uub), plus contributions from s-channel gluon
c---  diagrams only (uub->ssb, uub->ssb), plus a contribution from
c---  s-channel gluons and diagrams with a W on each line (uub->ddb)

c--- debug : to check 4-quarks only
c      do j=-nf,nf
c      do k=-nf,nf
c      msq(j,k)=0d0
c      enddo
c      enddo
c--- debug

c--- fill 4-quark matrix elements
       do j=-4,4
       do k=-4,4

       msq4q=0d0

       if         ((j  == 4) .and. (k  ==  1))  then
         msq4q=facqq*cdabs(cd_su)**2+us_ussq
       elseif     ((j  == 4) .and. (k  ==  2))  then
         msq4q=uc_ucsq
       elseif     ((j  == 4) .and. (k  ==  3))  then
         msq4q=ud_udsq
       elseif     ((j  == 4) .and. (k  ==  4))  then
         msq4q=uu_uusq*half
       elseif     ((j  == 4) .and. (k  ==  -1))  then
         msq4q=usb_usbsq
       elseif     ((j  == 4) .and. (k  ==  -2))  then
         msq4q=facqq*cdabs(cub_sdb)**2+ucb_ucbsq
       elseif     ((j  == 4) .and. (k  ==  -3))  then
         msq4q=udb_udbsq+udb_csbsq
       elseif     ((j  == 4) .and. (k  ==  -4))  then
         msq4q=uub_ssbsq+uub_ccbsq+uub_ddbsq+uub_uubsq

       elseif     ((j  == 3) .and. (k  ==  1))  then
         msq4q=ds_dssq
       elseif     ((j  == 3) .and. (k  ==  2))  then
         msq4q=facqq*cdabs(su_cd)**2+su_susq
       elseif     ((j  == 3) .and. (k  ==  3))  then
         msq4q=dd_ddsq*half
       elseif     ((j  == 3) .and. (k  ==  4))  then
         msq4q=du_dusq
       elseif     ((j  == 3) .and. (k  ==  -1))  then
         msq4q=facqq*cdabs(sdb_cub)**2+dsb_dsbsq
       elseif     ((j  == 3) .and. (k  ==  -2))  then
         msq4q=sub_subsq
       elseif     ((j  == 3) .and. (k  ==  -3))  then
         msq4q=ddb_ssbsq+ddb_ccbsq+ddb_ddbsq+ddb_uubsq
       elseif     ((j  == 3) .and. (k  ==  -4))  then
         msq4q=dub_dubsq+dub_scbsq

       elseif     ((j  == 2) .and. (k  ==  1))  then
         msq4q=ud_udsq
       elseif     ((j  == 2) .and. (k  ==  2))  then
         msq4q=uu_uusq*half
       elseif     ((j  == 2) .and. (k  ==  3))  then
         msq4q=facqq*cdabs(us_dc)**2+us_ussq
       elseif     ((j  == 2) .and. (k  ==  4))  then
         msq4q=uc_ucsq
       elseif     ((j  == 2) .and. (k  ==  -1))  then
         msq4q=udb_udbsq+udb_csbsq
       elseif     ((j  == 2) .and. (k  ==  -2))  then
         msq4q=uub_ssbsq+uub_ccbsq+uub_ddbsq+uub_uubsq
       elseif     ((j  == 2) .and. (k  ==  -3))  then
         msq4q=usb_usbsq
       elseif     ((j  == 2) .and. (k  ==  -4))  then
         msq4q=facqq*cdabs(ucb_dsb)**2+ucb_ucbsq

       elseif     ((j  == 1) .and. (k  ==  1))  then
         msq4q=dd_ddsq*half
       elseif     ((j  == 1) .and. (k  ==  2))  then
         msq4q=du_dusq
       elseif     ((j  == 1) .and. (k  ==  3))  then
         msq4q=ds_dssq
       elseif     ((j  == 1) .and. (k  ==  4))  then
         msq4q=facqq*cdabs(dc_us)**2+su_susq
       elseif     ((j  == 1) .and. (k  ==  -1))  then
         msq4q=ddb_ssbsq+ddb_ccbsq+ddb_ddbsq+ddb_uubsq
       elseif     ((j  == 1) .and. (k  ==  -2))  then
         msq4q=dub_dubsq+dub_scbsq
       elseif     ((j  == 1) .and. (k  ==  -3))  then
         msq4q=facqq*cdabs(dsb_ucb)**2+dsb_dsbsq
       elseif     ((j  == 1) .and. (k  ==  -4))  then
         msq4q=sub_subsq

       elseif     ((j  == -1) .and. (k  ==  1))  then
         msq4q=dbd_sbssq+dbd_cbcsq+dbd_dbdsq+dbd_ubusq
       elseif     ((j  == -1) .and. (k  ==  2))  then
         msq4q=dbu_dbusq+dbu_sbcsq
       elseif     ((j  == -1) .and. (k  ==  3))  then
         msq4q=facqq*cdabs(dbs_ubc)**2+dbs_dbssq
       elseif     ((j  == -1) .and. (k  ==  4))  then
         msq4q=sbu_sbusq
       elseif     ((j  == -1) .and. (k  ==  -1))  then
         msq4q=dbdb_dbdbsq*half
       elseif     ((j  == -1) .and. (k  ==  -2))  then
         msq4q=dbub_dbubsq
       elseif     ((j  == -1) .and. (k  ==  -3))  then
         msq4q=dbsb_dbsbsq
       elseif     ((j  == -1) .and. (k  ==  -4))  then
         msq4q=facqq*cdabs(dbcb_ubsb)**2+sbub_sbubsq

       elseif     ((j  == -2) .and. (k  ==  1))  then
         msq4q=ubd_ubdsq+ubd_cbssq
       elseif     ((j  == -2) .and. (k  ==  2))  then
         msq4q=ubu_sbssq+ubu_cbcsq+ubu_dbdsq+ubu_ubusq
       elseif     ((j  == -2) .and. (k  ==  3))  then
         msq4q=ubs_ubssq
       elseif     ((j  == -2) .and. (k  ==  4))  then
         msq4q=facqq*cdabs(ubc_dbs)**2+ubc_ubcsq
       elseif     ((j  == -2) .and. (k  ==  -1))  then
         msq4q=ubdb_ubdbsq
       elseif     ((j  == -2) .and. (k  ==  -2))  then
         msq4q=ubub_ububsq*half
       elseif     ((j  == -2) .and. (k  ==  -3))  then
         msq4q=facqq*cdabs(ubsb_dbcb)**2+ubsb_ubsbsq
       elseif     ((j  == -2) .and. (k  ==  -4))  then
         msq4q=ubcb_ubcbsq

       elseif     ((j  == -3) .and. (k  ==  1))  then
         msq4q=facqq*cdabs(sbd_cbu)**2+dbs_dbssq
       elseif     ((j  == -3) .and. (k  ==  2))  then
         msq4q=sbu_sbusq
       elseif     ((j  == -3) .and. (k  ==  3))  then
         msq4q=dbd_sbssq+dbd_cbcsq+dbd_dbdsq+dbd_ubusq
       elseif     ((j  == -3) .and. (k  ==  4))  then
         msq4q=dbu_dbusq+dbu_sbcsq
       elseif     ((j  == -3) .and. (k  ==  -1))  then
         msq4q=dbsb_dbsbsq
       elseif     ((j  == -3) .and. (k  ==  -2))  then
         msq4q=facqq*cdabs(sbub_cbdb)**2+sbub_sbubsq
       elseif     ((j  == -3) .and. (k  ==  -3))  then
         msq4q=dbdb_dbdbsq*half
       elseif     ((j  == -3) .and. (k  ==  -4))  then
         msq4q=dbub_dbubsq

       elseif     ((j  == -4) .and. (k  ==  -1))  then
         msq4q=facqq*cdabs(cbdb_sbub)**2+ubsb_ubsbsq
       elseif     ((j  == -4) .and. (k  ==  -2))  then
         msq4q=ubcb_ubcbsq
       elseif     ((j  == -4) .and. (k  ==  -3))  then
         msq4q=ubdb_ubdbsq
       elseif     ((j  == -4) .and. (k  ==  -4))  then
         msq4q=ubub_ububsq*half
       elseif     ((j  == -4) .and. (k  ==  +1))  then
         msq4q=ubs_ubssq
       elseif     ((j  == -4) .and. (k  ==  +2))  then
         msq4q=facqq*cdabs(cbu_sbd)**2+ubc_ubcsq
       elseif     ((j  == -4) .and. (k  ==  +3))  then
         msq4q=ubd_ubdsq+ubd_cbssq
       elseif     ((j  == -4) .and. (k  ==  +4))  then
         msq4q=ubu_sbssq+ubu_cbcsq+ubu_dbdsq+ubu_ubusq

       endif

c--- overall factor from Z propagator if we're not generating on-shell
c       if (zerowidth .eqv. .false.) then
c         msq4q=msq4q*cdabs(propzg)**2
c       endif

       msq(j,k)=msq(j,k)+msq4q

       enddo
       enddo

      return
      end

