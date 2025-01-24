!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module Multichannel
    implicit none

    public :: multichan,dipoleconfig

    private

    contains

    function multichan(x1_in,x2_in,x3,x4,pin,pout,wtdip, dipconfig_in)
    use types
    implicit none
!***********************************************************************
!                                                                      *
!     Given a n-parton phase space point, generate an (n+1)-parton     *
!     configuration using a multi-channel procedure based on the       *
!     dipole structure of the process at hand                          *
!                                                                      *
!     x1,x2,x3: random numbers uniform on [0,1]                        *
!     x4: random number uniform on [0,1], used for choosing channel    *
!                                                                      *
!     Author: J.M.Campbell                                             *
!       Date: 19th March 2009                                          *
!                                                                      *
!     Based on routine originally used in:                             *
!                                                                      *
!     "Generalized unitarity at work: first NLO QCD results            *
!      for hadronic W+3 jet production"                                *
!      R. K. Ellis, K. Melnikov and G. Zanderighi                      *
!      arXiv:0901.4101 (2009)                                          *
!                                                                      *
!***********************************************************************
          
    include 'constants.f'
    include 'nf.f'
    include 'mxpart.f'
    include 'cplx.h'
    include 'npart.f'
    include 'maxd.f'
    include 'debug.f'

    logical :: multichan
    logical :: swapjpkp
    real(dp), intent(in) :: x1_in,x2_in,x3,x4, pin(mxpart,4)
    real(dp), intent(out) :: pout(mxpart,4), wtdip
    real(dp), intent(in), optional :: dipconfig_in(:,:)

    integer::maxdip,dipconfig(maxd,3),ichan,ip,jp,kp,nu,m,iskip,ipp, &
    kpp,mp
    real(dp) :: tiny,ranchan,vtilde,x,u,z,y,phi,sik,pt,v1(4),v2(4),vperp(4), &
    qi(4),qj(4),qk(4),dot,k(4),kt(4),ks(4),kDk,ksDks,ktDp,ksDp, &
    xjac,sij,sjk,sim,sjm,qm(4),ptsq,ptmp(4)
    real(dp) :: x1,x2
    real(dp) :: dotpr
    integer :: diptype
    integer:: maxperms,perm2(2),perm3(2)
    integer, parameter:: inin=1, infi=2, fiin=3, fifi=4, phot=5
    parameter(tiny=1.e-14_dp)

    x1 = x1_in
    x2 = x2_in
         
    if (dosingcheck) then
! For approaching collinear region
      x1=1.e-6_dp*x1+1.e-5_dp*(one-x1)
      x2=0.5_dp-(1.e-4_dp*x2+1.e-3_dp*(one-x2))
!      write(6,*) 'WARNING: multichan set to check collinear limits'
! For approaching soft region for initial-initial dipoles
!       x1=1.e-6_dp*x1
!       x2=1.e-6_dp*x2

    endif
!      write(6,*) 'WARNING: multichan set to check soft limits'

!--- set dipole weight and momenta to zero, in case of early return
    wtdip=zip
    pout(:,:)=zip

    if (present(dipconfig_in)) then
        perm2 = 1
        perm3 = 1
        maxperms = 1
        maxdip = size(dipconfig_in, 1)
        dipconfig(1:maxdip,1:3) = dipconfig_in(1:maxdip,1:3)
    else
!---    Step 1: obtain number of dipoles and their configurations
        call dipoleconfig(maxdip,dipconfig,maxperms,perm2,perm3)
    endif

!--- catch processes for which dipole configurations are not specified
    if (maxperms < 0) then
      write(6,*) 'Dipole configurations for this process not properly specified in dipoleconfig.f'
      stop
    endif

                
!--- Step 2: choose dipole channel
!---  (note: this is uniform and not chosen adaptively)
    ranchan=x4*real(maxdip*maxperms,dp)-tiny
    ichan=int(ranchan)+1

! Relevant channels to check singletop+jet collinear limits
!      ichan=4      ! s78
!      ichan=3      ! s67
!      ichan=15     ! s68
!      ichan=1      ! s17
!      ichan=13     ! s18
!      ichan=11     ! s16
!      ichan=7      ! s27
!      ichan=19     ! s28
!      ichan=12     ! s26

! Relevant channels to check singletop+jet soft limits
!      ichan=3 ! DEBUG -- p7 soft
!      ichan=15 ! DEBUG -- p8 soft

    if ((ichan < 1) .OR. (ichan > maxdip*maxperms)) then
        write(6,*) 'Chosen channel is out of bounds: ',ichan
        write(6,*) ' (should be between 1 and',maxdip*maxperms,')'
        write(6,*) 'setting ichan to ',maxdip*maxperms
        ichan = maxdip*maxperms
    endif
               
!--- Step 3: determine type of dipole (ii/if/fi/ff)
    if (ichan <= maxdip) then
        ip=dipconfig(ichan,1)
        jp=dipconfig(ichan,2)
        kp=dipconfig(ichan,3)
    elseif ((ichan > maxdip) .AND. (ichan <= 2*maxdip)) then
        ip=dipconfig(ichan-maxdip,1)
        jp=dipconfig(ichan-maxdip,2)
        kp=dipconfig(ichan-maxdip,3)
    elseif ((ichan > 2*maxdip) .AND. (ichan <= 3*maxdip)) then
        ip=dipconfig(ichan-2*maxdip,1)
        jp=dipconfig(ichan-2*maxdip,2)
        kp=dipconfig(ichan-2*maxdip,3)
    endif
!      write(6,*) 'Chosen dipole: ',ip,jp,kp
!      write(6,*) 'Chosen channel ',ichan,': ',ip,jp,kp
          
    swapjpkp=.false.
! easier to assume jp > kp (can only apply in final state)
    if (kp > jp) then
      kpp=kp
      kp=jp
      jp=kpp
      swapjpkp = .true.
    endif

    if     ((ip <= 2) .AND. (kp <= 2)) then
        diptype=inin
    elseif ((ip <= 2) .AND. (kp > 2)) then
        diptype=infi
    elseif ((ip > 2) .AND. (kp <= 2)) then
        diptype=fiin
    elseif ((ip > 2) .AND. (kp > 2)) then
        diptype=fifi
    endif

!--- special notation for identified photon dipole
    if (kp == 0) then
        diptype=phot
    endif

!--- Step 4: generate dipole phase-space according to type

!--- NOTE: the section below still implicitly assumes that the
!---       emitted parton is in the last position, since the
!---       other momenta are assumed to exist in pin.
!---       Need to modify so that this is not the case.
!---       IN FACT OKAY FOR inin AND FIXED NOW FOR infi

    if    (diptype == inin) then
    !---    initial-initial dipoles
        vtilde=half*(one-(one-two*x1)*sqrt(one+four*x1*(one-x1)))
        x=one-x2**2
        phi=two*pi*x3
        if (x == 0._dp) then
            multichan = .false.
            return
        endif
                
        if (vtilde > one-x) then
            multichan = .false.
            return    ! c.f. (5.151) of C.-S.
        endif
                
        sik=two*dot(pin,ip,kp)
        if ((pin(ip,4) == 0._dp) .or. (pin(kp,4) == 0._dp)) then
            multichan = .false.
            return
        endif
    !--- catch generation of exceptionally collinear points
        if (abs(sik/(two*pin(ip,4)*pin(kp,4))) < tiny) then
            multichan = .false.
            return
        endif
        ptsq=sik*vtilde*(one-vtilde-x)/(x+tiny)
        if (ptsq < 0) then
            multichan = .false.
            return
        endif
        pt=sqrt(ptsq)
                
        if (.not. getperp(pin,ip,kp,v1,v2)) then
            multichan = .false.
            return
        endif
                
        do nu=1,4
            vperp(nu)=cos(phi)*v1(nu)+sin(phi)*v2(nu)
                    
            qi(nu)=pin(ip,nu)/(x+tiny)
            qj(nu)=-(one-vtilde-x)/(x+tiny)*pin(ip,nu) &
            -vtilde*pin(kp,nu) &
            +pt*vperp(nu)
            qk(nu)=pin(kp,nu)
        enddo
                
        wtdip=8._dp*sik/x**2/16._dp/pi**2
                
    !--- note: we may end up with initial state momentum qi with
    !---       energy larger than that of the beam; that should be
    !---       rejected outside of this routine
                
    elseif (diptype == infi) then
    !---    initial-final dipoles
        u=half*(one-(one-two*x1)*sqrt(one+four*x1*(one-x1)))
        x=one-x2**2
        phi=two*pi*x3
        if (x == 0._dp) then
            multichan = .false.
            return
        endif
                
    !--- account for dipoles with emitted partons not in the last position
        kpp=kp
        if (kp > jp) kpp=kp-1
                
        sik=two*dot(pin,ip,kpp)
        if ((pin(ip,4) == 0._dp) .or. (pin(kp,4) == 0._dp)) then
            multichan = .false.
            return
        endif
    !--- catch generation of exceptionally collinear points
        if (abs(sik/pin(ip,4)/pin(kpp,4)) < tiny) then
            multichan = .false.
            return
        endif
        ptsq=-sik*u*(one-u)*(one-x)/(x+tiny)
        if (ptsq < 0) then
            multichan = .false.
            return
        endif
        pt=sqrt(ptsq)
                
        if (.not. getperp(pin,ip,kpp,v1,v2)) then
            multichan = .false.
            return
        endif
                
        do nu=1,4
            vperp(nu)=cos(phi)*v1(nu)+sin(phi)*v2(nu)

            qi(nu)=pin(ip,nu)/(x+tiny)
            qj(nu)=-(one-u)*(one-x)/(x+tiny)*pin(ip,nu) &
            +u*pin(kpp,nu)+pt*vperp(nu)
            qk(nu)=-u*(one-x)/(x+tiny)*pin(ip,nu) &
            +(one-u)*pin(kpp,nu)-pt*vperp(nu)
        enddo
                
        wtdip=-8._dp*sik/x**2/16._dp/pi**2

    !--- note: we may end up with initial state momentum qi with
    !---       energy larger than that of the beam; that should be
    !---       rejected outside of this routine
                
    elseif (diptype == fiin) then
    !---    final-initial dipoles
        z=half*(one-(one-two*x2)*sqrt(one+four*x2*(one-x2)))
        x=one-x1**2
        phi=two*pi*x3
        if (x == 0._dp) then
            multichan = .false.
            return
        endif
                
        sik=two*dot(pin,ip,kp)
        if ((pin(ip,4) == 0._dp) .or. (pin(kp,4) == 0._dp)) then
            multichan = .false.
            return
        endif
    !--- catch generation of exceptionally collinear points
        if (abs(sik/(two*pin(ip,4)*pin(kp,4))) < tiny) then
            multichan = .false.
            return
        endif
        ptsq=-sik*z*(one-z)*(one-x)/(x+tiny)
        if (ptsq < 0) then
            multichan = .false.
            return
        endif
        pt=sqrt(ptsq)
                
        if (.not. getperp(pin,ip,kp,v1,v2)) then
            multichan = .false.
            return
        endif
                
        do nu=1,4
            vperp(nu)=cos(phi)*v1(nu)+sin(phi)*v2(nu)
                    
            qk(nu)=pin(kp,nu)/(x+tiny)
            qi(nu)=-(one-z)*(one-x)/(x+tiny)*pin(kp,nu) &
            +z*pin(ip,nu)+pt*vperp(nu)
            qj(nu)=-z*(one-x)/(x+tiny)*pin(kp,nu) &
            +(one-z)*pin(ip,nu)-pt*vperp(nu)
        enddo
                
        wtdip=-8._dp*sik/x**2/16._dp/pi**2

    !--- note: we may end up with initial state momentum qi with
    !---       energy larger than that of the beam; that should be
    !---       rejected outside of this routine
                
    elseif (diptype == fifi) then
    !---    final-final dipoles
        z=half*(one-(one-two*x2)*sqrt(one+four*x2*(one-x2)))
        y=x1**2
        phi=two*pi*x3
                
    !--- account for dipoles with emitted partons not in the last position
        ipp=ip
        kpp=kp
        if (jp < ip) ipp=jp
        if (jp < kp) kpp=jp
        if (ipp == kpp) then
            write(6,*) 'unanticipated dipole in multichan.f'
            write(6,*) 'ip,jp,kp = ',ip,jp,kp
            stop
        endif

        sik=two*dot(pin,ipp,kpp)
        if ((pin(ip,4) == 0._dp) .or. (pin(kp,4) == 0._dp)) then
            multichan = .false.
            return
        endif
    !--- catch generation of exceptionally collinear points
        if (abs(sik/(two*pin(ipp,4)*pin(kpp,4))) < tiny) then
            multichan = .false.
            return
        endif
        ptsq=sik*z*(one-z)*y
        if (ptsq < 0) then
            multichan = .false.
            return
        endif
        pt=sqrt(ptsq)
                
        if (.not. getperp(pin,ipp,kpp,v1,v2)) then
            multichan = .false.
            return
        endif
                
        do nu=1,4
            vperp(nu)=cos(phi)*v1(nu)+sin(phi)*v2(nu)
                    
            qi(nu)=z*pin(ipp,nu)+y*(one-z)*pin(kpp,nu)+pt*vperp(nu)
            qj(nu)=(one-z)*pin(ipp,nu)+y*z*pin(kpp,nu)-pt*vperp(nu)
            qk(nu)=(one-y)*pin(kpp,nu)
        enddo
                
        wtdip=8._dp*sik*(one-y)/16._dp/pi**2

    elseif (diptype == phot) then
    !---    identified photon dipoles
        z=x2**2
        y=x1
        phi=two*pi*x3
                
        if (npart == 3) then

        !--- work out other spectator
            mp=3+4+5-ip-jp
                    
            sim=two*dot(pin,3,4)
            ptsq=sim*y*(one-y)*(1-z)
            if (ptsq < 0) then
                multichan = .false.
                return
            endif
            pt=sqrt(ptsq)
                    
            if (.not. getperp(pin,ip,mp,v1,v2)) then
                multichan = .false.
                return
            endif
                    
            do nu=1,4
                vperp(nu)=cos(phi)*v1(nu)+sin(phi)*v2(nu)
                        
                qi(nu)=z*pin(ip,nu)
                qj(nu)=(one-z)*(one-y)*pin(ip,nu)+y*pin(mp,nu)+pt*vperp(nu)
                qm(nu)=(one-z)*y*pin(ip,nu)+(one-y)*pin(mp,nu)-pt*vperp(nu)
            enddo

            wtdip=two*z*sim/16._dp/pi**2
                    
        else
            write(6,*) 'Identified photon dipole not written for npart>3'
            stop
        endif

    else
    !---    error or unimplemented type
        write(6,*) 'Unrecognized dipole type in multichan.f: ',diptype
        stop
    endif

!--- check for accidental NaN
    if (qi(4) /= qi(4)) then
        write(6,*) 'multichan: qi NaN, diptype = ',diptype
    !        write(6,*) 'z,pt,v1,v2',z,pt,v1,v2
        multichan = .false.
        return
    endif
    if (qj(4) /= qj(4)) then
        write(6,*) 'multichan: qj NaN, diptype = ',diptype
        multichan = .false.
        return
    endif
    if (qk(4) /= qk(4)) then
        write(6,*) 'multichan: qk NaN, diptype = ',diptype
        multichan = .false.
        return
    endif
          
!--- Step 5: reinstate all momenta
          
!c---    note: For now we assume that emitted partons in the
!c---          dipoles should always be added at the end of the list of
!c---          n-parton momenta.
!c---          This is true for simple cases (e.g. Drell-Yan, Wbb), but
!c---          not for more complicated ones.
!        if (jp .ne. npart+2) then
!          write(6,*) 'Dipole multichannel phase space, multichan.f,'
!          write(6,*) 'assumes that emitted radiation is the last'
!          write(6,*) 'momentum, but this is not the case: jp=',jp
!          stop
!        endif

    if    (diptype == inin) then
    !---    initial-initial dipoles
        do nu=1,4
            pout(ip,nu)=qi(nu)
            pout(kp,nu)=qk(nu)
            pout(jp,nu)=qj(nu)
        !---       must now perform LT on remaining final state partons
        !---        (c.f. relevant section of transform.f)
            k(nu)=-qi(nu)-qk(nu)-qj(nu)
            kt(nu)=-x*qi(nu)-qk(nu)
            ks(nu)=k(nu)+kt(nu)
        enddo
        kDk=dotpr(k,k)
        ksDks=dotpr(ks,ks)
        if ((kDk == 0._dp) .or. (ksDks == 0._dp)) then
            multichan = .false.
            return
        endif
        iskip=0
        do m=3,npart+2
            if (m == jp) then
                iskip=1
            else
!                ktDp=kt(4)*pin(m-iskip,4)-kt(1)*pin(m-iskip,1) &
!                -kt(2)*pin(m-iskip,2)-kt(3)*pin(m-iskip,3)
!                ksDp=ks(4)*pin(m-iskip,4)-ks(1)*pin(m-iskip,1) &
!                -ks(2)*pin(m-iskip,2)-ks(3)*pin(m-iskip,3)
                ktDp=dotpr(kt,pin(m-iskip,:))
                ksDp=dotpr(ks,pin(m-iskip,:))
                do nu=1,4
                    pout(m,nu)=pin(m-iskip,nu)-two*ksDp*ks(nu)/ksDks &
                    +two*ktDp*k(nu)/kDk
                enddo
            endif
        enddo
                
    elseif (( diptype == infi) &
         .OR. (diptype == fiin) &
         .OR. (diptype == fifi) ) then
    !---    all other dipoles are straightforward
        do nu=1,4
            pout(ip,nu)=qi(nu)
            pout(kp,nu)=qk(nu)
            pout(jp,nu)=qj(nu)
        enddo
        if (swapjpkp) then
           pout(kp,:)=qj(:)
           pout(jp,:)=qk(:)
        endif
        iskip=0
        do m=1,npart+2
            if     (m == jp) then
                iskip=1
            elseif ((m /= ip) .AND. (m /= kp)) then
                do nu=1,4
                    pout(m,nu)=pin(m-iskip,nu)
                enddo
            endif
        enddo

    elseif (diptype == phot) then
    !---    identified photon dipole
        if (npart == 3) then
            do nu=1,4
                pout(1,nu)=pin(1,nu)
                pout(2,nu)=pin(2,nu)
                pout(ip,nu)=qi(nu)
                pout(jp,nu)=qj(nu)
                pout(mp,nu)=qm(nu)
            enddo
        else
            write(6,*) 'Should not get here, npart=',npart
            stop
        endif

    else
    !---    error or unimplemented type
        write(6,*) 'Unrecognized dipole type in multichan.f: ',diptype
        stop
    endif

!--- now perform permutations if necessary
    if     ((ichan > maxdip) .AND. (ichan <= 2*maxdip)) then
        ptmp(:)=pout(perm2(1),:)
        pout(perm2(1),:)=pout(perm2(2),:)
        pout(perm2(2),:)=ptmp(:)
    elseif ((ichan > 2*maxdip) .AND. (ichan <= 3*maxdip)) then
        ptmp(:)=pout(perm3(1),:)
        pout(perm3(1),:)=pout(perm3(2),:)
        pout(perm3(2),:)=ptmp(:)
    endif

!--- Step 6: compute the sum of inverse Jacobians for all channels
    xjac=zip
          
    do m=1,maxdip*maxperms

        if (m <= maxdip) then
            ip=dipconfig(m,1)
            jp=dipconfig(m,2)
            kp=dipconfig(m,3)
        elseif ((m > maxdip) .AND. (m <= 2*maxdip)) then
            ip=dipconfig(m-maxdip,1)
            jp=dipconfig(m-maxdip,2)
            kp=dipconfig(m-maxdip,3)
            if     (ip == perm2(1)) then
                ip=perm2(2)
            elseif (ip == perm2(2)) then
                ip=perm2(1)
            endif
            if     (jp == perm2(1)) then
                jp=perm2(2)
            elseif (jp == perm2(2)) then
                jp=perm2(1)
            endif
            if     (kp == perm2(1)) then
                kp=perm2(2)
            elseif (kp == perm2(2)) then
                kp=perm2(1)
            endif
        elseif ((m > 2*maxdip) .AND. (m <= 3*maxdip)) then
            ip=dipconfig(m-2*maxdip,1)
            jp=dipconfig(m-2*maxdip,2)
            kp=dipconfig(m-2*maxdip,3)
            if     (ip == perm3(1)) then
                ip=perm3(2)
            elseif (ip == perm3(2)) then
                ip=perm3(1)
            endif
            if     (jp == perm3(1)) then
                jp=perm3(2)
            elseif (jp == perm3(2)) then
                jp=perm3(1)
            endif
            if     (kp == perm3(1)) then
                kp=perm3(2)
            elseif (kp == perm3(2)) then
                kp=perm3(1)
            endif
        endif

        if     ((ip <= 2) .AND. (kp <= 2)) then
            diptype=inin
        elseif ((ip <= 2) .AND. (kp > 2)) then
            diptype=infi
        elseif ((ip > 2) .AND. (kp <= 2)) then
            diptype=fiin
        elseif ((ip > 2) .AND. (kp > 2)) then
            diptype=fifi
        endif
    !--- special notation for identified photon dipole
        if (kp == 0) then
            diptype=phot
        endif

    !--- compute relevant invariants
        if (diptype /= phot) then
            sij=two*dot(pout,ip,jp)
            sik=two*dot(pout,ip,kp)
            sjk=two*dot(pout,jp,kp)
        else
            mp=3+4+5-ip-jp
            sij=two*dot(pout,ip,jp)
            sim=two*dot(pout,ip,mp)
            sjm=two*dot(pout,jp,mp)
        endif

        if     (diptype == inin) then
        !---    initial-initial dipoles
            if (sik == 0._dp) then
              multichan = .false.
              return
            endif
            x=one+(sij+sjk)/sik
            vtilde=-sij/sik
            if ((min(x,vtilde) < 0._dp) .or. (max(x,vtilde) > 1._dp)) then
              multichan = .false.
              return
            endif
            xjac=xjac+(sqrt(vtilde)+sqrt(one-vtilde)) &
            /(sqrt(vtilde*(one-vtilde)*(one-x))+tiny)
        elseif (diptype == infi) then
        !---    initial-final dipoles
            if (sij+sik == 0._dp) then
              multichan = .false.
              return
            endif
            x=one+sjk/(sij+sik)
            u=sij/(sij+sik)
            if ((min(x,u) < 0._dp) .or. (max(x,u) > 1._dp)) then
              multichan = .false.
              return
            endif
            xjac=xjac+(sqrt(u)+sqrt(one-u)) &
            /(sqrt(u*(one-u)*(one-x))+tiny)

        elseif (diptype == fiin) then
        !---    final-initial dipoles
            if (sjk+sik == 0._dp) then
              multichan = .false.
              return
            endif
            x=one+sij/(sjk+sik)
            z=sik/(sik+sjk)
            if ((min(x,z) < 0._dp) .or. (max(x,z) > 1._dp)) then
              multichan = .false.
              return
            endif
            xjac=xjac+(sqrt(z)+sqrt(one-z)) &
            /(sqrt(z*(one-z)*(one-x))+tiny)

        elseif (diptype == fifi) then
        !---    final-final dipoles
            if (sjk+sik == 0._dp) then
              multichan = .false.
              return
            endif
            y=sij/(sij+sjk+sik)
            z=sik/(sjk+sik)
            if ((min(y,z) < 0._dp) .or. (max(y,z) > 1._dp)) then
              multichan = .false.
              return
            endif
            xjac=xjac+(sqrt(z)+sqrt(one-z)) &
            /(sqrt(z*(one-z)*y)+tiny)

        elseif (diptype == phot) then
        !---    identified photon dipoles
            z=(sij+sim)/(sij+sim+sjm)
            xjac=xjac+one/(sqrt(z)+tiny)

        else
        !---    error or unimplemented type
            write(6,*) 'Unrecognized dipole type in multichan.f: ',diptype
            stop
        endif

    enddo
          
!--- Step 7: compute final weight
    wtdip=wtdip*real(maxdip*maxperms,dp)/xjac

!--- catch NaN, typically generated by rounding errors in generation
!--- of momenta causing p^2<0
    if (wtdip /= wtdip) then
        wtdip=zip
        multichan = .false.
        return
    endif

!--- for checking
!      if (ichan == 8) then
!      write(6,*) 'wtdip = ',wtdip
!      write(6,*) 'PIN'
!      call writeout(pin)
!      write(6,*)
!      write(6,*) 'POUT'
!      call writeout(pout)
!      pause
!      endif
                
    multichan = .true.
    return
    end function multichan

    subroutine dipoleconfig(maxdip,dipconfig,maxperms,perm2,perm3)
    use VVconfig_m
    use types
    implicit none
!***********************************************************************
!                                                                      *
!     Given a process, return an array that specifies all the          *
!     dipole configurations that are present in the real subtractions  *
!                                                                      *
!          maxdip: maximum no. of dipoles for this process             *
!       dipconfig: array of dipole configs (i,j,k) in standard MCFM    *
!                  notation, i.e. (emitter, emitted, spectator)        *
!                                                                      *
!     Author: J.M.Campbell                                             *
!       Date: 19th March 2009                                          *
!                                                                      *
!***********************************************************************
          
    include 'constants.f'
    include 'nf.f'
    include 'mxpart.f'
    include 'cplx.h'
    include 'kprocess.f'
    include 'ptilde.f'
    include 'frag.f'
    include 'taucut.f'
    include 'kpart.f'
    include 'npart.f'
    include 'hdecaymode.f'
    include 'ipsgen.f'
    include 'nwz.f'
    integer:: maxdip,dipconfig(maxd,3),i7,i8
    integer:: maxperms,perm2(2),perm3(2)

    integer :: dipstage
    common/dipstage/dipstage
!$omp threadprivate(/dipstage/)

! default: do not use any permutations to generate dipoles
    maxperms=1
          
    if    ((kcase==kW_only) .OR. (kcase==kZ_only) &
      .OR. (kcase==kggfus0) .OR. (kcase==kHgagaI)) then
        maxdip=2
        dipconfig(1,:)= (/ 1,5,2 /)
        dipconfig(2,:)= (/ 2,5,1 /)
    elseif ((kcase==kWln_ew) .or. (kcase==kWln_aq)) then
      if (ipsgen == 1) then
        maxdip=2
        dipconfig(1,:)= (/ 1,5,2 /)
        dipconfig(2,:)= (/ 2,5,1 /)
      elseif (ipsgen == 2) then
        maxdip=1
        if (nwz == +1) then
          dipconfig(1,:)= (/ 4,5,3 /)
        else
          dipconfig(1,:)= (/ 3,5,4 /)
        endif
      endif
    elseif (kcase==kHi_Zga) then
        maxdip=2
        dipconfig(1,:)= (/ 1,6,2 /)
        dipconfig(2,:)= (/ 2,6,1 /)
    elseif (kcase==kWgaj_a) then
        maxdip=2
        dipconfig(1,:)= (/ 1,6,2 /)
        dipconfig(2,:)= (/ 2,6,1 /)
    elseif (kcase==kWgajja) then
        maxdip=2
        dipconfig(1,:)= (/ 1,7,2 /)
        dipconfig(2,:)= (/ 2,7,1 /)
        maxperms=2
        perm2(1)=6
        perm2(2)=7
    elseif (kcase==kWga_ew) then
      if (ipsgen == 1) then
        maxdip=2
        dipconfig(1,:)= (/ 1,6,2 /)
        dipconfig(2,:)= (/ 2,6,1 /)
        maxperms=2
        perm2(1)=5
        perm2(2)=6
      elseif (ipsgen == 2) then
        maxdip=1
        if (nwz == +1) then
          dipconfig(1,:)= (/ 4,6,3 /)
        else
          dipconfig(1,:)= (/ 3,6,4 /)
        endif
        maxperms=2
        perm2(1)=5
        perm2(2)=6
      else
        maxdip=1
        if (nwz == +1) then
          dipconfig(1,:)= (/ 4,6,3 /)
        else
          dipconfig(1,:)= (/ 3,6,4 /)
        endif
      endif
    elseif (kcase==kWgajew) then
      if (ipsgen == 1) then
        maxdip=4
        dipconfig(1,:)= (/ 1,7,2 /)
        dipconfig(2,:)= (/ 2,7,1 /)
        dipconfig(3,:)= (/ 6,7,1 /)
        dipconfig(4,:)= (/ 6,7,2 /)
        maxperms=2
        perm2(1)=5
        perm2(2)=7
      elseif (ipsgen == 2) then
        maxdip=1
        if (nwz == +1) then
          dipconfig(1,:)= (/ 4,7,3 /)
        else
          dipconfig(1,:)= (/ 3,7,4 /)
        endif
        maxperms=2
        perm2(1)=5
        perm2(2)=7
      else
        maxdip=1
        if (nwz == +1) then
          dipconfig(1,:)= (/ 4,7,3 /)
        else
          dipconfig(1,:)= (/ 3,7,4 /)
        endif
      endif
    elseif ((kcase==kWWqqbr) .OR. (kcase==kWWnpol) &
         .OR. (kcase==kWZbbar) .OR. (kcase==kZZlept) &
         .OR. (kcase==kWHbbar) .OR. (kcase==kZHbbar) &
         .OR. (kcase==kWHgaga) .OR. (kcase==kZHgaga) &
         .OR. (kcase==kHZZ_4l) .OR. (kcase==kHmZZ4l) &
         .OR. (kcase==kZHlept) .OR. (kcase==kVVlept)) then
        maxdip=2
        dipconfig(1,:)= (/ 1,7,2 /)
        dipconfig(2,:)= (/ 2,7,1 /)
    elseif (kcase==kqq_Hqq) then
        maxdip=8
        dipconfig(1,:)= (/ 1,7,5 /)
        dipconfig(2,:)= (/ 5,7,1 /)
        dipconfig(3,:)= (/ 2,7,6 /)
        dipconfig(4,:)= (/ 6,7,2 /)
        dipconfig(5,:)= (/ 1,5,2 /)
        dipconfig(6,:)= (/ 2,6,1 /)
        dipconfig(7,:)= (/ 1,6,2 /)
        dipconfig(8,:)= (/ 2,7,1 /)
    elseif (kcase==kHi_Zaj) then
        maxdip=10
        dipconfig(1,:)= (/ 1,6,2 /)
        dipconfig(2,:)= (/ 2,6,1 /)
        dipconfig(3,:)= (/ 1,7,2 /)
        dipconfig(4,:)= (/ 2,7,1 /)
        dipconfig(5,:)= (/ 1,6,7 /)
        dipconfig(6,:)= (/ 6,7,1 /)
        dipconfig(7,:)= (/ 1,7,6 /)
        dipconfig(8,:)= (/ 2,7,6 /)
        dipconfig(9,:)= (/ 6,7,2 /)
        dipconfig(10,:)=(/ 2,6,7 /)
    elseif ((kcase==kW_1jet) .OR. (kcase==kZ_1jet) &
         .OR. (kcase==kggfus1) .OR. (kcase==khjetma) &
         .OR. (kcase==kgmgmjt) .OR. (kcase==kHgagaj)) then
        maxdip=6
        dipconfig(1,:)= (/ 1,6,2 /)
        dipconfig(2,:)= (/ 2,6,1 /)
        dipconfig(3,:)= (/ 1,6,5 /)
        dipconfig(4,:)= (/ 2,6,5 /)
        dipconfig(5,:)= (/ 5,6,1 /)
        dipconfig(6,:)= (/ 5,6,2 /)

        if ((frag) .AND. (kcase==kgmgmjt)) then
            maxdip=8
            dipconfig(7,:)=(/ 3,6,1 /)
            dipconfig(8,:)=(/ 4,6,1 /)
        endif
    ! generate remaining combinations by permuting
    !  perms2: exchange 5 and 6
        maxperms=2
        perm2(1)=5
        perm2(2)=6

    elseif (kcase==kbq_tpq_jet) then
        maxperms = 1
        if (kpart == klord .OR. kpart == kvirt) then
            maxdip=2
        ! for corr_on_beam == 1
        !             if (corr_on_beam == 1) then
            dipconfig(1,:) = [1,7,6]
            dipconfig(2,:) = [6,7,1]
        !                 !dipconfig(1,:) = [6,7,1]
        !             else
        !                 dipconfig(1,:) = [2,7,6]
        !                 dipconfig(2,:) = [6,7,2]
        !             endif
        elseif (kpart == kreal) then
            maxdip = 1
            if (dipstage == 1) then
            ! ipconfig(1,:) = [1,7,6]
                dipconfig(1,:) = [6,7,1]
            elseif (dipstage == 2) then
                dipconfig(1,:) = [6,8,1]
            ! ipconfig(1,:) = [6,8,7]
            else
                dipconfig(1,:) = [7,8,6]
            endif
        else
            error stop "dipoleconfig"

        endif
    !       maxdip=6
    !       if (corr_on_beam == 1) then
    !           dipconfig(1,:)=  (/ 1,7,6 /)
    !           dipconfig(2,:)=  (/ 1,7,8 /)
    !           dipconfig(3,:)=  (/ 6,7,1 /)
    !           dipconfig(4,:)=  (/ 7,8,1 /)
    !           dipconfig(5,:)=  (/ 6,7,8 /)
    !           dipconfig(6,:)=  (/ 7,8,6 /)
    !       else
    !           dipconfig(1,:)=  (/ 2,7,6 /)
    !           dipconfig(2,:)=  (/ 2,7,8 /)
    !           dipconfig(3,:)=  (/ 6,7,2 /)
    !           dipconfig(4,:)= (/ 7,8,2 /)
    !           dipconfig(5,:)= (/ 1,6,7 /)
    !           dipconfig(6,:)= (/ 2,6,7 /)
    !       endif

    ! generate remaining combinations by permuting
    !  perms2: exchange 7 and 8
    ! axperms=2
    ! erm2(1)=7
    ! erm2(2)=8

    elseif ((kcase==kW_2jet) .OR. (kcase==kZ_2jet) &
         .OR. (kcase==kggfus2) .OR. (kcase==kgagajj)) then
        maxdip=12
        dipconfig(1,:)= (/ 1,7,2 /)
        dipconfig(2,:)= (/ 2,7,1 /)
        dipconfig(3,:)= (/ 1,7,5 /)
        dipconfig(4,:)= (/ 2,7,6 /)
        dipconfig(5,:)= (/ 1,7,6 /)
        dipconfig(6,:)= (/ 2,7,5 /)
        dipconfig(7,:)= (/ 5,7,1 /)
        dipconfig(8,:)= (/ 6,7,2 /)
        dipconfig(9,:)= (/ 6,7,1 /)
        dipconfig(10,:)=(/ 5,7,2 /)
        dipconfig(11,:)=(/ 5,7,6 /)
        dipconfig(12,:)=(/ 6,7,5 /)

    ! generate remaining combinations by permuting
    !  perms2: exchange 6 and 7, perm3: exchange 5 and 7
        maxperms=3
        perm2(1)=6
        perm2(2)=7
        perm3(1)=5
        perm3(2)=7

    !        dipconfig(13,:)=(/ 1,6,2 /)
    !        dipconfig(14,:)=(/ 2,6,1 /)
    !        dipconfig(15,:)=(/ 1,6,5 /)
    !        dipconfig(16,:)=(/ 2,6,7 /)
    !        dipconfig(17,:)=(/ 1,6,7 /)
    !        dipconfig(18,:)=(/ 2,6,5 /)
    !        dipconfig(19,:)=(/ 5,6,1 /)
    !c        dipconfig(20,:)=(/ 7,6,2 /)
    !c        dipconfig(21,:)=(/ 7,6,1 /)
    !        dipconfig(20,:)=(/ 5,6,2 /)
    !        dipconfig(21,:)=(/ 5,6,7 /)
    !c        dipconfig(24,:)=(/ 7,6,5 /)
    
    !        dipconfig(22,:)=(/ 1,5,2 /)
    !        dipconfig(23,:)=(/ 2,5,1 /)
    !        dipconfig(24,:)=(/ 1,5,7 /)
    !        dipconfig(25,:)=(/ 2,5,6 /)
    !        dipconfig(26,:)=(/ 1,5,6 /)
    !        dipconfig(27,:)=(/ 2,5,7 /)
    !!        dipconfig(28,:)=(/ 7,5,1 /)
    !c        dipconfig(32,:)=(/ 6,5,2 /)
    !c        dipconfig(33,:)=(/ 6,5,1 /)
    !c        dipconfig(34,:)=(/ 7,5,2 /)
    !c        dipconfig(35,:)=(/ 7,5,6 /)
    !c        dipconfig(36,:)=(/ 6,5,7 /)
    elseif (kcase==kgam_2j) then
        maxdip=12
        dipconfig(1,:)= (/ 1,6,2 /)
        dipconfig(2,:)= (/ 2,6,1 /)
        dipconfig(3,:)= (/ 1,6,4 /)
        dipconfig(4,:)= (/ 2,6,5 /)
        dipconfig(5,:)= (/ 1,6,5 /)
        dipconfig(6,:)= (/ 2,6,4 /)
        dipconfig(7,:)= (/ 4,6,1 /)
        dipconfig(8,:)= (/ 5,6,2 /)
        dipconfig(9,:)= (/ 5,6,1 /)
        dipconfig(10,:)=(/ 4,6,2 /)
        dipconfig(11,:)=(/ 4,6,5 /)
        dipconfig(12,:)=(/ 5,6,4 /)

! generate remaining combinations by permuting
!  perms2: exchange 5 and 6, perm3: exchange 4 and 6
        maxperms=3
        perm2(1)=5
        perm2(2)=6
        perm3(1)=4
        perm3(2)=6

    elseif ((kcase==kH_tjet) .OR. (kcase==kZ_tjet)) then
        maxdip=8
        dipconfig(1,:)= (/ 1,7,6 /)
        dipconfig(2,:)= (/ 6,7,1 /)
        dipconfig(3,:)= (/ 2,7,6 /)
        dipconfig(4,:)= (/ 6,7,2 /)
        dipconfig(5,:)= (/ 1,7,2 /)
        dipconfig(6,:)= (/ 2,7,1 /)
        dipconfig(7,:)= (/ 1,6,2 /)
        dipconfig(8,:)= (/ 2,6,1 /)
    elseif ((kcase==kH_tdkj) .OR. (kcase==kZ_tdkj)) then
        maxdip=8
        dipconfig(1,:)= (/ 1,9,8 /)
        dipconfig(2,:)= (/ 8,9,1 /)
        dipconfig(3,:)= (/ 2,9,8 /)
        dipconfig(4,:)= (/ 8,9,2 /)
        dipconfig(5,:)= (/ 1,9,2 /)
        dipconfig(6,:)= (/ 2,9,1 /)
        dipconfig(7,:)= (/ 1,8,2 /)
        dipconfig(8,:)= (/ 2,8,1 /)
    elseif (kcase==khflgam) then
        maxdip=6
        dipconfig(1,:)= (/ 1,5,2 /)
        dipconfig(2,:)= (/ 2,5,1 /)
        dipconfig(3,:)= (/ 4,5,1 /)
        dipconfig(4,:)= (/ 1,5,4 /)
        dipconfig(5,:)= (/ 4,5,2 /)
        dipconfig(6,:)= (/ 2,5,4 /)
    elseif ((kcase==kgamgam) .OR. (kcase==kgg2gam)) then
        dipconfig(1,:)= (/ 1,5,2 /)
        dipconfig(2,:)= (/ 2,5,1 /)
        if (frag) then
            dipconfig(3,:)= (/ 3,5,0 /)
            dipconfig(4,:)= (/ 4,5,0 /)
            maxdip=4
        else
            maxdip=2
        endif
    elseif (kcase==ktwo_ew) then
        dipconfig(1,:)= (/ 1,5,2 /)
        dipconfig(2,:)= (/ 2,5,1 /)
        dipconfig(3,:)= (/ 3,5,1 /)
        dipconfig(4,:)= (/ 4,5,1 /)
        dipconfig(5,:)= (/ 1,3,2 /)
        dipconfig(6,:)= (/ 2,4,1 /)
        dipconfig(7,:)= (/ 1,4,2 /)
        dipconfig(8,:)= (/ 2,3,1 /)
        maxdip=8
    elseif ((kcase==kdirgam) .OR. (kcase==kgamjet)) then
        maxdip=6
        dipconfig(1,:)= (/ 1,5,2 /)
        dipconfig(2,:)= (/ 2,5,1 /)
        dipconfig(3,:)= (/ 1,5,4 /)
        dipconfig(4,:)= (/ 2,5,4 /)
        dipconfig(5,:)= (/ 4,5,1 /)
        dipconfig(6,:)= (/ 4,5,2 /)
        if (frag) then
            maxdip=7
            dipconfig(7,:)= (/ 3,5,1 /)
        endif

    ! generate remaining combinations by permuting
    !  perms2: exchange 4 and 5
        maxperms=2
        perm2(1)=4
        perm2(2)=5

    elseif (kcase==kWgamma) then
        dipconfig(1,:)= (/ 1,6,2 /)
        dipconfig(2,:)= (/ 2,6,1 /)
        if (frag) then
            dipconfig(3,:)= (/ 5,6,2 /)
            maxdip=3
        else
            maxdip=2
        endif
    elseif (kcase==kZgamma) then
        if (decayChannel() == decayQuarks) then
            dipconfig(1,:) = (/ 3,6,4 /)
            dipconfig(2,:) = (/ 4,6,3 /)
            maxdip = 2
        else
            dipconfig(1,:) = (/ 1,6,2 /)
            dipconfig(2,:) = (/ 2,6,1 /)
            if (frag) then
                dipconfig(3,:) = (/ 5,6,2 /)
                maxdip=3
            else
                maxdip=2
            endif
        endif
    elseif (kcase==kZgajet) then
        if (decayChannel() == decayQuarks) then
            dipconfig(1,:)  = (/ 3,7,4 /)
            dipconfig(2,:)  = (/ 3,7,6 /)
            dipconfig(3,:)  = (/ 4,7,3 /)
            dipconfig(4,:)  = (/ 4,7,6 /)
            dipconfig(5,:)  = (/ 6,7,3 /)
            dipconfig(6,:)  = (/ 6,7,4 /)
            maxdip = 6
        else
            dipconfig(1,:)= (/ 1,7,2 /)
            dipconfig(2,:)= (/ 2,7,1 /)
            dipconfig(3,:)= (/ 6,7,1 /)
            dipconfig(4,:)= (/ 6,7,2 /)
            if (frag) then
                dipconfig(5,:)= (/ 5,7,2 /)
                maxdip=5
            else
                maxdip=4
            endif
        endif
    elseif ((kcase==kW_2gam) .OR. (kcase==kZ_2gam)) then
        dipconfig(1,:)= (/ 1,7,2 /)
        dipconfig(2,:)= (/ 2,7,1 /)
        if (frag) then
            dipconfig(3,:)= (/ 5,7,2 /)
            dipconfig(4,:)= (/ 6,7,2 /)
            maxdip=4
        else
            maxdip=2
        endif
    elseif (kcase==kdm_gam) then
        dipconfig(1,:)= (/ 1,6,2 /)
        dipconfig(2,:)= (/ 2,6,1 /)
        if (frag) then
            dipconfig(3,:)= (/ 5,6,2 /)
            maxdip=3
        else
            maxdip=2
        endif
    elseif (kcase==ktrigam) then
        dipconfig(1,:)= (/ 1,6,2 /)
        dipconfig(2,:)= (/ 2,6,1 /)
        if (frag) then
            dipconfig(3,:)= (/ 3,6,2 /)
            dipconfig(4,:)= (/ 4,6,2 /)
            dipconfig(5,:)= (/ 5,6,2 /)
            maxdip=5
        else
            maxdip=2
        endif
    elseif (kcase==ktottth) then
        dipconfig(1,:)= (/ 1,6,2 /)
        dipconfig(2,:)= (/ 2,6,1 /)
        maxdip=2
    elseif (kcase==ktt_mix) then
        dipconfig(1,:)= (/ 1,5,2 /)
        dipconfig(2,:)= (/ 2,5,1 /)
        maxdip=2
    elseif (kcase==kW_cjet) then
        dipconfig(1,:)= (/ 1,6,2 /)
        dipconfig(2,:)= (/ 2,6,1 /)
        maxdip=2
    elseif ((kcase==kWH1jet) .OR. (kcase==kZH1jet) &
        .OR.(kcase==kWW_jet) .OR. (kcase==kWZ_jet) &
        .OR.(kcase==kZZ_jet)) then
        if (hdecaymode == 'wpwm') then
            i7=9
            i8=10
        else
            i7=7
            i8=8
        endif
        maxdip=6
        dipconfig(1,:)= (/ 1,i8,2 /)
        dipconfig(2,:)= (/ 2,i8,1 /)
        dipconfig(3,:)= (/ 1,i8,i7 /)
        dipconfig(4,:)= (/ 2,i8,i7 /)
        dipconfig(5,:)= (/ i7,i8,1 /)
        dipconfig(6,:)= (/ i7,i8,2 /)
    ! generate remaining combinations by permuting
    !  perms2: exchange 5 and 6
        maxperms=2
        perm2(1)=i7
        perm2(2)=i8
    elseif (kcase==kWgajet) then
        maxdip=6
        dipconfig(1,:)= (/ 1,7,2 /)
        dipconfig(2,:)= (/ 2,7,1 /)
        dipconfig(3,:)= (/ 1,7,6 /)
        dipconfig(4,:)= (/ 2,7,6 /)
        dipconfig(5,:)= (/ 6,7,1 /)
        dipconfig(6,:)= (/ 6,7,2 /)
! generate remaining combinations by permuting
!  perms2: exchange 6 and 7
        maxperms=2
        perm2(1)=6
        perm2(2)=7
    elseif(kcase==kfourga) then
        dipconfig(1,:)= (/ 1,7,2 /)
        dipconfig(2,:)= (/ 2,7,1 /)
        if (frag) then
            dipconfig(3,:)= (/ 3,7,1 /)
            dipconfig(4,:)= (/ 4,7,1 /)
            dipconfig(5,:)= (/ 5,7,1 /)
            dipconfig(6,:)= (/ 6,7,1 /)
            maxdip=6
        else
            maxdip=2
        endif
    elseif(kcase==kqg_tbq) then
        maxdip=4
        dipconfig(1,:)= (/ 1,6,5 /)
        dipconfig(2,:)= (/ 2,6,5 /)
        dipconfig(3,:)= (/ 5,6,1 /)
        dipconfig(4,:)= (/ 5,6,2 /)
    elseif (kcase==ktopanom) then
        maxdip = 14
        dipconfig(1,:)  = (/ 1,7,6 /)
        dipconfig(2,:)  = (/ 6,7,1 /)
        dipconfig(3,:)  = (/ 2,7,6 /)
        dipconfig(4,:)  = (/ 6,7,2 /)
        dipconfig(5,:)  = (/ 1,7,2 /)
        dipconfig(6,:)  = (/ 1,6,2 /)
        dipconfig(7,:)  = (/ 2,7,1 /)
        dipconfig(8,:)  = (/ 2,6,1 /)
        dipconfig(9,:)  = (/ 2,7,5 /)
        dipconfig(10,:) = (/ 5,7,2 /)
        dipconfig(11,:) = (/ 1,7,5 /)
        dipconfig(12,:) = (/ 5,7,1 /)
        dipconfig(13,:) = (/ 1,7,2 /)
        dipconfig(14,:) = (/ 2,7,1 /)

    elseif (kcase==kbq_tpq) then
        if (ipsgen == 1) then
            maxdip = 1
            dipconfig(1,:)  = (/ 6,7,1 /)
        ! ipconfig(2,:)  = (/ 1,7,6 /)
        !           dipconfig(3,:)  = (/ 2,7,6 /)
        !           dipconfig(4,:)  = (/ 6,7,2 /)
        !           dipconfig(5,:)  = (/ 1,7,2 /)
        !           dipconfig(6,:)  = (/ 1,6,2 /)
        !           dipconfig(7,:)  = (/ 2,7,1 /)
        !           dipconfig(8,:)  = (/ 2,6,1 /)
        elseif (ipsgen == 2) then
            call abort
        elseif (ipsgen == 3) then
            call abort
        else
            error stop "unknown ipsgen in dipoleconfig.f kbq_tpq"
        endif
    !       dipconfig(9,:)  = (/ 2,7,5 /)
    !       dipconfig(10,:) = (/ 5,7,2 /)
    !       dipconfig(11,:) = (/ 1,7,5 /)
    !       dipconfig(12,:) = (/ 5,7,1 /)
    !       dipconfig(13,:) = (/ 1,7,2 /)
    !       dipconfig(14,:) = (/ 2,7,1 /)
    else
! Dipole configurations not yet supplied; this value may be trapped outside this routine
        maxperms=-1
    endif
          
    return
    end subroutine dipoleconfig
          

    function getperp(p,i1,i2,v1,v2)
    use types
    implicit none
!***********************************************************************
!                                                                      *
!     Given two momenta, p(i1,..) and p(i2, ..), construct two         *
!     vectors v1 and v2 that satisfy:                                  *
!                                                                      *
!      v1.p1 = 0, v1.p2 = 0                                            *
!      v2.p1 = 0, v2.p2 = 0                                            *
!      v1.v2 = 0                                                       *
!      v1.v1 = -1, v2.v2 = -1                                          *
!                                                                      *
!     Author: J.M.Campbell                                             *
!       Date: 19th March 2009                                          *
!                                                                      *
!***********************************************************************
          
    include 'constants.f'
    include 'nf.f'
    include 'mxpart.f'
    include 'cplx.h'
    real(dp):: p(mxpart,4),p1(4),p2(4),v1(4),v2(4),p1Dp2, &
    offset,vsq,vDp1,vDp2,v1Dv2,tiny
    real(dp) :: dotpr
    integer:: i1,i2,nu
    parameter(tiny=1.e-15_dp)

    logical :: getperp
          
!      p(i1,1)=   0.0000000000000000._dp
!      p(i1,2)=  0.0000000000000000._dp
!      p(i1,3)=624.49264065711168._dp
!      p(i1,4)= -624.49264065711168._dp
!      p(i2,1)= -8.42527279039759410e-006_dp
!      p(i2,2)= -1.62606796446384395e-006_dp
!      p(i2,3)= -421.23910530482470._dp
!      p(i2,4)=  421.23910530482482._dp

!--- set up vectors p1,p2 and normalize
    do nu=1,4
        p1(nu)=p(i1,nu)/p(i1,4)
        p2(nu)=p(i2,nu)/p(i2,4)
    enddo
          
    p1Dp2=dotpr(p1,p2)
    if (p1Dp2 == 0._dp) then
        getperp = .false.
        return
    endif

!--- construct momentum v1 orthogonal to p1 and p2
          
!--- offset is to avoid exceptional configurations, usually it is zero
    offset=0._dp
          
    22 continue
       
    v1(4)=0._dp
    v1(1)=5._dp+offset
    v1(2)=3._dp-offset
    v1(3)=4._dp-offset
          
    vsq=dotpr(v1,v1)
    vDp1=dotpr(v1,p1)
    vDp2=dotpr(v1,p2)

!--- avoid the case where we construct a momentum v1 with v1.v1=0
    if (abs(vsq-two*vDp1*vDp2/p1Dp2) < tiny) then
        offset=offset+0.7_dp
        if (offset > 3._dp) then
            getperp = .false.
            return      ! avoid looping forever
        endif
        goto 22
    endif
          
    do nu=1,4
        v1(nu)=v1(nu)-(vDp1*p2(nu)+vDp2*p1(nu))/p1Dp2
    enddo

!--- normalize so that v1.v1=-1
    vsq=dotpr(v1,v1)
    if (vsq >= 0._dp) then
        getperp = .false.
        return
    endif
    do nu=1,4
        v1(nu)=v1(nu)/sqrt(-vsq)
    enddo

!--- construct momentum v2 orthogonal to p1 and p2
          
!--- offset is to avoid exceptional configurations, usually it is zero
    offset=0._dp
          
    33 continue
       
    v2(4)=0._dp
    v2(1)=1._dp-offset
    v2(2)=2._dp+offset
    v2(3)=3._dp+offset
          
    vsq=dotpr(v2,v2)
    vDp1=dotpr(v2,p1)
    vDp2=dotpr(v2,p2)

!--- avoid the case where we construct a momentum v2 with v2.v2=0
    if (abs(vsq-two*vDp1*vDp2/p1Dp2) < tiny) then
        offset=offset+0.6_dp
        if (offset > 3._dp) then
            getperp = .false.
            return      ! avoid looping forever
        endif
        goto 33
    endif
          
    do nu=1,4
        v2(nu)=v2(nu)-(vDp1*p2(nu)+vDp2*p1(nu))/p1Dp2
    enddo

!--- now ensure that v1.v2=0
    v1Dv2=dotpr(v1,v2)
    vsq=dotpr(v1,v1)
    if (vsq == 0._dp) then
        getperp = .false.
        return
    endif
    do nu=1,4
        v2(nu)=v2(nu)-v1(nu)*v1Dv2/vsq
    enddo
          
!--- normalize so that v2.v2=-1
    vsq=dotpr(v2,v2)
!--- avoid the case where we construct a momentum v2 with v2.v2=0
!--- ( -vsq should be positive)
    if (-vsq < tiny) then
        offset=offset+0.6_dp
        if (offset > 3._dp) then
            getperp = .false.
            return      ! avoid looping forever
        endif
        goto 33
    endif
          
    do nu=1,4
        v2(nu)=v2(nu)/sqrt(-vsq)
    enddo

!--- for checking
!      write(6,*) 'v1.v1 = ',dotpr(v1,v1)
!      write(6,*) 'v2.v2 = ',dotpr(v2,v2)
!      write(6,*) 'v1.v2 = ',dotpr(v1,v2)
!      write(6,*) 'v1.p1 = ',dotpr(v1,p1)
!      write(6,*) 'v1.p2 = ',dotpr(v1,p2)
!      write(6,*) 'v2.p1 = ',dotpr(v2,p1)
!      write(6,*) 'v2.p2 = ',dotpr(v2,p2)

    if ((v1(4) /= v1(4)) .OR. (v2(4) /= v2(4))) then
    !        write(6,*) 'warning: could not generate momenta in getperp.f'
        getperp = .false.
        return
    !        write(6,*) 'i1,i2',i1,i2
    !        write(6,*) 'p(i1,:)',p(i1,:)
    !        write(6,*) 'p(i2,:)',p(i2,:)
    !        write(6,*) 'p1',p1
    !        write(6,*) 'p2',p2
    !        write(6,*) 'p1Dp2',p1Dp2
    !        write(6,*) 'vDp1',vDp1
    !        write(6,*) 'vDp2',vDp2
    !        write(6,*) 'v1',v1
    !        write(6,*) 'v2',v2
    !        write(6,*) 'offset',offset
    !        stop
    endif

    getperp = .true.
    return
    end function getperp
          
end module
