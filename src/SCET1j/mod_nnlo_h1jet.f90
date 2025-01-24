!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

! note that mod_nnlo_z1jet needs to be kept in sync, since the
! passed_taucut_1jet routine is there, which is reused for H+jet

module nnlo_h1jet
    use types
    use SCET
    use SCET_Jet
    implicit none

    private

    public :: lumxmsq_h1jet

    contains

      subroutine lumxmsq_h1jet(p,xx,z1,z2,QB,order,xmsq_central,central)
          use SCET
          use LHAPDF
      implicit none
!---- Matrix element for H production
!---- in the heavy quark (mt=Infinity) limit.
!---- averaged over initial colours and spins
!     g(-p1) + g(-p2) --> H(-->b(p3)+b~(p4)) + g(p5)
!---
      include 'src/Inc/nf.f'
      include 'src/Inc/mxpart.f'
      include 'src/Inc/masses.f'
      include 'src/Inc/ewcouple.f'
      include 'src/Inc/qcdcouple.f'
      include 'src/Inc/scale.f'
      include 'src/Inc/facscale.f'
      include 'src/Inc/hdecaymode.f'
      include 'src/Inc/taucut.f'
      include 'src/Inc/constants.f'
      include 'src/Inc/beamtype.f'
      integer:: j,k,m,n,order
      real(dp),intent(out):: xmsq_central
      real(dp) :: xmsq
      real(dp):: p(mxpart,4),fac,Asq,xx(2),hard(2)
      real(dp) :: soft1gg(-1:1),soft2gg(-1:3),soft2gg_nab(-1:3)
      real(dp) :: soft1qg(-1:1),soft2qg(-1:3),soft2qg_nab(-1:3)
      real(dp) :: soft1gq(-1:1),soft2gq(-1:3),soft2gq_nab(-1:3)
      real(dp) :: soft1qa(-1:1),soft2qa(-1:3),soft2qa_nab(-1:3)
      real(dp) :: beama0(-5:5),beamb0(-5:5),tauc
      real(dp) :: beama1(-5:5,-1:1),beamb1(-5:5,-1:1)
      real(dp) :: beama2(-5:5,-1:3),beamb2(-5:5,-1:3)
      real(dp) :: jet1q(-1:1),jet2q(-1:3),jet1g(-1:1),jet2g(-1:3)
      real(dp) :: z1,z2,QB(2),lum0,lum1(-1:1),lum2(-1:3),bit
      real(dp) :: msqgg,msqqa,msqqg,msqgq,assemblejet,y12,y15,y25,Iijm(5)
      real(dp) :: hdecay,msqgamgam,ss,tt,uu,mhsq,s(mxpart,mxpart),s34
      real(dp):: gg,qqb,qg,gqb
      real(dp):: Hgg(2),Hqqb(2),Hqg(2),Hgqb(2)
      integer, parameter:: iglue=5

      logical, intent(in) :: central
      character(len=1) num

      real(dp) :: ptpure, puremass

      call dotem(iglue,p,s)

!--- note: important to use this statement (rather than dot product) in
!--- order to retain the possibility of massive decay products
      s34=(p(3,4)+p(4,4))**2 &
         -(p(3,1)+p(4,1))**2-(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2

      if (hdecaymode == 'tlta') then
          call htautaudecay(p,3,4,hdecay)
      elseif (hdecaymode == 'bqba') then
          call hbbdecay(p,3,4,hdecay)
      elseif (hdecaymode == 'gaga') then
          hdecay=msqgamgam(sqrt(s34))
      else
      write(6,*) 'Unimplemented process in gg_hg'
      stop
      endif
      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)

      call Hampgggsq(order,5,1,2,3,4,gg,Hgg(1),Hgg(2),qg,Hqg(1),Hqg(2))
      call Hampgggsq(order,2,5,1,3,4,gg,Hgg(1),Hgg(2),gqb,Hgqb(1),Hgqb(2))
      call Hampgggsq(order,2,1,5,3,4,gg,Hgg(1),Hgg(2),qqb,Hqqb(1),Hqqb(2))
      
      Asq=(as/(3._dp*pi))**2/vevsq
      fac=Asq*gsq*hdecay

      msqgg=avegg*fac*V*xn*gg
      msqqa=+aveqq*fac*V/two*qqb
      msqgq=-aveqg*fac*V/two*gqb
      msqqg=-aveqg*fac*V/two*qg

      if (order >= 0) then
      call fdist(ih1,xx(1),facscale,beama0,1)
      call fdist(ih2,xx(2),facscale,beamb0,2)
      endif
      if (order >= 1) then
      call xbeam1bis(ih1,z1,xx(1),QB(1),beama1,1)
      call xbeam1bis(ih2,z2,xx(2),QB(2),beamb1,2)
      endif
      if (order >= 2) then
      call xbeam2bis(ih1,z1,xx(1),QB(1),beama2,1)
      call xbeam2bis(ih2,z2,xx(2),QB(2),beamb2,2)
      endif

      call jetq(order,two*p(iglue,4),jet1q,jet2q, disttau=.true.)
      call jetg(order,two*p(iglue,4),jet1g,jet2g)
      
      y12=s(1,2)/p(1,4)/p(2,4)/four
      y15=s(1,iglue)/p(1,4)/p(iglue,4)/four
      y25=s(2,iglue)/p(2,4)/p(iglue,4)/four

      call computeIijm(y12,y25,y15,Iijm)
      call soft_ab_ggg(order,y12,y25,y15,Iijm,1,2,3,4,5,3,soft1gg,soft2gg)
      call soft_ab_qgq(order,y12,y25,y15,Iijm,1,2,3,4,5,3,soft1qg,soft2qg)
      call soft_ab_qgq(order,y12,y15,y25,Iijm,2,1,3,5,4,3,soft1gq,soft2gq)
      call soft_ab_qag(order,y12,y25,y15,Iijm,1,2,3,4,5,3,soft1qa,soft2qa)
      if (order > 1) then
        call soft_nab_ggg(order,y12,y25,y15,Iijm,1,2,3,4,5,3,soft2gg_nab)
        call soft_nab_qgq(order,y12,y25,y15,Iijm,1,2,3,4,5,3,soft2qg_nab)
        call soft_nab_qgq(order,y12,y15,y25,Iijm,2,1,3,5,4,3,soft2gq_nab)
        call soft_nab_qag(order,y12,y25,y15,Iijm,1,2,3,4,5,3,soft2qa_nab)
        soft2gg(:)=soft2gg(:)+soft2gg_nab(:)
        soft2qg(:)=soft2qg(:)+soft2qg_nab(:)
        soft2gq(:)=soft2gq(:)+soft2gq_nab(:)
        soft2qa(:)=soft2qa(:)+soft2qa_nab(:)
      endif

      if (dynamictau) then
          tauc = getdynamictau(p, taucut)
      else
          tauc = taucut
      endif
        
      xmsq=zip
      do j=-nf,nf
      do k=-nf,nf
      
      if ((j .ne. 0) .and. (k .ne. 0) &
       .and. (j .ne. -k)) cycle
      
      if (((j == 0) .and. (k == 0)) &
       .or. (j*k .ne. 0)) then !  gg -> g or qa -> g => jet is a gluon
        if ((j == 0) .and. (k == 0)) then
          bit=assemblejet(order,tauc, &
           beama0(j),beamb0(k),beama1(j,:),beamb1(k,:), &
           beama2(j,:),beamb2(k,:),soft1gg,soft2gg,jet1g,jet2g,Hgg)
        else
          bit=assemblejet(order,tauc, &
           beama0(j),beamb0(k),beama1(j,:),beamb1(k,:), &
           beama2(j,:),beamb2(k,:),soft1qa,soft2qa,jet1g,jet2g,Hqqb)
        endif
      else                     ! qg -> q  => jet is a quark
        if (j == 0) then
          bit=assemblejet(order,tauc, &
           beama0(j),beamb0(k),beama1(j,:),beamb1(k,:), &
           beama2(j,:),beamb2(k,:),soft1gq,soft2gq,jet1q,jet2q,Hgqb)
        elseif (k == 0) then
          bit=assemblejet(order,tauc, &
           beama0(j),beamb0(k),beama1(j,:),beamb1(k,:), &
           beama2(j,:),beamb2(k,:),soft1qg,soft2qg,jet1q,jet2q,Hqg)
        endif
      endif
      
      if     ((j == 0) .and. (k == 0)) then
        bit=bit*msqgg
      elseif (j == 0) then
        bit=bit*msqgq
      elseif (k == 0) then
        bit=bit*msqqg
      elseif (j == -k) then
        bit=bit*msqqa
      else
        bit=zip
      endif

      xmsq=xmsq+bit

      enddo
      enddo

! SNAP

      xmsq_central = xmsq

      if (central .and. doMultiTaucut) then
          scetreweight(:) = 0._dp
          if (xmsq_central /= 0._dp) then
              do m=1,size(tcutarray)
                  if (dynamictau) then
                      tauc = getdynamictau(p, tcutarray(m))
                  else
                      tauc = tcutarray(m)
                  endif

! SNIP FROM ABOVE
      xmsq=zip
      do j=-nf,nf
      do k=-nf,nf
      
      if ((j .ne. 0) .and. (k .ne. 0) &
       .and. (j .ne. -k)) cycle
      
      if (((j == 0) .and. (k == 0)) &
       .or. (j*k .ne. 0)) then !  gg -> g or qa -> g => jet is a gluon
        if ((j == 0) .and. (k == 0)) then
          bit=assemblejet(order,tauc, &
           beama0(j),beamb0(k),beama1(j,:),beamb1(k,:), &
           beama2(j,:),beamb2(k,:),soft1gg,soft2gg,jet1g,jet2g,Hgg)
        else
          bit=assemblejet(order,tauc, &
           beama0(j),beamb0(k),beama1(j,:),beamb1(k,:), &
           beama2(j,:),beamb2(k,:),soft1qa,soft2qa,jet1g,jet2g,Hqqb)
        endif
      else                     ! qg -> q  => jet is a quark
        if (j == 0) then
          bit=assemblejet(order,tauc, &
           beama0(j),beamb0(k),beama1(j,:),beamb1(k,:), &
           beama2(j,:),beamb2(k,:),soft1gq,soft2gq,jet1q,jet2q,Hgqb)
        elseif (k == 0) then
          bit=assemblejet(order,tauc, &
           beama0(j),beamb0(k),beama1(j,:),beamb1(k,:), &
           beama2(j,:),beamb2(k,:),soft1qg,soft2qg,jet1q,jet2q,Hqg)
        endif
      endif
      
      if     ((j == 0) .and. (k == 0)) then
        bit=bit*msqgg
      elseif (j == 0) then
        bit=bit*msqgq
      elseif (k == 0) then
        bit=bit*msqqg
      elseif (j == -k) then
        bit=bit*msqqa
      else
        bit=zip
      endif

      xmsq=xmsq+bit

      enddo
      enddo
! SNAP

                  scetreweight(m) = xmsq
              enddo
              scetreweight(:) = scetreweight(:) / xmsq_central
          endif
      endif

      end subroutine

    subroutine Hampgggsq(order,p1,p2,p3,p4,p5,msq0,hard1,hard2, &
    msq0qq,hard1qq,hard2qq)
!--- returns msq0, hard1, hard2 which can be used to reconstruct the
!--- complete hard function for ab -> H+c using:
!---  H = msq0 * (1 + [as/2/pi]*hard1 + [as/2/pi]^2*hard2)
      use constants
    implicit none
    include 'src/Inc/nf.f'
    include 'src/Inc/mxpart.f'
    include 'src/Inc/sprods_com.f'
    include 'src/Inc/hpls.f'
    include 'src/Inc/scale.f'
    include 'src/Inc/masses.f'
!     This routine calculates the amplitude for a right-handed quark
!     0--> q^+(1)+qb^1(2)+g^+(3)+l^-(4)+lb^+(5)
!     according to Eq.22 of 1309.3245v3
!     index is helicity of outgoing gluon line
    logical, parameter:: checkBBLM= .FALSE. ! use to check result in Eq. (73)
    real(dp):: uu,vv,qsq,xlf,hard1,hard2
    integer:: order,p1,p2,p3,p4,p5,q1,q2,q3,iperm,region,k
    complex(dp):: al,be,ga,de,alC,beC,gaC,deC,sq0,sq1
    complex(dp):: amp(0:2,2),ampqqbgll,C0,lnrat
    real(dp):: hnppp,hnppm,hnqqg, &
    s12,s13,s23,msq0,msq1,msq2,res0,res1,res2, &
    evolve1Lgg,evolve2Lgg,evolve1Lqqb,evolve2Lqqb, &
    c1,res0qq,res1qq,res2qq,msq0qq,msq1qq,msq2qq,hard1qq,hard2qq,Ct1,Ct2
    complex(dp):: alp0,bet0,gam0,alp1,bet1,gam1,alp2,bet2,gam2
    real(dp):: alR(0:2),alI(0:2), beR(0:2),beI(0:2)
    real(dp):: gaR(0:2),gaI(0:2)
    real(dp):: LgamHgg0,LgamHgg1,LgamHqq0,LgamHqq1,Lrat

    real(dp), parameter :: zip = 0._dp
    real(dp), parameter :: three = 3._dp, four = 4._dp, eight = 8._dp
    real(dp), parameter :: one = 1._dp, two = 2._dp, six = 6._dp
    real(dp), parameter :: be0 = 11/three*CA-4/three*TR*NF
    real(dp), parameter :: be1= 34/three*CA**2-(20/three*CA+4*CF)*TR*NF

    real(dp), parameter:: &
    UGamHgg0= 6*CA, &
    UGamHgg1= &
    - 40/3._dp*CA*TR*nf &
    + 134/3._dp*CA**2 &
    - 2*CA**2*pisq, &
    UGamHqq0= + 2*CA + 4*CF, &
    UGamHqq1= &
    - 40/9._dp*CA*TR*nf &
    + 134/9._dp*CA**2 &
    - 2/3._dp*CA**2*pisq &
    - 80/9._dp*CF*TR*nf &
    + 268/9._dp*CF*CA &
    - 4/3._dp*CF*CA*pisq
!--- statement functions for hn of Eqs.(50), (51)
    hnppp(s12,s23,s13,qsq,alp0,alp1) &
    =qsq**4/(two*s12*s23*s13)*real(alp0*conjg(alp1),dp)
    hnppm(s12,s23,s13,qsq,bet0,bet1) &
    =s12**3/(two*s23*s13)*real(bet0*conjg(bet1),dp)
!--- statement function for qqg hn of Eq.(54)
    hnqqg(s12,s23,s13,gam0,gam1) &
    =s23**2/(two*s12)*real(gam0*conjg(gam1),dp)
!--- statement function for one-loop evolution
    evolve1Lgg(musq)= &
    -UGamHgg0*log(abs(s(p1,p2)/musq))**2/two &
    -LgamHgg0*log(abs(s(p1,p2)/musq))
!--- statement function for two-loop evolution
    evolve2Lgg(musq,c1)= &
    + log(abs(s(p1,p2)/musq))**4/eight * ( &
    + UGamHgg0**2 &
    ) &
    + log(abs(s(p1,p2)/musq))**3/six * ( &
    + 3*UGamHgg0*LgamHgg0 &
    + be0*UGamHgg0 &
    ) &
    + log(abs(s(p1,p2)/musq))**2/two * ( &
    + LgamHgg0**2 &
    - UGamHgg1 &
    - c1*UGamHgg0 &
    + be0*LgamHgg0 &
    ) &
    + log(abs(s(p1,p2)/musq)) * ( &
    - LgamHgg1 &
    - c1*LgamHgg0 &
    - be0*c1 &
    )
!--- statement function for one-loop evolution
    evolve1Lqqb(musq)= &
    -UGamHqq0*log(abs(s(p1,p2)/musq))**2/two &
    -LgamHqq0*log(abs(s(p1,p2)/musq))
!--- statement function for two-loop evolution
    evolve2Lqqb(musq,c1)= &
    + log(abs(s(p1,p2)/musq))**4/eight * ( &
    + UGamHqq0**2 &
    ) &
    + log(abs(s(p1,p2)/musq))**3/six * ( &
    + 3*UGamHqq0*LgamHqq0 &
    + be0*UGamHqq0 &
    ) &
    + log(abs(s(p1,p2)/musq))**2/two * ( &
    + LgamHqq0**2 &
    - UGamHqq1 &
    - c1*UGamHqq0 &
    + be0*LgamHqq0 &
    ) &
    + log(abs(s(p1,p2)/musq)) * ( &
    - LgamHqq1 &
    - c1*LgamHqq0 &
    - be0*c1 &
    )

    xlf=real(nf,dp)
          
    if (checkBBLM) then
    !---- ! Check Becher et al. result at a single PS point
        if ((p1 == 2) .AND. (p2 == 1)) then
            s(p1,p2)=1._dp
            s(p1,p3)=-0.4_dp
            s(p4,p5)=0.1_dp**2
            s(p2,p3)=s(p4,p5)-s(p1,p2)-s(p1,p3)
            s(p2,p1)=s(p1,p2)
            s(p3,p1)=s(p1,p3)
            s(p3,p2)=s(p2,p3)
        endif
        if ((p1 == 5) .AND. (p2 == 1)) then
            s(p2,p3)=1._dp
            s(p1,p3)=-0.4_dp
            s(p4,p5)=0.1_dp**2
            s(p1,p2)=s(p4,p5)-s(p1,p3)-s(p2,p3)
            s(p2,p1)=s(p1,p2)
            s(p3,p1)=s(p1,p3)
            s(p3,p2)=s(p2,p3)
        endif
        if ((p1 == 2) .AND. (p2 == 5)) then
            s(p1,p3)=1._dp
            s(p2,p3)=-0.4_dp
            s(p4,p5)=0.1_dp**2
            s(p1,p2)=s(p4,p5)-s(p1,p3)-s(p2,p3)
            s(p2,p1)=s(p1,p2)
            s(p3,p1)=s(p1,p3)
            s(p3,p2)=s(p2,p3)
        endif
        scale=0.6_dp
        musq=scale**2
    endif
          
!      write(6,*) 'region',region
!      write(6,*) 's(p1,p2)',s(p1,p2)
!      write(6,*) 's(p1,p3)',s(p1,p3)
!      write(6,*) 's(p2,p3)',s(p2,p3)

    qsq=s(p1,p2)+s(p1,p3)+s(p2,p3)

!--- logarithm appearing in evolution expressions
!--- note: take absolute value in order to work for all crossings
    Lrat=log(abs(s(p1,p2)**2/(s(p1,p3)*s(p2,p3))))
!--- constants for evolution components
    LgamHgg0 = - 2*CA*Lrat
    LgamHgg1 = &
    + 256/9._dp*CA*TR*nf &
    - 2/3._dp*CA*TR*nf*pisq &
    + 40/9._dp*CA*Lrat*TR*nf &
    - 692/9._dp*CA**2 &
    + 6*CA**2*zeta3 &
    + 11/6._dp*CA**2*pisq &
    - 134/9._dp*CA**2*Lrat &
    + 2/3._dp*CA**2*Lrat*pisq &
    + 12*CF*TR*nf &
    + 3*be1 &
    + 2*be1 ! Ct correction term, cf == (B.5)

!--- note: different from functions in W/Z routines due to additional
!--- beta-function factors arising from more powers of alpha-s at LO
    LgamHqq0 = &
    - 2*CA*Lrat &
    - 6*CF &
    + 2*be0
    LgamHqq1 = &
    + 256/27._dp*CA*TR*nf &
    - 2/9._dp*CA*TR*nf*pisq &
    + 40/9._dp*CA*Lrat*TR*nf &
    - 692/27._dp*CA**2 &
    + 2*CA**2*zeta3 &
    + 11/18._dp*CA**2*pisq &
    - 134/9._dp*CA**2*Lrat &
    + 2/3._dp*CA**2*Lrat*pisq &
    + 368/27._dp*CF*TR*nf &
    + 4/3._dp*CF*TR*nf*pisq &
    - 961/27._dp*CF*CA &
    + 52*CF*CA*zeta3 &
    - 11/3._dp*CF*CA*pisq &
    - 3*CF**2 &
    - 48*CF**2*zeta3 &
    + 4*CF**2*pisq &
    + 3*be1 &
    + 2*be1 ! Ct correction term, cf == (B.5)
        
    res0=zip
    res1=zip
    res2=zip
    res0qq=zip
    res1qq=zip
    res2qq=zip
    do iperm=1,4
        if     (iperm == 1) then
            q1=p1
            q2=p2
            q3=p3
        elseif (iperm == 2) then
            q1=p3
            q2=p2
            q3=p1
        elseif (iperm == 3) then
            q1=p3
            q2=p1
            q3=p2
        elseif (iperm == 4) then
            q1=p2
            q2=p1
            q3=p3
        endif

    !--- determine correct region from invariants
        if (s(q1,q2) > 0) then
            region=2
        elseif (s(q1,q3) > 0) then
            region=3
        elseif (s(q2,q3) > 0) then
            region=4
        else
            write(6,*) 'region not found:'
            write(6,*) 's12,s13,s23',s(q1,q2),s(q1,q3),s(q2,q3)
            stop
        endif

        if (region == 2) then
            uu=-s(q1,q3)/s(q1,q2)
            vv=+qsq/s(q1,q2)
        elseif (region == 3) then
            uu=-s(q2,q3)/s(q1,q3)
            vv=+qsq/s(q1,q3)
        elseif (region == 4) then
            uu=-s(q1,q3)/s(q2,q3)
            vv=+qsq/s(q2,q3)
        endif
    !      write(6,*) 'iperm,region',iperm,region
              
        if (order > 0) then
        !--- fill arrays for 2DHPLs
            call tdhpl(uu,vv,2*order,G1,G2,G3,G4,H1,H2,H3,H4)
        else
        !--- fill with dummy values
            G1=zip
            G2=zip
            G3=zip
            G4=zip
            H1=zip
            H2=zip
            H3=zip
            H4=zip
        endif

    !--- Coefficients for Region 2
        if (region == 2) then
        !--- region 2: relevant for gg->Hg (alpha,beta) and qa->Hg (gamma)
            alR(0)=Halpha_2a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
            alI(0)=Halpha_2a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
            beR(0)=Hbeta_2a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
            beI(0)=Hbeta_2a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
            gaR(0)=Hgamma_2a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
            gaI(0)=Hgamma_2a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
            if (order >= 1) then
                alR(1)=Halpha_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
                alI(1)=Halpha_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
                beR(1)=Hbeta_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
                beI(1)=Hbeta_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
                gaR(1)=Hgamma_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
                gaI(1)=Hgamma_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
            endif
            if (order >= 2) then
                alR(2)=Halpha_2a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
                alI(2)=Halpha_2a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
                beR(2)=Hbeta_2a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
                beI(2)=Hbeta_2a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
                gaR(2)=Hgamma_2a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
                gaI(2)=Hgamma_2a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
            endif

        !--- Coefficients for Region 3
        elseif (region == 3) then
        !--- region 3: relevant for qq->Hg (gamma)
            gaR(0)=Hgamma_3a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
            gaI(0)=Hgamma_3a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
            if (order >= 1) then
                gaR(1)=Hgamma_3a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
                gaI(1)=Hgamma_3a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
            endif
            if (order >= 2) then
                gaR(2)=Hgamma_3a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
                gaI(2)=Hgamma_3a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
            endif

        !--- Coefficients for Region 4
        elseif (region == 4) then
        !--- region 4: relevant for gg->Hg (beta) and qq->Hg (gamma)
            beR(0)=Hbeta_4a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
            beI(0)=Hbeta_4a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
            gaR(0)=Hgamma_4a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
            gaI(0)=Hgamma_4a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
            if (order >= 1) then
                beR(1)=Hbeta_4a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
                beI(1)=Hbeta_4a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
                gaR(1)=Hgamma_4a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
                gaI(1)=Hgamma_4a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
            endif
            if (order >= 2) then
                beR(2)=Hbeta_4a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
                beI(2)=Hbeta_4a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
                gaR(2)=Hgamma_4a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
                gaI(2)=Hgamma_4a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
            endif
                  
        else
            write(6,*) 'Hampgggsq: region should be 2,3 or 4: ',region
            stop

        endif

        alp0=cmplx(alR(0),alI(0),dp)
        bet0=cmplx(beR(0),beI(0),dp)
        gam0=cmplx(gaR(0),gaI(0),dp)
        alp1=cmplx(alR(1),alI(1),dp)
        bet1=cmplx(beR(1),beI(1),dp)
        gam1=cmplx(gaR(1),gaI(1),dp)
        alp2=cmplx(alR(2),alI(2),dp)
        bet2=cmplx(beR(2),beI(2),dp)
        gam2=cmplx(gaR(2),gaI(2),dp)
              
        if (iperm == 1) then
            res0=res0 &
            +two*hnppp(s(q1,q2),s(q2,q3),s(q1,q3),qsq,alp0,alp0)
            res1=res1 &
            +two*hnppp(s(q1,q2),s(q2,q3),s(q1,q3),qsq,alp0,alp1)
            res2=res2 &
            +two*hnppp(s(q1,q2),s(q2,q3),s(q1,q3),qsq,alp0,alp2) &
            +hnppp(s(q1,q2),s(q2,q3),s(q1,q3),qsq,alp1,alp1)
        endif
        if (iperm <= 3) then
            res0=res0 &
            +two*hnppm(s(q1,q2),s(q2,q3),s(q1,q3),qsq,bet0,bet0)
            res1=res1 &
            +two*hnppm(s(q1,q2),s(q2,q3),s(q1,q3),qsq,bet0,bet1)
            res2=res2 &
            +two*hnppm(s(q1,q2),s(q2,q3),s(q1,q3),qsq,bet0,bet2) &
            +hnppm(s(q1,q2),s(q2,q3),s(q1,q3),qsq,bet1,bet1)
        endif
              
        if ((iperm == 1) .OR. (iperm == 4)) then
        !--- note interchange of q1 and q2 in this call wrt. above,
        !--- per footnote 3 of 1309.3245
        !--- seems to be necessary for qqb~ channel but not for qg
            if ((p1 == 2) .AND. (p2 == 1)) then
                res0qq=res0qq &
                +two*hnqqg(s(q1,q2),s(q1,q3),s(q2,q3),gam0,gam0)
                res1qq=res1qq &
                +two*hnqqg(s(q1,q2),s(q1,q3),s(q2,q3),gam0,gam1)
                res2qq=res2qq &
                +two*hnqqg(s(q1,q2),s(q1,q3),s(q2,q3),gam0,gam2) &
                +hnqqg(s(q1,q2),s(q1,q3),s(q2,q3),gam1,gam1)
            else
                res0qq=res0qq &
                +two*hnqqg(s(q1,q2),s(q2,q3),s(q1,q3),gam0,gam0)
                res1qq=res1qq &
                +two*hnqqg(s(q1,q2),s(q2,q3),s(q1,q3),gam0,gam1)
                res2qq=res2qq &
                +two*hnqqg(s(q1,q2),s(q2,q3),s(q1,q3),gam0,gam2) &
                +hnqqg(s(q1,q2),s(q2,q3),s(q1,q3),gam1,gam1)
            endif
        endif
              
    !      write(6,*) 'alp0,alp1',alp0,alp1
    !      write(6,*) 'bet0,bet1',bet0,bet1
    !      write(6,*) 'gam0,gam1',gam0,gam1
                  
    enddo ! end of iperm
          
    msq0=res0
    msq0qq=res0qq
    if (order > 0) then
        msq1=two*res1
    !--- RG running to restore scale dependence from musq=qsq
        msq1=msq1 &
        +(evolve1Lgg(musq)-evolve1Lgg(qsq))*res0
    !        write(6,*) 'gg 1L/LO/fourpi',msq1/msq0/fourpi
        msq1qq=two*res1qq
    !--- RG running to restore scale dependence from musq=qsq
        msq1qq=msq1qq &
        +(evolve1Lqqb(musq)-evolve1Lqqb(qsq))*res0qq
    !        write(6,*) 'qq 1L/LO/fourpi',msq1qq/msq0qq/fourpi
        if (order == 1) then
        !--- if just doing NLO calculation, add correction from
        !--- Wilson coefficient |Ct|^2 here, in units of (as/4/pi)
            msq1=msq1+(5*CA-3*CF)*res0*two
            msq1qq=msq1qq+(5*CA-3*CF)*res0qq*two
        endif
    else
        msq1=zip
        msq1qq=zip
    endif
    if (order > 1) then
        msq2=two*res2
    !--- RG running to restore scale dependence from musq=qsq
        c1=two*res1/res0-evolve1Lgg(qsq)
        msq2=msq2 &
        +(evolve2Lgg(musq,c1) &
        -evolve2Lgg(qsq,c1))*res0
        if (checkBBLM) then
            write(6,*) 'gg 2L/LO/fourpi**2',msq2/msq0/fourpi**2
            pause
        endif
        msq2qq=two*res2qq
    !--- RG running to restore scale dependence from musq=qsq
        c1=two*res1qq/res0qq-evolve1Lqqb(qsq)
        msq2qq=msq2qq &
        +(evolve2Lqqb(musq,c1) &
        -evolve2Lqqb(qsq,c1))*res0qq
    !        write(6,*) 'qq 2L/LO/fourpi**2',msq2qq/msq0qq/fourpi**2
    !        pause
    else
        msq2=zip
        msq2qq=zip
    endif
          
!--- convert from coefficient of [as/4/pi]^n to coefficient of [as/2/pi]^n
    msq1=msq1/two
    msq1qq=msq1qq/two
    msq2=msq2/four
    msq2qq=msq2qq/four

!--- normalize coefficients to LO
    hard1=msq1/msq0
    hard1qq=msq1qq/msq0qq
    hard2=msq2/msq0
    hard2qq=msq2qq/msq0qq
          
!      hard1=hard1-pisq/twelve*3._dp*CA
!      hard1qq=hard1qq-pisq/twelve*(2._dp*CF+CA)
          
    return
    end subroutine Hampgggsq

!*
! elicity_amplitudes_H._dpcpp

! reatedon:Dec3,2013
! uthor:lorentzen
!/
          
          
!#      include "peter/helicity_amplitudes_H._dphpp"
          
          
! eter::math_constants
          
          
!**********************************************************************
! egion2a:g+g-->H+g*
!*********************************************************************/
          
    function Halpha_2a0re(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=1
    return
    end function Halpha_2a0re
          
          
    function Halpha_2a0im(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=0
    return
    end function Halpha_2a0im
          
          
    function Halpha_2a1re(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    use constants
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=-3*G2(0,2)+6*G2(1,2)-3*G2(2,0)+6*G2(3,0)-6*G1(3)*H1(0)+ &
    G1(1)*(-6*H1(0)-6*H1(1))-3*H1(0)*H1(1)+G1(0)*(9*H1(0)+3*H1(1))+ &
    G1(2)*(9*H1(0)+3*H1(1))-6*H2(0,0)+ &
    (4*omv*(nf+3*u)+4*nf*u*(-1+u+v)+ &
    &  3*(-4+4*v-4*(u*u)+114*zeta2*(v*v)))/(12._dp*(v*v))- &
    (3*(G1(0)*G1(0)))/2._dp-(3*(G1(2)*G1(2)))/2._dp-(9*(H1(0)*H1(0)))/2._dp- &
    (3*(H1(1)*H1(1)))/2._dp
    return
    end function Halpha_2a1re
          
          
    function Halpha_2a1im(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    use constants
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=-3*pi*G1(0)+6*pi*G1(1)-3*pi*G1(2)+6*pi*G1(3)+3*pi*H1(0)+3*pi*H1(1)
    return
    end function Halpha_2a1im
          
          
    function Halpha_2a2re(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    use constants
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=108*zeta2*G2(3,2)+(66-4*nf)*G3(2,0,2)+(66-4*nf)*G3(2,2,0)+ &
    &  36*G4(0,0,1,2)+18*G4(0,0,2,2)-36*G4(0,0,3,0)+36*G4(0,1,0,2)+ &
    &  36*G4(0,1,2,0)-36*G4(0,1,2,2)+18*G4(0,2,0,2)-18*G4(0,2,1,2)+ &
    &  18*G4(0,2,2,0)-18*G4(0,2,3,0)-18*G4(0,3,0,2)-18*G4(0,3,2,0)+ &
    &  72*G4(0,3,3,0)-36*G4(1,0,2,2)-72*G4(1,1,1,2)+72*G4(1,1,2,2)- &
    &  36*G4(1,2,0,2)+72*G4(1,2,1,2)-36*G4(1,2,2,0)+18*G4(2,0,0,2)- &
    &  18*G4(2,0,1,2)+18*G4(2,0,2,0)-18*G4(2,0,3,0)-18*G4(2,1,0,2)+ &
    &  72*G4(2,1,1,2)-18*G4(2,1,2,0)+18*G4(2,2,0,0)-36*G4(2,2,1,2)+ &
    &  36*G4(2,2,3,0)-36*G4(2,3,0,0)+36*G4(2,3,0,2)+36*G4(2,3,2,0)- &
    &  36*G4(3,0,0,2)-36*G4(3,0,2,0)+72*G4(3,0,3,0)-36*G4(3,2,0,0)+ &
    &  72*G4(3,3,0,0)-72*G4(3,3,3,0)+G3(3,3,0)*(-66+4*nf-72*H1(0))+ &
    G3(0,2,2)*(66-4*nf-36*H1(0))+36*G3(0,0,3)*H1(0)+18*G3(0,2,3)*H1(0)+ &
    &  18*G3(0,3,2)*H1(0)-72*G3(0,3,3)*H1(0)+18*G3(2,0,3)*H1(0)- &
    &  36*G3(2,2,3)*H1(0)-36*G3(2,3,2)*H1(0)+18*G3(3,0,2)*H1(0)- &
    &  72*G3(3,0,3)*H1(0)+18*G3(3,2,0)*H1(0)+72*G3(3,3,3)*H1(0)+ &
    G3(2,1,2)*(5*(-33+2*nf)+36*H1(0))+G3(1,2,2)*(-132+8*nf+72*H1(0))+ &
    G3(1,2,1)*(-72*H1(0)-72*H1(1))+G3(2,1,1)*(-72*H1(0)-72*H1(1))+ &
    G3(1,1,2)*(-66+4*nf-72*H1(0)-72*H1(1))+ &
    G3(0,0,1)*(-36*H1(0)-36*H1(1))+G3(0,1,0)*(-36*H1(0)-36*H1(1))+ &
    G3(2,3,0)*(18*H1(0)-36*H1(1))+G3(0,0,2)*(66-4*nf-18*H1(1))+ &
    G3(0,2,0)*(66-4*nf-18*H1(1))+ &
    G3(2,0,0)*(66-4*nf-36*H1(0)-18*H1(1))+G3(0,2,1)*(18*H1(0)+18*H1(1))+ &
    G3(2,0,1)*(18*H1(0)+18*H1(1))+G3(2,1,0)*(18*H1(0)+18*H1(1))+ &
    G3(0,3,0)*(5*(-33+2*nf)+36*H1(0)+18*H1(1))+ &
    G3(0,1,2)*(18*H1(0)+36*H1(1))+G3(1,0,2)*(18*H1(0)+36*H1(1))+ &
    G3(1,2,0)*(18*H1(0)+36*H1(1))+G3(2,2,1)*(36*H1(0)+36*H1(1))+ &
    G3(3,0,0)*(-132+8*nf+72*H1(0)+36*H1(1))+ &
    G3(1,1,1)*(72*H1(0)+72*H1(1))+ &
    G2(3,3)*(-144*zeta2+(66-4*nf)*H1(0)+72*H2(0,0))+ &
    G2(2,1)*(-36*zeta2+(165-10*nf)*H1(0)+(165-10*nf)*H1(1)+36*H2(0,1)- &
    &  36*H2(1,0))+G2(0,3)*(36*zeta2+(165-10*nf)*H1(0)-72*H2(0,0)- &
    &  18*H2(0,1)-18*H2(1,0))+G2(2,3)*(-108*zeta2+36*H2(0,1)+36*H2(1,0))+ &
    G2(0,1)*(-108*zeta2-18*H2(0,1)-18*H2(1,0)-36*H2(1,1))+ &
    G2(1,0)*(108*zeta2-18*H2(0,1)-18*H2(1,0)-36*H2(1,1))+ &
    G2(1,1)*(-72*zeta2+(66-4*nf)*H1(0)+(66-4*nf)*H1(1)+72*H2(1,0)+ &
    &  72*H2(1,1))+(-363+22*nf)*H3(0,0,0)+(33-2*nf)*H3(1,0,1)+ &
    (-66+4*nf)*H3(1,1,0)+108*H4(0,0,0,0)+ &
    H2(1,1)*(363._dp/4._dp-11*nf-54*zeta2+(nf*nf)/3._dp)+ &
    G2(2,2)*(363._dp/4._dp-11*nf+18*zeta2+(-132+8*nf)*H1(0)-36*H2(0,1)+ &
    &  36*H2(1,0)+(nf*nf)/3._dp)+ &
    G2(0,0)*(363._dp/4._dp-11*nf-54*zeta2+(-132+8*nf)*H1(0)+ &
    (-66+4*nf)*H1(1)+72*H2(0,0)+36*H2(0,1)+36*H2(1,0)+18*H2(1,1)+ &
    (nf*nf)/3._dp)+H2(0,1)*(605._dp/4._dp-(55*nf)/3._dp-108*zeta2+ &
    (5*(nf*nf))/9._dp)+H1(1)*((33-2*nf)*H2(0,0)+ &
    ((270-100*nf)*(omv*omv)+ &
    u*(-2*u*(nf*nf)+nf*(150-233*v-50*(-3+4*v)*(u*u)- &
    &  100*(u*u*u)+u*(-61+233*v-100*(v*v))+102*(v*v))+ &
    &  3*(-135+324*v+45*(-3+4*v)*(u*u)+90*(u*u*u)+ &
    (-382+27*zeta3)*(v*v)+u*(57-324*v+90*(v*v)))))/ &
    (9._dp*u*(v*v)))+(H2(1,0)*(omv*omv* &
    (5445-660*nf-3888*zeta2+20*(nf*nf))*(v*v)+ &
    &  2*omv*u*(-20*(nf*nf)*(v*v)+ &
    &  4*nf*(-2+3*v+200*(u*u*u*u)+288*(v*v)-62*(v*v*v))+ &
    &  9*(-24+36*v-240*(u*u*u*u)+(-917+432*zeta2)*(v*v)+ &
    &  150*(v*v*v)))+u*u* &
    (20*(nf*nf)*(v*v)-9*(-120+1272*v-120*(u*u*u*u)+ &
    (-3197+432*zeta2)*(v*v)-12*(u*u)*(62-145*v+60*(v*v))+ &
    &  1380*(v*v*v)-12*u*(-48+201*v-195*(v*v)+40*(v*v*v))- &
    &  120*(v*v*v*v))-4*nf* &
    (-10-632*v+100*(u*u*u*u)+1503*(v*v)+ &
    &  2*(u*u)*(299-662*v+300*(v*v))-772*(v*v*v)+ &
    &  2*u*(-196+783*v-786*(v*v)+200*(v*v*v))+100*(v*v*v*v))))) &
    /(36._dp*(v*v)*(w*w))+(H2(0,0)* &
    (-3*(omv*omv)*(-72+72*v+(-3837+4212*zeta2)*(v*v)+ &
    &  4*nf*(6-6*v+145*(v*v)))+ &
    &  2*omv*u*(80*(-27+10*nf)*(u*u*u*u)+ &
    &  9*(-36+360*v+(-1603+1404*zeta2)*(v*v)+150*(v*v*v)))+ &
    &  9*(u*u)*(360-1944*v+120*(u*u*u*u)+(3943-1404*zeta2)*(v*v)+ &
    &  12*(u*u)*(64-145*v+60*(v*v))-1380*(v*v*v)+ &
    &  12*u*(-50+207*v-195*(v*v)+40*(v*v*v))+120*(v*v*v*v))- &
    &  4*nf*u*(100*(u*u*u*u*u)+4*(u*u*u)*(154-331*v+150*(v*v))+ &
    &  2*(u*u)*(-225+810*v-786*(v*v)+200*(v*v*v))- &
    &  2*(27-186*v+726*(v*v)-629*(v*v*v)+62*(v*v*v*v))+ &
    u*(270-1014*v+1827*(v*v)-772*(v*v*v)+100*(v*v*v*v)))+ &
    &  60*(nf*nf)*(v*v)*(w*w)))/(36._dp*(v*v)*(w*w))+ &
    (-360*zeta2*(omv*omv)*(342-342*v-8241*(v*v)+ &
    &  2*nf*(-57+57*v+205*(v*v)))+ &
    &  120*u*zeta2*(-9*(480*(u*u*u*u*u)+ &
    &  30*(u*u*u)*(103-232*v+96*(v*v))+ &
    &  6*(u*u)*(-421+1665*v-1560*(v*v)+320*(v*v*v))+ &
    u*(1296-7146*v+7963*(v*v)-5520*(v*v*v)+480*(v*v*v*v))- &
    &  2*(207-1269*v-380*(v*v)+842*(v*v*v)+600*(v*v*v*v)))+ &
    &  2*nf*(-501+2397*v+800*(u*u*u*u*u)-2805*(v*v)+ &
    u*u*u*(4955-10592*v+4800*(v*v))+1901*(v*v*v)+ &
    u*u*(-3669+13041*v-12576*(v*v)+3200*(v*v*v))- &
    &  992*(v*v*v*v)+u*(1644-7533*v+10602*(v*v)-6176*(v*v*v)+ &
    &  800*(v*v*v*v))))+2753190*zeta4*(v*v)*(w*w)- &
    &  80*(20*u*(-1+u+v)*(nf*nf)- &
    &  3*nf*(554-554*v+554*(u*u)-(857+330*zeta3)*(v*v))+ &
    &  9*(494-494*v+494*(u*u)-(1321+495*zeta3)*(v*v)))*(w*w)- &
    &  160*omv*(480*(-27+10*nf)*zeta2*(u*u*u*u*u)+ &
    (-2223*u+831*nf*u+10*(nf*nf))*(w*w)))/(4320._dp*(v*v)*(w*w))+ &
    ((33-2*nf)/6._dp-(27*H1(0))/2._dp-(9*H1(1))/2._dp)*(G1(0)*G1(0)*G1(0))+ &
    (9*(G1(0)*G1(0)*G1(0)*G1(0)))/8._dp+ &
    ((33-2*nf)/6._dp-(27*H1(0))/2._dp-(9*H1(1))/2._dp)*(G1(2)*G1(2)*G1(2))+ &
    (9*(G1(2)*G1(2)*G1(2)*G1(2)))/8._dp+ &
    (-33._dp/2._dp+nf+(27*H1(1))/2._dp)*(H1(0)*H1(0)*H1(0))+ &
    (81*(H1(0)*H1(0)*H1(0)*H1(0)))/8._dp+ &
    G2(1,2)*(H1(0)*(33-2*nf-18*H1(1))+(99-6*nf)*H1(1)-36*H2(0,0)- &
    &  72*H2(1,0)+(20*(-27+10*nf)*(omv*omv)- &
    &  40*omv*(-27+10*nf)*(u*u*u*u*u)+ &
    u*(4*nf*(50*(u*u*u*u*u)+u*u*(-7+70*v-62*(v*v))+ &
    u*(58-70*v-30*(v*v))+2*(u*u*u)*(29-81*v+25*(v*v))- &
    &  2*(50-81*v+31*(v*v)))+ &
    &  9*(-60*(u*u*u*u*u)+30*(4-9*v+5*(v*v))- &
    &  6*(u*u*u)*(14-45*v+10*(v*v))+ &
    &  6*(u*u)*(6-29*v+25*(v*v))+ &
    u*(-84+174*v+(268+342*zeta2)*(v*v)))))/(18._dp*(u*u)*(v*v)) &
    -27*(H1(0)*H1(0))-9*(H1(1)*H1(1)))+ &
    G2(3,0)*(H1(0)*(33-2*nf-18*H1(1))+(-33+2*nf)*H1(1)-108*H2(0,0)- &
    &  18*H2(0,1)-18*H2(1,0)+(-3*(omv*omv)* &
    (4*nf*(-3+3*v+10*(v*v))- &
    &  3*(-12+12*v+(268+486*zeta2)*(v*v)))- &
    &  2*omv*u*(40*(-27+10*nf)*(u*u*u*u)+ &
    &  9*(-18+180*v+(106+486*zeta2)*(v*v)+75*(v*v*v)))+ &
    u*(-9*u*(180-972*v+60*(u*u*u*u)+(1064-486*zeta2)*(v*v)+ &
    &  6*(u*u)*(64-145*v+60*(v*v))-690*(v*v*v)+ &
    &  6*u*(-50+207*v-195*(v*v)+40*(v*v*v))+60*(v*v*v*v))+ &
    &  4*nf*(-27+186*v+50*(u*u*u*u*u)-231*(v*v)+ &
    u*u*u*(308-662*v+300*(v*v))+134*(v*v*v)+ &
    u*u*(-225+810*v-786*(v*v)+200*(v*v*v))-62*(v*v*v*v)+ &
    u*(135-507*v+666*(v*v)-386*(v*v*v)+50*(v*v*v*v)))))/ &
    (18._dp*(v*v)*(w*w))-27*(H1(0)*H1(0))-9*(H1(1)*H1(1)))+ &
    (9*H2(0,0)+(-4*(nf*nf)*(v*v)+ &
    &  4*nf*(-3+3*omv*u+3*v-3*(u*u)+43*(v*v))- &
    &  3*(-12+12*omv*u+12*v-12*(u*u)+631*(v*v)+270*zeta2*(v*v)))/ &
    (24._dp*(v*v)))*(H1(1)*H1(1))+ &
    G1(0)*G1(0)*((9*G2(0,2))/2._dp-9*G2(1,2)+(9*G2(2,0))/2._dp-9*G2(3,0)+ &
    &  9*G1(3)*H1(0)+G1(2)*((-27*H1(0))/2._dp-(9*H1(1))/2._dp)+(33._dp/2._dp-nf)*H1(1)+ &
    G1(1)*(9*H1(0)+9*H1(1))+H1(0)*(33._dp/2._dp-nf+(27*H1(1))/2._dp)+9*H2(0,0)+ &
    (-4*(nf*nf)*(v*v)+4*nf*(-3+3*omv*u+3*v-3*(u*u)+43*(v*v))- &
    &  3*(-12+12*omv*u+12*v-12*(u*u)+631*(v*v)+270*zeta2*(v*v)))/ &
    (24._dp*(v*v))+(9*(G1(2)*G1(2)))/4._dp+(117*(H1(0)*H1(0)))/4._dp+ &
    (9*(H1(1)*H1(1)))/4._dp)+G2(2,0)* &
    ((-99._dp/2._dp+3*nf)*H1(1)+H1(0)*(-66+4*nf+9*H1(1))+18*H2(0,0)- &
    &  18*H2(0,1)-18*H2(1,0)+(4*(nf*nf)*(v*v)- &
    &  27*(-4+4*v+80*omv*(u*u*u)-40*(u*u*u*u)+ &
    (49+162*zeta2)*(v*v)-4*(u*u)*(13-45*v+10*(v*v))+ &
    &  4*u*(3-28*v+25*(v*v)))- &
    &  4*nf*(-200*omv*(u*u*u)+100*(u*u*u*u)+ &
    u*(-7+131*v-124*(v*v))+3*(3-3*v+v*v)+ &
    u*u*(107-324*v+100*(v*v))))/(36._dp*(v*v))+ &
    (27*(H1(0)*H1(0)))/2._dp+(9*(H1(1)*H1(1)))/2._dp)+ &
    G2(0,2)*((-99._dp/2._dp+3*nf)*H1(1)+H1(0)*(-66+4*nf+9*H1(1))+18*H2(0,0)+ &
    &  18*H2(0,1)+18*H2(1,0)+(4*(nf*nf)*(v*v)- &
    &  27*(-4+4*v+80*omv*(u*u*u)-40*(u*u*u*u)+ &
    (49+162*zeta2)*(v*v)-4*(u*u)*(13-45*v+10*(v*v))+ &
    &  4*u*(3-28*v+25*(v*v)))- &
    &  4*nf*(-200*omv*(u*u*u)+100*(u*u*u*u)+ &
    u*(-7+131*v-124*(v*v))+3*(3-3*v+v*v)+ &
    u*u*(107-324*v+100*(v*v))))/(36._dp*(v*v))+ &
    (27*(H1(0)*H1(0)))/2._dp+(9*(H1(1)*H1(1)))/2._dp)+ &
    G1(2)*G1(2)*((9*G2(0,2))/2._dp-9*G2(1,2)+(9*G2(2,0))/2._dp-9*G2(3,0)+ &
    &  9*G1(3)*H1(0)+(-33._dp/2._dp+nf)*H1(1)+H1(0)*(33._dp/2._dp-nf+(63*H1(1))/2._dp)+ &
    &  9*H2(0,0)+(-4*(nf*nf)*(v*v)+ &
    &  4*nf*(-3+3*omv*u+3*v-3*(u*u)+43*(v*v))- &
    &  3*(-12+12*omv*u+12*v-12*(u*u)+631*(v*v)+270*zeta2*(v*v)))/ &
    (24._dp*(v*v))+(117*(H1(0)*H1(0)))/4._dp+(27*(H1(1)*H1(1)))/4._dp)+ &
    H1(0)*H1(0)*((-33._dp/2._dp+nf)*H1(1)+27*H2(0,0)+ &
    (-20*(nf*nf)*(v*v)+12*nf*(-3+3*omv*u+3*v-3*(u*u)+65*(v*v))+ &
    &  27*(4-4*omv*u-4*v+4*(u*u)-3*(97+6*zeta2)*(v*v)))/ &
    (24._dp*(v*v))+(45*(H1(1)*H1(1)))/4._dp)+ &
    G1(3)*(10*(-33+2*nf)*zeta2-108*zeta2*H1(1)+(-132+8*nf)*H2(0,0)+ &
    &  108*H3(0,0,0)+(99-6*nf+18*H1(1))*(H1(0)*H1(0))+ &
    &  27*(H1(0)*H1(0)*H1(0))+H1(0)* &
    ((33-2*nf)*H1(1)+(3*(omv*omv)* &
    (4*nf*(-3+3*v+10*(v*v))- &
    &  3*(-12+12*v+(268+486*zeta2)*(v*v)))+ &
    &  2*omv*u*(40*(-27+10*nf)*(u*u*u*u)+ &
    &  9*(-18+180*v+(106+486*zeta2)*(v*v)+75*(v*v*v)))+ &
    u*(9*u*(180-972*v+60*(u*u*u*u)+(1064-486*zeta2)*(v*v)+ &
    &  6*(u*u)*(64-145*v+60*(v*v))-690*(v*v*v)+ &
    &  6*u*(-50+207*v-195*(v*v)+40*(v*v*v))+60*(v*v*v*v)) &
    -4*nf*(-27+186*v+50*(u*u*u*u*u)-231*(v*v)+ &
    u*u*u*(308-662*v+300*(v*v))+134*(v*v*v)+ &
    u*u*(-225+810*v-786*(v*v)+200*(v*v*v))- &
    &  62*(v*v*v*v)+ &
    u*(135-507*v+666*(v*v)-386*(v*v*v)+50*(v*v*v*v)))))/ &
    (18._dp*(v*v)*(w*w))+9*(H1(1)*H1(1))))+ &
    G1(2)*(G2(0,2)*(-33._dp/2._dp+nf-9*H1(0)-9*H1(1))+ &
    G2(2,0)*(-33._dp/2._dp+nf-9*H1(0)-9*H1(1))+ &
    G2(1,2)*(33-2*nf+18*H1(0)+18*H1(1))+ &
    G2(3,0)*(33-2*nf+18*H1(0)+18*H1(1))+(165-10*nf)*H2(0,0)+ &
    (-33+2*nf)*H2(0,1)+(132-8*nf)*H2(1,0)-36*H3(0,0,0)-72*H3(1,1,0)+ &
    H1(1)*(-18*H2(0,0)+(-36+36*omv*u+36*v-36*(u*u)+804*(v*v)+ &
    &  594*zeta2*(v*v)-4*nf*(-3+3*omv*u+3*v-3*(u*u)+10*(v*v))) &
    /(12._dp*(v*v)))+(10*(-27+10*nf)*(omv*omv)+ &
    u*(2*u*(nf*nf)+3*(135-324*v-45*(-3+4*v)*(u*u)-90*(u*u*u)+ &
    u*(-57+324*v-90*(v*v))+ &
    (382+198*zeta2+189*zeta3)*(v*v))+ &
    nf*(-150+233*v+50*(-3+4*v)*(u*u)+100*(u*u*u)- &
    &  6*(17+6*zeta2)*(v*v)+u*(61-233*v+100*(v*v)))))/ &
    (9._dp*u*(v*v))+G1(3)*(H1(0)*(-33+2*nf-18*H1(1))-18*(H1(0)*H1(0)))+ &
    (-165._dp/2._dp+5*nf-(81*H1(1))/2._dp)*(H1(0)*H1(0))- &
    (81*(H1(0)*H1(0)*H1(0)))/2._dp+ &
    H1(0)*(-18*H2(0,0)+(40*(-27+10*nf)*(omv*omv)- &
    &  80*omv*(-27+10*nf)*(u*u*u*u*u)+ &
    u*(4*nf*(100*(u*u*u*u*u)+u*u*(-23+149*v-124*(v*v))+ &
    u*(125-149*v-90*(v*v))-4*(50-81*v+31*(v*v))+ &
    u*u*u*(125-324*v+100*(v*v)))+ &
    &  27*(-40*(u*u*u*u*u)-20*(u*u*u)*(3-9*v+2*(v*v))+ &
    &  20*(4-9*v+5*(v*v))+4*(u*u)*(7-30*v+25*(v*v))+ &
    u*(-60+120*v+(268+150*zeta2)*(v*v)))))/ &
    (36._dp*(u*u)*(v*v))-(45*(H1(1)*H1(1)))/2._dp)+ &
    (33._dp/2._dp-nf)*(H1(1)*H1(1))-(9*(H1(1)*H1(1)*H1(1)))/2._dp)+ &
    G1(0)*(G2(0,2)*(-33._dp/2._dp+nf-9*H1(0))+G2(2,0)*(-33._dp/2._dp+nf-9*H1(0))+ &
    G2(1,2)*(33-2*nf+18*H1(0))+G2(3,0)*(33-2*nf+18*H1(0))+ &
    (330-20*nf)*H2(0,0)+(231._dp/2._dp-7*nf)*H2(0,1)+ &
    (231._dp/2._dp-7*nf)*H2(1,0)+(66-4*nf)*H2(1,1)-108*H3(0,0,0)- &
    &  18*H3(0,0,1)-18*H3(0,1,0)-18*H3(1,0,0)+ &
    (H1(1)*(9*(-12+12*v+240*omv*(u*u*u)-120*(u*u*u*u)+ &
    (268+486*zeta2)*(v*v)-12*(u*u)*(13-45*v+10*(v*v))+ &
    &  12*u*(3-28*v+25*(v*v)))+ &
    &  4*nf*(9-9*v-200*omv*(u*u*u)+100*(u*u*u*u)+ &
    u*(-7+131*v-124*(v*v))-30*(v*v)+ &
    u*u*(107-324*v+100*(v*v)))))/(36._dp*(v*v))+ &
    (10*(-27+10*nf)*(-5+4*v)*(u*u*u)+20*(-27+10*nf)*(u*u*u*u)+ &
    &  2*(u*u)*(-9*(64-153*v+30*(v*v))+ &
    nf*(211-383*v+100*(v*v)))+ &
    omv*(2*nf*(-39+39*v+(-19+54*zeta2)*(v*v))- &
    &  3*(-66+66*v+(-386+594*zeta2+54*zeta3)*(v*v)))+ &
    u*(-2*nf*(-78-55*v+3*(38+18*zeta2)*(v*v))+ &
    &  3*(-132-336*v+(82+594*zeta2+54*zeta3)*(v*v)))+ &
    &  4*(nf*nf)*(w*w))/(18._dp*w*(v*v))+ &
    ((-27*H1(0))/2._dp-(9*H1(1))/2._dp)*(G1(2)*G1(2))+ &
    G1(3)*((-33+2*nf)*H1(0)-18*(H1(0)*H1(0)))+ &
    G1(1)*(H1(0)*(-33+2*nf-18*H1(1))+(-33+2*nf)*H1(1)- &
    &  18*(H1(0)*H1(0)))+(-165._dp/2._dp+5*nf-(81*H1(1))/2._dp)*(H1(0)*H1(0))- &
    (81*(H1(0)*H1(0)*H1(0)))/2._dp+ &
    H1(0)*((-165._dp/2._dp+5*nf)*H1(1)-18*H2(0,0)+ &
    (-9*(omv*omv)*(4*nf*(-3+3*v+10*(v*v))- &
    &  3*(-12+12*v+(268+342*zeta2)*(v*v)))- &
    &  2*omv*u*(80*(-27+10*nf)*(u*u*u*u)+ &
    &  27*(-18+128*v+(158+342*zeta2)*(v*v)+50*(v*v*v)))+ &
    u*(-27*u*(136-676*v+40*(u*u*u*u)+(632-342*zeta2)*(v*v)+ &
    &  20*(u*u)*(13-29*v+12*(v*v))-460*(v*v*v)+ &
    &  4*u*(-53+210*v-195*(v*v)+40*(v*v*v))+40*(v*v*v*v)) &
    +4*nf*(-81+435*v+100*(u*u*u*u*u)-447*(v*v)+ &
    u*u*u*(625-1324*v+600*(v*v))+217*(v*v*v)+ &
    u*u*(-477+1647*v-1572*(v*v)+400*(v*v*v))- &
    &  124*(v*v*v*v)+ &
    u*(306-1077*v+1329*(v*v)-772*(v*v*v)+100*(v*v*v*v))) &
    ))/(36._dp*(v*v)*(w*w))-(45*(H1(1)*H1(1)))/2._dp)+ &
    (-33._dp/2._dp+nf)*(H1(1)*H1(1))+ &
    G1(2)*((33._dp/2._dp-nf)*H1(1)+H1(0)*(66-4*nf+36*H1(1))- &
    ((33-2*nf)*(33-2*nf))/36._dp+45*(H1(0)*H1(0))+9*(H1(1)*H1(1)))- &
    (9*(H1(1)*H1(1)*H1(1)))/2._dp)+((-33+2*nf)*(H1(1)*H1(1)*H1(1)))/6._dp+ &
    H1(0)*((99-6*nf)*H2(0,0)+H1(1)* &
    (18*H2(0,0)+(-20*(nf*nf)*(v*v)+ &
    &  12*nf*(-3+3*omv*u+3*v-3*(u*u)+65*(v*v))+ &
    &  27*(4-4*omv*u-4*v+4*(u*u)-3*(97+6*zeta2)*(v*v)))/ &
    (36._dp*(v*v)))+(-2*(omv*omv)*(135*v-50*nf*v+u*(nf*nf))+ &
    omv*u*(-4*(nf*nf)*(u*u)+ &
    nf*(39-239*v+(78-200*v)*(u*u)-2*(-70+54*zeta2)*(v*v))+ &
    &  9*(-11+71*v+(-22+60*v)*(u*u)+ &
    (-256+198*zeta2+27*zeta3)*(v*v)))+ &
    u*u*(2*(nf*nf)*(2-3*v+u*u+v*v)- &
    &  9*(-22+123*v+(-11+30*v)*(u*u)+ &
    &  3*(-109+66*zeta2+9*zeta3)*(v*v)+30*(v*v*v))+ &
    nf*(-78+417*v+(-39+100*v)*(u*u)+ &
    (-379+108*zeta2)*(v*v)+100*(v*v*v))))/(9._dp*u*w*(v*v))+ &
    (-33._dp/2._dp+nf)*(H1(1)*H1(1))+(9*(H1(1)*H1(1)*H1(1)))/2._dp)+ &
    G1(1)*(6*(-44+(8*nf)/3._dp)*zeta2-72*zeta3+6*(-33+2*nf)*H2(0,0)+ &
    &  6*(-33+2*nf)*H2(0,1)+(-132+8*nf)*H2(1,0)+(-132+8*nf)*H2(1,1)+ &
    &  36*H3(0,0,0)+36*H3(0,0,1)+36*H3(0,1,0)+36*H3(1,0,0)+72*H3(1,1,0)+ &
    (H1(1)*((540-200*nf)*(omv*omv)+40*omv*(-27+10*nf)*(u*u*u*u*u)+ &
    u*(-4*nf*(50*(u*u*u*u*u)+u*u*(-7+70*v-62*(v*v))+ &
    u*(58-70*v-30*(v*v))+ &
    &  2*(u*u*u)*(29-81*v+25*(v*v))-2*(50-81*v+31*(v*v))) &
    -9*(-60*(u*u*u*u*u)+30*(4-9*v+5*(v*v))- &
    &  6*(u*u*u)*(14-45*v+10*(v*v))+ &
    &  6*(u*u)*(6-29*v+25*(v*v))+ &
    u*(-84+174*v+(268+342*zeta2)*(v*v))))))/ &
    (18._dp*(u*u)*(v*v))+(9*H1(0)+9*H1(1))*(G1(2)*G1(2))+ &
    (99-6*nf+45*H1(1))*(H1(0)*H1(0))+27*(H1(0)*H1(0)*H1(0))+ &
    G1(2)*(H1(0)*(-33+2*nf-36*H1(1))+(-33+2*nf)*H1(1)- &
    &  18*(H1(0)*H1(0))-18*(H1(1)*H1(1)))+(33-2*nf)*(H1(1)*H1(1))+ &
    H1(0)*((132-8*nf)*H1(1)+ &
    ((540-200*nf)*(omv*omv)+40*omv*(-27+10*nf)*(u*u*u*u*u)+ &
    u*(-4*nf*(50*(u*u*u*u*u)+u*u*(-7+70*v-62*(v*v))+ &
    u*(58-70*v-30*(v*v))+ &
    &  2*(u*u*u)*(29-81*v+25*(v*v))-2*(50-81*v+31*(v*v))) &
    -9*(-60*(u*u*u*u*u)+30*(4-9*v+5*(v*v))- &
    &  6*(u*u*u)*(14-45*v+10*(v*v))+ &
    &  6*(u*u)*(6-29*v+25*(v*v))+ &
    u*(-84+174*v+(268+198*zeta2)*(v*v)))))/ &
    (18._dp*(u*u)*(v*v))+27*(H1(1)*H1(1)))+9*(H1(1)*H1(1)*H1(1)))+ &
    (9*(H1(1)*H1(1)*H1(1)*H1(1)))/8._dp
    return
    end function Halpha_2a2re
          
          
    function Halpha_2a2im(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    use constants
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=36*pi*G3(0,0,1)+18*pi*G3(0,0,2)-36*pi*G3(0,0,3)+36*pi*G3(0,1,0)+ &
    &  18*pi*G3(0,2,0)-18*pi*G3(0,2,1)+18*pi*G3(0,2,2)-18*pi*G3(0,2,3)- &
    &  18*pi*G3(0,3,0)-18*pi*G3(0,3,2)+72*pi*G3(0,3,3)-36*pi*G3(1,0,2)- &
    &  72*pi*G3(1,1,1)+72*pi*G3(1,1,2)-36*pi*G3(1,2,0)+72*pi*G3(1,2,1)- &
    &  36*pi*G3(1,2,2)+18*pi*G3(2,0,0)-18*pi*G3(2,0,1)+18*pi*G3(2,0,2)- &
    &  18*pi*G3(2,0,3)-18*pi*G3(2,1,0)+72*pi*G3(2,1,1)-18*pi*G3(2,1,2)+ &
    &  18*pi*G3(2,2,0)-36*pi*G3(2,2,1)+36*pi*G3(2,2,3)+36*pi*G3(2,3,2)- &
    &  36*pi*G3(3,0,0)-36*pi*G3(3,0,2)+72*pi*G3(3,0,3)-36*pi*G3(3,2,0)+ &
    &  72*pi*G3(3,3,0)-72*pi*G3(3,3,3)-18*pi*G2(0,1)*H1(0)+18*pi*G2(3,2)*H1(0)+ &
    G2(3,3)*(2*(-33+2*nf)*pi-72*pi*H1(0))+ &
    G2(2,2)*((66-4*nf)*pi+36*pi*H1(0))+36*pi*G2(1,2)*H1(1)+ &
    G2(1,1)*(2*(-33+2*nf)*pi-72*pi*H1(1))+ &
    G2(2,3)*(-18*pi*H1(0)-36*pi*H1(1))+ &
    G2(0,0)*((66-4*nf)*pi-36*pi*H1(0)-18*pi*H1(1))+ &
    G2(0,2)*((66-4*nf)*pi-9*pi*H1(0)-18*pi*H1(1))+ &
    G2(2,0)*((66-4*nf)*pi-9*pi*H1(0)-18*pi*H1(1))+ &
    G2(2,1)*(5*(-33+2*nf)*pi-18*pi*H1(0)+18*pi*H1(1))+ &
    G2(0,3)*(5*(-33+2*nf)*pi+54*pi*H1(0)+18*pi*H1(1))+ &
    G2(1,0)*(18*pi*H1(0)+36*pi*H1(1))+G2(3,0)*(72*pi*H1(0)+36*pi*H1(1))+ &
    &  7*(33-2*nf)*pi*H2(0,0)+(7*(33-2*nf)*pi*H2(0,1))/2._dp+ &
    (5*(33-2*nf)*pi*H2(1,0))/2._dp+(66-4*nf)*pi*H2(1,1)-108*pi*H3(0,0,0)- &
    &  18*pi*H3(0,0,1)-18*pi*H3(0,1,0)-18*pi*H3(1,0,0)+ &
    (pi*H1(1)*(9*(-12+12*v+240*omv*(u*u*u)-120*(u*u*u*u)+ &
    (268+270*zeta2)*(v*v)-12*(u*u)*(13-45*v+10*(v*v))+ &
    &  12*u*(3-28*v+25*(v*v)))+ &
    &  4*nf*(9-9*v-200*omv*(u*u*u)+100*(u*u*u*u)+ &
    u*(-7+131*v-124*(v*v))-30*(v*v)+ &
    u*u*(107-324*v+100*(v*v)))))/(36._dp*(v*v))+ &
    (pi*(10*(-27+10*nf)*(omv*omv*omv)+ &
    omv*u*(675-1242*v+(906+594*zeta2+81*zeta3)*(v*v)- &
    nf*(250-333*v+2*(89+18*zeta2)*(v*v)))+ &
    u*u*(nf*(211-355*v-(39+100*v)*(u*u)+ &
    u*(78+122*v-200*(v*v))+(339+36*zeta2)*(v*v)- &
    &  100*(v*v*v))-3*(192-438*v-3*(11+30*v)*(u*u)+ &
    u*(66+114*v-180*(v*v))+ &
    (449+198*zeta2+27*zeta3)*(v*v)-90*(v*v*v))+ &
    &  2*(nf*nf)*(w*w))))/(9._dp*u*w*(v*v))+ &
    ((-33._dp/2._dp+nf)*pi-9*pi*G1(1)+(9*pi*G1(2))/2._dp-9*pi*G1(3)- &
    (27*pi*H1(0))/2._dp-(9*pi*H1(1))/2._dp)*(G1(0)*G1(0))+ &
    (9*pi*(G1(0)*G1(0)*G1(0)))/2._dp+ &
    ((-33._dp/2._dp+nf)*pi-9*pi*G1(3)-(27*pi*H1(0))/2._dp-(27*pi*H1(1))/2._dp)* &
    (G1(2)*G1(2))+(9*pi*(G1(2)*G1(2)*G1(2)))/2._dp+ &
    ((7*(-33+2*nf)*pi)/2._dp-(45*pi*H1(1))/2._dp)*(H1(0)*H1(0))- &
    (27*pi*(H1(0)*H1(0)*H1(0)))/2._dp+(-33._dp/2._dp+nf)*pi*(H1(1)*H1(1))+ &
    H1(0)*((5*(-33+2*nf)*pi*H1(1))/2._dp+18*pi*H2(0,0)+ &
    (pi*(-3*(omv*omv)*(4*nf*(-3+3*v+10*(v*v))- &
    &  3*(-12+12*v+(268+270*zeta2)*(v*v)))- &
    &  2*omv*u*(80*(-27+10*nf)*(u*u*u*u)+ &
    &  9*(-18+336*v+5*(-10+54*zeta2)*(v*v)+150*(v*v*v)))+ &
    u*(-9*u*(312-1860*v+120*(u*u*u*u)- &
    &  5*(-472+54*zeta2)*(v*v)+ &
    &  12*(u*u)*(63-145*v+60*(v*v))-1380*(v*v*v)+ &
    &  12*u*(-47+204*v-195*(v*v)+40*(v*v*v))+120*(v*v*v*v)) &
    +4*nf*(-27+309*v+100*(u*u*u*u*u)-477*(v*v)+ &
    u*u*u*(607-1324*v+600*(v*v))+319*(v*v*v)+ &
    u*u*(-423+1593*v-1572*(v*v)+400*(v*v*v))- &
    &  124*(v*v*v*v)+ &
    u*(234-951*v+1335*(v*v)-772*(v*v*v)+100*(v*v*v*v))))) &
    )/(36._dp*(v*v)*(w*w))-(27*pi*(H1(1)*H1(1)))/2._dp)+ &
    G1(3)*((-33+2*nf)*pi*H1(1)+H1(0)*(3*(-33+2*nf)*pi-18*pi*H1(1))- &
    &  108*pi*H2(0,0)-18*pi*H2(0,1)-18*pi*H2(1,0)+ &
    (pi*(-3*(omv*omv)*(4*nf*(-3+3*v+10*(v*v))- &
    &  3*(-12+12*v+(268+270*zeta2)*(v*v)))- &
    &  2*omv*u*(40*(-27+10*nf)*(u*u*u*u)+ &
    &  9*(-18+180*v+(106+270*zeta2)*(v*v)+75*(v*v*v)))+ &
    u*(-9*u*(180-972*v+60*(u*u*u*u)+(1064-270*zeta2)*(v*v)+ &
    &  6*(u*u)*(64-145*v+60*(v*v))-690*(v*v*v)+ &
    &  6*u*(-50+207*v-195*(v*v)+40*(v*v*v))+60*(v*v*v*v))+ &
    &  4*nf*(-27+186*v+50*(u*u*u*u*u)-231*(v*v)+ &
    u*u*u*(308-662*v+300*(v*v))+134*(v*v*v)+ &
    u*u*(-225+810*v-786*(v*v)+200*(v*v*v))- &
    &  62*(v*v*v*v)+ &
    u*(135-507*v+666*(v*v)-386*(v*v*v)+50*(v*v*v*v))))))/ &
    (18._dp*(v*v)*(w*w))-9*pi*(H1(0)*H1(0))-9*pi*(H1(1)*H1(1)))+ &
    G1(1)*((-33+2*nf)*pi*H1(0)+(-33+2*nf)*pi*H1(1)+ &
    G1(2)*((33-2*nf)*pi+18*pi*H1(0)+18*pi*H1(1))-36*pi*H2(0,0)- &
    &  18*pi*H2(0,1)-18*pi*H2(1,0)-36*pi*H2(1,1)+ &
    (pi*(20*(-27+10*nf)*(omv*omv)-40*omv*(-27+10*nf)*(u*u*u*u*u)+ &
    u*(4*nf*(50*(u*u*u*u*u)+u*u*(-7+70*v-62*(v*v))+ &
    u*(58-70*v-30*(v*v))+ &
    &  2*(u*u*u)*(29-81*v+25*(v*v))-2*(50-81*v+31*(v*v))) &
    +9*(-60*(u*u*u*u*u)+30*(4-9*v+5*(v*v))- &
    &  6*(u*u*u)*(14-45*v+10*(v*v))+ &
    &  6*(u*u)*(6-29*v+25*(v*v))+ &
    u*(-84+174*v+(268+270*zeta2)*(v*v))))))/ &
    (18._dp*(u*u)*(v*v))-9*pi*(G1(2)*G1(2))-9*pi*(H1(0)*H1(0))- &
    &  9*pi*(H1(1)*H1(1)))+G1(0)*(G1(1)*((33-2*nf)*pi+18*pi*H1(0))+ &
    G1(3)*((33-2*nf)*pi+18*pi*H1(0))+(-33+2*nf)*pi*H1(1)+ &
    G1(2)*((-33+2*nf)*pi-18*pi*H1(0)-9*pi*H1(1))+ &
    H1(0)*((-33+2*nf)*pi+9*pi*H1(1))+90*pi*H2(0,0)+36*pi*H2(0,1)+ &
    &  36*pi*H2(1,0)+18*pi*H2(1,1)- &
    (pi*(9*(-12+12*v+240*omv*(u*u*u)-120*(u*u*u*u)+ &
    (268+270*zeta2)*(v*v)-12*(u*u)*(13-45*v+10*(v*v))+ &
    &  12*u*(3-28*v+25*(v*v)))+ &
    &  4*nf*(9-9*v-200*omv*(u*u*u)+100*(u*u*u*u)+ &
    u*(-7+131*v-124*(v*v))-30*(v*v)+ &
    u*u*(107-324*v+100*(v*v)))))/(36._dp*(v*v))+ &
    (9*pi*(G1(2)*G1(2)))/2._dp+(9*pi*(H1(0)*H1(0)))/2._dp+(9*pi*(H1(1)*H1(1)))/2._dp) &
    +G1(2)*((-33+2*nf)*pi*H1(1)+ &
    G1(3)*((33-2*nf)*pi+18*pi*H1(0)+18*pi*H1(1))+ &
    H1(0)*(4*(33-2*nf)*pi+27*pi*H1(1))+18*pi*H2(0,0)-18*pi*H2(0,1)- &
    &  18*pi*H2(1,0)-(pi*(9*(-12+12*v+240*omv*(u*u*u)-120*(u*u*u*u)+ &
    (268+270*zeta2)*(v*v)-12*(u*u)*(13-45*v+10*(v*v))+ &
    &  12*u*(3-28*v+25*(v*v)))+ &
    &  4*nf*(9-9*v-200*omv*(u*u*u)+100*(u*u*u*u)+ &
    u*(-7+131*v-124*(v*v))-30*(v*v)+ &
    u*u*(107-324*v+100*(v*v)))))/(36._dp*(v*v))+ &
    (9*pi*(H1(0)*H1(0)))/2._dp+(27*pi*(H1(1)*H1(1)))/2._dp)- &
    (9*pi*(H1(1)*H1(1)*H1(1)))/2._dp
    return
    end function Halpha_2a2im
          
          
    function Hbeta_2a0re(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=1
    return
    end function Hbeta_2a0re
          
          
    function Hbeta_2a0im(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=0
    return
    end function Hbeta_2a0im
          
          
    function Hbeta_2a1re(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    use constants
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=-((-3+nf)*u*w)/3._dp+(57*zeta2)/2._dp-3*G2(0,2)+6*G2(1,2)-3*G2(2,0)+ &
    &  6*G2(3,0)-6*G1(3)*H1(0)+G1(1)*(-6*H1(0)-6*H1(1))-3*H1(0)*H1(1)+ &
    G1(0)*(9*H1(0)+3*H1(1))+G1(2)*(9*H1(0)+3*H1(1))-6*H2(0,0)- &
    (3*(G1(0)*G1(0)))/2._dp-(3*(G1(2)*G1(2)))/2._dp-(9*(H1(0)*H1(0)))/2._dp- &
    (3*(H1(1)*H1(1)))/2._dp
    return
    end function Hbeta_2a1re
          
          
    function Hbeta_2a1im(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    use constants
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=-3*pi*G1(0)+6*pi*G1(1)-3*pi*G1(2)+6*pi*G1(3)+3*pi*H1(0)+3*pi*H1(1)
    return
    end function Hbeta_2a1im
          
          
    function Hbeta_2a2re(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    use constants
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=943._dp/6._dp+(337*u)/3._dp-33*v+(48*v)/omu+(48*v)/omw-(337*u*v)/3._dp- &
    (99*zeta3)/2._dp+72*u*zeta3+72*v*zeta3-36*u*v*zeta3+ &
    (10197*zeta4)/16._dp+(66-4*nf)*G3(2,0,2)+(66-4*nf)*G3(2,2,0)+ &
    &  36*G4(0,0,1,2)+18*G4(0,0,2,2)-36*G4(0,0,3,0)+36*G4(0,1,0,2)+ &
    &  36*G4(0,1,2,0)-36*G4(0,1,2,2)+18*G4(0,2,0,2)-18*G4(0,2,1,2)+ &
    &  18*G4(0,2,2,0)-18*G4(0,2,3,0)-18*G4(0,3,0,2)-18*G4(0,3,2,0)+ &
    &  72*G4(0,3,3,0)-36*G4(1,0,2,2)-72*G4(1,1,1,2)+72*G4(1,1,2,2)- &
    &  36*G4(1,2,0,2)+72*G4(1,2,1,2)-36*G4(1,2,2,0)+18*G4(2,0,0,2)- &
    &  18*G4(2,0,1,2)+18*G4(2,0,2,0)-18*G4(2,0,3,0)-18*G4(2,1,0,2)+ &
    &  72*G4(2,1,1,2)-18*G4(2,1,2,0)+18*G4(2,2,0,0)-36*G4(2,2,1,2)+ &
    &  36*G4(2,2,3,0)-36*G4(2,3,0,0)+36*G4(2,3,0,2)+36*G4(2,3,2,0)- &
    &  36*G4(3,0,0,2)-36*G4(3,0,2,0)+72*G4(3,0,3,0)-36*G4(3,2,0,0)+ &
    &  72*G4(3,3,0,0)-72*G4(3,3,3,0)+G3(0,2,2)*(66-4*nf-36*H1(0))+ &
    &  36*G3(0,0,3)*H1(0)+18*G3(0,2,3)*H1(0)+18*G3(0,3,2)*H1(0)- &
    &  72*G3(0,3,3)*H1(0)+18*G3(2,0,3)*H1(0)-36*G3(2,2,3)*H1(0)- &
    &  36*G3(2,3,2)*H1(0)-72*G3(3,0,3)*H1(0)+72*G3(3,3,3)*H1(0)+ &
    G3(1,2,2)*(-132+8*nf+72*H1(0))+G3(1,2,1)*(-72*H1(0)-72*H1(1))+ &
    G3(2,1,1)*(-72*H1(0)-72*H1(1))+G3(0,0,1)*(-36*H1(0)-36*H1(1))+ &
    G3(0,1,0)*(-36*H1(0)-36*H1(1))+G3(0,0,2)*(66-4*nf-18*H1(1))+ &
    G3(0,2,0)*(66-4*nf-18*H1(1))+ &
    G3(2,0,0)*(66-4*nf-36*H1(0)-18*H1(1))+G3(0,2,1)*(18*H1(0)+18*H1(1))+ &
    G3(2,0,1)*(18*H1(0)+18*H1(1))+G3(2,1,0)*(18*H1(0)+18*H1(1))+ &
    G3(2,2,1)*(36*H1(0)+36*H1(1))+ &
    G3(3,0,0)*(-132+8*nf+72*H1(0)+36*H1(1))+ &
    G3(1,1,1)*(72*H1(0)+72*H1(1))+(-66+4*nf)*H3(1,1,0)+108*H4(0,0,0,0)+ &
    H2(1,1)*(363._dp/4._dp-11*nf-54*zeta2+(nf*nf)/3._dp)+ &
    G2(2,2)*(363._dp/4._dp-11*nf+18*zeta2+(-132+8*nf)*H1(0)-36*H2(0,1)+ &
    &  36*H2(1,0)+(nf*nf)/3._dp)+ &
    G2(0,0)*(363._dp/4._dp-11*nf-54*zeta2+(-132+8*nf)*H1(0)+ &
    (-66+4*nf)*H1(1)+72*H2(0,0)+36*H2(0,1)+36*H2(1,0)+18*H2(1,1)+ &
    (nf*nf)/3._dp)+(10*u*w*(nf*nf))/27._dp-(337*(u*u))/3._dp-18*zeta3*(u*u)+ &
    &  36*v*zeta3*(u*u)+omu*H3(0,1,0)* &
    (-33+3*u-6*(u*u)+nf*(2-u+2*(u*u)))+ &
    G3(1,1,2)*(-72*H1(0)-72*H1(1)+ &
    &  2*u*(36-9*u+nf*(-3+3*u-2*(u*u))+6*(u*u)))+ &
    G2(0,1)*(-108*zeta2-18*H2(0,1)-18*H2(1,0)-36*H2(1,1)+ &
    u*H1(0)*(36-9*u+nf*(-3+3*u-2*(u*u))+6*(u*u))+ &
    u*H1(1)*(36-9*u+nf*(-3+3*u-2*(u*u))+6*(u*u)))+ &
    G2(1,0)*(108*zeta2-18*H2(0,1)-18*H2(1,0)-36*H2(1,1)+ &
    u*H1(0)*(36-9*u+nf*(-3+3*u-2*(u*u))+6*(u*u))+ &
    u*H1(1)*(36-9*u+nf*(-3+3*u-2*(u*u))+6*(u*u)))+ &
    G3(0,1,2)*(18*H1(0)+36*H1(1)+ &
    u*(-36+9*u-6*(u*u)+nf*(3-3*u+2*(u*u))))+ &
    G3(1,0,2)*(18*H1(0)+36*H1(1)+ &
    u*(-36+9*u-6*(u*u)+nf*(3-3*u+2*(u*u))))+ &
    G3(1,2,0)*(18*H1(0)+36*H1(1)+ &
    u*(-36+9*u-6*(u*u)+nf*(3-3*u+2*(u*u))))+ &
    G2(1,1)*(-72*zeta2+72*H2(1,0)+72*H2(1,1)+ &
    &  2*u*H1(0)*(-36+9*u-6*(u*u)+nf*(3-3*u+2*(u*u)))+ &
    &  2*u*H1(1)*(-36+9*u-6*(u*u)+nf*(3-3*u+2*(u*u))))+ &
    &  12*zeta3*(u*u*u)+G3(0,3,0)*(-165+36*u+36*H1(0)+18*H1(1)-9*(u*u)+ &
    nf*(10-3*u+3*(u*u)-2*(u*u*u))+6*(u*u*u))+ &
    G2(0,3)*(36*zeta2-72*H2(0,0)-18*H2(0,1)-18*H2(1,0)+ &
    H1(0)*(165-36*u+9*(u*u)-6*(u*u*u)+ &
    nf*(-10+3*u-3*(u*u)+2*(u*u*u))))-18*zeta3*(v*v)+ &
    &  36*u*zeta3*(v*v)+H2(0,1)*(605._dp/4._dp+6*v-12*u*v-108*zeta2+ &
    (5*(nf*nf))/9._dp-6*(v*v)+nf*(-55._dp/3._dp+(-2+4*u)*v+2*(v*v)))+ &
    w*H3(0,0,1)*(3*omw-3*(11+4*u*v+2*(u*u)+2*(v*v))+ &
    nf*(2-u-v+4*u*v+2*(u*u)+2*(v*v)))+ &
    w*H3(1,0,0)*(3*omw-3*(11+4*u*v+2*(u*u)+2*(v*v))+ &
    nf*(2-u-v+4*u*v+2*(u*u)+2*(v*v)))- &
    omw*H3(1,0,1)*(nf*(3-3*v+u*(-3+4*v)+2*(u*u)+2*(v*v))- &
    &  3*(12-3*v+u*(-3+4*v)+2*(u*u)+2*(v*v)))+ &
    G3(3,3,0)*(-72*H1(0)-2*w*(3*omw-3*(11+4*u*v+2*(u*u)+2*(v*v))+ &
    nf*(2-u-v+4*u*v+2*(u*u)+2*(v*v))))+ &
    G3(3,0,2)*(18*H1(0)+w*(3*omw-3*(11+4*u*v+2*(u*u)+2*(v*v))+ &
    nf*(2-u-v+4*u*v+2*(u*u)+2*(v*v))))+ &
    G3(3,2,0)*(18*H1(0)+w*(3*omw-3*(11+4*u*v+2*(u*u)+2*(v*v))+ &
    nf*(2-u-v+4*u*v+2*(u*u)+2*(v*v))))+ &
    G3(2,3,0)*(18*H1(0)-36*H1(1)+ &
    w*(3*omw-3*(11+4*u*v+2*(u*u)+2*(v*v))+ &
    nf*(2-u-v+4*u*v+2*(u*u)+2*(v*v))))+ &
    G2(3,2)*(108*zeta2-w*H1(0)*(3*omw-3*(11+4*u*v+2*(u*u)+2*(v*v))+ &
    nf*(2-u-v+4*u*v+2*(u*u)+2*(v*v))))+ &
    G2(2,3)*(-108*zeta2+36*H2(0,1)+36*H2(1,0)- &
    w*H1(0)*(3*omw-3*(11+4*u*v+2*(u*u)+2*(v*v))+ &
    nf*(2-u-v+4*u*v+2*(u*u)+2*(v*v))))+ &
    G2(3,3)*(-144*zeta2+72*H2(0,0)+ &
    &  2*w*H1(0)*(3*omw-3*(11+4*u*v+2*(u*u)+2*(v*v))+ &
    nf*(2-u-v+4*u*v+2*(u*u)+2*(v*v))))+ &
    H1(1)*((33-2*nf)*H2(0,0)+((270-100*nf)*(omv*omv)+ &
    u*(3*(-292+234*v+27*zeta3+135*(-7+4*v)*(u*u*u*u)+ &
    &  270*(u*u*u*u*u)+ &
    u*u*(-1093+1026*v+27*zeta3-351*(v*v))-135*(v*v)- &
    &  2*u*(-301+162*v+27*zeta3+9*(v*v))+ &
    &  18*(u*u*u)*(76-72*v+15*(v*v)))+ &
    nf*(254-385*v-150*(-7+4*v)*(u*u*u*u)-300*(u*u*u*u*u)+ &
    u*u*u*(-1417+1486*v-300*(v*v))+ &
    u*(-425+474*v-72*(v*v))+150*(v*v)+ &
    u*u*(938-1175*v+436*(v*v)))))/(9._dp*u*(omu*omu)))+ &
    &  12*zeta3*(v*v*v)+H3(0,0,0)*(nf* &
    (26+6*omv*(-1+u)*u-3*v-4*(u*u*u)+3*(v*v)-6*u*(v*v)- &
    &  2*(v*v*v))+3*(-143+12*v-6*omv*(u*u)+4*(u*u*u)-3*(v*v)+ &
    &  6*u*(4-v+v*v)+2*(v*v*v)))+ &
    G3(2,1,2)*(36*H1(0)-3*(44+12*v+(-3+6*v)*(u*u)+2*(u*u*u)- &
    &  3*(v*v)+6*u*(2-v+v*v)+2*(v*v*v))+ &
    nf*(8+3*v+(-3+6*v)*(u*u)+2*(u*u*u)-3*(v*v)+ &
    u*(3-6*v+6*(v*v))+2*(v*v*v)))+ &
    G2(2,1)*(-36*zeta2+36*H2(0,1)-36*H2(1,0)+ &
    H1(0)*(3*(44+12*v+(-3+6*v)*(u*u)+2*(u*u*u)-3*(v*v)+ &
    &  6*u*(2-v+v*v)+2*(v*v*v))- &
    nf*(8+3*v+(-3+6*v)*(u*u)+2*(u*u*u)-3*(v*v)+ &
    u*(3-6*v+6*(v*v))+2*(v*v*v)))+ &
    H1(1)*(3*(44+12*v+(-3+6*v)*(u*u)+2*(u*u*u)-3*(v*v)+ &
    &  6*u*(2-v+v*v)+2*(v*v*v))- &
    nf*(8+3*v+(-3+6*v)*(u*u)+2*(u*u*u)-3*(v*v)+ &
    u*(3-6*v+6*(v*v))+2*(v*v*v))))+ &
    (nf*(6*(-529+4*v*(32-57/omu-57/omw-27*zeta3)-186*zeta3+ &
    (754+108*zeta3-216*v*zeta3)*(u*u)-72*zeta3*(u*u*u)+ &
    &  108*zeta3*(v*v)-2*u* &
    (377+54*zeta3-v*(377+108*zeta3)+108*zeta3*(v*v))- &
    &  72*zeta3*(v*v*v))+6*zeta2* &
    (-4800*omv*(u*u*u)+2400*(u*u*u*u)+ &
    u*(-1091+763*v+256*(v*v))+u*u*(3491-4544*v+2400*(v*v))+ &
    (1478*omv*u-615*(omv*omv)-663*(u*u))/(w*w))))/108._dp+ &
    &  6*zeta2*(120*omv*(u*u*u)-60*(u*u*u*u)+ &
    u*u*(-451._dp/4._dp+104*v-60*(v*v))-(u*(-211+139*v+64*(v*v)))/4._dp+ &
    (-5794*omv*u+2747*(omv*omv)+2927*(u*u))/(24._dp*(w*w)))- &
    (H2(0,0)*(-480*omv*(-27+10*nf)*(u*u*u*u*u)+ &
    &  120*(-27+10*nf)*(u*u*u*u*u*u)- &
    &  3*(omv*omv)*(3837+72*v-72*(v*v)+ &
    &  4*nf*(-145-6*v+6*(v*v)))+ &
    &  4*(u*u*u*u)*(-27*(205-352*v+180*(v*v))+ &
    &  4*nf*(481-892*v+450*(v*v)))+ &
    &  2*u*(-1+v)*(4*nf*(435-72*v-6*(v*v)+16*(v*v*v))- &
    &  9*(1279-204*v+6*(v*v)+48*(v*v*v)))+ &
    &  12*(u*u*u)*(4*nf*(-131+320*v-292*(v*v)+100*(v*v*v))- &
    &  9*(-195+407*v-336*(v*v)+120*(v*v*v)))+ &
    &  3*(u*u)*(-7077+8496*v-7092*(v*v)+3456*(v*v*v)- &
    &  1080*(v*v*v*v)+4*nf* &
    (361-586*v+630*(v*v)-368*(v*v*v)+100*(v*v*v*v)))+ &
    &  6*(2106*zeta2-10*(nf*nf))*(w*w)))/(36._dp*(w*w))- &
    (H2(1,0)*(165*(-33+4*nf)*(omv*omv)- &
    &  480*omv*(-27+10*nf)*(u*u*u*u*u)+ &
    &  120*(-27+10*nf)*(u*u*u*u*u*u)+ &
    &  12*(-1+v)*(u*u*u)*(-27*(63-72*v+40*(v*v))+ &
    nf*(506-768*v+400*(v*v)))+ &
    &  4*(u*u*u*u)*(-27*(203-352*v+180*(v*v))+ &
    &  2*nf*(953-1784*v+900*(v*v)))+ &
    &  2*u*(-1+v)*(4*nf*(280-90*v+21*(v*v)+16*(v*v*v))- &
    &  9*(893-228*v+42*(v*v)+48*(v*v*v)))+ &
    &  3*(u*u)*(-5919+8424*v-7236*(v*v)+3456*(v*v*v)- &
    &  1080*(v*v*v*v)+4*nf* &
    (269-580*v+642*(v*v)-368*(v*v*v)+100*(v*v*v*v)))+ &
    &  4*(972*zeta2-5*(nf*nf))*(w*w)))/(36._dp*(w*w))+ &
    ((33-2*nf)/6._dp-(27*H1(0))/2._dp-(9*H1(1))/2._dp)*(G1(0)*G1(0)*G1(0))+ &
    (9*(G1(0)*G1(0)*G1(0)*G1(0)))/8._dp+ &
    ((33-2*nf)/6._dp-(27*H1(0))/2._dp-(9*H1(1))/2._dp)*(G1(2)*G1(2)*G1(2))+ &
    (9*(G1(2)*G1(2)*G1(2)*G1(2)))/8._dp+ &
    (-33._dp/2._dp+nf+(27*H1(1))/2._dp)*(H1(0)*H1(0)*H1(0))+ &
    (81*(H1(0)*H1(0)*H1(0)*H1(0)))/8._dp+ &
    G2(1,2)*((99-6*nf)*H1(1)-36*H2(0,0)-72*H2(1,0)+ &
    H1(0)*(66+36*u-18*H1(1)-9*(u*u)+6*(u*u*u)- &
    nf*(4+3*u-3*(u*u)+2*(u*u*u)))+ &
    (20*(-27+10*nf)*(omv*omv)- &
    &  2*omv*u*(135-1620*(u*u*u*u)+nf*(76+600*(u*u*u*u)))+ &
    u*u*(4*nf*(-42-9*v+150*(u*u*u*u)+9*(v*v)+ &
    &  2*u*(-31+32*v+8*(v*v))+2*(u*u)*(106-142*v+75*(v*v))) &
    +9*(342*zeta2-2* &
    (-179-6*v+90*(u*u*u*u)+6*(v*v)+ &
    &  3*u*(-25+21*v+8*(v*v))+3*(u*u)*(55-52*v+30*(v*v))) &
    )))/(18._dp*(u*u))-27*(H1(0)*H1(0))-9*(H1(1)*H1(1)))+ &
    G2(3,0)*(-108*H2(0,0)-18*H2(0,1)-18*H2(1,0)+ &
    omw*H1(1)*(nf*(3-3*v+u*(-3+4*v)+2*(u*u)+2*(v*v))- &
    &  3*(12-3*v+u*(-3+4*v)+2*(u*u)+2*(v*v)))+ &
    H1(0)*(-18*H1(1)-3*(-33+12*v+(-3+6*v)*(u*u)+2*(u*u*u)- &
    &  3*(v*v)+6*u*(2-v+v*v)+2*(v*v*v))+ &
    nf*(-6+3*v+(-3+6*v)*(u*u)+2*(u*u*u)-3*(v*v)+ &
    u*(3-6*v+6*(v*v))+2*(v*v*v)))+ &
    (-240*omv*(-27+10*nf)*(u*u*u*u*u)+60*(-27+10*nf)*(u*u*u*u*u*u)- &
    &  12*(omv*omv)*(-201+9*v-9*(v*v)+nf*(10-3*v+3*(v*v)))+ &
    &  2*(u*u*u*u)*(-27*(205-352*v+180*(v*v))+ &
    &  4*nf*(481-892*v+450*(v*v)))+ &
    &  2*u*(-1+v)*(9*(268+102*v-3*(v*v)-24*(v*v*v))+ &
    &  4*nf*(-30-36*v-3*(v*v)+8*(v*v*v)))+ &
    &  6*(u*u*u)*(4*nf*(-131+320*v-292*(v*v)+100*(v*v*v))- &
    &  9*(-195+407*v-336*(v*v)+120*(v*v*v)))+ &
    &  6*(u*u)*(2*nf*(98-293*v+315*(v*v)-184*(v*v*v)+ &
    &  50*(v*v*v*v))-3* &
    (136-708*v+591*(v*v)-288*(v*v*v)+90*(v*v*v*v)))+ &
    &  4374*zeta2*(w*w))/(18._dp*(w*w))-27*(H1(0)*H1(0))-9*(H1(1)*H1(1)))+ &
    (9*H2(0,0)+(-4*(nf*nf)-3*(631+12*omv*u+270*zeta2-12*(u*u))+ &
    &  4*nf*(43+3*omv*u-3*(u*u)))/24._dp)*(H1(1)*H1(1))+ &
    G1(0)*G1(0)*((9*G2(0,2))/2._dp-9*G2(1,2)+(9*G2(2,0))/2._dp-9*G2(3,0)+ &
    &  9*G1(3)*H1(0)+G1(2)*((-27*H1(0))/2._dp-(9*H1(1))/2._dp)+(33._dp/2._dp-nf)*H1(1)+ &
    G1(1)*(9*H1(0)+9*H1(1))+H1(0)*(33._dp/2._dp-nf+(27*H1(1))/2._dp)+9*H2(0,0)+ &
    (-4*(nf*nf)-3*(631+12*omv*u+270*zeta2-12*(u*u))+ &
    &  4*nf*(43+3*omv*u-3*(u*u)))/24._dp+(9*(G1(2)*G1(2)))/4._dp+ &
    (117*(H1(0)*H1(0)))/4._dp+(9*(H1(1)*H1(1)))/4._dp)+ &
    G2(2,0)*(-147._dp/4._dp-72*u+48*u*v-(243*zeta2)/2._dp+ &
    (-99._dp/2._dp+3*nf)*H1(1)+H1(0)*(-66+4*nf+9*H1(1))+18*H2(0,0)- &
    &  18*H2(0,1)-18*H2(1,0)+(nf*nf)/9._dp+162*(u*u)-156*v*(u*u)- &
    &  180*(u*u*u)+180*v*(u*u*u)+90*(u*u*u*u)+24*u*(v*v)+ &
    &  90*(u*u)*(v*v)-(nf*(3-600*omv*(u*u*u)+300*(u*u*u*u)+ &
    u*(-115+83*v+32*(v*v))+u*u*(415-568*v+300*(v*v))))/9._dp+ &
    (27*(H1(0)*H1(0)))/2._dp+(9*(H1(1)*H1(1)))/2._dp)+ &
    G2(0,2)*(-147._dp/4._dp-72*u+48*u*v-(243*zeta2)/2._dp+ &
    (-99._dp/2._dp+3*nf)*H1(1)+H1(0)*(-66+4*nf+9*H1(1))+18*H2(0,0)+ &
    &  18*H2(0,1)+18*H2(1,0)+(nf*nf)/9._dp+162*(u*u)-156*v*(u*u)- &
    &  180*(u*u*u)+180*v*(u*u*u)+90*(u*u*u*u)+24*u*(v*v)+ &
    &  90*(u*u)*(v*v)-(nf*(3-600*omv*(u*u*u)+300*(u*u*u*u)+ &
    u*(-115+83*v+32*(v*v))+u*u*(415-568*v+300*(v*v))))/9._dp+ &
    (27*(H1(0)*H1(0)))/2._dp+(9*(H1(1)*H1(1)))/2._dp)+ &
    G1(2)*G1(2)*((9*G2(0,2))/2._dp-9*G2(1,2)+(9*G2(2,0))/2._dp-9*G2(3,0)+ &
    &  9*G1(3)*H1(0)+(-33._dp/2._dp+nf)*H1(1)+H1(0)*(33._dp/2._dp-nf+(63*H1(1))/2._dp)+ &
    &  9*H2(0,0)+(-4*(nf*nf)-3*(631+12*omv*u+270*zeta2-12*(u*u))+ &
    &  4*nf*(43+3*omv*u-3*(u*u)))/24._dp+(117*(H1(0)*H1(0)))/4._dp+ &
    (27*(H1(1)*H1(1)))/4._dp)+H1(0)*H1(0)* &
    ((-33._dp/2._dp+nf)*H1(1)+27*H2(0,0)+ &
    (-20*(nf*nf)-27*(291+4*u-4*u*v+18*zeta2-4*(u*u))+ &
    &  12*nf*(65+3*omv*u-3*(u*u)))/24._dp+(45*(H1(1)*H1(1)))/4._dp)+ &
    G1(3)*(-108*zeta2*H1(1)+108*H3(0,0,0)+ &
    w*H2(0,1)*(3*omw-3*(11+4*u*v+2*(u*u)+2*(v*v))+ &
    nf*(2-u-v+4*u*v+2*(u*u)+2*(v*v)))+ &
    w*H2(1,0)*(3*omw-3*(11+4*u*v+2*(u*u)+2*(v*v))+ &
    nf*(2-u-v+4*u*v+2*(u*u)+2*(v*v)))- &
    &  2*H2(0,0)*(-3*(-44+12*v+(-3+6*v)*(u*u)+2*(u*u*u)-3*(v*v)+ &
    &  6*u*(2-v+v*v)+2*(v*v*v))+ &
    nf*(-8+3*v+(-3+6*v)*(u*u)+2*(u*u*u)-3*(v*v)+ &
    u*(3-6*v+6*(v*v))+2*(v*v*v)))+ &
    &  4*zeta2*(-3*(22+12*v+(-3+6*v)*(u*u)+2*(u*u*u)-3*(v*v)+ &
    &  6*u*(2-v+v*v)+2*(v*v*v))+ &
    nf*(4+3*v+(-3+6*v)*(u*u)+2*(u*u*u)-3*(v*v)+ &
    u*(3-6*v+6*(v*v))+2*(v*v*v)))+ &
    (99-6*nf+18*H1(1))*(H1(0)*H1(0))+27*(H1(0)*H1(0)*H1(0))+ &
    H1(0)*((33-2*nf)*H1(1)-(-240*omv*(-27+10*nf)*(u*u*u*u*u)+ &
    &  60*(-27+10*nf)*(u*u*u*u*u*u)- &
    &  12*(omv*omv)*(-201+9*v-9*(v*v)+nf*(10-3*v+3*(v*v)))+ &
    &  2*(u*u*u*u)*(-27*(205-352*v+180*(v*v))+ &
    &  4*nf*(481-892*v+450*(v*v)))+ &
    &  2*u*(-1+v)*(9*(268+102*v-3*(v*v)-24*(v*v*v))+ &
    &  4*nf*(-30-36*v-3*(v*v)+8*(v*v*v)))+ &
    &  6*(u*u*u)*(4*nf*(-131+320*v-292*(v*v)+100*(v*v*v))- &
    &  9*(-195+407*v-336*(v*v)+120*(v*v*v)))+ &
    &  6*(u*u)*(2*nf*(98-293*v+315*(v*v)-184*(v*v*v)+ &
    &  50*(v*v*v*v))- &
    &  3*(136-708*v+591*(v*v)-288*(v*v*v)+90*(v*v*v*v)))+ &
    &  4374*zeta2*(w*w))/(18._dp*(w*w))+9*(H1(1)*H1(1))))+ &
    G1(2)*(G2(0,2)*(-33._dp/2._dp+nf-9*H1(0)-9*H1(1))+ &
    G2(2,0)*(-33._dp/2._dp+nf-9*H1(0)-9*H1(1))+ &
    G2(1,2)*(33-2*nf+18*H1(0)+18*H1(1))+ &
    G2(3,0)*(33-2*nf+18*H1(0)+18*H1(1))+(132-8*nf)*H2(1,0)- &
    &  36*H3(0,0,0)-72*H3(1,1,0)+ &
    H1(1)*(67+3*u-3*u*v+(99*zeta2)/2._dp-18*H2(0,0)-3*(u*u)+ &
    nf*(-10._dp/3._dp-omv*u+u*u))+ &
    omw*H2(0,1)*(nf*(3-3*v+u*(-3+4*v)+2*(u*u)+2*(v*v))- &
    &  3*(12-3*v+u*(-3+4*v)+2*(u*u)+2*(v*v)))+ &
    (10*(-27+10*nf)*(omv*omv)+ &
    u*(nf*(-254+385*v-36*zeta2*((-1+u)*(-1+u))+ &
    &  150*(-7+4*v)*(u*u*u*u)+300*(u*u*u*u*u)+ &
    u*u*(-938+1175*v-436*(v*v))-150*(v*v)+ &
    u*(425-474*v+72*(v*v))+u*u*u*(1417-1486*v+300*(v*v))) &
    +3*(292-234*v+189*zeta3+198*zeta2*((-1+u)*(-1+u))- &
    &  135*(-7+4*v)*(u*u*u*u)-270*(u*u*u*u*u)+135*(v*v)+ &
    &  2*u*(-301+162*v-189*zeta3+9*(v*v))- &
    &  18*(u*u*u)*(76-72*v+15*(v*v))+ &
    u*u*(1093-1026*v+189*zeta3+351*(v*v)))))/ &
    (9._dp*u*(omu*omu))+H2(0,0)* &
    (-3*(-66+12*v+(-3+6*v)*(u*u)+2*(u*u*u)-3*(v*v)+ &
    &  6*u*(2-v+v*v)+2*(v*v*v))+ &
    nf*(-12+3*v+(-3+6*v)*(u*u)+2*(u*u*u)-3*(v*v)+ &
    u*(3-6*v+6*(v*v))+2*(v*v*v)))+ &
    G1(3)*(H1(0)*(-33+2*nf-18*H1(1))-18*(H1(0)*H1(0)))+ &
    (-165._dp/2._dp+5*nf-(81*H1(1))/2._dp)*(H1(0)*H1(0))- &
    (81*(H1(0)*H1(0)*H1(0)))/2._dp+ &
    H1(0)*(-18*H2(0,0)+(40*(-27+10*nf)*(omv*omv)- &
    &  4*omv*u*(135-1620*(u*u*u*u)+nf*(76+600*(u*u*u*u)))+ &
    u*u*(4*nf*(-114+300*(u*u*u*u)+u*(-133+101*v+32*(v*v))+ &
    u*u*(433-568*v+300*(v*v)))+ &
    &  27*(150*zeta2-8* &
    (-41+15*(u*u*u*u)+u*(-13+9*v+4*(v*v))+ &
    u*u*(28-26*v+15*(v*v))))))/(36._dp*(u*u))- &
    (45*(H1(1)*H1(1)))/2._dp)+(33._dp/2._dp-nf)*(H1(1)*H1(1))- &
    (9*(H1(1)*H1(1)*H1(1)))/2._dp)+ &
    G1(0)*(G2(0,2)*(-33._dp/2._dp+nf-9*H1(0))+G2(2,0)*(-33._dp/2._dp+nf-9*H1(0))+ &
    G2(1,2)*(33-2*nf+18*H1(0))+G2(3,0)*(33-2*nf+18*H1(0))+ &
    (231._dp/2._dp-7*nf)*H2(1,0)+(66-4*nf)*H2(1,1)-108*H3(0,0,0)- &
    &  18*H3(0,0,1)-18*H3(0,1,0)-18*H3(1,0,0)+ &
    H2(0,0)*(330-36*u+9*(u*u)-6*(u*u*u)+ &
    nf*(-20+3*u-3*(u*u)+2*(u*u*u)))+ &
    H2(0,1)*(231._dp/2._dp-36*u+9*(u*u)-6*(u*u*u)+ &
    nf*(-7+3*u-3*(u*u)+2*(u*u*u)))+ &
    (-36*omw*(-30+nf)*v*w+12*(72-19*nf)*w*(v*v)+ &
    omw*omw*(3*(-180+36*v*(10+w*(-5+3*u+15*(u*u)))+ &
    w*(-82+576*u-54*zeta3-810*(u*u)+540*(u*u*u)+ &
    &  18*zeta2*(-33+24*u-6*(u*u)+4*(u*u*u)))-180*(v*v))+ &
    &  2*nf*(100-2*v*(100+w*(-50+7*u+150*(u*u)))+ &
    w*(13-217*u+450*(u*u)+ &
    &  6*zeta2*(9-9*u+9*(u*u)-6*(u*u*u))-300*(u*u*u))+ &
    &  100*(v*v))))/(18._dp*w*(omw*omw))+ &
    H1(1)*(67+72*u-48*u*v+(243*zeta2)/2._dp-162*(u*u)+156*v*(u*u)+ &
    &  180*(u*u*u)-180*v*(u*u*u)-90*(u*u*u*u)-24*u*(v*v)- &
    &  90*(u*u)*(v*v)+(nf*(-30-600*omv*(u*u*u)+300*(u*u*u*u)+ &
    u*(-115+83*v+32*(v*v))+u*u*(415-568*v+300*(v*v))))/9._dp) &
    +((-27*H1(0))/2._dp-(9*H1(1))/2._dp)*(G1(2)*G1(2))+ &
    G1(3)*((-33+2*nf)*H1(0)-18*(H1(0)*H1(0)))+ &
    G1(1)*(H1(0)*(-33+2*nf-18*H1(1))+(-33+2*nf)*H1(1)- &
    &  18*(H1(0)*H1(0)))+(-165._dp/2._dp+5*nf-(81*H1(1))/2._dp)*(H1(0)*H1(0))- &
    (81*(H1(0)*H1(0)*H1(0)))/2._dp+ &
    H1(0)*((-165._dp/2._dp+5*nf)*H1(1)-18*H2(0,0)+ &
    ((7236-360*nf)*(omv*omv)-480*omv*(-27+10*nf)*(u*u*u*u*u)+ &
    &  120*(-27+10*nf)*(u*u*u*u*u*u)+ &
    &  12*(-1+v)*(u*u*u)* &
    (-54*(33-36*v+20*(v*v))+nf*(533-768*v+400*(v*v)))+ &
    &  4*(u*u*u*u)*(-54*(103-176*v+90*(v*v))+ &
    nf*(1933-3568*v+1800*(v*v)))+ &
    &  4*u*(-1+v)*(-27*(-133-44*v+10*(v*v)+8*(v*v*v))+ &
    nf*(-171-234*v+69*(v*v)+32*(v*v*v)))+ &
    &  12*(u*u)*(-18*(13-126*v+105*(v*v)-48*(v*v*v)+ &
    &  15*(v*v*v*v))+ &
    nf*(195-634*v+669*(v*v)-368*(v*v*v)+100*(v*v*v*v)))+ &
    &  9234*zeta2*(w*w))/(36._dp*(w*w))-(45*(H1(1)*H1(1)))/2._dp)+ &
    (-33._dp/2._dp+nf)*(H1(1)*H1(1))+ &
    G1(2)*((33._dp/2._dp-nf)*H1(1)+H1(0)*(66-4*nf+36*H1(1))- &
    ((33-2*nf)*(33-2*nf))/36._dp+45*(H1(0)*H1(0))+9*(H1(1)*H1(1)))- &
    (9*(H1(1)*H1(1)*H1(1)))/2._dp)+((-33+2*nf)*(H1(1)*H1(1)*H1(1)))/6._dp+ &
    H1(0)*(-91+11*u+48*v-(60*v)/omu-(60*v)/omw-(30*v)/u-101*u*v- &
    (30*v)/w+27*zeta3+(99-6*nf)*H2(0,0)+(2*u*w*(nf*nf))/9._dp+ &
    H1(1)*(18*H2(0,0)-(5*(nf*nf))/9._dp- &
    (3*(291+4*u-4*u*v+18*zeta2-4*(u*u)))/4._dp+ &
    nf*(65._dp/3._dp+u-u*v-u*u))-11*(u*u)+90*v*(u*u)+ &
    &  6*zeta2*(33._dp/2._dp-24*u+6*(u*u)-4*(u*u*u))+63*(v*v)+ &
    (30*(v*v))/u+90*u*(v*v)+(30*(v*v))/w-(48*(v*v))/(omu*omu)- &
    (48*(v*v))/(omw*omw)+(nf* &
    (-57-103*v+(18*v)/omu+(18*v)/omw+(100*omv*v)/u+ &
    (100*v)/w+(39-300*v)*(u*u)+ &
    &  18*zeta2*(-3+6*u-6*(u*u)+4*(u*u*u))+ &
    u*(-39+339*v-300*(v*v))-164*(v*v)-(100*(v*v))/w+ &
    (114*(v*v))/(omu*omu)+(114*(v*v))/(omw*omw)))/9._dp+ &
    (-33._dp/2._dp+nf)*(H1(1)*H1(1))+(9*(H1(1)*H1(1)*H1(1)))/2._dp)+ &
    G1(1)*(-72*zeta3+6*(-33+2*nf)*H2(0,0)+(-132+8*nf)*H2(1,1)+ &
    &  36*H3(0,0,0)+36*H3(0,0,1)+36*H3(0,1,0)+36*H3(1,0,0)+72*H3(1,1,0)+ &
    &  6*zeta2*(-66+12*u-3*(u*u)+nf*(4-u+u*u-(2*(u*u*u))/3._dp)+ &
    &  2*(u*u*u))+H2(0,1)*(-165+36*u-9*(u*u)+ &
    nf*(10-3*u+3*(u*u)-2*(u*u*u))+6*(u*u*u))+ &
    H2(1,0)*(nf*(10+3*u-3*(u*u)+2*(u*u*u))- &
    &  3*(55+12*u-3*(u*u)+2*(u*u*u)))+ &
    (H1(1)*((540-200*nf)*(omv*omv)+ &
    &  2*omv*u*(135-1620*(u*u*u*u)+nf*(76+600*(u*u*u*u)))+ &
    u*u*(-4*nf*(-42-9*v+150*(u*u*u*u)+9*(v*v)+ &
    &  2*u*(-31+32*v+8*(v*v))+2*(u*u)*(106-142*v+75*(v*v)) &
    )-9*(342*zeta2- &
    &  2*(-179-6*v+90*(u*u*u*u)+6*(v*v)+ &
    &  3*u*(-25+21*v+8*(v*v))+ &
    &  3*(u*u)*(55-52*v+30*(v*v)))))))/(18._dp*(u*u))+ &
    (9*H1(0)+9*H1(1))*(G1(2)*G1(2))+(99-6*nf+45*H1(1))*(H1(0)*H1(0))+ &
    &  27*(H1(0)*H1(0)*H1(0))+G1(2)* &
    (H1(0)*(-33+2*nf-36*H1(1))+(-33+2*nf)*H1(1)-18*(H1(0)*H1(0))- &
    &  18*(H1(1)*H1(1)))+(33-2*nf)*(H1(1)*H1(1))+ &
    H1(0)*((132-8*nf)*H1(1)+ &
    ((540-200*nf)*(omv*omv)+ &
    &  2*omv*u*(135-1620*(u*u*u*u)+nf*(76+600*(u*u*u*u)))+ &
    u*u*(-4*nf*(-42-9*v+150*(u*u*u*u)+9*(v*v)+ &
    &  2*u*(-31+32*v+8*(v*v))+ &
    &  2*(u*u)*(106-142*v+75*(v*v)))- &
    &  9*(198*zeta2-2* &
    (-179-6*v+90*(u*u*u*u)+6*(v*v)+ &
    &  3*u*(-25+21*v+8*(v*v))+ &
    &  3*(u*u)*(55-52*v+30*(v*v))))))/(18._dp*(u*u))+ &
    &  27*(H1(1)*H1(1)))+9*(H1(1)*H1(1)*H1(1)))+(9*(H1(1)*H1(1)*H1(1)*H1(1)))/8._dp
    return
    end function Hbeta_2a2re
          
          
    function Hbeta_2a2im(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    use constants
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=36*pi*G3(0,0,1)+18*pi*G3(0,0,2)-36*pi*G3(0,0,3)+36*pi*G3(0,1,0)+ &
    &  18*pi*G3(0,2,0)-18*pi*G3(0,2,1)+18*pi*G3(0,2,2)-18*pi*G3(0,2,3)- &
    &  18*pi*G3(0,3,0)-18*pi*G3(0,3,2)+72*pi*G3(0,3,3)-36*pi*G3(1,0,2)- &
    &  72*pi*G3(1,1,1)+72*pi*G3(1,1,2)-36*pi*G3(1,2,0)+72*pi*G3(1,2,1)- &
    &  36*pi*G3(1,2,2)+18*pi*G3(2,0,0)-18*pi*G3(2,0,1)+18*pi*G3(2,0,2)- &
    &  18*pi*G3(2,0,3)-18*pi*G3(2,1,0)+72*pi*G3(2,1,1)-18*pi*G3(2,1,2)+ &
    &  18*pi*G3(2,2,0)-36*pi*G3(2,2,1)+36*pi*G3(2,2,3)+36*pi*G3(2,3,2)- &
    &  36*pi*G3(3,0,0)-36*pi*G3(3,0,2)+72*pi*G3(3,0,3)-36*pi*G3(3,2,0)+ &
    &  72*pi*G3(3,3,0)-72*pi*G3(3,3,3)+G2(2,2)*((66-4*nf)*pi+36*pi*H1(0))+ &
    G2(0,0)*((66-4*nf)*pi-36*pi*H1(0)-18*pi*H1(1))+ &
    G2(0,2)*((66-4*nf)*pi-9*pi*H1(0)-18*pi*H1(1))+ &
    G2(2,0)*((66-4*nf)*pi-9*pi*H1(0)-18*pi*H1(1))+(66-4*nf)*pi*H2(1,1)- &
    &  108*pi*H3(0,0,0)-18*pi*H3(0,0,1)-18*pi*H3(0,1,0)-18*pi*H3(1,0,0)+ &
    G2(1,1)*(-72*pi*H1(1)+2*pi*u* &
    (36-9*u+nf*(-3+3*u-2*(u*u))+6*(u*u)))+ &
    G2(0,1)*(-18*pi*H1(0)+pi*u*(-36+9*u-6*(u*u)+ &
    nf*(3-3*u+2*(u*u))))+ &
    G2(1,2)*(36*pi*H1(1)+pi*u*(-36+9*u-6*(u*u)+ &
    nf*(3-3*u+2*(u*u))))+ &
    G2(1,0)*(18*pi*H1(0)+36*pi*H1(1)+ &
    pi*u*(-36+9*u-6*(u*u)+nf*(3-3*u+2*(u*u))))+ &
    (pi*H2(0,1)*(3*(77-24*u+6*(u*u)-4*(u*u*u))+ &
    &  2*nf*(-7+3*u-3*(u*u)+2*(u*u*u))))/2._dp+ &
    G2(0,3)*(54*pi*H1(0)+18*pi*H1(1)+ &
    pi*(-165+36*u-9*(u*u)+nf*(10-3*u+3*(u*u)-2*(u*u*u))+ &
    &  6*(u*u*u)))+G2(3,3)*(-72*pi*H1(0)- &
    &  2*pi*w*(3*omw-3*(11+4*u*v+2*(u*u)+2*(v*v))+ &
    nf*(2-u-v+4*u*v+2*(u*u)+2*(v*v))))+ &
    G2(3,2)*(18*pi*H1(0)+pi*w*(3*omw-3*(11+4*u*v+2*(u*u)+2*(v*v))+ &
    nf*(2-u-v+4*u*v+2*(u*u)+2*(v*v))))+ &
    G2(2,3)*(-18*pi*H1(0)-36*pi*H1(1)+ &
    pi*w*(3*omw-3*(11+4*u*v+2*(u*u)+2*(v*v))+ &
    nf*(2-u-v+4*u*v+2*(u*u)+2*(v*v))))+ &
    G2(3,0)*(72*pi*H1(0)+36*pi*H1(1)+ &
    pi*w*(3*omw-3*(11+4*u*v+2*(u*u)+2*(v*v))+ &
    nf*(2-u-v+4*u*v+2*(u*u)+2*(v*v))))+ &
    (pi*H1(1)*(4*nf*(-30-600*omv*(u*u*u)+300*(u*u*u*u)+ &
    u*(-115+83*v+32*(v*v))+u*u*(415-568*v+300*(v*v)))+ &
    &  9*(270*zeta2+4*(67+180*omv*(u*u*u)-90*(u*u*u*u)- &
    &  24*u*(-3+2*v+v*v)-6*(u*u)*(27-26*v+15*(v*v))))))/36._dp &
    +pi*H2(0,0)*(3*(88-12*v+6*omv*(u*u)-4*(u*u*u)+3*(v*v)- &
    &  6*u*(4-v+v*v)-2*(v*v*v))+ &
    nf*(-16-6*omv*(-1+u)*u+3*v+4*(u*u*u)-3*(v*v)+6*u*(v*v)+ &
    &  2*(v*v*v)))+(pi*H2(1,0)* &
    (2*nf*(-7+3*v+(-3+6*v)*(u*u)+2*(u*u*u)-3*(v*v)+ &
    u*(3-6*v+6*(v*v))+2*(v*v*v))- &
    &  3*(-77+24*v+6*(-1+2*v)*(u*u)+4*(u*u*u)-6*(v*v)+ &
    &  12*u*(2-v+v*v)+4*(v*v*v))))/2._dp+ &
    (pi*(2*u*w*(nf*nf)+3*(-202-144*v+(180*v)/omu+(180*v)/omw- &
    &  90/w+(180*v)/w+27*zeta3-(90*(omv*omv))/u- &
    &  3*(11+90*v)*(u*u)+u*(33+237*v-270*(v*v))-189*(v*v)- &
    (90*(v*v))/w+(144*(v*v))/(omu*omu)+(144*(v*v))/(omw*omw)- &
    &  18*zeta2*(-22+6*v*(2-u+u*u)+(-3+6*u)*(v*v)+2*(v*v*v))) &
    +nf*(-98+103*v-(18*v)/omu-(18*v)/omw+100/w-(200*v)/w+ &
    (100*(omv*omv))/u+(39+300*v)*(u*u)+164*(v*v)+ &
    (100*(v*v))/w-(114*(v*v))/(omu*omu)- &
    (114*(v*v))/(omw*omw)+3*u*(-13-87*v+100*(v*v))+ &
    &  18*zeta2*(-4+v*(3-6*u+6*(u*u))+(-3+6*u)*(v*v)+ &
    &  2*(v*v*v)))))/9._dp+ &
    G2(2,1)*(-18*pi*H1(0)+18*pi*H1(1)+ &
    pi*(-3*(44+12*v+(-3+6*v)*(u*u)+2*(u*u*u)-3*(v*v)+ &
    &  6*u*(2-v+v*v)+2*(v*v*v))+ &
    nf*(8+3*v+(-3+6*v)*(u*u)+2*(u*u*u)-3*(v*v)+ &
    u*(3-6*v+6*(v*v))+2*(v*v*v))))+ &
    ((-33._dp/2._dp+nf)*pi-9*pi*G1(1)+(9*pi*G1(2))/2._dp-9*pi*G1(3)- &
    (27*pi*H1(0))/2._dp-(9*pi*H1(1))/2._dp)*(G1(0)*G1(0))+ &
    (9*pi*(G1(0)*G1(0)*G1(0)))/2._dp+ &
    ((-33._dp/2._dp+nf)*pi-9*pi*G1(3)-(27*pi*H1(0))/2._dp-(27*pi*H1(1))/2._dp)* &
    (G1(2)*G1(2))+(9*pi*(G1(2)*G1(2)*G1(2)))/2._dp+ &
    ((7*(-33+2*nf)*pi)/2._dp-(45*pi*H1(1))/2._dp)*(H1(0)*H1(0))- &
    (27*pi*(H1(0)*H1(0)*H1(0)))/2._dp+(-33._dp/2._dp+nf)*pi*(H1(1)*H1(1))+ &
    H1(0)*((5*(-33+2*nf)*pi*H1(1))/2._dp+18*pi*H2(0,0)+ &
    (pi*(-480*omv*(-27+10*nf)*(u*u*u*u*u)+ &
    &  120*(-27+10*nf)*(u*u*u*u*u*u)- &
    &  12*(omv*omv)*(2*nf*(5-3*v+3*(v*v))- &
    &  3*(67-6*v+6*(v*v)))+ &
    &  4*(u*u*u*u)*(-108*(51-88*v+45*(v*v))+ &
    nf*(1915-3568*v+1800*(v*v)))+ &
    &  4*u*(-1+v)*(9*(137+96*v-24*(v*v*v))+ &
    nf*(-69-126*v-21*(v*v)+32*(v*v*v)))+ &
    &  12*(u*u*u)*(-36*(-48+101*v-84*(v*v)+30*(v*v*v))+ &
    nf*(-515+1271*v-1168*(v*v)+400*(v*v*v)))+ &
    &  12*(u*u)*(-6*(97-345*v+291*(v*v)-144*(v*v*v)+ &
    &  45*(v*v*v*v))+ &
    nf*(197-568*v+621*(v*v)-368*(v*v*v)+100*(v*v*v*v)))+ &
    &  2430*zeta2*(w*w)))/(36._dp*(w*w))-(27*pi*(H1(1)*H1(1)))/2._dp)+ &
    G1(3)*(-108*pi*H2(0,0)-18*pi*H2(0,1)-18*pi*H2(1,0)+ &
    omw*pi*H1(1)*(nf*(3-3*v+u*(-3+4*v)+2*(u*u)+2*(v*v))- &
    &  3*(12-3*v+u*(-3+4*v)+2*(u*u)+2*(v*v)))+ &
    H1(0)*(-18*pi*H1(1)+2*omw*pi* &
    (nf*(3-3*v+u*(-3+4*v)+2*(u*u)+2*(v*v))- &
    &  3*(12-3*v+u*(-3+4*v)+2*(u*u)+2*(v*v))))+ &
    (pi*(-240*omv*(-27+10*nf)*(u*u*u*u*u)+ &
    &  60*(-27+10*nf)*(u*u*u*u*u*u)- &
    &  12*(omv*omv)*(-201+9*v-9*(v*v)+nf*(10-3*v+3*(v*v)))+ &
    &  2*(u*u*u*u)*(-27*(205-352*v+180*(v*v))+ &
    &  4*nf*(481-892*v+450*(v*v)))+ &
    &  2*u*(-1+v)*(9*(268+102*v-3*(v*v)-24*(v*v*v))+ &
    &  4*nf*(-30-36*v-3*(v*v)+8*(v*v*v)))+ &
    &  6*(u*u*u)*(4*nf*(-131+320*v-292*(v*v)+100*(v*v*v))- &
    &  9*(-195+407*v-336*(v*v)+120*(v*v*v)))+ &
    &  6*(u*u)*(2*nf*(98-293*v+315*(v*v)-184*(v*v*v)+ &
    &  50*(v*v*v*v))- &
    &  3*(136-708*v+591*(v*v)-288*(v*v*v)+90*(v*v*v*v)))+ &
    &  2430*zeta2*(w*w)))/(18._dp*(w*w))-9*pi*(H1(0)*H1(0))- &
    &  9*pi*(H1(1)*H1(1)))+G1(1)*(2*(-33+2*nf)*pi*H1(0)+ &
    G1(2)*((33-2*nf)*pi+18*pi*H1(0)+18*pi*H1(1))-36*pi*H2(0,0)- &
    &  18*pi*H2(0,1)-18*pi*H2(1,0)-36*pi*H2(1,1)+ &
    omu*pi*H1(1)*(-33+3*u-6*(u*u)+nf*(2-u+2*(u*u)))+ &
    (pi*(20*(-27+10*nf)*(omv*omv)- &
    &  2*omv*u*(135-1620*(u*u*u*u)+nf*(76+600*(u*u*u*u)))+ &
    u*u*(4*nf*(-42-9*v+150*(u*u*u*u)+9*(v*v)+ &
    &  2*u*(-31+32*v+8*(v*v))+2*(u*u)*(106-142*v+75*(v*v)) &
    )+9*(270*zeta2- &
    &  2*(-179-6*v+90*(u*u*u*u)+6*(v*v)+ &
    &  3*u*(-25+21*v+8*(v*v))+ &
    &  3*(u*u)*(55-52*v+30*(v*v)))))))/(18._dp*(u*u))- &
    &  9*pi*(G1(2)*G1(2))-9*pi*(H1(0)*H1(0))-9*pi*(H1(1)*H1(1)))+ &
    G1(0)*(G1(1)*((33-2*nf)*pi+18*pi*H1(0))+ &
    G1(3)*((33-2*nf)*pi+18*pi*H1(0))+(-33+2*nf)*pi*H1(1)+ &
    G1(2)*((-33+2*nf)*pi-18*pi*H1(0)-9*pi*H1(1))+90*pi*H2(0,0)+ &
    &  36*pi*H2(0,1)+36*pi*H2(1,0)+18*pi*H2(1,1)+ &
    H1(0)*(9*pi*H1(1)+omu*pi* &
    (-33+3*u-6*(u*u)+nf*(2-u+2*(u*u))))- &
    (pi*(4*nf*(-30-600*omv*(u*u*u)+300*(u*u*u*u)+ &
    u*(-115+83*v+32*(v*v))+u*u*(415-568*v+300*(v*v)))+ &
    &  9*(270*zeta2+4*(67+180*omv*(u*u*u)-90*(u*u*u*u)- &
    &  24*u*(-3+2*v+v*v)-6*(u*u)*(27-26*v+15*(v*v))))))/ &
    &  36._dp+(9*pi*(G1(2)*G1(2)))/2._dp+(9*pi*(H1(0)*H1(0)))/2._dp+ &
    (9*pi*(H1(1)*H1(1)))/2._dp)+G1(2)* &
    ((-33+2*nf)*pi*H1(1)+G1(3)*((33-2*nf)*pi+18*pi*H1(0)+18*pi*H1(1))+ &
    &  18*pi*H2(0,0)-18*pi*H2(0,1)-18*pi*H2(1,0)- &
    (pi*(4*nf*(-30-600*omv*(u*u*u)+300*(u*u*u*u)+ &
    u*(-115+83*v+32*(v*v))+u*u*(415-568*v+300*(v*v)))+ &
    &  9*(270*zeta2+4*(67+180*omv*(u*u*u)-90*(u*u*u*u)- &
    &  24*u*(-3+2*v+v*v)-6*(u*u)*(27-26*v+15*(v*v))))))/ &
    &  36._dp+H1(0)*(27*pi*H1(1)+ &
    pi*(3*(33+12*v+(-3+6*v)*(u*u)+2*(u*u*u)-3*(v*v)+ &
    &  6*u*(2-v+v*v)+2*(v*v*v))- &
    nf*(6+3*v+(-3+6*v)*(u*u)+2*(u*u*u)-3*(v*v)+ &
    u*(3-6*v+6*(v*v))+2*(v*v*v))))+(9*pi*(H1(0)*H1(0)))/2._dp+ &
    (27*pi*(H1(1)*H1(1)))/2._dp)-(9*pi*(H1(1)*H1(1)*H1(1)))/2._dp
    return
    end function Hbeta_2a2im
          
          
!**********************************************************************
! egion2a:q+qbar-->V+g*
!*********************************************************************/
          
    function Hgamma_2a0re(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=1
    return
    end function Hgamma_2a0re
          
          
    function Hgamma_2a0im(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=0
    return
    end function Hgamma_2a0im
          
          
    function Hgamma_2a1re(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    use constants
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=(528-40*nf+60/u+126*zeta2)/36._dp+G2(0,2)/3._dp+(8*G2(1,2))/3._dp+ &
    G2(2,0)/3._dp+(8*G2(3,0))/3._dp-(8*G1(3)*H1(0))/3._dp+ &
    H1(0)*(7-(2*nf)/3._dp-3*H1(1))+G1(1)*((-8*H1(0))/3._dp-(8*H1(1))/3._dp)+ &
    G1(0)*((17*H1(0))/3._dp-H1(1)/3._dp)+G1(2)*((17*H1(0))/3._dp+3*H1(1))- &
    (8*H2(0,0))/3._dp+(10*H2(1,0))/3._dp-(3*(G1(0)*G1(0)))/2._dp-(3*(G1(2)*G1(2)))/2._dp- &
    (17*(H1(0)*H1(0)))/6._dp-(3*(H1(1)*H1(1)))/2._dp
    return
    end function Hgamma_2a1re
          
          
    function Hgamma_2a1im(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    use constants
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=(7-(2*nf)/3._dp)*pi+(pi*G1(0))/3._dp+(8*pi*G1(1))/3._dp+(pi*G1(2))/3._dp+ &
    (8*pi*G1(3))/3._dp+3*pi*H1(0)-(pi*H1(1))/3._dp
    return
    end function Hgamma_2a1im
          
          
    function Hgamma_2a2re(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    use constants
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=(64*G4(0,0,1,2))/9._dp-(118*G4(0,0,2,2))/9._dp+(16*G4(0,0,3,0))/9._dp+ &
    (154*G4(0,1,0,2))/9._dp-(80*G4(0,1,1,2))/9._dp+(154*G4(0,1,2,0))/9._dp+ &
    (116*G4(0,1,2,2))/9._dp-(118*G4(0,2,0,2))/9._dp+(88*G4(0,2,1,2))/9._dp- &
    (118*G4(0,2,2,0))/9._dp+(88*G4(0,2,3,0))/9._dp-8*G4(0,3,0,2)-8*G4(0,3,2,0)+ &
    (128*G4(0,3,3,0))/9._dp+(80*G4(1,0,1,2))/9._dp+4*G4(1,0,2,2)+ &
    (10*G4(1,0,3,0))/9._dp-(128*G4(1,1,1,2))/9._dp+(128*G4(1,1,2,2))/9._dp+ &
    &  4*G4(1,2,0,2)+(16*G4(1,2,1,2))/3._dp+4*G4(1,2,2,0)+(10*G4(1,2,3,0))/9._dp- &
    (70*G4(1,3,0,1))/9._dp-(70*G4(1,3,2,1))/9._dp-(118*G4(2,0,0,2))/9._dp+ &
    (88*G4(2,0,1,2))/9._dp-(118*G4(2,0,2,0))/9._dp+(88*G4(2,0,3,0))/9._dp- &
    &  8*G4(2,1,0,2)+(128*G4(2,1,1,2))/9._dp+(70*G4(2,1,1,3))/3._dp-8*G4(2,1,2,0)- &
    (140*G4(2,1,2,1))/9._dp-(140*G4(2,1,2,3))/9._dp+(70*G4(2,1,3,1))/3._dp- &
    (118*G4(2,2,0,0))/9._dp-(280*G4(2,2,1,1))/9._dp+(16*G4(2,2,1,2))/9._dp- &
    (70*G4(2,2,1,3))/9._dp+(64*G4(2,2,3,0))/9._dp-(70*G4(2,2,3,1))/9._dp+ &
    (116*G4(2,3,0,0))/9._dp-(70*G4(2,3,0,1))/9._dp+(154*G4(2,3,0,2))/9._dp+ &
    (70*G4(2,3,1,1))/3._dp+(154*G4(2,3,2,0))/9._dp-(80*G4(2,3,3,0))/9._dp+ &
    &  4*G4(3,0,0,2)-(140*G4(3,0,1,1))/9._dp-(20*G4(3,0,1,2))/3._dp+4*G4(3,0,2,0)- &
    (70*G4(3,0,2,1))/9._dp+(16*G4(3,0,3,0))/3._dp-(70*G4(3,1,0,1))/9._dp- &
    (70*G4(3,1,2,1))/9._dp+4*G4(3,2,0,0)-(70*G4(3,2,0,1))/9._dp+ &
    (70*G4(3,2,1,1))/9._dp-(20*G4(3,2,1,2))/3._dp+(70*G4(3,2,2,1))/9._dp+ &
    (80*G4(3,2,3,0))/9._dp+(128*G4(3,3,0,0))/9._dp-(128*G4(3,3,3,0))/9._dp+ &
    G3(3,3,0)*((8*(-33+2*nf))/9._dp-(128*H1(0))/9._dp)+ &
    G3(3,0,2)*(-49._dp/3._dp+(4*nf)/9._dp-(28*H1(0))/9._dp)+ &
    G3(3,2,0)*(-49._dp/3._dp+(4*nf)/9._dp-(28*H1(0))/9._dp)-(16*G3(0,0,3)*H1(0))/9._dp- &
    (88*G3(0,2,3)*H1(0))/9._dp+8*G3(0,3,2)*H1(0)-(128*G3(0,3,3)*H1(0))/9._dp- &
    (10*G3(1,0,3)*H1(0))/9._dp-(10*G3(1,2,3)*H1(0))/9._dp-(88*G3(2,0,3)*H1(0))/9._dp- &
    (70*G3(2,1,3)*H1(0))/9._dp+(2*G3(2,2,3)*H1(0))/3._dp-(70*G3(2,3,1)*H1(0))/9._dp- &
    (154*G3(2,3,2)*H1(0))/9._dp+(80*G3(2,3,3)*H1(0))/9._dp-(16*G3(3,0,3)*H1(0))/3._dp- &
    (70*G3(3,2,2)*H1(0))/9._dp-(80*G3(3,2,3)*H1(0))/9._dp+(128*G3(3,3,3)*H1(0))/9._dp+ &
    G3(2,2,0)*((29+nf)/3._dp+(100*H1(0))/9._dp)+ &
    G3(0,2,2)*((29+nf)/3._dp+(116*H1(0))/9._dp)+ &
    G3(2,0,2)*((29+nf)/3._dp+(350*H1(0))/9._dp)+ &
    G3(0,1,0)*((-154*H1(0))/9._dp-(154*H1(1))/9._dp)+ &
    G3(2,3,0)*(39-(32*nf)/9._dp-(98*H1(0))/9._dp-(154*H1(1))/9._dp)+ &
    G3(2,1,1)*((-128*H1(0))/9._dp-(128*H1(1))/9._dp)+ &
    G3(0,2,1)*((-88*H1(0))/9._dp-(88*H1(1))/9._dp)+ &
    G3(2,0,1)*((-88*H1(0))/9._dp-(88*H1(1))/9._dp)+ &
    G3(1,0,1)*((-80*H1(0))/9._dp-(80*H1(1))/9._dp)+ &
    G3(0,0,1)*((-64*H1(0))/9._dp-(64*H1(1))/9._dp)+ &
    G3(1,2,1)*((-16*H1(0))/3._dp-(16*H1(1))/3._dp)+ &
    G3(1,2,0)*((4*(-30+nf))/9._dp-(98*H1(0))/9._dp-4*H1(1))+ &
    G3(1,0,2)*((4*(-30+nf))/9._dp-(28*H1(0))/9._dp-4*H1(1))+ &
    G3(3,0,0)*((8*(-21+nf))/3._dp+(128*H1(0))/9._dp-4*H1(1))+ &
    G3(2,2,1)*((124*H1(0))/9._dp-(16*H1(1))/9._dp)+ &
    G3(3,2,1)*((-10*H1(0))/9._dp-(10*H1(1))/9._dp)+ &
    G3(3,0,1)*((20*H1(0))/3._dp-(10*H1(1))/9._dp)+ &
    G3(0,3,0)*((4*(-96+5*nf))/9._dp+(64*H1(0))/9._dp+8*H1(1))+ &
    G3(2,1,0)*(8*H1(0)+8*H1(1))+G3(0,1,1)*((80*H1(0))/9._dp+(80*H1(1))/9._dp)+ &
    G3(0,0,2)*(-7+nf/3._dp+(100*H1(0))/9._dp+(118*H1(1))/9._dp)+ &
    G3(2,0,0)*(-7+nf/3._dp+(116*H1(0))/9._dp+(118*H1(1))/9._dp)+ &
    G3(0,2,0)*(-7+nf/3._dp+(350*H1(0))/9._dp+(118*H1(1))/9._dp)+ &
    G3(1,1,1)*((128*H1(0))/9._dp+(128*H1(1))/9._dp)+ &
    G2(3,3)*((-256*zeta2)/9._dp+(8*(33-2*nf)*H1(0))/9._dp+(128*H2(0,0))/9._dp)+ &
    G2(0,3)*((304*zeta2)/9._dp-(4*(-96+5*nf)*H1(0))/9._dp-(128*H2(0,0))/9._dp- &
    &  8*H2(0,1)-8*H2(1,0))+G2(3,2)* &
    ((-58*zeta2)/9._dp+(49._dp/3._dp-(4*nf)/9._dp)*H1(0)+(100*H2(0,0))/9._dp+ &
    (10*H2(0,1))/9._dp+(70*H2(1,0))/9._dp)+ &
    G2(2,3)*((-1192*zeta2)/9._dp+(-39+(32*nf)/9._dp)*H1(0)+(80*H2(0,0))/9._dp+ &
    (154*H2(0,1))/9._dp+(154*H2(1,0))/9._dp)+ &
    G2(1,0)*((-128*zeta2)/9._dp-(4*(-30+nf)*H1(0))/9._dp- &
    (4*(-30+nf)*H1(1))/9._dp+(170*H2(0,0))/9._dp+12*H2(0,1)+(98*H2(1,0))/9._dp+ &
    &  4*H2(1,1))+(-132+8*nf)*H3(0,0,0)+(-49._dp/3._dp+(4*nf)/9._dp)*H3(0,0,1)+ &
    (226._dp/3._dp-(46*nf)/9._dp)*H3(0,1,0)+(149-(92*nf)/9._dp)*H3(1,0,0)+ &
    (64*H4(0,0,0,0))/3._dp-20*H4(0,0,0,1)-(100*H4(0,0,1,0))/9._dp- &
    (100*H4(0,1,0,0))/9._dp-(10*H4(0,1,0,1))/9._dp-(70*H4(0,1,1,0))/9._dp- &
    (250*H4(1,0,0,0))/9._dp-(10*H4(1,0,0,1))/9._dp-(20*H4(1,0,1,0))/9._dp+ &
    (80*H4(1,1,0,0))/9._dp+(250*H4(1,1,0,1))/9._dp-(430*H4(1,1,1,0))/9._dp+ &
    H2(1,1)*(80-(41*nf)/6._dp-40/(9._dp*u)-(406*zeta2)/9._dp+(5*(nf*nf))/36._dp)+ &
    G2(0,0)*(80-(41*nf)/6._dp+(194*zeta2)/9._dp+(8*(-21+nf)*H1(0))/3._dp+ &
    (7-nf/3._dp)*H1(1)+(128*H2(0,0))/9._dp-4*H2(0,1)-(116*H2(1,0))/9._dp- &
    (118*H2(1,1))/9._dp+(5*(nf*nf))/36._dp)-(4*(-25+54*zeta2)*(nf*nf))/81._dp+ &
    (H2(0,1)*(-1395*nf+30*(495-(2*w)/omw+(2*(-8*v+w))/u)+696*zeta2+ &
    &  30*(nf*nf)))/108._dp+G3(1,1,2)* &
    ((-128*H1(0))/9._dp-(128*H1(1))/9._dp+ &
    (8*(-5+17*u+(-21+2*nf)*(u*u)))/(9._dp*(u*u)))+ &
    G2(1,1)*((-128*zeta2)/9._dp+(128*H2(1,0))/9._dp+(128*H2(1,1))/9._dp- &
    (8*H1(0)*(-5+17*u+(-21+2*nf)*(u*u)))/(9._dp*(u*u))- &
    (8*H1(1)*(-5+17*u+(-21+2*nf)*(u*u)))/(9._dp*(u*u)))+ &
    G3(1,2,2)*((128*H1(0))/9._dp+(-40+190*u+6*(-109+4*nf)*(u*u))/ &
    (9._dp*(u*u)))+G2(2,2)*(80-(41*nf)/6._dp-40/(9._dp*u)+(242*zeta2)/9._dp+ &
    (80*H2(0,0))/9._dp+(16*H2(0,1))/9._dp-(16*H2(1,0))/9._dp+(5*(nf*nf))/36._dp+ &
    (H1(0)*(-40+190*u+6*(-109+4*nf)*(u*u)))/(9._dp*(u*u)))+ &
    G2(2,1)*((176*zeta2)/9._dp+(64*H2(0,1))/9._dp-(64*H2(1,0))/9._dp+ &
    (H1(0)*(2*v*(13+(78-10*nf)*v)*(u*u)+ &
    (-27+(861-40*nf)*v)*(u*u*u)+(507-20*nf)*(u*u*u*u)+ &
    u*v*(-10+143*v-243*(v*v))-5*(v*v)*(1-18*v+9*(v*v))))/ &
    (9._dp*(omw*omw)*(u*u))+ &
    (H1(1)*(2*v*(13+(78-10*nf)*v)*(u*u)+ &
    (-27+(861-40*nf)*v)*(u*u*u)+(507-20*nf)*(u*u*u*u)+ &
    u*v*(-10+143*v-243*(v*v))-5*(v*v)*(1-18*v+9*(v*v))))/ &
    (9._dp*(omw*omw)*(u*u)))+ &
    G3(2,1,2)*((64*H1(0))/9._dp+(2*v*(-13+2*(-39+5*nf)*v)*(u*u)+ &
    (27+(-861+40*nf)*v)*(u*u*u)+(-507+20*nf)*(u*u*u*u)+ &
    &  5*(v*v)*(1-18*v+9*(v*v))+u*v*(10-143*v+243*(v*v)))/ &
    (9._dp*(omw*omw)*(u*u)))+ &
    G3(0,1,2)*((-98*H1(0))/9._dp-(116*H1(1))/9._dp+ &
    (10*u*(8-19*v)*v-4*(50+(-189+16*nf)*v)*(u*u*u)+ &
    (378-32*nf)*(u*u*u*u)+40*(v*v)+ &
    u*u*(35-390*v+(378-32*nf)*(v*v)))/(9._dp*(omw*omw)*(u*u)))+ &
    G2(0,1)*((-1112*zeta2)/9._dp+2*H2(0,1)+(98*H2(1,0))/9._dp+(116*H2(1,1))/9._dp+ &
    (H1(0)*(10*u*v*(-8+19*v)+4*(50+(-189+16*nf)*v)*(u*u*u)+ &
    (-378+32*nf)*(u*u*u*u)-40*(v*v)+ &
    u*u*(-35+390*v+(-378+32*nf)*(v*v))))/(9._dp*(omw*omw)*(u*u)) &
    +(H1(1)*(10*u*v*(-8+19*v)+4*(50+(-189+16*nf)*v)*(u*u*u)+ &
    (-378+32*nf)*(u*u*u*u)-40*(v*v)+ &
    u*u*(-35+390*v+(-378+32*nf)*(v*v))))/(9._dp*(omw*omw)*(u*u))) &
    +H1(1)*((-4*(-30+nf)*H2(0,0))/9._dp+(5*(-30+nf)*H2(1,0))/9._dp- &
    ((4*nf*(-138+84*v+42*(-1+2*v)*zeta2)+ &
    &  3*(1981+60*(-5+8*v)*zeta2-1702*zeta3+ &
    &  4*v*(-624+851*zeta3)))*(u*u*u*u)+ &
    &  2*(14*nf*(6+6*zeta2)+3*(-618+120*zeta2+851*zeta3))* &
    (u*u*u*u*u)-6*omv*(-407+54*nf)*(v*v)- &
    &  3*u*v*(-1628+3187*v-542*(v*v)+4*nf*(54-113*v+23*(v*v)))+ &
    u*u*u*(4*nf*(177-v*(345+84*zeta2)+7*(6+6*zeta2)*(v*v))+ &
    &  3*(-1559+v*(4504-3404*zeta3)+2*(-642+851*zeta3)*(v*v)+ &
    &  30*zeta2*(1-18*v+8*(v*v))))- &
    u*u*(4*nf*(81-435*v+(276+42*zeta2)*(v*v))+ &
    &  3*(-814+3932*v+(-3065+1702*zeta3)*(v*v)+ &
    &  30*zeta2*(-1-2*v+8*(v*v))+24*(v*v*v))))/ &
    (54._dp*omu*(omw*omw)*(u*u)))+ &
    (H2(0,0)*((39639-5336*nf-432*v+4272*zeta2+180*(nf*nf))*(u*u*u)- &
    &  144*(u*u*u*u)-2*(u*u)* &
    (35613-39399*v+nf*(-5012+5336*v)+4272*omv*zeta2+ &
    &  180*omv*(nf*nf)+216*(v*v))+ &
    u*(34653-70266*v+4272*zeta2*(omv*omv)+ &
    &  180*(omv*omv)*(nf*nf)+38679*(v*v)- &
    &  4*nf*(1253-2506*v+1334*(v*v))-144*(v*v*v))- &
    &  480*(1-2*(v*v)+v*v*v)))/(108._dp*u*(w*w))+ &
    (H3(1,1,0)*(-723+42*nf-(20*w)/omw-(5*(w*w))/(omw*omw)))/9._dp+ &
    (H3(1,0,1)*(138-5*nf+(20*w)/omw+(5*(w*w))/(omw*omw)))/9._dp+ &
    nf*(-53237._dp/486._dp-352/(27._dp*u)+(38*v)/(3-3*v)-(589*zeta3)/27._dp+ &
    (zeta2*(108*u*v+54*(u*u)+54*(v*v)+1055*(w*w)))/(18._dp*(w*w)))+ &
    (H2(1,0)*((29460-576*v)*(u*u*u*u)-144*(u*u*u*u*u)+ &
    &  12*(u*u*u)*(-5559+27*omw*nf+7370*v-72*(v*v))+ &
    &  6*(u*u)*(5845+2*(-11183+54*omw*nf)*v+14760*(v*v)- &
    &  96*(v*v*v))+60*v*(-3+15*v-13*(v*v)+v*v*v)+ &
    u*(-120+35970*v+12*(-5689+27*omw*nf)*(v*v)+29640*(v*v*v)- &
    &  144*(v*v*v*v)+omw*(-2555*nf+12036*zeta2+30*(nf*nf))*(w*w)))) &
    /(108._dp*omw*u*(w*w))+(80* &
    (134459+24336/u-(7776*v)/omv-11133*zeta3)+2582370*zeta4- &
    (360*zeta2*(60*u*(-1+u+v)*(w*w)+ &
    omw*(4*(407*v-4*w*(67+12*w))*(u*u)+814*(u*u*u)- &
    &  10*w*(77*w+116*v*w+16*(v*v))+ &
    u*(-16*v*w*(77+12*w)+814*(v*v)+9933*(w*w)))))/ &
    (omw*u*(w*w)))/12960._dp+ &
    ((33-2*nf)/6._dp-(17*H1(0))/2._dp+H1(1)/2._dp)*(G1(0)*G1(0)*G1(0))+ &
    (9*(G1(0)*G1(0)*G1(0)*G1(0)))/8._dp+ &
    ((33-2*nf)/6._dp-(17*H1(0))/2._dp-(9*H1(1))/2._dp)*(G1(2)*G1(2)*G1(2))+ &
    (9*(G1(2)*G1(2)*G1(2)*G1(2)))/8._dp+ &
    ((68*(-12+nf))/27._dp+(17*H1(1))/2._dp)*(H1(0)*H1(0)*H1(0))+ &
    (289*(H1(0)*H1(0)*H1(0)*H1(0)))/72._dp+ &
    G2(1,2)*(H1(0)*(137._dp/3._dp-(20*nf)/9._dp-3/u-8*H1(1))-(154*H2(0,0))/9._dp- &
    (80*H2(0,1))/9._dp-22*H2(1,0)+ &
    (H1(1)*(40-190*u+(534-20*nf)*(u*u)))/(9._dp*(u*u))+ &
    (6*omv*(407-54*nf)*v+ &
    (4827-302*nf+144*v-2076*zeta2)*(u*u*u*u)+72*(u*u*u*u*u)+ &
    u*u*(-4023+4254*v-54*nf*(-10+9*v)-30*(v*v))+ &
    u*u*u*(4002+4797*v-2*nf*(189+151*v)-2076*v*zeta2+ &
    &  72*(v*v))-3*u*(-814+2165*v-74*(v*v)+36*nf*(3-8*v+v*v)) &
    )/(54._dp*omw*(u*u*u))-(68*(H1(0)*H1(0)))/9._dp-4*(H1(1)*H1(1)))+ &
    G2(3,0)*(H1(0)*(32-(20*nf)/9._dp-8*H1(1))+3*H1(1)-(64*H2(0,0))/3._dp+ &
    &  12*H2(0,1)+(28*H2(1,0))/9._dp+ &
    ((4683-302*nf+216*v-1848*zeta2)*(u*u*u)+72*(u*u*u*u)+ &
    &  2*(u*u)*(-6696+nf*(464-302*v)+4803*v+1848*omv*zeta2+ &
    &  108*(v*v))+240*(1-2*(v*v)+v*v*v)+ &
    u*(7176-13872*v-1848*zeta2*(omv*omv)+ &
    nf*(-464+928*v-302*(v*v))+5163*(v*v)+72*(v*v*v)))/ &
    (54._dp*u*(w*w))-(68*(H1(0)*H1(0)))/9._dp-4*(H1(1)*H1(1)))+ &
    G2(2,0)*(H1(0)*(95._dp/3._dp-(32*nf)/9._dp-H1(1))+(-34._dp/3._dp-(5*nf)/18._dp)*H1(1)- &
    (268*H2(0,0))/9._dp-(262*H2(0,1))/9._dp-(172*H2(1,0))/9._dp+ &
    (-60*(1+8*v)+u*(5034-296*nf-144*v+7518*zeta2+3*(nf*nf))- &
    &  144*(u*u))/(108._dp*u)-(17*(H1(0)*H1(0)))/18._dp-(H1(1)*H1(1))/2._dp)+ &
    G2(0,2)*((-34._dp/3._dp-(5*nf)/18._dp)*H1(1)-(268*H2(0,0))/9._dp-(28*H2(0,1))/9._dp+ &
    (62*H2(1,0))/9._dp+(-60*(1+8*v)+ &
    u*(5034-296*nf-144*v+7518*zeta2+3*(nf*nf))-144*(u*u))/ &
    (108._dp*u)+H1(0)*(-H1(1)+ &
    (10*u*(8-19*v)*v+(-200+(570-64*nf)*v)*(u*u*u)+ &
    (285-32*nf)*(u*u*u*u)+40*(v*v)+ &
    u*u*(35-390*v+(285-32*nf)*(v*v)))/(9._dp*(omw*omw)*(u*u)))- &
    (17*(H1(0)*H1(0)))/18._dp-(H1(1)*H1(1))/2._dp)+ &
    (4*H2(0,0)-5*H2(1,0)+(486*nf+9*(-764-20/u+30*zeta2)-5*(nf*nf))/ &
    &  72._dp)*(H1(1)*H1(1))+G1(0)*G1(0)* &
    (-G2(0,2)/2._dp-4*G2(1,2)-G2(2,0)/2._dp-4*G2(3,0)+4*G1(3)*H1(0)+ &
    G1(2)*((-17*H1(0))/2._dp-(9*H1(1))/2._dp)+((-30+nf)*H1(1))/18._dp+ &
    H1(0)*(-41._dp/3._dp+(14*nf)/9._dp+(7*H1(1))/2._dp)+G1(1)*(4*H1(0)+4*H1(1))+ &
    &  4*H2(0,0)-5*H2(1,0)+(486*nf+9*(-764-20/u+30*zeta2)-5*(nf*nf))/ &
    &  72._dp+(9*(G1(2)*G1(2)))/4._dp+(67*(H1(0)*H1(0)))/4._dp+(9*(H1(1)*H1(1)))/4._dp)+ &
    G1(2)*G1(2)*(-G2(0,2)/2._dp-4*G2(1,2)-G2(2,0)/2._dp-4*G2(3,0)+4*G1(3)*H1(0)+ &
    (-33._dp/2._dp+nf)*H1(1)+H1(0)*(-41._dp/3._dp+(14*nf)/9._dp+(43*H1(1))/2._dp)+ &
    &  4*H2(0,0)-5*H2(1,0)+(486*nf+9*(-764-20/u+30*zeta2)-5*(nf*nf))/ &
    &  72._dp+(67*(H1(0)*H1(0)))/4._dp+(27*(H1(1)*H1(1)))/4._dp)+ &
    H1(0)*H1(0)*((-75._dp/2._dp+3*nf)*H1(1)+(68*H2(0,0))/9._dp-(85*H2(1,0))/9._dp+ &
    (-1020+u*(-58041+4996*nf+306*zeta2-84*(nf*nf)))/(216._dp*u)+ &
    (35*(H1(1)*H1(1)))/4._dp)+G1(3)* &
    ((6*(-687+64*nf)*zeta2+480*zeta3)/27._dp+(58*zeta2*H1(1))/9._dp+ &
    (16*(-33+2*nf)*H2(0,0))/9._dp+(-49._dp/3._dp+(4*nf)/9._dp)*H2(0,1)+ &
    (-49._dp/3._dp+(4*nf)/9._dp)*H2(1,0)+(64*H3(0,0,0))/3._dp-20*H3(0,0,1)- &
    (100*H3(0,1,0))/9._dp-(100*H3(1,0,0))/9._dp-(10*H3(1,0,1))/9._dp- &
    (70*H3(1,1,0))/9._dp+(76._dp/3._dp-(8*nf)/9._dp+8*H1(1))*(H1(0)*H1(0))+ &
    (68*(H1(0)*H1(0)*H1(0)))/9._dp+ &
    H1(0)*((-4*(-30+nf)*H1(1))/9._dp+ &
    ((-4683+302*nf-216*v+1848*zeta2)*(u*u*u)-72*(u*u*u*u)- &
    &  2*(u*u)*(-6696+nf*(464-302*v)+4803*v+1848*omv*zeta2+ &
    &  108*(v*v))+u*(-7176+13872*v+1848*zeta2*(omv*omv)- &
    &  5163*(v*v)+nf*(464-928*v+302*(v*v))-72*(v*v*v))- &
    &  240*(1-2*(v*v)+v*v*v))/(54._dp*u*(w*w))+4*(H1(1)*H1(1))))+ &
    G1(2)*(G2(0,2)*(5._dp/3._dp-nf/18._dp+H1(0)+H1(1))+ &
    G2(2,0)*(5._dp/3._dp-nf/18._dp+H1(0)+H1(1))+ &
    G2(1,2)*((-4*(-30+nf))/9._dp+8*H1(0)+8*H1(1))+ &
    G2(3,0)*((-4*(-30+nf))/9._dp+8*H1(0)+8*H1(1))+ &
    ((107-4*nf)*H2(0,0))/3._dp+ &
    H1(1)*(111-(20*nf)/3._dp+85/(9._dp*u)-(619*zeta2)/18._dp-8*H2(0,0)+ &
    &  10*H2(1,0))-(74*H3(0,0,0))/9._dp+(10*H3(0,0,1))/9._dp+ &
    (20*H3(0,1,0))/9._dp-(80*H3(1,0,0))/9._dp-(250*H3(1,0,1))/9._dp- &
    (218*H3(1,1,0))/9._dp+(H2(1,0)* &
    ((163+(915-40*nf)*v)*(u*u*u)+(534-20*nf)*(u*u*u*u)- &
    &  45*(omv*omv)*(v*v)-9*u*v*(10-37*v+27*(v*v))+ &
    u*u*(-40+406*v+(183-20*nf)*(v*v))))/(9._dp*(omw*omw)*(u*u))+ &
    ((-8*nf*(69-42*v+6*(-1+2*v)*zeta2)+ &
    &  3*(1981-2496*v+6*(-350+578*v)*zeta2-406*zeta3+ &
    &  812*v*zeta3))*(u*u*u*u)+ &
    (-8*nf*(-21+6*zeta2)+6*(-618+714*zeta2+203*zeta3))* &
    (u*u*u*u*u)-6*omv*(-407+54*nf+90*omv*zeta2)*(v*v)+ &
    u*u*u*(-4*nf*(-177+v*(345-24*zeta2)+2*(-21+6*zeta2)*(v*v))+ &
    &  3*(-1559+v*(4504-812*zeta3)+2*(-642+203*zeta3)*(v*v)+ &
    &  6*zeta2*(137-852*v+472*(v*v))))+ &
    &  3*u*v*(1628-3187*v+542*(v*v)-4*nf*(54-113*v+23*(v*v))+ &
    &  36*zeta2*(-10+42*v-37*(v*v)+5*(v*v*v)))+ &
    u*u*(4*nf*(-81+435*v+2*(-138+6*zeta2)*(v*v))+ &
    &  3*(814-3932*v+(3065-406*zeta3)*(v*v)-24*(v*v*v)+ &
    &  6*zeta2*(-25+334*v-694*(v*v)+162*(v*v*v)))))/ &
    (54._dp*omu*(omw*omw)*(u*u))+ &
    (H2(0,1)*(24-4*nf-(20*w)/omw-(63*w)/u-(5*(w*w))/(omw*omw)+ &
    (45*(w*w))/(u*u)))/9._dp+ &
    G1(3)*(H1(0)*((4*(-30+nf))/9._dp-8*H1(1))-8*(H1(0)*H1(0)))+ &
    (73._dp/6._dp-(19*nf)/9._dp-(51*H1(1))/2._dp)*(H1(0)*H1(0))- &
    (289*(H1(0)*H1(0)*H1(0)))/18._dp+ &
    H1(0)*((122._dp/3._dp-(32*nf)/9._dp)*H1(1)-8*H2(0,0)+10*H2(1,0)+ &
    (12*omv*(407-54*nf)*v+ &
    (21642-1324*nf+288*v-21450*zeta2)*(u*u*u*u)+ &
    &  144*(u*u*u*u*u)+6*(u*u)* &
    (-1341+1498*v-18*nf*(-10+9*v)+80*(v*v))+ &
    u*u*u*(8544+22122*v-4*nf*(189+331*v)-21450*v*zeta2+ &
    &  144*(v*v))-6*u* &
    (-814+2165*v-74*(v*v)+36*nf*(3-8*v+v*v)))/ &
    (108._dp*omw*(u*u*u))-(35*(H1(1)*H1(1)))/2._dp)+ &
    (33._dp/2._dp-nf)*(H1(1)*H1(1))-(9*(H1(1)*H1(1)*H1(1)))/2._dp)+ &
    G1(0)*((omv*(-330-2352*u*zeta2+20*nf*(1+6*u*zeta2))+ &
    u*(1111-1188*v+nf*(-22+56*v+u*(56-120*zeta2))- &
    &  322*zeta3+322*v*zeta3+2*u*(-594+1176*zeta2+161*zeta3)))/ &
    (18._dp*u*w)+G2(0,2)*(5._dp/3._dp-nf/18._dp+H1(0))+ &
    G2(2,0)*(5._dp/3._dp-nf/18._dp+H1(0))+G2(1,2)*((-4*(-30+nf))/9._dp+8*H1(0))+ &
    G2(3,0)*((-4*(-30+nf))/9._dp+8*H1(0))-(32*(-33+2*nf)*H2(0,0))/9._dp+ &
    (43._dp/6._dp+nf/9._dp)*H2(0,1)+(-289._dp/6._dp+(37*nf)/9._dp)*H2(1,0)+ &
    ((29+nf)*H2(1,1))/3._dp-(64*H3(0,0,0))/3._dp+12*H3(0,0,1)+ &
    (28*H3(0,1,0))/9._dp+(268*H3(1,0,0))/9._dp+(100*H3(1,0,1))/9._dp+ &
    (10*H3(1,1,0))/9._dp+(H1(1)*(60*(1+8*v)+ &
    u*(-2334+116*nf+144*v-7518*zeta2)+144*(u*u)))/(108._dp*u)+ &
    ((-17*H1(0))/2._dp+H1(1)/2._dp)*(G1(2)*G1(2))+ &
    G1(3)*((4*(-30+nf)*H1(0))/9._dp-8*(H1(0)*H1(0)))+ &
    G1(1)*(H1(0)*((4*(-30+nf))/9._dp-8*H1(1))+(4*(-30+nf)*H1(1))/9._dp- &
    &  8*(H1(0)*H1(0)))+(73._dp/6._dp-(19*nf)/9._dp-(289*H1(1))/18._dp)* &
    (H1(0)*H1(0))-(289*(H1(0)*H1(0)*H1(0)))/18._dp+ &
    H1(0)*(((-61+2*nf)*H1(1))/6._dp-8*H2(0,0)+10*H2(1,0)+ &
    ((21354-1324*nf+432*v-4578*zeta2)*(u*u*u)+144*(u*u*u*u)+ &
    &  2*(u*u)*(-4*nf*(-412+331*v)+4578*omv*zeta2+ &
    &  6*(-4185+3599*v+36*(v*v)))+ &
    &  60*(17-18*v-7*(v*v)+8*(v*v*v))+ &
    u*(-4578*zeta2*(omv*omv)-4*nf*(412-824*v+331*(v*v))+ &
    &  6*(4210-8440*v+3719*(v*v)+24*(v*v*v))))/(108._dp*u*(w*w))- &
    (15*(H1(1)*H1(1)))/2._dp)+ &
    G1(2)*(((-30+nf)*H1(1))/18._dp+H1(0)*((-8*(-30+nf))/9._dp+16*H1(1))- &
    ((-30+nf)*(-30+nf))/36._dp+25*(H1(0)*H1(0))-H1(1)*H1(1))+ &
    (5._dp/3._dp-nf/18._dp)*(H1(1)*H1(1))+(H1(1)*H1(1)*H1(1))/2._dp)+ &
    ((-33+2*nf)*(H1(1)*H1(1)*H1(1)))/6._dp+ &
    G1(1)*((-88+(16*nf)/3._dp)*H2(0,0)+(-71+(28*nf)/9._dp+3/u)*H2(1,0)+ &
    (154*H3(0,0,0))/9._dp+16*H3(0,0,1)+(134*H3(0,1,0))/9._dp+ &
    (154*H3(1,0,0))/9._dp+(80*H3(1,0,1))/9._dp+22*H3(1,1,0)+ &
    (H2(1,1)*(-40+190*u+6*(-109+4*nf)*(u*u)))/(9._dp*(u*u))+ &
    (H2(0,1)*(-40+163*u+(-807+44*nf)*(u*u)))/(9._dp*(u*u))+ &
    (-114*zeta3*(u*u)+6*zeta2*(100-421*u+(-507+56*nf)*(u*u)))/ &
    (27._dp*(u*u))+(H1(1)*(6*omv*(-407+54*nf)*v+ &
    (-4827+302*nf-144*v+2076*zeta2)*(u*u*u*u)-72*(u*u*u*u*u)+ &
    u*u*u*(-4002-4797*v+nf*(378+302*v)+2076*v*zeta2- &
    &  72*(v*v))+u*u*(4023-4254*v+54*nf*(-10+9*v)+ &
    &  30*(v*v))+3*u*(-814+2165*v-74*(v*v)+ &
    &  36*nf*(3-8*v+v*v))))/(54._dp*omw*(u*u*u))+ &
    (4*H1(0)+4*H1(1))*(G1(2)*G1(2))+ &
    (76._dp/3._dp-(8*nf)/9._dp+(140*H1(1))/9._dp)*(H1(0)*H1(0))+ &
    (68*(H1(0)*H1(0)*H1(0)))/9._dp+ &
    G1(2)*(H1(0)*((4*(-30+nf))/9._dp-16*H1(1))+(4*(-30+nf)*H1(1))/9._dp- &
    &  8*(H1(0)*H1(0))-8*(H1(1)*H1(1)))-(4*(-30+nf)*(H1(1)*H1(1)))/9._dp+ &
    H1(0)*((-4*(-29+nf)*H1(1))/3._dp+ &
    (6*omv*(-407+54*nf)*v+ &
    (-4827+302*nf-144*v+2844*zeta2)*(u*u*u*u)- &
    &  72*(u*u*u*u*u)+u*u*u* &
    (-4002-4797*v+nf*(378+302*v)+2844*v*zeta2-72*(v*v))+ &
    u*u*(4023-4254*v+54*nf*(-10+9*v)+30*(v*v))+ &
    &  3*u*(-814+2165*v-74*(v*v)+36*nf*(3-8*v+v*v)))/ &
    (54._dp*omw*(u*u*u))+12*(H1(1)*H1(1)))+4*(H1(1)*H1(1)*H1(1)))+ &
    H1(0)*((76._dp/3._dp-(8*nf)/9._dp)*H2(0,0)+(5*(-57+2*nf)*H2(1,0))/9._dp+ &
    H1(1)*(8*H2(0,0)-10*H2(1,0)+ &
    (-180+u*(-8946+705*nf+294*zeta2-10*(nf*nf)))/(36._dp*u))+ &
    (3*(omv*omv*omv)*(4884*v-648*nf*v-11862*zeta2*(u*u)+ &
    &  1356*nf*zeta2*(u*u)+160*(nf*nf)*(u*u))+ &
    &  3*u*(omv*omv)*(160*(-2+u+v)*(nf*nf)*(u*u)+ &
    &  2*nf*(678*(-2+u+v)*zeta2*(u*u)-6*(20-147*v+46*(v*v)))+ &
    &  3*(-3954*(-2+u+v)*zeta2*(u*u)+2*(540-2303*v+542*(v*v)))) &
    +2*omv*(u*u)*(-2*nf*(4465-8291*v+2800*(v*v))+ &
    &  3*(14570+4758*zeta3-v*(24601+9516*zeta3)+ &
    (7367+4758*zeta3)*(v*v)+72*(v*v*v)))- &
    &  4*(u*u*u)*(nf*(-9470+23441*v-16582*(v*v)+ &
    u*(4825-9488*v+3637*(v*v))+3637*(v*v*v))- &
    &  3*(-2*(8500+2379*zeta3)+v*(39629+11895*zeta3)- &
    &  2*(13643+4758*zeta3)*(v*v)+ &
    u*(8905+2379*zeta3-2*v*(8095+2379*zeta3)+ &
    (5989+2379*zeta3)*(v*v))+(5953+2379*zeta3)*(v*v*v))))/ &
    (324._dp*omu*w*(omv*omv)*(u*u))+(-27+2*nf)*(H1(1)*H1(1))+ &
    (9*(H1(1)*H1(1)*H1(1)))/2._dp)+(9*(H1(1)*H1(1)*H1(1)*H1(1)))/8._dp
    return
    end function Hgamma_2a2re
          
          
    function Hgamma_2a2im(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    use constants
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=(64*pi*G3(0,0,1))/9._dp-(118*pi*G3(0,0,2))/9._dp+(16*pi*G3(0,0,3))/9._dp+ &
    (154*pi*G3(0,1,0))/9._dp-(80*pi*G3(0,1,1))/9._dp+30*pi*G3(0,1,2)- &
    (118*pi*G3(0,2,0))/9._dp+(88*pi*G3(0,2,1))/9._dp-(118*pi*G3(0,2,2))/9._dp+ &
    (88*pi*G3(0,2,3))/9._dp-8*pi*G3(0,3,0)-8*pi*G3(0,3,2)+ &
    (128*pi*G3(0,3,3))/9._dp+(80*pi*G3(1,0,1))/9._dp+4*pi*G3(1,0,2)+ &
    (10*pi*G3(1,0,3))/9._dp-(128*pi*G3(1,1,1))/9._dp+(128*pi*G3(1,1,2))/9._dp+ &
    &  4*pi*G3(1,2,0)+(16*pi*G3(1,2,1))/3._dp+4*pi*G3(1,2,2)+ &
    (10*pi*G3(1,2,3))/9._dp-(118*pi*G3(2,0,0))/9._dp+(88*pi*G3(2,0,1))/9._dp- &
    (118*pi*G3(2,0,2))/9._dp+(88*pi*G3(2,0,3))/9._dp-8*pi*G3(2,1,0)+ &
    (128*pi*G3(2,1,1))/9._dp-8*pi*G3(2,1,2)-(118*pi*G3(2,2,0))/9._dp+ &
    (16*pi*G3(2,2,1))/9._dp+(64*pi*G3(2,2,3))/9._dp+30*pi*G3(2,3,0)+ &
    (154*pi*G3(2,3,2))/9._dp-(80*pi*G3(2,3,3))/9._dp+4*pi*G3(3,0,0)+ &
    (10*pi*G3(3,0,1))/9._dp+4*pi*G3(3,0,2)+(16*pi*G3(3,0,3))/3._dp+ &
    &  4*pi*G3(3,2,0)+(10*pi*G3(3,2,1))/9._dp+(80*pi*G3(3,2,3))/9._dp+ &
    (128*pi*G3(3,3,0))/9._dp-(128*pi*G3(3,3,3))/9._dp+ &
    G2(3,3)*((8*(-33+2*nf)*pi)/9._dp-(128*pi*H1(0))/9._dp)+ &
    G2(3,2)*(((-147+4*nf)*pi)/9._dp-(38*pi*H1(0))/9._dp)+ &
    G2(2,2)*(((29+nf)*pi)/3._dp+(28*pi*H1(0))/3._dp)+ &
    G2(2,3)*((39-(32*nf)/9._dp)*pi-28*pi*H1(0)-(154*pi*H1(1))/9._dp)+ &
    G2(1,0)*((4*(-30+nf)*pi)/9._dp-12*pi*H1(0)-4*pi*H1(1))+ &
    G2(3,0)*(((7-4*nf)*pi)/3._dp+(28*pi*H1(0))/9._dp-4*pi*H1(1))+ &
    G2(0,3)*((4*(-96+5*nf)*pi)/9._dp+(136*pi*H1(0))/9._dp+8*pi*H1(1))+ &
    G2(0,0)*(((-21+nf)*pi)/3._dp+4*pi*H1(0)+(118*pi*H1(1))/9._dp)+ &
    G2(0,2)*(((108+nf)*pi)/9._dp+(379*pi*H1(0))/9._dp+(118*pi*H1(1))/9._dp)+ &
    G2(2,0)*(((108+nf)*pi)/9._dp+(379*pi*H1(0))/9._dp+(118*pi*H1(1))/9._dp)+ &
    ((651-32*nf)*pi*H2(0,0))/9._dp+((129+2*nf)*pi*H2(0,1))/18._dp+ &
    ((29+nf)*pi*H2(1,1))/3._dp-(4*pi*H3(0,0,0))/3._dp+12*pi*H3(0,0,1)+ &
    (38*pi*H3(0,1,0))/9._dp+(278*pi*H3(1,0,0))/9._dp+(100*pi*H3(1,0,1))/9._dp- &
    (80*pi*H3(1,1,0))/3._dp+(pi*H1(1)* &
    (60*(9+8*v)+u*(-2334+116*nf+144*v-5310*zeta2)+144*(u*u)))/ &
    (108._dp*u)+G2(1,1)*((-128*pi*H1(1))/9._dp+ &
    (8*pi*(-5+17*u+(-21+2*nf)*(u*u)))/(9._dp*(u*u)))+ &
    G2(1,2)*((-10*pi*H1(0))/9._dp-4*pi*H1(1)- &
    (2*pi*(20-95*u+(51+6*nf)*(u*u)))/(9._dp*(u*u)))- &
    (pi*H2(1,0)*(10+20*v+2*u*(10+3*(-89+2*nf)*v)+ &
    &  3*(-89+2*nf)*(u*u)+3*(-89+2*nf)*(v*v)))/(18._dp*(omw*omw))+ &
    G2(2,1)*((8*pi*H1(0))/9._dp+8*pi*H1(1)+ &
    (pi*(2*v*(-13+2*(-39+5*nf)*v)*(u*u)+ &
    (27+(-861+40*nf)*v)*(u*u*u)+(-507+20*nf)*(u*u*u*u)+ &
    &  5*(v*v)*(1-18*v+9*(v*v))+u*v*(10-143*v+243*(v*v))))/ &
    (9._dp*(omw*omw)*(u*u)))+ &
    G2(0,1)*((-172*pi*H1(0))/9._dp-30*pi*H1(1)+ &
    (pi*(10*u*(8-19*v)*v-4*(50+(-189+16*nf)*v)*(u*u*u)+ &
    (378-32*nf)*(u*u*u*u)+40*(v*v)+ &
    u*u*(35-390*v+(378-32*nf)*(v*v))))/(9._dp*(omw*omw)*(u*u)))+ &
    (pi*(-36*(-407+54*nf)*(omv*omv)+ &
    &  18*omv*(-1041+58*nf)*zeta2*(u*u)+ &
    u*(480*omu*u*w*(nf*nf)+ &
    &  2*nf*(2*u*(-6544+6229*v)+ &
    (21802-10658*v+522*(-2+v)*zeta2)*(u*u)+ &
    (-10658+522*zeta2)*(u*u*u)+36*(81-104*v+23*(v*v)))- &
    &  3*((101096-50308*v+6246*(-2+v)*zeta2+16824*zeta3- &
    &  8412*v*zeta3)*(u*u)+ &
    (-50164+6246*zeta2-8412*zeta3)*(u*u*u)- &
    &  2*u*(29513+4206*zeta3-v*(28385+4206*zeta3)+72*(v*v))+ &
    &  6*(2163-2705*v+542*(v*v))))))/(324._dp*omu*w*(u*u))+ &
    (((-159+17*nf)*pi)/18._dp-4*pi*G1(1)-(pi*G1(2))/2._dp-4*pi*G1(3)- &
    (7*pi*H1(0))/2._dp+(pi*H1(1))/2._dp)*(G1(0)*G1(0))-(pi*(G1(0)*G1(0)*G1(0)))/2._dp+ &
    (((-159+17*nf)*pi)/18._dp-4*pi*G1(3)-(7*pi*H1(0))/2._dp+(3*pi*H1(1))/2._dp)* &
    (G1(2)*G1(2))-(pi*(G1(2)*G1(2)*G1(2)))/2._dp+ &
    ((-41+(22*nf)/9._dp)*pi-(145*pi*H1(1))/18._dp)*(H1(0)*H1(0))- &
    (17*pi*(H1(0)*H1(0)*H1(0)))/2._dp+((-159+17*nf)*pi*(H1(1)*H1(1)))/18._dp+ &
    G1(3)*(3*pi*H1(1)+H1(0)*(((-237+16*nf)*pi)/9._dp-8*pi*H1(1))- &
    (4*pi*H2(0,0))/3._dp+12*pi*H2(0,1)+(38*pi*H2(1,0))/9._dp+ &
    (pi*((4683-302*nf+216*v-1200*zeta2)*(u*u*u)+72*(u*u*u*u)+ &
    &  2*(u*u)*(-6696+nf*(464-302*v)+4803*v+1200*omv*zeta2+ &
    &  108*(v*v))+240*(1-2*(v*v)+v*v*v)+ &
    u*(7176-13872*v-1200*zeta2*(omv*omv)+ &
    nf*(-464+928*v-302*(v*v))+5163*(v*v)+72*(v*v*v))))/ &
    (54._dp*u*(w*w))-(76*pi*(H1(0)*H1(0)))/9._dp-4*pi*(H1(1)*H1(1)))+ &
    G1(1)*(G1(2)*((-4*(-30+nf)*pi)/9._dp+8*pi*H1(0)+8*pi*H1(1))+ &
    (26*pi*H2(0,0))/9._dp+12*pi*H2(0,1)+2*pi*H2(1,0)+4*pi*H2(1,1)+ &
    (2*pi*H1(1)*(20-95*u+(-9+8*nf)*(u*u)))/(9._dp*(u*u))+ &
    H1(0)*((-80*pi*H1(1))/9._dp+(pi*(40-163*u+27*(u*u)))/(9._dp*(u*u)))- &
    (pi*(6*omv*(-407+54*nf)*v+ &
    (-4827+302*nf-144*v+1320*zeta2)*(u*u*u*u)-72*(u*u*u*u*u)+ &
    u*u*u*(-4002-4797*v+nf*(378+302*v)+1320*v*zeta2- &
    &  72*(v*v))+u*u*(4023-4254*v+54*nf*(-10+9*v)+ &
    &  30*(v*v))+3*u*(-814+2165*v-74*(v*v)+ &
    &  36*nf*(3-8*v+v*v))))/(54._dp*omw*(u*u*u))- &
    &  4*pi*(G1(2)*G1(2))-(76*pi*(H1(0)*H1(0)))/9._dp-4*pi*(H1(1)*H1(1)))+ &
    H1(0)*(((-187+14*nf)*pi*H1(1))/6._dp-(8*pi*H2(0,0))/9._dp+(10*pi*H2(1,0))/9._dp+ &
    (pi*((24306-2988*nf+576*v-2490*zeta2+96*(nf*nf))*(u*u*u*u)+ &
    &  144*(u*u*u*u*u)+u*u*u* &
    (-56604+nf*(6624-8964*v)+73938*v-2490*(-2+3*v)*zeta2+ &
    &  96*(-2+3*v)*(nf*nf)+864*(v*v))+ &
    &  60*v*(8+9*v-34*(v*v)+17*(v*v*v))+ &
    u*u*(-2490*zeta2*(1-4*v+3*(v*v))+ &
    &  96*(nf*nf)*(1-4*v+3*(v*v))- &
    &  36*nf*(92-368*v+249*(v*v))+ &
    &  6*(4862-19208*v+12663*(v*v)+96*(v*v*v)))+ &
    u*(540+v*(29712-3312*nf-2490*zeta2+96*(nf*nf))+ &
    (-60684+6624*nf+4980*zeta2-192*(nf*nf))*(v*v)+ &
    (27366-2988*nf-2490*zeta2+96*(nf*nf))*(v*v*v)+ &
    &  144*(v*v*v*v))))/(108._dp*omw*u*(w*w))-(7*pi*(H1(1)*H1(1)))/2._dp)+ &
    G1(2)*((8-(20*nf)/9._dp)*pi*H1(1)+ &
    G1(3)*((-4*(-30+nf)*pi)/9._dp+8*pi*H1(0)+8*pi*H1(1))- &
    (278*pi*H2(0,0))/9._dp-(262*pi*H2(0,1))/9._dp+(26*pi*H2(1,0))/3._dp+ &
    (pi*(-60*(9+8*v)+u*(2334-116*nf-144*v+5310*zeta2)- &
    &  144*(u*u)))/(108._dp*u)+ &
    H1(0)*(7*pi*H1(1)+(pi*((163+(2085-180*nf)*v)*(u*u*u)+ &
    (1119-90*nf)*(u*u*u*u)-45*(omv*omv)*(v*v)- &
    &  9*u*v*(10-37*v+27*(v*v))+ &
    u*u*(-40+406*v+(768-90*nf)*(v*v))))/(9._dp*(omw*omw)*(u*u)) &
    )+(161*pi*(H1(0)*H1(0)))/18._dp-(3*pi*(H1(1)*H1(1)))/2._dp)+ &
    G1(0)*(G1(1)*((-4*(-30+nf)*pi)/9._dp+8*pi*H1(0))+ &
    G1(3)*((-4*(-30+nf)*pi)/9._dp+8*pi*H1(0))-(46*pi*H1(1))/3._dp+ &
    H1(0)*((2*(192-19*nf)*pi)/9._dp-(19*pi*H1(1))/9._dp)+ &
    G1(2)*(-((-30+nf)*pi)/9._dp+2*pi*H1(0)+pi*H1(1))+(20*pi*H2(0,0))/9._dp- &
    &  4*pi*H2(0,1)-24*pi*H2(1,0)-(118*pi*H2(1,1))/9._dp- &
    (pi*(60*(1+8*v)+u*(-2334+116*nf+144*v+930*zeta2)+ &
    &  144*(u*u)))/(108._dp*u)-(pi*(G1(2)*G1(2)))/2._dp+ &
    (161*pi*(H1(0)*H1(0)))/18._dp-(pi*(H1(1)*H1(1)))/2._dp)+ &
    (pi*(H1(1)*H1(1)*H1(1)))/2._dp
    return
    end function Hgamma_2a2im
          
          
!**********************************************************************
! egion3a:qbar+g-->V+qbar*
!*********************************************************************/
          
    function Hgamma_3a0re(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=1
    return
    end function Hgamma_3a0re
          
          
    function Hgamma_3a0im(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=0
    return
    end function Hgamma_3a0im
          
          
    function Hgamma_3a1re(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    use constants
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=(588-40*nf-(60*omv)/u+846*zeta2)/36._dp-3*G2(0,2)+6*G2(1,2)- &
    &  3*G2(2,0)+(8*G2(3,0))/3._dp-(8*G1(3)*H1(0))/3._dp+G1(1)*(-6*H1(0)-6*H1(1))+ &
    G1(2)*(-7+(2*nf)/3._dp+(17*H1(0))/3._dp-H1(1)/3._dp)+ &
    H1(0)*(7-(2*nf)/3._dp+H1(1)/3._dp)+(7-(2*nf)/3._dp)*H1(1)+ &
    G1(0)*((17*H1(0))/3._dp+3*H1(1))-(8*H2(0,0))/3._dp-(10*H2(1,0))/3._dp- &
    (3*(G1(0)*G1(0)))/2._dp+(G1(2)*G1(2))/6._dp-(17*(H1(0)*H1(0)))/6._dp+(H1(1)*H1(1))/6._dp
    return
    end function Hgamma_3a1re
          
          
    function Hgamma_3a1im(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    use constants
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=-3*pi*G1(0)+6*pi*G1(1)-3*pi*G1(2)+(8*pi*G1(3))/3._dp-(pi*H1(0))/3._dp+3*pi*H1(1)
    return
    end function Hgamma_3a1im
          
          
    function Hgamma_3a2re(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    use constants
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=158795._dp/162._dp-1352/(9._dp*u)+(48*v)/omu+(1352*v)/(9._dp*u)- &
    (1237*zeta3)/18._dp+(75493*zeta4)/144._dp+(72-6*nf)*G3(2,2,0)- &
    (16*G4(0,0,1,2))/9._dp+18*G4(0,0,2,2)+(16*G4(0,0,3,0))/9._dp+26*G4(0,1,0,2)+ &
    &  26*G4(0,1,2,0)-36*G4(0,1,2,2)+18*G4(0,2,0,2)-(92*G4(0,2,1,2))/9._dp+ &
    &  18*G4(0,2,2,0)-(152*G4(0,2,3,0))/9._dp+(8*G4(0,3,0,2))/9._dp+ &
    (8*G4(0,3,2,0))/9._dp+(128*G4(0,3,3,0))/9._dp-36*G4(1,0,2,2)-72*G4(1,1,1,2)+ &
    &  72*G4(1,1,2,2)-36*G4(1,2,0,2)+72*G4(1,2,1,2)-36*G4(1,2,2,0)+ &
    &  18*G4(2,0,0,2)-18*G4(2,0,1,2)+18*G4(2,0,2,0)-18*G4(2,0,3,0)- &
    &  18*G4(2,1,0,2)+72*G4(2,1,1,2)-18*G4(2,1,2,0)+18*G4(2,2,0,0)- &
    &  36*G4(2,2,1,2)+36*G4(2,2,3,0)-(244*G4(2,3,0,0))/9._dp+ &
    (164*G4(2,3,0,2))/9._dp+(164*G4(2,3,2,0))/9._dp+(80*G4(2,3,3,0))/9._dp- &
    (164*G4(3,0,0,2))/9._dp-(70*G4(3,0,1,2))/9._dp-(164*G4(3,0,2,0))/9._dp+ &
    (16*G4(3,0,3,0))/3._dp-(164*G4(3,2,0,0))/9._dp-(70*G4(3,2,1,2))/9._dp- &
    (80*G4(3,2,3,0))/9._dp+(128*G4(3,3,0,0))/9._dp-(128*G4(3,3,3,0))/9._dp+ &
    G3(0,2,2)*(72-6*nf-36*H1(0))+G3(2,0,2)*(72-6*nf-(160*H1(0))/9._dp)+ &
    G3(3,3,0)*((8*(-33+2*nf))/9._dp-(128*H1(0))/9._dp)-(16*G3(0,0,3)*H1(0))/9._dp+ &
    (152*G3(0,2,3)*H1(0))/9._dp-(8*G3(0,3,2)*H1(0))/9._dp-(128*G3(0,3,3)*H1(0))/9._dp+ &
    &  18*G3(2,0,3)*H1(0)-36*G3(2,2,3)*H1(0)-(164*G3(2,3,2)*H1(0))/9._dp- &
    (80*G3(2,3,3)*H1(0))/9._dp-(16*G3(3,0,3)*H1(0))/3._dp+(80*G3(3,2,3)*H1(0))/9._dp+ &
    (128*G3(3,3,3)*H1(0))/9._dp+G3(3,2,0)*(15+(92*H1(0))/9._dp)+ &
    G3(3,0,2)*(15+(34*H1(0))/3._dp)+G3(1,2,2)*(12*(-12+nf)+72*H1(0))+ &
    G3(1,2,1)*(-72*H1(0)-72*H1(1))+G3(2,1,1)*(-72*H1(0)-72*H1(1))+ &
    G3(0,1,0)*(-26*H1(0)-26*H1(1))+ &
    G3(2,3,0)*(-121._dp/3._dp+4*nf+18*H1(0)-(164*H1(1))/9._dp)+ &
    G3(2,0,0)*(63-3*nf-(244*H1(0))/9._dp-18*H1(1))+ &
    G3(0,0,2)*(63-3*nf-20*H1(0)-18*H1(1))+ &
    G3(0,2,0)*(63-3*nf-10*H1(0)-18*H1(1))+ &
    G3(0,3,0)*((4*(-96+5*nf))/9._dp+(64*H1(0))/9._dp-(8*H1(1))/9._dp)+ &
    G3(0,0,1)*((16*H1(0))/9._dp+(16*H1(1))/9._dp)+ &
    G3(3,0,1)*((70*H1(0))/9._dp+(70*H1(1))/9._dp)+ &
    G3(3,2,1)*((70*H1(0))/9._dp+(70*H1(1))/9._dp)+ &
    G3(0,2,1)*((92*H1(0))/9._dp+(92*H1(1))/9._dp)+G3(2,0,1)*(18*H1(0)+18*H1(1))+ &
    G3(2,1,0)*(18*H1(0)+18*H1(1))+ &
    G3(3,0,0)*((8*(-21+nf))/3._dp+(128*H1(0))/9._dp+(164*H1(1))/9._dp)+ &
    G3(1,0,2)*(-30+nf+18*H1(0)+36*H1(1))+ &
    G3(1,2,0)*(-30+nf+18*H1(0)+36*H1(1))+G3(2,2,1)*(36*H1(0)+36*H1(1))+ &
    G3(1,1,1)*(72*H1(0)+72*H1(1))+ &
    G2(3,3)*((-256*zeta2)/9._dp+(8*(33-2*nf)*H1(0))/9._dp+(128*H2(0,0))/9._dp)+ &
    G2(3,2)*((622*zeta2)/9._dp-15*H1(0)-(100*H2(0,0))/9._dp-(70*H2(0,1))/9._dp- &
    (10*H2(1,0))/9._dp)+G2(0,3)*((-176*zeta2)/9._dp-(4*(-96+5*nf)*H1(0))/9._dp- &
    (128*H2(0,0))/9._dp+(8*H2(0,1))/9._dp+(8*H2(1,0))/9._dp)+ &
    G2(2,3)*((-332*zeta2)/9._dp+(121._dp/3._dp-4*nf)*H1(0)-(80*H2(0,0))/9._dp+ &
    (164*H2(0,1))/9._dp+(164*H2(1,0))/9._dp)+ &
    G2(1,0)*(108*zeta2+(30-nf)*H1(0)+(30-nf)*H1(1)-18*H2(0,1)- &
    &  18*H2(1,0)-36*H2(1,1))+(-132+8*nf)*H3(0,0,0)+15*H3(0,0,1)+ &
    (10*(-69+5*nf)*H3(0,1,0))/9._dp+((-451+32*nf)*H3(1,0,0))/3._dp+ &
    (64*H4(0,0,0,0))/3._dp+20*H4(0,0,0,1)+(100*H4(0,0,1,0))/9._dp+ &
    (100*H4(0,1,0,0))/9._dp+(70*H4(0,1,0,1))/9._dp+(10*H4(0,1,1,0))/9._dp+ &
    (250*H4(1,0,0,0))/9._dp-(20*H4(1,0,0,1))/9._dp-(10*H4(1,0,1,0))/3._dp- &
    (170*H4(1,1,0,0))/9._dp-(230*H4(1,1,0,1))/9._dp+40*H4(1,1,1,0)+ &
    (H3(1,1,0)*(-915+20/omw+56*nf-5/(omw*omw)))/9._dp+ &
    (H3(1,0,1)*(5-20*omw+123*(omw*omw)))/(9._dp*(omw*omw))+ &
    G2(0,0)*(80-(41*nf)/6._dp-(406*zeta2)/9._dp+(8*(-21+nf)*H1(0))/3._dp+ &
    &  3*(-21+nf)*H1(1)+(128*H2(0,0))/9._dp+(164*H2(0,1))/9._dp+ &
    (244*H2(1,0))/9._dp+18*H2(1,1)+(5*(nf*nf))/36._dp)+ &
    H2(1,1)*(455._dp/4._dp-(61*nf)/3._dp+(104*zeta2)/9._dp+(8*(nf*nf))/9._dp)+ &
    G2(2,2)*(455._dp/4._dp-(61*nf)/3._dp+18*zeta2+12*(-12+nf)*H1(0)-36*H2(0,1)+ &
    &  36*H2(1,0)+(8*(nf*nf))/9._dp)+(100*(nf*nf))/81._dp+ &
    (H2(0,1)*(-3150*nf+15*(1287-4/omw+(4+32*v)/u)-7464*zeta2+ &
    &  120*(nf*nf)))/108._dp+G2(0,1)* &
    (-48*zeta2-18*H2(0,1)-18*H2(1,0)-36*H2(1,1)+ &
    (H1(0)*(-9*nf+5*(87-4/omw-(22*omv)/u+1/(omw*omw)- &
    (8*(omv*omv))/(u*u))))/9._dp+ &
    (H1(1)*(-9*nf+5*(87-4/omw-(22*omv)/u+1/(omw*omw)- &
    (8*(omv*omv))/(u*u))))/9._dp)+ &
    G3(0,1,2)*(18*H1(0)+36*H1(1)+ &
    (9*nf+5*(-87+4/omw+(22*omv)/u-1/(omw*omw)+ &
    (8*(omv*omv))/(u*u)))/9._dp)+ &
    G3(2,1,2)*(36*H1(0)+(-1431+20/omw+108*nf-(83*omv)/u- &
    &  5/(omw*omw)-(40*(omv*omv))/(u*u))/9._dp)+ &
    G2(2,1)*(-36*zeta2+36*H2(0,1)-36*H2(1,0)+ &
    (H1(0)*(1431-20/omw-108*nf+(83*omv)/u+5/(omw*omw)+ &
    (40*(omv*omv))/(u*u)))/9._dp+ &
    (H1(1)*(1431-20/omw-108*nf+(83*omv)/u+5/(omw*omw)+ &
    (40*(omv*omv))/(u*u)))/9._dp)+ &
    G3(1,1,2)*(42-4*nf+14/u-72*H1(0)-72*H1(1)+10/(u*u))+ &
    G2(1,1)*(-72*zeta2+72*H2(1,0)+72*H2(1,1)+ &
    H1(0)*(4*nf-(2*(5+7*u+21*(u*u)))/(u*u))+ &
    H1(1)*(4*nf-(2*(5+7*u+21*(u*u)))/(u*u)))+ &
    nf*((-59573+(6336*omv)/u-(6156*v)/omu-10602*zeta3)/486._dp+ &
    &  6*zeta2*(-359._dp/36._dp+4*u*v+2*(u*u)+2*(v*v)))+ &
    (zeta2*(38753+60/omw-4928*v-16*u*(268+407*v)-48/w- &
    &  3256*(u*u)-3256*(v*v)-(10*(185-193*v+64*(v*v)))/u))/36._dp+ &
    (H2(0,0)*(-21648*zeta2+180*(nf*nf)- &
    &  4*nf*(1253+162*u*v+81*(u*u)+81*(v*v))+ &
    &  3*(11023+1232*v+4*u*(268+407*v)+48/w+814*(u*u)+ &
    &  814*(v*v)+(160*(1-3*v+v*v))/u)))/108._dp+ &
    (H2(1,0)*(-17724*zeta2+120*(nf*nf)- &
    &  2*nf*(995+324*u*v+162*(u*u)+162*(v*v))+ &
    &  3*(-137+20/omw+1232*v+4*u*(268+407*v)-48/w+814*(u*u)+ &
    &  814*(v*v)+(20*(-3-6*v+8*(v*v)))/u)))/108._dp+ &
    H1(1)*((-4*H2(0,0))/3._dp-(5*H2(1,0))/3._dp+(40*(nf*nf))/27._dp+ &
    nf*(-5086._dp/81._dp-3*u-6*v+(2*v)/omu-(107*zeta2)/9._dp- &
    (6*(omv*omv))/(u*u)+(38*(v*v))/(3._dp*(omu*omu))- &
    (3*(-2+v+v*v))/u)+ &
    (49396+2442*u+4884*v-(6480*v)/omu-2628*zeta3+ &
    (90*zeta2*(-2+8*omw+113*(omw*omw)))/(omw*omw)+ &
    (4884*(omv*omv))/(u*u)-(5184*(v*v))/(omu*omu)+ &
    (6*(-1093+686*v+407*(v*v)))/u)/108._dp)+ &
    ((33-2*nf)/6._dp-(17*H1(0))/2._dp-(9*H1(1))/2._dp)*(G1(0)*G1(0)*G1(0))+ &
    (9*(G1(0)*G1(0)*G1(0)*G1(0)))/8._dp+ &
    ((4*(-12+nf))/27._dp+(17*H1(0))/18._dp-H1(1)/18._dp)*(G1(2)*G1(2)*G1(2))+ &
    (G1(2)*G1(2)*G1(2)*G1(2))/72._dp+((68*(-12+nf))/27._dp-(17*H1(1))/18._dp)* &
    (H1(0)*H1(0)*H1(0))+(289*(H1(0)*H1(0)*H1(0)*H1(0)))/72._dp+ &
    G1(3)*((6*(-1137+64*nf)*zeta2+480*zeta3)/27._dp-(622*zeta2*H1(1))/9._dp+ &
    (16*(-33+2*nf)*H2(0,0))/9._dp+15*H2(0,1)+15*H2(1,0)+ &
    (64*H3(0,0,0))/3._dp+20*H3(0,0,1)+(100*H3(0,1,0))/9._dp+ &
    (100*H3(1,0,0))/9._dp+(70*H3(1,0,1))/9._dp+(10*H3(1,1,0))/9._dp+ &
    (76._dp/3._dp-(8*nf)/9._dp-(8*H1(1))/9._dp)*(H1(0)*H1(0))+ &
    (68*(H1(0)*H1(0)*H1(0)))/9._dp+ &
    H1(0)*((-4*H1(1))/3._dp+(-8232*zeta2- &
    &  2*nf*(-232+162*u*v+81*(u*u)+81*(v*v))+ &
    &  3*(-2656+616*v+u*(536+814*v)+24/w+407*(u*u)+ &
    &  407*(v*v)+(80*(1-3*v+v*v))/u))/54._dp-(4*(H1(1)*H1(1)))/9._dp)) &
    +G2(2,0)*(H1(0)*(-373._dp/3._dp+9*nf-H1(1))+(-147._dp/2._dp+6*nf)*H1(1)+ &
    &  18*H2(0,0)-(2*H2(0,1))/9._dp-(82*H2(1,0))/9._dp+ &
    ((3216+(4884-648*nf)*v)*(u*u)-6*(-407+54*nf)*(u*u*u)+ &
    &  60*(7-15*v+8*(v*v))+ &
    u*(-14556+3696*v-13542*zeta2+12*(nf*nf)+2442*(v*v)- &
    &  81*nf*(-7+4*(v*v))))/(108._dp*u)+(17*(H1(0)*H1(0)))/2._dp- &
    (H1(1)*H1(1))/2._dp)+G2(0,2)*((-147._dp/2._dp+6*nf)*H1(1)+(332*H2(0,0))/9._dp+ &
    (232*H2(0,1))/9._dp+(242*H2(1,0))/9._dp+ &
    ((3216+(4884-648*nf)*v)*(u*u)-6*(-407+54*nf)*(u*u*u)+ &
    &  60*(7-15*v+8*(v*v))+ &
    u*(-14556+3696*v-12462*zeta2+12*(nf*nf)+2442*(v*v)- &
    &  81*nf*(-7+4*(v*v))))/(108._dp*u)+ &
    H1(0)*(-H1(1)+(2*(65+(-841+45*nf)*v)*(u*u*u)+ &
    (-786+45*nf)*(u*u*u*u)+40*(omv*omv)*(v*v)- &
    &  10*u*v*(-8+5*v+3*(v*v))+ &
    u*u*(35+160*v+(-966+45*nf)*(v*v)))/(9._dp*(omw*omw)*(u*u))) &
    +(17*(H1(0)*H1(0)))/2._dp-(H1(1)*H1(1))/2._dp)+ &
    H1(0)*H1(0)*((-47._dp/3._dp+(14*nf)/9._dp)*H1(1)+(68*H2(0,0))/9._dp+(85*H2(1,0))/9._dp- &
    (-1020*omv+u*(59061-4996*nf+6174*zeta2+84*(nf*nf)))/(216._dp*u)- &
    (5*(H1(1)*H1(1)))/12._dp)+G1(0)*G1(0)* &
    (-98+(27*nf)/4._dp+(5*omv)/(2._dp*u)-(105*zeta2)/4._dp+(9*G2(0,2))/2._dp- &
    &  9*G2(1,2)+(9*G2(2,0))/2._dp-4*G2(3,0)+4*G1(3)*H1(0)+ &
    G1(2)*(21._dp/2._dp-nf-(17*H1(0))/2._dp+H1(1)/2._dp)+((9+nf)*H1(1))/2._dp+ &
    H1(0)*(-41._dp/3._dp+(14*nf)/9._dp+(17*H1(1))/2._dp)+G1(1)*(9*H1(0)+9*H1(1))+ &
    &  4*H2(0,0)+5*H2(1,0)-(5*(nf*nf))/72._dp-(G1(2)*G1(2))/4._dp+ &
    (67*(H1(0)*H1(0)))/4._dp-(H1(1)*H1(1))/4._dp)+ &
    G1(2)*G1(2)*((-188*nf+15*(181-(4*omv)/u+42*zeta2))/216._dp-G2(0,2)/2._dp+ &
    G2(1,2)-G2(2,0)/2._dp+(4*G2(3,0))/9._dp-(4*G1(3)*H1(0))/9._dp+ &
    H1(0)*(7._dp/3._dp-(4*nf)/9._dp-(11*H1(1))/6._dp)-(4*(-12+nf)*H1(1))/9._dp- &
    (4*H2(0,0))/9._dp-(5*H2(1,0))/9._dp-(29*(H1(0)*H1(0)))/12._dp+(H1(1)*H1(1))/12._dp) &
    +G2(3,0)*(H1(0)*(32-(20*nf)/9._dp+(8*H1(1))/9._dp)-(41*H1(1))/3._dp- &
    (64*H2(0,0))/3._dp-(172*H2(0,1))/9._dp-(92*H2(1,0))/9._dp+ &
    (8232*zeta2+2*nf*(-232+162*u*v+81*(u*u)+81*(v*v))+ &
    &  3*(2656-616*v-2*u*(268+407*v)-24/w-407*(u*u)-407*(v*v)- &
    (80*(1-3*v+v*v))/u))/54._dp-(68*(H1(0)*H1(0)))/9._dp+ &
    (4*(H1(1)*H1(1)))/9._dp)+((-188*nf+15*(181-(4*omv)/u+42*zeta2))/ &
    &  216._dp-(4*H2(0,0))/9._dp-(5*H2(1,0))/9._dp)*(H1(1)*H1(1))+ &
    G2(1,2)*((147-12*nf)*H1(1)-36*H2(0,0)-72*H2(1,0)+ &
    H1(0)*(135-12*nf+7/u+2*H1(1)+5/(u*u))+ &
    (2*(-407+54*nf)*v*(omv*omv)+ &
    (-536+3*(-407+54*nf)*v)*(u*u*u*u*u)+ &
    (-407+54*nf)*(u*u*u*u*u*u)- &
    omv*u*(814-1905*v-463*(v*v)+36*nf*(-3+7*v+v*v))+ &
    u*u*u*u*(5726-1152*v+2898*zeta2-1221*(v*v)+ &
    &  6*nf*(-58+27*(v*v)))+ &
    u*u*(1101-1692*v+21*(v*v)+18*nf*(-8+11*v+v*v)- &
    &  80*(v*v*v))+u*u*u* &
    (-1074+69*v*(90+42*zeta2)-696*(v*v)-407*(v*v*v)+ &
    &  6*nf*(15-61*v+9*(v*v*v))))/(18._dp*omw*(u*u*u))- &
    &  17*(H1(0)*H1(0))+H1(1)*H1(1))+ &
    G1(1)*(6*(-33+2*nf)*H2(0,0)+12*(-12+nf)*H2(1,1)+36*H3(0,0,0)+ &
    &  36*H3(0,0,1)+36*H3(0,1,0)+36*H3(1,0,0)+72*H3(1,1,0)+ &
    H2(1,0)*(-192+14*nf-7/u-5/(u*u))+ &
    H2(0,1)*(-150+10*nf+7/u+5/(u*u))+ &
    (-216*zeta3*(u*u)+6*zeta2*(-10-14*u+(-141+10*nf)*(u*u)))/ &
    (3._dp*(u*u))+(H1(1)*(2*(407-54*nf)*v*(omv*omv)+ &
    omv*u*(814-1905*v-463*(v*v)+36*nf*(-3+7*v+v*v))+ &
    u*u*(-1101+1692*v+(536+1221*v)*(u*u*u)+407*(u*u*u*u)- &
    &  21*(v*v)+u*u*(-5726+1152*v-2898*zeta2+1221*(v*v))+ &
    &  80*(v*v*v)+u*(1074-69*v*(90+42*zeta2)+696*(v*v)+ &
    &  407*(v*v*v))- &
    &  6*nf*(27*v*(u*u*u)+9*(u*u*u*u)+3*(-8+11*v+v*v)+ &
    u*u*(-58+27*(v*v))+u*(15-61*v+9*(v*v*v))))))/ &
    (18._dp*omw*(u*u*u))+(-H1(0)-H1(1))*(G1(2)*G1(2))+ &
    (57-2*nf+15*H1(1))*(H1(0)*H1(0))+17*(H1(0)*H1(0)*H1(0))+ &
    H1(0)*((54-2*nf)*H1(1)+(2*(407-54*nf)*v*(omv*omv)+ &
    omv*u*(814-1905*v-463*(v*v)+36*nf*(-3+7*v+v*v))+ &
    u*u*(-1101+1692*v+(536+1221*v)*(u*u*u)+407*(u*u*u*u)- &
    &  21*(v*v)+u*u*(-5726+1152*v-1602*zeta2+1221*(v*v))+ &
    &  80*(v*v*v)+u*(1074-3*v*(2070+534*zeta2)+696*(v*v)+ &
    &  407*(v*v*v))- &
    &  6*nf*(27*v*(u*u*u)+9*(u*u*u*u)+3*(-8+11*v+v*v)+ &
    u*u*(-58+27*(v*v))+u*(15-61*v+9*(v*v*v)))))/ &
    (18._dp*omw*(u*u*u))-3*(H1(1)*H1(1)))-3*(H1(1)*H1(1))+ &
    G1(2)*(3*H1(1)+H1(0)*(3+4*H1(1))+2*(H1(0)*H1(0))+2*(H1(1)*H1(1)))- &
    H1(1)*H1(1)*H1(1))+G1(2)*(G2(1,2)*(-3-2*H1(0)-2*H1(1))+ &
    G2(3,0)*(-4._dp/3._dp-(8*H1(0))/9._dp-(8*H1(1))/9._dp)+ &
    G2(0,2)*(3._dp/2._dp+H1(0)+H1(1))+G2(2,0)*(3._dp/2._dp+H1(0)+H1(1))+ &
    (719._dp/3._dp-16*nf)*H2(0,0)-(404*H3(0,0,0))/9._dp-(160*H3(0,0,1))/9._dp- &
    (160*H3(0,1,0))/9._dp+(160*H3(1,0,1))/9._dp-(568*H3(1,1,0))/9._dp+ &
    H2(1,0)*(-8*nf-(3*omv)/u+(5*(198-4/omw+1/(omw*omw)))/9._dp)- &
    (40*(nf*nf))/27._dp+H1(1)*((8*H2(0,0))/9._dp+(10*H2(1,0))/9._dp- &
    (-60*omv+u*(15000-2384*nf+3534*zeta2+96*(nf*nf)))/(108._dp*u)) &
    +(H2(0,1)*(-135+20/omw-(83*omv)/u-5/(omw*omw)- &
    (40*(omv*omv))/(u*u)))/9._dp+ &
    nf*(5086._dp/81._dp+3*u+6*v-(2*v)/omu+(91*zeta2)/9._dp+ &
    (6*(omv*omv))/(u*u)-(38*(v*v))/(3._dp*(omu*omu))+ &
    (3*(-2+v+v*v))/u)+ &
    (-49396-2442*u-4884*v+(6480*v)/omu+3084*zeta3+ &
    &  6*zeta2*(-1551-120/omw+(112*omv)/u+30/(omw*omw)+ &
    (80*(omv*omv))/(u*u))-(4884*(omv*omv))/(u*u)+ &
    (5184*(v*v))/(omu*omu)-(6*(-1093+686*v+407*(v*v)))/u)/108._dp+ &
    G1(3)*(H1(0)*(4._dp/3._dp+(8*H1(1))/9._dp)+(8*(H1(0)*H1(0)))/9._dp)+ &
    ((4*(-93+nf))/9._dp+(17*H1(1))/6._dp)*(H1(0)*H1(0))- &
    (289*(H1(0)*H1(0)*H1(0)))/18._dp+ &
    H1(0)*(((-69+8*nf)*H1(1))/9._dp+(8*H2(0,0))/9._dp+(10*H2(1,0))/9._dp+ &
    (12*(-407+54*nf)*v*(omv*omv)+ &
    &  6*(-536+3*(-407+54*nf)*v)*(u*u*u*u*u)+ &
    &  6*(-407+54*nf)*(u*u*u*u*u*u)- &
    &  6*omv*u*(814-1905*v-463*(v*v)+36*nf*(-3+7*v+v*v))+ &
    u*u*u*u*(19356-6912*v+10830*zeta2-96*(nf*nf)-7326*(v*v)+ &
    &  4*nf*(74+243*(v*v)))+ &
    &  6*(u*u)*(1101-1692*v-69*(v*v)+18*nf*(-8+11*v+v*v)- &
    &  80*(v*v*v))+u*u*u* &
    (-6384+5*v*(4344+2166*zeta2)-96*v*(nf*nf)-4176*(v*v)- &
    &  2442*(v*v*v)+4*nf*(135+47*v+81*(v*v*v))))/ &
    (108._dp*omw*(u*u*u))+(5*(H1(1)*H1(1)))/6._dp)+ &
    (4*(-12+nf)*(H1(1)*H1(1)))/9._dp-(H1(1)*H1(1)*H1(1))/18._dp)+ &
    H1(0)*((76._dp/3._dp-(8*nf)/9._dp)*H2(0,0)-(5*(-57+2*nf)*H2(1,0))/9._dp+ &
    H1(1)*((-8*H2(0,0))/9._dp-(10*H2(1,0))/9._dp+ &
    (-60*omv+u*(-4305+766*nf+6462*zeta2-24*(nf*nf)))/(108._dp*u)) &
    +(36*(-407+54*nf)*v*(omv*omv)+ &
    &  3*omv*u*(160*u*(nf*nf)*((-1+u)*(-1+u))+ &
    nf*(-1668*u*zeta2*((-1+u)*(-1+u))+ &
    &  12*(20-163*v+81*(v*v)))+ &
    &  3*(5910*u*zeta2*((-1+u)*(-1+u))-6*(180-949*v+407*(v*v))) &
    )+2*(u*u)*(2*nf*(-5365+6832*v-369*(v*v)+ &
    u*u*(-5005+4834*v+243*(v*v))-1026*(v*v*v)+ &
    u*(10190-10505*v-72*(v*v)+243*(v*v*v)))- &
    &  3*(-22670+35333*v-4758*zeta3+4758*v*zeta3-9183*(v*v)+ &
    u*u*(-19430+19097*v-4758*zeta3+4758*v*zeta3+ &
    &  1221*(v*v))-2592*(v*v*v)+ &
    u*(40480-46711*v+9516*zeta3-9516*v*zeta3+3234*(v*v)+ &
    &  1221*(v*v*v)))))/(324._dp*omv*(omu*omu)*(u*u))- &
    (4*(-12+nf)*(H1(1)*H1(1)))/9._dp+(H1(1)*H1(1)*H1(1))/18._dp)+ &
    G1(0)*(G2(0,2)*((-30+nf)/2._dp-9*H1(0))+G2(2,0)*((-30+nf)/2._dp-9*H1(0))+ &
    G2(3,0)*((-4*(-30+nf))/9._dp+8*H1(0))+G2(1,2)*(30-nf+18*H1(0))- &
    (32*(-33+2*nf)*H2(0,0))/9._dp+(135._dp/2._dp-5*nf)*H2(0,1)+ &
    (737._dp/6._dp-9*nf)*H2(1,0)+(72-6*nf)*H2(1,1)-(64*H3(0,0,0))/3._dp- &
    (172*H3(0,0,1))/9._dp-(92*H3(0,1,0))/9._dp-(332*H3(1,0,0))/9._dp- &
    (20*H3(1,0,1))/3._dp+(10*H3(1,1,0))/9._dp+ &
    (330*omv+u*(451+407*u+407*v-192*zeta2-322*zeta3)- &
    &  2*nf*(10-9*u-10*v+27*u*v+27*(u*u)))/(18._dp*u)+ &
    (H1(1)*(12*(-268+(-407+54*nf)*v)*(u*u)+6*(-407+54*nf)*(u*u*u)- &
    &  60*(7-15*v+8*(v*v))+ &
    u*(18066-3696*v+13422*zeta2-2442*(v*v)+ &
    &  36*nf*(-29+9*(v*v)))))/(108._dp*u)+ &
    ((17*H1(0))/18._dp+H1(1)/2._dp)*(G1(2)*G1(2))+ &
    G1(1)*(H1(0)*(-30+nf-18*H1(1))+(-30+nf)*H1(1)-18*(H1(0)*H1(0)))+ &
    G1(3)*((4*(-30+nf)*H1(0))/9._dp-8*(H1(0)*H1(0)))+ &
    (73._dp/6._dp-(19*nf)/9._dp-(119*H1(1))/18._dp)*(H1(0)*H1(0))- &
    (289*(H1(0)*H1(0)*H1(0)))/18._dp+ &
    G1(2)*(-((-30+nf)*(-39+4*nf))/36._dp+ &
    H1(0)*(23._dp/3._dp+nf-(26*H1(1))/9._dp)-(3*H1(1))/2._dp+ &
    (145*(H1(0)*H1(0)))/9._dp-H1(1)*H1(1))+(3*(H1(1)*H1(1)))/2._dp+ &
    H1(0)*((-37._dp/6._dp-nf)*H1(1)-8*H2(0,0)-10*H2(1,0)+ &
    (18*(43+nf*(18-54*v)+407*v)*(u*u*u)- &
    &  6*(-407+54*nf)*(u*u*u*u)+ &
    u*u*(-31680+2028*v-16302*zeta2+ &
    nf*(1648+648*v-972*(v*v))+7326*(v*v))+ &
    &  60*(-17+50*v-41*(v*v)+8*(v*v*v))+ &
    u*(16302*omv*zeta2- &
    &  4*nf*(412-412*v-81*(v*v)+81*(v*v*v))+ &
    &  6*(4890-5690*v+289*(v*v)+407*(v*v*v))))/(108._dp*u*w)+ &
    (35*(H1(1)*H1(1)))/18._dp)+(H1(1)*H1(1)*H1(1))/2._dp)- &
    (4*(-12+nf)*(H1(1)*H1(1)*H1(1)))/27._dp+(H1(1)*H1(1)*H1(1)*H1(1))/72._dp
    return
    end function Hgamma_3a2re
          
          
    function Hgamma_3a2im(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    use constants
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=(-16*pi*G3(0,0,1))/9._dp+18*pi*G3(0,0,2)+(16*pi*G3(0,0,3))/9._dp+ &
    &  26*pi*G3(0,1,0)-10*pi*G3(0,1,2)+18*pi*G3(0,2,0)-(92*pi*G3(0,2,1))/9._dp+ &
    &  18*pi*G3(0,2,2)-(152*pi*G3(0,2,3))/9._dp+(8*pi*G3(0,3,0))/9._dp+ &
    (8*pi*G3(0,3,2))/9._dp+(128*pi*G3(0,3,3))/9._dp-36*pi*G3(1,0,2)- &
    &  72*pi*G3(1,1,1)+72*pi*G3(1,1,2)-36*pi*G3(1,2,0)+72*pi*G3(1,2,1)- &
    &  36*pi*G3(1,2,2)+18*pi*G3(2,0,0)-18*pi*G3(2,0,1)+18*pi*G3(2,0,2)- &
    &  18*pi*G3(2,0,3)-18*pi*G3(2,1,0)+72*pi*G3(2,1,1)-18*pi*G3(2,1,2)+ &
    &  18*pi*G3(2,2,0)-36*pi*G3(2,2,1)+36*pi*G3(2,2,3)-(80*pi*G3(2,3,0))/9._dp+ &
    (164*pi*G3(2,3,2))/9._dp+(80*pi*G3(2,3,3))/9._dp-(164*pi*G3(3,0,0))/9._dp- &
    (70*pi*G3(3,0,1))/9._dp-(164*pi*G3(3,0,2))/9._dp+(16*pi*G3(3,0,3))/3._dp- &
    (164*pi*G3(3,2,0))/9._dp-(70*pi*G3(3,2,1))/9._dp-(80*pi*G3(3,2,3))/9._dp+ &
    (128*pi*G3(3,3,0))/9._dp-(128*pi*G3(3,3,3))/9._dp+ &
    G2(3,3)*((8*(-33+2*nf)*pi)/9._dp-(128*pi*H1(0))/9._dp)+ &
    G2(3,2)*(15*pi+18*pi*H1(0))+G2(2,2)*(-6*(-12+nf)*pi+36*pi*H1(0))+ &
    G2(2,3)*((-121._dp/3._dp+4*nf)*pi-(2*pi*H1(0))/9._dp-(164*pi*H1(1))/9._dp)+ &
    G2(0,2)*((-9*(-15+nf)*pi)/2._dp-(241*pi*H1(0))/9._dp-18*pi*H1(1))+ &
    G2(0,0)*(-3*(-21+nf)*pi-(164*pi*H1(0))/9._dp-18*pi*H1(1))+ &
    G2(2,0)*((-9*(-15+nf)*pi)/2._dp-(161*pi*H1(0))/9._dp-18*pi*H1(1))+ &
    G2(0,3)*((4*(-96+5*nf)*pi)/9._dp+(56*pi*H1(0))/9._dp-(8*pi*H1(1))/9._dp)+ &
    G2(3,0)*((19-(4*nf)/3._dp)*pi+(76*pi*H1(0))/3._dp+(164*pi*H1(1))/9._dp)+ &
    G2(1,2)*((-39+4*nf)*pi+36*pi*H1(1))+ &
    G2(1,0)*((-30+nf)*pi+18*pi*H1(0)+36*pi*H1(1))+ &
    ((501-32*nf)*pi*H2(0,0))/9._dp+(5*(27-2*nf)*pi*H2(0,1))/2._dp- &
    &  6*(-12+nf)*pi*H2(1,1)-(124*pi*H3(0,0,0))/3._dp-(172*pi*H3(0,0,1))/9._dp- &
    &  18*pi*H3(0,1,0)-(104*pi*H3(1,0,0))/3._dp-(20*pi*H3(1,0,1))/3._dp+ &
    (80*pi*H3(1,1,0))/3._dp+G2(0,1)* &
    (-8*pi*H1(0)+10*pi*H1(1)+(pi* &
    (9*nf+5*(-87+4/omw+(22*omv)/u-1/(omw*omw)+ &
    (8*(omv*omv))/(u*u))))/9._dp)+ &
    G2(2,1)*(-18*pi*H1(0)+18*pi*H1(1)+ &
    (pi*(-1431+20/omw+108*nf-(83*omv)/u-5/(omw*omw)- &
    (40*(omv*omv))/(u*u)))/9._dp)+ &
    G2(1,1)*(-72*pi*H1(1)+(2*pi*(5+7*u+(21-2*nf)*(u*u)))/(u*u))+ &
    (pi*H2(1,0)*(-10+40*v+u*(40-2*(-915+82*nf)*v)+ &
    (915-82*nf)*(u*u)+(915-82*nf)*(v*v)))/(18._dp*(omw*omw))+ &
    (pi*H1(1)*(4*(-268+(-407+54*nf)*v)*(u*u)+2*(-407+54*nf)*(u*u*u)- &
    &  20*(7-15*v+8*(v*v))+ &
    u*(6022-1232*v+3450*zeta2-814*(v*v)+12*nf*(-29+9*(v*v)))))/ &
    (36._dp*u)-(pi*(-6*(-407+54*nf)*(omv*omv)+ &
    u*(2*nf*(456-699*v+(-54-81*v+204*zeta2)*(u*u*u)+ &
    u*u*(240+57*v-408*zeta2-81*(v*v))+243*(v*v)+ &
    u*(-480+399*v+204*zeta2+342*(v*v)))- &
    &  3*((-1305-407*v+2106*zeta2-778*zeta3)*(u*u*u)+ &
    u*u*(3493-742*v-4212*zeta2+1556*zeta3-407*(v*v))+ &
    &  3*(837-1244*v+407*(v*v))+ &
    u*(-3885+3253*v+2106*zeta2-778*zeta3+864*(v*v))))))/ &
    (54._dp*(omu*omu)*(u*u))+(((-30+nf)*pi)/2._dp-9*pi*G1(1)+ &
    (9*pi*G1(2))/2._dp-4*pi*G1(3)-(17*pi*H1(0))/2._dp-(9*pi*H1(1))/2._dp)* &
    (G1(0)*G1(0))+(9*pi*(G1(0)*G1(0)*G1(0)))/2._dp+ &
    ((3*pi)/2._dp+(4*pi*G1(3))/9._dp+(17*pi*H1(0))/18._dp+(3*pi*H1(1))/2._dp)* &
    (G1(2)*G1(2))-(pi*(G1(2)*G1(2)*G1(2)))/2._dp+ &
    (((-1131+70*nf)*pi)/18._dp-(155*pi*H1(1))/18._dp)*(H1(0)*H1(0))+ &
    (17*pi*(H1(0)*H1(0)*H1(0)))/18._dp+(3*pi*(H1(1)*H1(1)))/2._dp+ &
    G1(2)*(3*(-25+2*nf)*pi*H1(1)+ &
    G1(3)*((-4*pi)/3._dp-(8*pi*H1(0))/9._dp-(8*pi*H1(1))/9._dp)+ &
    (322*pi*H2(0,0))/9._dp-(2*pi*H2(0,1))/9._dp-(242*pi*H2(1,0))/9._dp+ &
    (pi*(-4*(-268+(-407+54*nf)*v)*(u*u)+(814-108*nf)*(u*u*u)+ &
    &  20*(7-15*v+8*(v*v))+ &
    u*(-6022+1232*v-3210*zeta2+814*(v*v)-12*nf*(-29+9*(v*v))) &
    ))/(36._dp*u)+H1(0)*((-17*pi*H1(1))/9._dp+ &
    (pi*((63+(733-36*nf)*v)*(u*u*u)-6*(-68+3*nf)*(u*u*u*u)+ &
    u*v*(80-77*v-3*(v*v))+40*(omv*omv)*(v*v)+ &
    u*u*(45+66*v-6*(-47+3*nf)*(v*v))))/(9._dp*(omw*omw)*(u*u))) &
    -(169*pi*(H1(0)*H1(0)))/18._dp-(3*pi*(H1(1)*H1(1)))/2._dp)+ &
    G1(0)*(G1(3)*((-4*(-30+nf)*pi)/9._dp+8*pi*H1(0))+ &
    G1(1)*(-((-30+nf)*pi)+18*pi*H1(0))+2*(-27+2*nf)*pi*H1(1)+ &
    H1(0)*((5*(-39+4*nf)*pi)/9._dp-pi*H1(1))+ &
    G1(2)*(((-27+nf)*pi)/2._dp-8*pi*H1(0)+pi*H1(1))+(100*pi*H2(0,0))/3._dp+ &
    (164*pi*H2(0,1))/9._dp+(304*pi*H2(1,0))/9._dp+18*pi*H2(1,1)+ &
    (pi*((3216+(4884-648*nf)*v)*(u*u)-6*(-407+54*nf)*(u*u*u)+ &
    &  60*(7-15*v+8*(v*v))+ &
    u*(-18066+3696*v-5070*zeta2+2442*(v*v)- &
    &  36*nf*(-29+9*(v*v)))))/(108._dp*u)-(pi*(G1(2)*G1(2)))/2._dp- &
    (pi*(H1(0)*H1(0)))/2._dp-(pi*(H1(1)*H1(1)))/2._dp)+ &
    G1(3)*((-41*pi*H1(1))/3._dp+H1(0)*((-43+(16*nf)/9._dp)*pi+(8*pi*H1(1))/9._dp)- &
    (124*pi*H2(0,0))/3._dp-(172*pi*H2(0,1))/9._dp-18*pi*H2(1,0)+ &
    (pi*(5280*zeta2+2*nf*(-232+162*u*v+81*(u*u)+81*(v*v))+ &
    &  3*(2656-616*v-2*u*(268+407*v)-24/w-407*(u*u)- &
    &  407*(v*v)-(80*(1-3*v+v*v))/u)))/54._dp+ &
    (4*pi*(H1(0)*H1(0)))/9._dp+(4*pi*(H1(1)*H1(1)))/9._dp)+ &
    H1(0)*((-289._dp/6._dp+3*nf)*pi*H1(1)+8*pi*H2(0,0)+10*pi*H2(1,0)- &
    (pi*(6*(-129-1628*v+54*nf*(-1+4*v))*(u*u*u*u)+ &
    &  6*(-407+54*nf)*(u*u*u*u*u)-60*(v*v)*(7-15*v+8*(v*v))+ &
    u*u*u*(6624-2802*v+4890*zeta2-14652*(v*v)+ &
    &  4*nf*(-52-243*v+486*(v*v)))+ &
    u*u*(4890*(-1+2*v)*zeta2+ &
    &  4*nf*(52-104*v-243*(v*v)+324*(v*v*v))- &
    &  6*(534-2358*v+627*(v*v)+1628*(v*v*v)))+ &
    u*(-60+v*(-3624+208*nf-4890*zeta2)+ &
    (8424-208*nf+4890*zeta2)*(v*v)-54*(41+6*nf)*(v*v*v)+ &
    &  6*(-407+54*nf)*(v*v*v*v))))/(108._dp*omw*u*w)+ &
    (17*pi*(H1(1)*H1(1)))/18._dp)+ &
    G1(1)*((42-4*nf)*pi*H1(1)+G1(2)*(-3*pi-2*pi*H1(0)-2*pi*H1(1))- &
    &  36*pi*H2(0,0)-18*pi*H2(0,1)-18*pi*H2(1,0)-36*pi*H2(1,1)+ &
    H1(0)*(20*pi*H1(1)-(pi*(5+7*u+12*(u*u)))/(u*u))+ &
    (pi*(2*(-407+54*nf)*v*(omv*omv)+ &
    (-536+3*(-407+54*nf)*v)*(u*u*u*u*u)+ &
    (-407+54*nf)*(u*u*u*u*u*u)- &
    omv*u*(814-1905*v-463*(v*v)+36*nf*(-3+7*v+v*v))+ &
    u*u*u*u*(5726-1152*v+2250*zeta2-1221*(v*v)+ &
    &  6*nf*(-58+27*(v*v)))+ &
    u*u*(1101-1692*v+21*(v*v)+18*nf*(-8+11*v+v*v)- &
    &  80*(v*v*v))+u*u*u* &
    (-1074+15*v*(414+150*zeta2)-696*(v*v)-407*(v*v*v)+ &
    &  6*nf*(15-61*v+9*(v*v*v)))))/(18._dp*omw*(u*u*u))+ &
    pi*(G1(2)*G1(2))+pi*(H1(0)*H1(0))+pi*(H1(1)*H1(1)))+ &
    (pi*(H1(1)*H1(1)*H1(1)))/2._dp
    return
    end function Hgamma_3a2im
          
          
!**********************************************************************
! egion4a:g+g-->V+g*
!*********************************************************************/
          
    function Hbeta_4a0re(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=1
    return
    end function Hbeta_4a0re
          
          
    function Hbeta_4a0im(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=0
    return
    end function Hbeta_4a0im
          
          
    function Hbeta_4a1re(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    use constants
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=-3*G2(0,2)+6*G2(1,2)-3*G2(2,0)+6*G2(3,0)-6*G1(3)*H1(0)+ &
    G1(1)*(-6*H1(0)-6*H1(1))-3*H1(0)*H1(1)+G1(0)*(9*H1(0)+3*H1(1))+ &
    G1(2)*(9*H1(0)+3*H1(1))-6*H2(0,0)+ &
    (4*(-3+nf)*u+342*zeta2*(w*w))/(12._dp*(w*w))-(3*(G1(0)*G1(0)))/2._dp- &
    (3*(G1(2)*G1(2)))/2._dp-(9*(H1(0)*H1(0)))/2._dp-(3*(H1(1)*H1(1)))/2._dp
    return
    end function Hbeta_4a1re
          
          
    function Hbeta_4a1im(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    use constants
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=-3*pi*G1(0)+6*pi*G1(1)-3*pi*G1(2)+6*pi*G1(3)+3*pi*H1(0)+3*pi*H1(1)
    return
    end function Hbeta_4a1im
          
          
    function Hbeta_4a2re(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    use constants
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=(66-4*nf)*G3(2,0,2)+(66-4*nf)*G3(2,2,0)+36*G4(0,0,1,2)+18*G4(0,0,2,2)- &
    &  36*G4(0,0,3,0)+36*G4(0,1,0,2)+36*G4(0,1,2,0)-36*G4(0,1,2,2)+ &
    &  18*G4(0,2,0,2)-18*G4(0,2,1,2)+18*G4(0,2,2,0)-18*G4(0,2,3,0)- &
    &  18*G4(0,3,0,2)-18*G4(0,3,2,0)+72*G4(0,3,3,0)-36*G4(1,0,2,2)- &
    &  72*G4(1,1,1,2)+72*G4(1,1,2,2)-36*G4(1,2,0,2)+72*G4(1,2,1,2)- &
    &  36*G4(1,2,2,0)+18*G4(2,0,0,2)-18*G4(2,0,1,2)+18*G4(2,0,2,0)- &
    &  18*G4(2,0,3,0)-18*G4(2,1,0,2)+72*G4(2,1,1,2)-18*G4(2,1,2,0)+ &
    &  18*G4(2,2,0,0)-36*G4(2,2,1,2)+36*G4(2,2,3,0)-36*G4(2,3,0,0)+ &
    &  36*G4(2,3,0,2)+36*G4(2,3,2,0)-36*G4(3,0,0,2)-36*G4(3,0,2,0)+ &
    &  72*G4(3,0,3,0)-36*G4(3,2,0,0)+72*G4(3,3,0,0)-72*G4(3,3,3,0)+ &
    G3(0,2,2)*(66-4*nf-36*H1(0))+36*G3(0,0,3)*H1(0)+18*G3(0,2,3)*H1(0)+ &
    &  18*G3(0,3,2)*H1(0)-72*G3(0,3,3)*H1(0)+18*G3(2,0,3)*H1(0)- &
    &  36*G3(2,2,3)*H1(0)-36*G3(2,3,2)*H1(0)-72*G3(3,0,3)*H1(0)+ &
    &  72*G3(3,3,3)*H1(0)+G3(3,0,2)*(-33+2*nf+18*H1(0))+ &
    G3(3,2,0)*(-33+2*nf+18*H1(0))+G3(2,1,2)*(-99+6*nf+36*H1(0))+ &
    G3(1,2,2)*(-132+8*nf+72*H1(0))+G2(3,2)*(108*zeta2+(33-2*nf)*H1(0))+ &
    G3(1,2,1)*(-72*H1(0)-72*H1(1))+G3(2,1,1)*(-72*H1(0)-72*H1(1))+ &
    G3(1,1,2)*(-66+4*nf-72*H1(0)-72*H1(1))+ &
    G3(0,0,1)*(-36*H1(0)-36*H1(1))+G3(0,1,0)*(-36*H1(0)-36*H1(1))+ &
    G3(2,3,0)*(-33+2*nf+18*H1(0)-36*H1(1))+ &
    G3(0,0,2)*(66-4*nf-18*H1(1))+G3(0,2,0)*(66-4*nf-18*H1(1))+ &
    G3(2,0,0)*(66-4*nf-36*H1(0)-18*H1(1))+G3(0,2,1)*(18*H1(0)+18*H1(1))+ &
    G3(2,0,1)*(18*H1(0)+18*H1(1))+G3(2,1,0)*(18*H1(0)+18*H1(1))+ &
    G3(0,1,2)*(18*H1(0)+36*H1(1))+G3(1,0,2)*(18*H1(0)+36*H1(1))+ &
    G3(1,2,0)*(18*H1(0)+36*H1(1))+G3(2,2,1)*(36*H1(0)+36*H1(1))+ &
    G3(3,0,0)*(-132+8*nf+72*H1(0)+36*H1(1))+ &
    G3(1,1,1)*(72*H1(0)+72*H1(1))+ &
    G2(2,1)*(-36*zeta2+(99-6*nf)*H1(0)+(99-6*nf)*H1(1)+36*H2(0,1)- &
    &  36*H2(1,0))+G2(2,3)*(-108*zeta2+(33-2*nf)*H1(0)+36*H2(0,1)+ &
    &  36*H2(1,0))+G2(0,1)*(-108*zeta2-18*H2(0,1)-18*H2(1,0)-36*H2(1,1))+ &
    G2(1,0)*(108*zeta2-18*H2(0,1)-18*H2(1,0)-36*H2(1,1))+ &
    G2(1,1)*(-72*zeta2+(66-4*nf)*H1(0)+(66-4*nf)*H1(1)+72*H2(1,0)+ &
    &  72*H2(1,1))+(-33+2*nf)*H3(0,0,1)+(-33+2*nf)*H3(1,0,0)+ &
    &  108*H4(0,0,0,0)+H2(1,1)*(363._dp/4._dp-11*nf-54*zeta2+(nf*nf)/3._dp)+ &
    G2(2,2)*(363._dp/4._dp-11*nf+18*zeta2+(-132+8*nf)*H1(0)-36*H2(0,1)+ &
    &  36*H2(1,0)+(nf*nf)/3._dp)+ &
    G2(0,0)*(363._dp/4._dp-11*nf-54*zeta2+(-132+8*nf)*H1(0)+ &
    (-66+4*nf)*H1(1)+72*H2(0,0)+36*H2(0,1)+36*H2(1,0)+18*H2(1,1)+ &
    (nf*nf)/3._dp)+H2(0,1)*(605._dp/4._dp-(55*nf)/3._dp-108*zeta2+ &
    (5*(nf*nf))/9._dp)+G3(0,3,0)* &
    (36*H1(0)+18*H1(1)+(-27*(-17+nf)*u*(omv*omv)+ &
    &  5*(-33+2*nf)*(omv*omv*omv)+27*omv*(-16+nf)*(u*u)+ &
    &  4*(33-2*nf)*(u*u*u))/(w*w*w))+ &
    G2(0,3)*(36*zeta2-72*H2(0,0)-18*H2(0,1)-18*H2(1,0)+ &
    (H1(0)*(27*(-17+nf)*u*(omv*omv)-5*(-33+2*nf)*(omv*omv*omv)- &
    &  27*omv*(-16+nf)*(u*u)+4*(-33+2*nf)*(u*u*u)))/(w*w*w))+ &
    G3(3,3,0)*(-72*H1(0)-(2*(-33+63*v-9*u*(-7+8*v)-36*(u*u)- &
    &  36*(v*v)+nf*(2-3*v+u*(-3+6*v)+3*(u*u)+3*(v*v))))/ &
    (w*w*w))+G2(3,3)*(-144*zeta2+72*H2(0,0)+ &
    (2*H1(0)*(-33+63*v-9*u*(-7+8*v)-36*(u*u)-36*(v*v)+ &
    nf*(2-3*v+u*(-3+6*v)+3*(u*u)+3*(v*v))))/(w*w*w))+ &
    H1(1)*((33-2*nf)*H2(0,0)+(-3*u* &
    (652-382*v-216*u*zeta2+nf*(-134+34*v+18*u*zeta2)- &
    &  27*zeta3+27*v*zeta3)*(omv*omv)+ &
    (270-100*nf)*(omv*omv*omv)+ &
    omv*(u*u)*(nf*(-667+506*v+54*u*zeta2-100*(v*v))- &
    &  9*(-491+442*v+126*u*zeta2+27*zeta3-27*v*zeta3-30*(v*v))) &
    +u*u*(-2*w*(nf*nf)- &
    &  3*(-652+270*v-198*zeta2+2*nf*(67-50*v+6*zeta2)+ &
    &  27*zeta3)*(u*u)+10*(-27+10*nf)*(u*u*u)+ &
    u*(-9*(491-562*v-27*zeta3+27*v*zeta3+90*(v*v))+ &
    nf*(667-906*v+300*(v*v)))))/(9._dp*u*(w*w*w)))+ &
    (2*H3(1,1,0)*(-9*(-15+nf)*u*(omv*omv)+(-33+2*nf)*(omv*omv*omv)+ &
    &  9*omv*(-18+nf)*(u*u)+2*(33-2*nf)*(u*u*u)))/(w*w*w)+ &
    (H3(0,1,0)*(-33+63*v+9*(-11+7*v)*(u*u)+33*(u*u*u)-36*(v*v)+ &
    &  9*u*(11-16*v+4*(v*v))+ &
    nf*(2-3*v-3*(-2+v)*(u*u)-2*(u*u*u)+3*(v*v)- &
    &  3*u*(2-4*v+v*v))))/(w*w*w)+ &
    (H3(0,0,0)*(-(nf*(-26+75*v+(-72+69*v)*(u*u)+22*(u*u*u)-75*(v*v)+ &
    &  3*u*(24-48*v+23*(v*v))+24*(v*v*v)))+ &
    &  3*(-143+417*v+(-387+375*v)*(u*u)+121*(u*u*u)-408*(v*v)+ &
    u*(405-792*v+384*(v*v))+132*(v*v*v))))/(w*w*w)+ &
    (H2(1,0)*(-4*(135*(13-8*v)+nf*(-524+400*v))*(u*u*u*u*u)- &
    &  40*(-27+10*nf)*(u*u*u*u*u*u)+ &
    u*u*u*u*(45*(505-528*v+144*(v*v))- &
    &  4*nf*(1261-1696*v+600*(v*v)))- &
    &  3*(omv*omv)*(4*nf*(55-116*v+55*(v*v))- &
    &  3*(605-1234*v+605*(v*v)))+ &
    &  4*u*(-1+v)*(2*nf*(-445+924*v-507*(v*v)+62*(v*v*v))- &
    &  9*(-749+1534*v-824*(v*v)+75*(v*v*v)))- &
    &  4*(u*u*u)*(-9*(-1244+1859*v-810*(v*v)+120*(v*v*v))+ &
    &  2*nf*(-955+1665*v-972*(v*v)+200*(v*v*v)))- &
    &  2*(u*u)*(-27*(857-1994*v+1317*(v*v)-280*(v*v*v)+ &
    &  20*(v*v*v*v))+4*nf* &
    (749-2047*v+1521*(v*v)-448*(v*v*v)+50*(v*v*v*v)))+ &
    &  4*(-972*zeta2+5*(nf*nf))*(w*w*w*w)))/(36._dp*(w*w*w*w))+ &
    (H2(0,0)*(-4*(135*(13-8*v)+nf*(-524+400*v))*(u*u*u*u*u)- &
    &  40*(-27+10*nf)*(u*u*u*u*u*u)- &
    &  3*(omv*omv)*(-3837+7746*v-3837*(v*v)+ &
    &  4*nf*(145-296*v+145*(v*v)))+ &
    u*u*u*u*(-4*nf*(1531-1696*v+600*(v*v))+ &
    &  9*(3199-2640*v+720*(v*v)))+ &
    &  4*u*(-1+v)*(2*nf*(-870+1857*v-1047*(v*v)+62*(v*v*v))- &
    &  9*(-1279+2690*v-1498*(v*v)+75*(v*v*v)))- &
    &  4*(u*u*u)*(-9*(-1774+2533*v-810*(v*v)+120*(v*v*v))+ &
    &  2*nf*(-1380+2205*v-972*(v*v)+200*(v*v*v)))- &
    &  2*(u*u)*(-27*(1459-3118*v+1991*(v*v)-280*(v*v*v)+ &
    &  20*(v*v*v*v))+4*nf* &
    (1629-3405*v+2331*(v*v)-448*(v*v*v)+50*(v*v*v*v)))+ &
    &  6*(-2106*zeta2+10*(nf*nf))*(w*w*w*w)))/(36._dp*(w*w*w*w))+ &
    (-360*omw*zeta2*(omv*omv*omv)* &
    (-8241+16698*v-8241*(v*v)+nf*(410-892*v+410*(v*v)))+ &
    &  2880*omw*(-72+19*nf)*(w*w*w*w)+ &
    omv*(-2880*(-72+19*nf)*v*(w*w*w*w)+ &
    omw*(-80*w*(-216*(-3+nf)*zeta3+ &
    w*(-18*(-337+337*v+54*zeta3)+ &
    &  6*nf*(-377+377*v+54*zeta3)+20*u*(nf*nf))+ &
    &  6*(-3*(337+99*v-216*zeta3)+nf*(377+64*v-54*zeta3))* &
    (w*w)+3*(-3693+757*nf-297*zeta3+258*nf*zeta3)*(w*w*w)) &
    +2753190*zeta4*(w*w*w*w)- &
    &  120*u*zeta2*(2*nf*(-2915+8534*v+(-2915+2424*v)*(u*u)+ &
    &  615*(u*u*u)-8043*(v*v)+2424*(v*v*v)+ &
    u*(-8390*v+3618*(v*v)-800*(-5+w*w*w*w))+ &
    &  992*(w*w*w*w))- &
    &  9*(2*(-5713+5482*v)*(u*u)+2747*(u*u*u)+ &
    &  6*u*(2833-5660*v+2739*(v*v)-80*(w*w*w*w))+ &
    &  2*(-5713+17028*v-16797*(v*v)+5482*(v*v*v)+ &
    &  600*(w*w*w*w)))))))/(4320._dp*omv*omw*(w*w*w*w))+ &
    ((33-2*nf)/6._dp-(27*H1(0))/2._dp-(9*H1(1))/2._dp)*(G1(0)*G1(0)*G1(0))+ &
    (9*(G1(0)*G1(0)*G1(0)*G1(0)))/8._dp+ &
    ((33-2*nf)/6._dp-(27*H1(0))/2._dp-(9*H1(1))/2._dp)*(G1(2)*G1(2)*G1(2))+ &
    (9*(G1(2)*G1(2)*G1(2)*G1(2)))/8._dp+ &
    (-33._dp/2._dp+nf+(27*H1(1))/2._dp)*(H1(0)*H1(0)*H1(0))+ &
    (81*(H1(0)*H1(0)*H1(0)*H1(0)))/8._dp+ &
    G2(1,2)*(134+75/u+75*u+6/w+171*zeta2+H1(0)*(33-2*nf-18*H1(1))+ &
    (99-6*nf)*H1(1)-36*H2(0,0)-72*H2(1,0)-30/(u*u)-30*(u*u)+ &
    (2*nf*(-62/u-62*u+50/(u*u)+50*(u*u)+ &
    (9-9*v-9*w-30*(w*w))/(w*w)))/9._dp-6/(w*w)+(6*v)/(w*w)- &
    &  27*(H1(0)*H1(0))-9*(H1(1)*H1(1)))+ &
    G2(3,0)*(-108*H2(0,0)-18*H2(0,1)-18*H2(1,0)+ &
    H1(0)*(-18*H1(1)+(6-9*w+36*(w*w)+66*(w*w*w)- &
    nf*(2-3*w+3*(w*w)+4*(w*w*w)))/(w*w*w))+ &
    (2*(135*(13-8*v)+nf*(-524+400*v))*(u*u*u*u*u)+ &
    &  20*(-27+10*nf)*(u*u*u*u*u*u)- &
    &  12*(omv*omv)*(nf*(10-17*v+10*(v*v))- &
    &  3*(67-131*v+67*(v*v)))+ &
    &  4*(u*u*u*u)*(-9*(173-330*v+90*(v*v))+ &
    nf*(518-848*v+300*(v*v)))- &
    &  2*u*(-1+v)*(2*nf*(120-123*v-57*(v*v)+62*(v*v*v))- &
    &  9*(536-940*v+317*(v*v)+75*(v*v*v)))+ &
    &  2*(u*u*u)*(-9*(41+718*v-810*(v*v)+120*(v*v*v))+ &
    &  2*nf*(-390+1215*v-972*(v*v)+200*(v*v*v)))+ &
    &  4*(u*u)*(-27*(-89+128*v+44*(v*v)-70*(v*v*v)+5*(v*v*v*v))+ &
    nf*(144-435*v+846*(v*v)-448*(v*v*v)+50*(v*v*v*v)))+ &
    &  4374*zeta2*(w*w*w*w))/(18._dp*(w*w*w*w))-27*(H1(0)*H1(0))- &
    &  9*(H1(1)*H1(1)))+(9*H2(0,0)+ &
    (u*(3822-3786*v+4*nf*(-89+86*v))+(-1893+172*nf)*(omv*omv)+ &
    (-1893+172*nf)*(u*u)-(810*zeta2+4*(nf*nf))*(w*w))/(24._dp*(w*w)))* &
    (H1(1)*H1(1))+G1(0)*G1(0)*((9*G2(0,2))/2._dp-9*G2(1,2)+(9*G2(2,0))/2._dp- &
    &  9*G2(3,0)+9*G1(3)*H1(0)+G1(2)*((-27*H1(0))/2._dp-(9*H1(1))/2._dp)+ &
    (33._dp/2._dp-nf)*H1(1)+G1(1)*(9*H1(0)+9*H1(1))+ &
    H1(0)*(33._dp/2._dp-nf+(27*H1(1))/2._dp)+9*H2(0,0)+ &
    (u*(3822-3786*v+4*nf*(-89+86*v))+(-1893+172*nf)*(omv*omv)+ &
    (-1893+172*nf)*(u*u)-(810*zeta2+4*(nf*nf))*(w*w))/(24._dp*(w*w)) &
    +(9*(G1(2)*G1(2)))/4._dp+(117*(H1(0)*H1(0)))/4._dp+(9*(H1(1)*H1(1)))/4._dp)+ &
    G2(2,0)*(-147._dp/4._dp-75*u-3/w-(243*zeta2)/2._dp+(-99._dp/2._dp+3*nf)*H1(1)+ &
    H1(0)*(-66+4*nf+9*H1(1))+18*H2(0,0)-18*H2(0,1)-18*H2(1,0)+ &
    (nf*nf)/9._dp+30*(u*u)+(nf* &
    (124*u-100*(u*u)-(3*(u*(1+u+2*v)+omv*omv))/(w*w)))/9._dp+ &
    &  3/(w*w)-(3*v)/(w*w)+(27*(H1(0)*H1(0)))/2._dp+(9*(H1(1)*H1(1)))/2._dp)+ &
    G2(0,2)*(-147._dp/4._dp-75*u-3/w-(243*zeta2)/2._dp+(-99._dp/2._dp+3*nf)*H1(1)+ &
    H1(0)*(-66+4*nf+9*H1(1))+18*H2(0,0)+18*H2(0,1)+18*H2(1,0)+ &
    (nf*nf)/9._dp+30*(u*u)+(nf* &
    (124*u-100*(u*u)-(3*(u*(1+u+2*v)+omv*omv))/(w*w)))/9._dp+ &
    &  3/(w*w)-(3*v)/(w*w)+(27*(H1(0)*H1(0)))/2._dp+(9*(H1(1)*H1(1)))/2._dp)+ &
    G1(2)*G1(2)*((9*G2(0,2))/2._dp-9*G2(1,2)+(9*G2(2,0))/2._dp-9*G2(3,0)+ &
    &  9*G1(3)*H1(0)+(-33._dp/2._dp+nf)*H1(1)+H1(0)*(33._dp/2._dp-nf+(63*H1(1))/2._dp)+ &
    &  9*H2(0,0)+(u*(3822-3786*v+4*nf*(-89+86*v))+ &
    (-1893+172*nf)*(omv*omv)+(-1893+172*nf)*(u*u)- &
    (810*zeta2+4*(nf*nf))*(w*w))/(24._dp*(w*w))+(117*(H1(0)*H1(0)))/4._dp+ &
    (27*(H1(1)*H1(1)))/4._dp)+H1(0)*H1(0)* &
    ((-33._dp/2._dp+nf)*H1(1)+27*H2(0,0)+ &
    (6*u*(2637-266*nf-2619*v+260*nf*v)+ &
    (-7857+780*nf)*(omv*omv)+(-7857+780*nf)*(u*u)- &
    (486*zeta2+20*(nf*nf))*(w*w))/(24._dp*(w*w))+(45*(H1(1)*H1(1)))/4._dp)+ &
    G1(3)*(-108*zeta2*H1(1)+(-33+2*nf)*H2(0,1)+(-33+2*nf)*H2(1,0)+ &
    &  108*H3(0,0,0)-(2*zeta2*(nf* &
    (-8+21*v+3*(-7+6*v)*(u*u)+6*(u*u*u)-21*(v*v)+ &
    &  3*u*(7-14*v+6*(v*v))+6*(v*v*v))- &
    &  3*(-44+120*v+3*(-37+33*v)*(u*u)+33*(u*u*u)-111*(v*v)+ &
    &  3*u*(40-74*v+33*(v*v))+33*(v*v*v))))/(w*w*w)+ &
    (2*H2(0,0)*(-6+9*w-36*(w*w)-99*(w*w*w)+ &
    nf*(2-3*w+3*(w*w)+6*(w*w*w))))/(w*w*w)+ &
    (99-6*nf+18*H1(1))*(H1(0)*H1(0))+27*(H1(0)*H1(0)*H1(0))+ &
    H1(0)*((33-2*nf)*H1(1)+((nf*(1048-800*v)+270*(-13+8*v))* &
    (u*u*u*u*u)+(540-200*nf)*(u*u*u*u*u*u)+ &
    &  12*(omv*omv)*(nf*(10-17*v+10*(v*v))- &
    &  3*(67-131*v+67*(v*v)))- &
    &  4*(u*u*u*u)*(-9*(173-330*v+90*(v*v))+ &
    nf*(518-848*v+300*(v*v)))+ &
    &  2*u*(-1+v)*(2*nf*(120-123*v-57*(v*v)+62*(v*v*v))- &
    &  9*(536-940*v+317*(v*v)+75*(v*v*v)))+ &
    &  2*(u*u*u)*(nf*(780-2430*v+1944*(v*v)-400*(v*v*v))+ &
    &  9*(41+718*v-810*(v*v)+120*(v*v*v)))- &
    &  4*(u*u)*(-27*(-89+128*v+44*(v*v)-70*(v*v*v)+ &
    &  5*(v*v*v*v))+ &
    nf*(144-435*v+846*(v*v)-448*(v*v*v)+50*(v*v*v*v)))- &
    &  4374*zeta2*(w*w*w*w))/(18._dp*(w*w*w*w))+9*(H1(1)*H1(1))))+ &
    G1(2)*(G2(0,2)*(-33._dp/2._dp+nf-9*H1(0)-9*H1(1))+ &
    G2(2,0)*(-33._dp/2._dp+nf-9*H1(0)-9*H1(1))+ &
    G2(1,2)*(33-2*nf+18*H1(0)+18*H1(1))+ &
    G2(3,0)*(33-2*nf+18*H1(0)+18*H1(1))+(198-12*nf)*H2(0,0)+ &
    (33-2*nf)*H2(0,1)+(99-6*nf)*H2(1,0)-36*H3(0,0,0)-72*H3(1,1,0)+ &
    (-3*omv*u*(-562-189*zeta3+v*(382+189*zeta3))+ &
    &  10*(-27+10*nf)*(omv*omv)+ &
    u*(2*u*(nf*nf)-3*u*(911-944*v+u*(-562+180*v-189*zeta3)+ &
    &  378*zeta3-378*v*zeta3+90*(u*u)+90*(v*v))+ &
    nf*(-302+404*v+(-302+200*v)*(u*u)+100*(u*u*u)- &
    &  102*(v*v)+u*(365-404*v+100*(v*v)))))/(9._dp*u*(w*w))+ &
    H1(1)*(-18*H2(0,0)+(-4*u*(411-23*nf-402*v+20*nf*v)+ &
    (804-40*nf)*(omv*omv)+(804-40*nf)*(u*u)+594*zeta2*(w*w)) &
    /(12._dp*(w*w)))+G1(3)* &
    (H1(0)*(-33+2*nf-18*H1(1))-18*(H1(0)*H1(0)))+ &
    (-165._dp/2._dp+5*nf-(81*H1(1))/2._dp)*(H1(0)*H1(0))- &
    (81*(H1(0)*H1(0)*H1(0)))/2._dp+ &
    H1(0)*(201+75/u+75*u+9/w+(225*zeta2)/2._dp-18*H2(0,0)-30/(u*u)- &
    &  30*(u*u)+(nf*(-124/u-124*u+100/(u*u)+100*(u*u)- &
    (9*(u*(-23+10*u+20*v)+10*(omv*omv)))/(w*w)))/9._dp- &
    &  9/(w*w)+(9*v)/(w*w)-(45*(H1(1)*H1(1)))/2._dp)+ &
    (33._dp/2._dp-nf)*(H1(1)*H1(1))-(9*(H1(1)*H1(1)*H1(1)))/2._dp)+ &
    G1(0)*(G2(0,2)*(-33._dp/2._dp+nf-9*H1(0))+G2(2,0)*(-33._dp/2._dp+nf-9*H1(0))+ &
    G2(1,2)*(33-2*nf+18*H1(0))+G2(3,0)*(33-2*nf+18*H1(0))+ &
    (231._dp/2._dp-7*nf)*H2(0,1)+(66-4*nf)*H2(1,1)-108*H3(0,0,0)- &
    &  18*H3(0,0,1)-18*H3(0,1,0)-18*H3(1,0,0)+ &
    (H1(1)*(4*nf*(-30*(omv*omv)+ &
    u*(-55+188*v+4*(-81+50*v)*(u*u)+100*(u*u*u)- &
    &  124*(v*v)+2*u*(159-224*v+50*(v*v))))+ &
    &  9*(-4*(-67*(omv*omv)+ &
    u*(62+16*v+15*(-9+4*v)*(u*u)+30*(u*u*u)-75*(v*v)+ &
    u*(113-210*v+30*(v*v))))+486*zeta2*(w*w))))/ &
    (36._dp*(w*w))+(H2(1,0)*(9*(-69+4*nf)*u*(omv*omv)- &
    &  7*(-33+2*nf)*(omv*omv*omv)-9*omv*(-63+4*nf)*(u*u)+ &
    &  5*(-33+2*nf)*(u*u*u)))/(2._dp*(w*w*w))+ &
    (H2(0,0)*(3*(-318+19*nf)*u*(omv*omv)+ &
    (330-20*nf)*(omv*omv*omv)+3*omv*(309-19*nf)*(u*u)+ &
    &  9*(-33+2*nf)*(u*u*u)))/(w*w*w)+ &
    ((-4*nf*(173-200*v+36*zeta2)+ &
    &  6*(248-360*v+396*zeta2+27*zeta3))*(u*u*u*u*u)+ &
    &  20*(-27+10*nf)*(u*u*u*u*u*u)+ &
    (2*nf*(-19+54*zeta2)-3*(-386+594*zeta2+54*zeta3))* &
    (omv*omv*omv)*(v*v)+ &
    u*u*u*u*(6*(-492+574*v+72*(-15+26*v)*zeta2-81*zeta3+ &
    &  135*v*zeta3-540*(v*v))+ &
    nf*(6*(63-111*v)*zeta2+10*(131-223*v+120*(v*v))))+ &
    omv*u*v*(nf*(416-484*v+40*(v*v)+ &
    &  54*zeta2*(4-15*v+11*(v*v)))- &
    &  3*(54*zeta2*(22-81*v+59*(v*v))+ &
    &  2*(82+v*(613-189*zeta3)+54*zeta3+ &
    (-749+135*zeta3)*(v*v))))+ &
    u*u*u*(3*(6*(65+27*zeta3)-24*v*(65+27*zeta3)+ &
    (82+540*zeta3)*(v*v)+18*zeta2*(111-462*v+395*(v*v))- &
    &  720*(v*v*v))+nf* &
    (-18*zeta2*(21-84*v+71*(v*v))+ &
    &  4*(-111+843*v-625*(v*v)+200*(v*v*v))))+ &
    u*u*(3*(54*zeta2*(-11+107*v-221*(v*v)+125*(v*v*v))- &
    &  2*(131+27*zeta3-243*v*(1+zeta3)+ &
    &  9*(-55+54*zeta3)*(v*v)+(841-270*zeta3)*(v*v*v)+ &
    &  90*(v*v*v*v)))+ &
    nf*(-54*zeta2*(-2+20*v-41*(v*v)+23*(v*v*v))+ &
    &  2*(113-729*v+1350*(v*v)-520*(v*v*v)+100*(v*v*v*v)))))/ &
    (18._dp*(omw*omw)*(w*w*w))+ &
    ((-27*H1(0))/2._dp-(9*H1(1))/2._dp)*(G1(2)*G1(2))+ &
    G1(3)*((-33+2*nf)*H1(0)-18*(H1(0)*H1(0)))+ &
    G1(1)*(H1(0)*(-33+2*nf-18*H1(1))+(-33+2*nf)*H1(1)- &
    &  18*(H1(0)*H1(0)))+(-165._dp/2._dp+5*nf-(81*H1(1))/2._dp)*(H1(0)*H1(0))- &
    (81*(H1(0)*H1(0)*H1(0)))/2._dp+ &
    H1(0)*((-165._dp/2._dp+5*nf)*H1(1)-18*H2(0,0)+ &
    ((7236-360*nf)*(omv*omv*omv*omv)+ &
    &  4*(135*(13-8*v)+nf*(-524+400*v))*(u*u*u*u*u)+ &
    &  40*(-27+10*nf)*(u*u*u*u*u*u)+ &
    &  4*(u*u*u*u)*(-27*(93-220*v+60*(v*v))+ &
    &  2*nf*(503-848*v+300*(v*v)))- &
    &  4*u*(-1+v)*(-27*(269-495*v+193*(v*v)+25*(v*v*v))+ &
    nf*(369-513*v-12*(v*v)+124*(v*v*v)))+ &
    &  4*(u*u*u)*(-54*(52+76*v-135*(v*v)+20*(v*v*v))+ &
    nf*(-651+2328*v-1944*(v*v)+400*(v*v*v)))+ &
    &  8*(u*u)*(-27*(-157+262*v-21*(v*v)-70*(v*v*v)+ &
    &  5*(v*v*v*v))+ &
    nf*(45-255*v+774*(v*v)-448*(v*v*v)+50*(v*v*v*v)))+ &
    &  9234*zeta2*(w*w*w*w))/(36._dp*(w*w*w*w))-(45*(H1(1)*H1(1)))/2._dp)+ &
    (-33._dp/2._dp+nf)*(H1(1)*H1(1))+ &
    G1(2)*((33._dp/2._dp-nf)*H1(1)+H1(0)*(66-4*nf+36*H1(1))- &
    ((33-2*nf)*(33-2*nf))/36._dp+45*(H1(0)*H1(0))+9*(H1(1)*H1(1)))- &
    (9*(H1(1)*H1(1)*H1(1)))/2._dp)+((-33+2*nf)*(H1(1)*H1(1)*H1(1)))/6._dp+ &
    H1(0)*((99-6*nf)*H2(0,0)+H1(1)* &
    (18*H2(0,0)+(6*u*(2637-266*nf-2619*v+260*nf*v)+ &
    (-7857+780*nf)*(omv*omv)+(-7857+780*nf)*(u*u)- &
    (486*zeta2+20*(nf*nf))*(w*w))/(36._dp*(w*w)))+ &
    (-2*u*(omv*omv*omv)*(v*v)* &
    (2*u*(nf*nf)-9*(-259+v*(590-54*zeta3)+27*zeta3+ &
    (-256+27*zeta3)*(v*v)))+ &
    &  20*(-27+10*nf)*(omv*omv*omv*omv)*(v*v*v)+ &
    u*(omv*omv)*(4*(-1+u+3*v)*(nf*nf)*(u*u*u)- &
    &  9*(6*zeta2*((-141+287*v)*(u*u*u*u)+55*(u*u*u*u*u)+ &
    u*u*u*(123-582*v+595*(v*v))+ &
    &  3*(v*v)*(-11+45*v-54*(v*v)+22*(v*v*v))+ &
    &  3*u*v*(-22+131*v-208*(v*v)+106*(v*v*v)))+ &
    &  2*u*v*(272-54*zeta3+v*(-1636+243*zeta3)+ &
    (2671-324*zeta3)*(v*v)+(-1298+135*zeta3)*(v*v*v)+ &
    &  30*(v*v*v*v)))+ &
    nf*(18*zeta2*((-24+53*v)*(u*u*u*u)+10*(u*u*u*u*u)+ &
    u*u*u*(24-96*v+109*(v*v))+ &
    &  3*(v*v)*(-2+9*v-9*(v*v)+4*(v*v*v))+ &
    &  3*u*v*(-4+26*v-34*(v*v)+19*(v*v*v)))+ &
    &  2*v*(5*v*(55-188*v+161*(v*v)-28*(v*v*v))+ &
    u*(4-728*v+1731*(v*v)-885*(v*v*v)+100*(v*v*v*v))))) &
    +2*(u*u*u*u)*(-9*(-262+v*(2117-486*zeta3)+81*zeta3+ &
    (-5495+999*zeta3)*(v*v)+(6215-864*zeta3)*(v*v*v)+ &
    u*u*(-91+v*(272-54*zeta3)+27*zeta3+ &
    (-259+27*zeta3)*(v*v)+30*(v*v*v))+ &
    (-2803+270*zeta3)*(v*v*v*v)+ &
    u*(262-81*zeta3+v*(-1291+297*zeta3)+ &
    (2180-351*zeta3)*(v*v)+(-1367+135*zeta3)*(v*v*v)+ &
    &  120*(v*v*v*v))+180*(v*v*v*v*v))+ &
    nf*(210-419*v-736*(v*v)+3191*(v*v*v)+ &
    u*u*(57+4*v-275*(v*v)+100*(v*v*v))-2960*(v*v*v*v)+ &
    &  2*u*(-105+168*v+368*(v*v)-745*(v*v*v)+200*(v*v*v*v))+ &
    &  600*(v*v*v*v*v)))+ &
    omv*(u*u*u)*(-4*v*(nf*nf)*(2-5*v+3*(v*v))+ &
    nf*(-114+558*v+2030*(v*v)-5570*(v*v*v)+4580*(v*v*v*v)- &
    &  54*zeta2*(2-27*v+74*(v*v)-86*(v*v*v)+37*(v*v*v*v))- &
    &  800*(v*v*v*v*v))+ &
    &  9*(18*zeta2*(11-138*v+428*(v*v)-506*(v*v*v)+ &
    &  205*(v*v*v*v))+ &
    &  2*(-91+27*zeta3-30*v*(-40+9*zeta3)+ &
    (-4295+729*zeta3)*(v*v)+(5599-756*zeta3)*(v*v*v)+ &
    (-2677+270*zeta3)*(v*v*v*v)+120*(v*v*v*v*v)))))/ &
    (18._dp*u*(omv*omv)*(omw*omw)*(w*w*w))+(-33._dp/2._dp+nf)*(H1(1)*H1(1))+ &
    (9*(H1(1)*H1(1)*H1(1)))/2._dp)+ &
    G1(1)*(6*(-44+(8*nf)/3._dp)*zeta2-72*zeta3+6*(-33+2*nf)*H2(0,0)+ &
    &  6*(-33+2*nf)*H2(0,1)+(-132+8*nf)*H2(1,0)+(-132+8*nf)*H2(1,1)+ &
    &  36*H3(0,0,0)+36*H3(0,0,1)+36*H3(0,1,0)+36*H3(1,0,0)+72*H3(1,1,0)+ &
    H1(1)*(-134-75/u-75*u-6/w-171*zeta2+30/(u*u)+30*(u*u)+ &
    &  6/(w*w)-(6*v)/(w*w)+ &
    (2*nf*(62/u+62*u-50/(u*u)-50*(u*u)+ &
    (-9+9*v+9*w+30*(w*w))/(w*w)))/9._dp)+ &
    (9*H1(0)+9*H1(1))*(G1(2)*G1(2))+(99-6*nf+45*H1(1))*(H1(0)*H1(0))+ &
    &  27*(H1(0)*H1(0)*H1(0))+G1(2)* &
    (H1(0)*(-33+2*nf-36*H1(1))+(-33+2*nf)*H1(1)-18*(H1(0)*H1(0))- &
    &  18*(H1(1)*H1(1)))+(33-2*nf)*(H1(1)*H1(1))+ &
    H1(0)*(-134-75/u-75*u-6/w-99*zeta2+(132-8*nf)*H1(1)+ &
    &  30/(u*u)+30*(u*u)+6/(w*w)-(6*v)/(w*w)+ &
    (2*nf*(62/u+62*u-50/(u*u)-50*(u*u)+ &
    (-9+9*v+9*w+30*(w*w))/(w*w)))/9._dp+27*(H1(1)*H1(1)))+ &
    &  9*(H1(1)*H1(1)*H1(1)))+(9*(H1(1)*H1(1)*H1(1)*H1(1)))/8._dp
    return
    end function Hbeta_4a2re
          
          
    function Hbeta_4a2im(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    use constants
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=36*pi*G3(0,0,1)+18*pi*G3(0,0,2)-36*pi*G3(0,0,3)+36*pi*G3(0,1,0)+ &
    &  18*pi*G3(0,2,0)-18*pi*G3(0,2,1)+18*pi*G3(0,2,2)-18*pi*G3(0,2,3)- &
    &  18*pi*G3(0,3,0)-18*pi*G3(0,3,2)+72*pi*G3(0,3,3)-36*pi*G3(1,0,2)- &
    &  72*pi*G3(1,1,1)+72*pi*G3(1,1,2)-36*pi*G3(1,2,0)+72*pi*G3(1,2,1)- &
    &  36*pi*G3(1,2,2)+18*pi*G3(2,0,0)-18*pi*G3(2,0,1)+18*pi*G3(2,0,2)- &
    &  18*pi*G3(2,0,3)-18*pi*G3(2,1,0)+72*pi*G3(2,1,1)-18*pi*G3(2,1,2)+ &
    &  18*pi*G3(2,2,0)-36*pi*G3(2,2,1)+36*pi*G3(2,2,3)+36*pi*G3(2,3,2)- &
    &  36*pi*G3(3,0,0)-36*pi*G3(3,0,2)+72*pi*G3(3,0,3)-36*pi*G3(3,2,0)+ &
    &  72*pi*G3(3,3,0)-72*pi*G3(3,3,3)-18*pi*G2(0,1)*H1(0)+ &
    G2(3,2)*((-33+2*nf)*pi+18*pi*H1(0))+ &
    G2(2,2)*((66-4*nf)*pi+36*pi*H1(0))+36*pi*G2(1,2)*H1(1)+ &
    G2(1,1)*(2*(-33+2*nf)*pi-72*pi*H1(1))+ &
    G2(2,3)*((-33+2*nf)*pi-18*pi*H1(0)-36*pi*H1(1))+ &
    G2(0,0)*((66-4*nf)*pi-36*pi*H1(0)-18*pi*H1(1))+ &
    G2(0,2)*((66-4*nf)*pi-9*pi*H1(0)-18*pi*H1(1))+ &
    G2(2,0)*((66-4*nf)*pi-9*pi*H1(0)-18*pi*H1(1))+ &
    G2(2,1)*(3*(-33+2*nf)*pi-18*pi*H1(0)+18*pi*H1(1))+ &
    G2(1,0)*(18*pi*H1(0)+36*pi*H1(1))+ &
    G2(3,0)*((-33+2*nf)*pi+72*pi*H1(0)+36*pi*H1(1))+ &
    (7*(33-2*nf)*pi*H2(0,1))/2._dp+(66-4*nf)*pi*H2(1,1)-108*pi*H3(0,0,0)- &
    &  18*pi*H3(0,0,1)-18*pi*H3(0,1,0)-18*pi*H3(1,0,0)+ &
    (pi*H1(1)*(4*nf*(-30*(omv*omv)+ &
    u*(-55+188*v+4*(-81+50*v)*(u*u)+100*(u*u*u)-124*(v*v)+ &
    &  2*u*(159-224*v+50*(v*v))))+ &
    &  9*(-4*(-67*(omv*omv)+ &
    u*(62+16*v+15*(-9+4*v)*(u*u)+30*(u*u*u)-75*(v*v)+ &
    u*(113-210*v+30*(v*v))))+270*zeta2*(w*w))))/(36._dp*(w*w)) &
    +G2(0,3)*(54*pi*H1(0)+18*pi*H1(1)+ &
    (pi*(-27*(-17+nf)*u*(omv*omv)+5*(-33+2*nf)*(omv*omv*omv)+ &
    &  27*omv*(-16+nf)*(u*u)+4*(33-2*nf)*(u*u*u)))/(w*w*w))+ &
    G2(3,3)*(-72*pi*H1(0)-(2*pi*(-33+63*v-9*u*(-7+8*v)-36*(u*u)- &
    &  36*(v*v)+nf*(2-3*v+u*(-3+6*v)+3*(u*u)+3*(v*v))))/ &
    (w*w*w))-(pi*H2(1,0)*(-9*(-69+4*nf)*u*(omv*omv)+ &
    &  7*(-33+2*nf)*(omv*omv*omv)+9*omv*(-63+4*nf)*(u*u)+ &
    &  5*(33-2*nf)*(u*u*u)))/(2._dp*(w*w*w))- &
    (pi*H2(0,0)*(-9*(-84+5*nf)*u*(omv*omv)+ &
    &  8*(-33+2*nf)*(omv*omv*omv)+9*omv*(-81+5*nf)*(u*u)+ &
    &  7*(33-2*nf)*(u*u*u)))/(w*w*w)+ &
    (pi*(20*(-27+10*nf)*(omv*omv*omv)*(v*v)+ &
    &  2*u*v*(omv*omv)*(2*nf* &
    (100-3*v*(113+12*zeta2)+(89+36*zeta2)*(v*v))- &
    &  3*(180-v*(752+396*zeta2+27*zeta3)+ &
    (302+396*zeta2+27*zeta3)*(v*v)))+ &
    omv*(u*u)*(6*(-90+2*v*(428+396*zeta2+27*zeta3)- &
    (2366+2664*zeta2+189*zeta3)*(v*v)+ &
    (1474+1872*zeta2+135*zeta3)*(v*v*v)-90*(v*v*v*v))+ &
    nf*(200-4*v*(455+72*zeta2)+3*(1314+318*zeta2)*(v*v)- &
    &  3*(650+222*zeta2)*(v*v*v)+200*(v*v*v*v)))+ &
    u*u*u*(-3*(18*zeta2*(-44+372*v+3*(-37+59*v)*(u*u)+ &
    &  33*(u*u*u)-723*(v*v)+3*u*(40-154*v+125*(v*v))+ &
    &  395*(v*v*v))- &
    &  2*(248+27*zeta3-27*v*(74+9*zeta3)+ &
    (-131+90*v-27*zeta3)*(u*u*u)+ &
    &  9*(463+54*zeta3)*(v*v)+ &
    u*u*(195+81*zeta3-v*(1051+135*zeta3)+360*(v*v))- &
    (3011+270*zeta3)*(v*v*v)+ &
    u*(-492-81*zeta3+12*v*(149+27*zeta3)- &
    (2669+270*zeta3)*(v*v)+540*(v*v*v))+360*(v*v*v*v))) &
    +nf*(18*zeta2*(-8+66*v+3*(-7+11*v)*(u*u)+6*(u*u*u)- &
    &  129*(v*v)+3*u*(7-28*v+23*(v*v))+71*(v*v*v))- &
    &  2*(346-2102*v+(-113+100*v)*(u*u*u)+3412*(v*v)+ &
    u*u*(222-880*v+400*(v*v))-2270*(v*v*v)+ &
    u*(-655+1622*v-2140*(v*v)+600*(v*v*v))+400*(v*v*v*v)) &
    ))))/(18._dp*u*(omw*omw)*(w*w*w))+ &
    ((-33._dp/2._dp+nf)*pi-9*pi*G1(1)+(9*pi*G1(2))/2._dp-9*pi*G1(3)- &
    (27*pi*H1(0))/2._dp-(9*pi*H1(1))/2._dp)*(G1(0)*G1(0))+ &
    (9*pi*(G1(0)*G1(0)*G1(0)))/2._dp+ &
    ((-33._dp/2._dp+nf)*pi-9*pi*G1(3)-(27*pi*H1(0))/2._dp-(27*pi*H1(1))/2._dp)* &
    (G1(2)*G1(2))+(9*pi*(G1(2)*G1(2)*G1(2)))/2._dp+ &
    ((7*(-33+2*nf)*pi)/2._dp-(45*pi*H1(1))/2._dp)*(H1(0)*H1(0))- &
    (27*pi*(H1(0)*H1(0)*H1(0)))/2._dp+(-33._dp/2._dp+nf)*pi*(H1(1)*H1(1))+ &
    H1(0)*((5*(-33+2*nf)*pi*H1(1))/2._dp+18*pi*H2(0,0)+ &
    (pi*((2412-120*nf)*(omv*omv*omv*omv)+ &
    &  4*(135*(13-8*v)+nf*(-524+400*v))*(u*u*u*u*u)+ &
    &  40*(-27+10*nf)*(u*u*u*u*u*u)+ &
    &  4*(u*u*u*u)*(-9*(413-660*v+180*(v*v))+ &
    &  2*nf*(533-848*v+300*(v*v)))- &
    &  4*u*(-1+v)*(-9*(265-407*v+43*(v*v)+75*(v*v*v))+ &
    nf*(111-15*v-252*(v*v)+124*(v*v*v)))+ &
    &  4*(u*u*u)*(-18*(-115+496*v-405*(v*v)+60*(v*v*v))+ &
    nf*(-909+2568*v-1944*(v*v)+400*(v*v*v)))+ &
    &  8*(u*u)*(-27*(-21-8*v+113*(v*v)-70*(v*v*v)+5*(v*v*v*v))+ &
    nf*(243-633*v+954*(v*v)-448*(v*v*v)+50*(v*v*v*v)))+ &
    &  2430*zeta2*(w*w*w*w)))/(36._dp*(w*w*w*w))-(27*pi*(H1(1)*H1(1)))/2._dp)+ &
    G1(3)*(-108*pi*H2(0,0)-18*pi*H2(0,1)-18*pi*H2(1,0)+ &
    H1(0)*(-18*pi*H1(1)-(omw*pi* &
    (-36+u*(63-66*v)+63*v-33*(u*u)-33*(v*v)+ &
    nf*(3-3*v+u*(-3+4*v)+2*(u*u)+2*(v*v))))/(w*w*w))+ &
    (pi*(2*(135*(13-8*v)+nf*(-524+400*v))*(u*u*u*u*u)+ &
    &  20*(-27+10*nf)*(u*u*u*u*u*u)- &
    &  12*(omv*omv)*(nf*(10-17*v+10*(v*v))- &
    &  3*(67-131*v+67*(v*v)))+ &
    &  4*(u*u*u*u)*(-9*(173-330*v+90*(v*v))+ &
    nf*(518-848*v+300*(v*v)))- &
    &  2*u*(-1+v)*(2*nf*(120-123*v-57*(v*v)+62*(v*v*v))- &
    &  9*(536-940*v+317*(v*v)+75*(v*v*v)))+ &
    &  2*(u*u*u)*(-9*(41+718*v-810*(v*v)+120*(v*v*v))+ &
    &  2*nf*(-390+1215*v-972*(v*v)+200*(v*v*v)))+ &
    &  4*(u*u)*(-27*(-89+128*v+44*(v*v)-70*(v*v*v)+5*(v*v*v*v))+ &
    nf*(144-435*v+846*(v*v)-448*(v*v*v)+50*(v*v*v*v)))+ &
    &  2430*zeta2*(w*w*w*w)))/(18._dp*(w*w*w*w))-9*pi*(H1(0)*H1(0))- &
    &  9*pi*(H1(1)*H1(1)))+G1(1)*((-33+2*nf)*pi*H1(0)+(-33+2*nf)*pi*H1(1)+ &
    G1(2)*((33-2*nf)*pi+18*pi*H1(0)+18*pi*H1(1))-36*pi*H2(0,0)- &
    &  18*pi*H2(0,1)-18*pi*H2(1,0)-36*pi*H2(1,1)+ &
    (pi*(4*nf*(-62/u-62*u+50/(u*u)+50*(u*u)+ &
    (9-9*v-9*w-30*(w*w))/(w*w))+ &
    &  9*(270*zeta2+2*(75/u+75*u-30/(u*u)-30*(u*u)+ &
    (2*(-3+3*v+3*w+67*(w*w)))/(w*w)))))/18._dp- &
    &  9*pi*(G1(2)*G1(2))-9*pi*(H1(0)*H1(0))-9*pi*(H1(1)*H1(1)))+ &
    G1(0)*(G1(1)*((33-2*nf)*pi+18*pi*H1(0))+ &
    G1(3)*((33-2*nf)*pi+18*pi*H1(0))+(-33+2*nf)*pi*H1(1)+ &
    G1(2)*((-33+2*nf)*pi-18*pi*H1(0)-9*pi*H1(1))+ &
    H1(0)*((-33+2*nf)*pi+9*pi*H1(1))+90*pi*H2(0,0)+36*pi*H2(0,1)+ &
    &  36*pi*H2(1,0)+18*pi*H2(1,1)- &
    (pi*(4*nf*(-30*(omv*omv)+ &
    u*(-55+188*v+4*(-81+50*v)*(u*u)+100*(u*u*u)- &
    &  124*(v*v)+2*u*(159-224*v+50*(v*v))))+ &
    &  9*(-4*(-67*(omv*omv)+ &
    u*(62+16*v+15*(-9+4*v)*(u*u)+30*(u*u*u)-75*(v*v)+ &
    u*(113-210*v+30*(v*v))))+270*zeta2*(w*w))))/ &
    (36._dp*(w*w))+(9*pi*(G1(2)*G1(2)))/2._dp+(9*pi*(H1(0)*H1(0)))/2._dp+ &
    (9*pi*(H1(1)*H1(1)))/2._dp)+G1(2)* &
    ((-33+2*nf)*pi*H1(1)+G1(3)*((33-2*nf)*pi+18*pi*H1(0)+18*pi*H1(1))+ &
    H1(0)*((66-4*nf)*pi+27*pi*H1(1))+18*pi*H2(0,0)-18*pi*H2(0,1)- &
    &  18*pi*H2(1,0)-(pi*(4*nf*(-30*(omv*omv)+ &
    u*(-55+188*v+4*(-81+50*v)*(u*u)+100*(u*u*u)- &
    &  124*(v*v)+2*u*(159-224*v+50*(v*v))))+ &
    &  9*(-4*(-67*(omv*omv)+ &
    u*(62+16*v+15*(-9+4*v)*(u*u)+30*(u*u*u)-75*(v*v)+ &
    u*(113-210*v+30*(v*v))))+270*zeta2*(w*w))))/ &
    (36._dp*(w*w))+(9*pi*(H1(0)*H1(0)))/2._dp+(27*pi*(H1(1)*H1(1)))/2._dp)- &
    (9*pi*(H1(1)*H1(1)*H1(1)))/2._dp
    return
    end function Hbeta_4a2im
          
          
!**********************************************************************
! egion4a:q+g-->V+q*
!*********************************************************************/
          
    function Hgamma_4a0re(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=1
    return
    end function Hgamma_4a0re
          
          
    function Hgamma_4a0im(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=0
    return
    end function Hgamma_4a0im
          
          
    function Hgamma_4a1re(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    use constants
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=(588-40*nf-60*u-60*v+846*zeta2)/36._dp-3*G2(0,2)+6*G2(1,2)- &
    &  3*G2(2,0)+(8*G2(3,0))/3._dp-(8*G1(3)*H1(0))/3._dp+G1(1)*(-6*H1(0)-6*H1(1))+ &
    G1(2)*(-7+(2*nf)/3._dp+(17*H1(0))/3._dp-H1(1)/3._dp)+ &
    H1(0)*(7-(2*nf)/3._dp+H1(1)/3._dp)+(7-(2*nf)/3._dp)*H1(1)+ &
    G1(0)*((17*H1(0))/3._dp+3*H1(1))-(8*H2(0,0))/3._dp-(10*H2(1,0))/3._dp- &
    (3*(G1(0)*G1(0)))/2._dp+(G1(2)*G1(2))/6._dp-(17*(H1(0)*H1(0)))/6._dp+(H1(1)*H1(1))/6._dp
    return
    end function Hgamma_4a1re
          
          
    function Hgamma_4a1im(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    use constants
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=-3*pi*G1(0)+6*pi*G1(1)-3*pi*G1(2)+(8*pi*G1(3))/3._dp-(pi*H1(0))/3._dp+3*pi*H1(1)
    return
    end function Hgamma_4a1im
          
          
    function Hgamma_4a2re(u,v, &
    H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    use constants
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=(72-6*nf)*G3(2,2,0)-(16*G4(0,0,1,2))/9._dp+18*G4(0,0,2,2)+ &
    (16*G4(0,0,3,0))/9._dp+26*G4(0,1,0,2)+26*G4(0,1,2,0)-36*G4(0,1,2,2)+ &
    &  18*G4(0,2,0,2)-(92*G4(0,2,1,2))/9._dp+18*G4(0,2,2,0)-(152*G4(0,2,3,0))/9._dp+ &
    (8*G4(0,3,0,2))/9._dp+(8*G4(0,3,2,0))/9._dp+(128*G4(0,3,3,0))/9._dp- &
    &  36*G4(1,0,2,2)-72*G4(1,1,1,2)+72*G4(1,1,2,2)-36*G4(1,2,0,2)+ &
    &  72*G4(1,2,1,2)-36*G4(1,2,2,0)+18*G4(2,0,0,2)-18*G4(2,0,1,2)+ &
    &  18*G4(2,0,2,0)-18*G4(2,0,3,0)-18*G4(2,1,0,2)+72*G4(2,1,1,2)- &
    &  18*G4(2,1,2,0)+18*G4(2,2,0,0)-36*G4(2,2,1,2)+36*G4(2,2,3,0)- &
    (244*G4(2,3,0,0))/9._dp+(164*G4(2,3,0,2))/9._dp+(164*G4(2,3,2,0))/9._dp+ &
    (80*G4(2,3,3,0))/9._dp-(164*G4(3,0,0,2))/9._dp-(70*G4(3,0,1,2))/9._dp- &
    (164*G4(3,0,2,0))/9._dp+(16*G4(3,0,3,0))/3._dp-(164*G4(3,2,0,0))/9._dp- &
    (70*G4(3,2,1,2))/9._dp-(80*G4(3,2,3,0))/9._dp+(128*G4(3,3,0,0))/9._dp- &
    (128*G4(3,3,3,0))/9._dp+G3(0,2,2)*(72-6*nf-36*H1(0))+ &
    G3(2,0,2)*(72-6*nf-(160*H1(0))/9._dp)-(16*G3(0,0,3)*H1(0))/9._dp+ &
    (152*G3(0,2,3)*H1(0))/9._dp-(8*G3(0,3,2)*H1(0))/9._dp-(128*G3(0,3,3)*H1(0))/9._dp+ &
    &  18*G3(2,0,3)*H1(0)-36*G3(2,2,3)*H1(0)-(164*G3(2,3,2)*H1(0))/9._dp- &
    (80*G3(2,3,3)*H1(0))/9._dp-(16*G3(3,0,3)*H1(0))/3._dp+(80*G3(3,2,3)*H1(0))/9._dp+ &
    (128*G3(3,3,3)*H1(0))/9._dp+G3(1,2,2)*(12*(-12+nf)+72*H1(0))+ &
    G3(1,2,1)*(-72*H1(0)-72*H1(1))+G3(2,1,1)*(-72*H1(0)-72*H1(1))+ &
    G3(0,1,0)*(-26*H1(0)-26*H1(1))+G3(0,0,1)*((16*H1(0))/9._dp+(16*H1(1))/9._dp)+ &
    G3(3,0,1)*((70*H1(0))/9._dp+(70*H1(1))/9._dp)+ &
    G3(3,2,1)*((70*H1(0))/9._dp+(70*H1(1))/9._dp)+ &
    G3(0,2,1)*((92*H1(0))/9._dp+(92*H1(1))/9._dp)+G3(2,0,1)*(18*H1(0)+18*H1(1))+ &
    G3(2,1,0)*(18*H1(0)+18*H1(1))+G3(2,2,1)*(36*H1(0)+36*H1(1))+ &
    G3(1,1,1)*(72*H1(0)+72*H1(1))+(4*(-201+14*nf)*H3(1,1,0))/9._dp+ &
    (64*H4(0,0,0,0))/3._dp+20*H4(0,0,0,1)+(100*H4(0,0,1,0))/9._dp+ &
    (100*H4(0,1,0,0))/9._dp+(70*H4(0,1,0,1))/9._dp+(10*H4(0,1,1,0))/9._dp+ &
    (250*H4(1,0,0,0))/9._dp-(20*H4(1,0,0,1))/9._dp-(10*H4(1,0,1,0))/3._dp- &
    (170*H4(1,1,0,0))/9._dp-(230*H4(1,1,0,1))/9._dp+40*H4(1,1,1,0)+ &
    H2(1,1)*(455._dp/4._dp-(61*nf)/3._dp+(104*zeta2)/9._dp+(8*(nf*nf))/9._dp)+ &
    G2(2,2)*(455._dp/4._dp-(61*nf)/3._dp+18*zeta2+12*(-12+nf)*H1(0)-36*H2(0,1)+ &
    &  36*H2(1,0)+(8*(nf*nf))/9._dp)+ &
    (H2(0,1)*(19305-3150*nf-480*v-(60*u*v)/omv-7464*zeta2+ &
    &  120*(nf*nf)))/108._dp+G2(1,1)* &
    (-72*zeta2+72*H2(1,0)+72*H2(1,1)+ &
    &  2*H1(0)*(-21+2*nf-7*u-5*(u*u))+ &
    &  2*H1(1)*(-21+2*nf-7*u-5*(u*u)))+ &
    G3(1,0,2)*(-18+nf-7*u+18*H1(0)+36*H1(1)-5*(u*u))+ &
    G3(1,2,0)*(-18+nf-7*u+18*H1(0)+36*H1(1)-5*(u*u))+ &
    G3(1,1,2)*(42-4*nf+14*u-72*H1(0)-72*H1(1)+10*(u*u))+ &
    (H3(0,1,0)*(-813+50*nf+(63+20/omv)*u+(45-5/(omv*omv))*(u*u)))/ &
    &  9._dp+H3(0,0,0)*(-460._dp/3._dp+8*nf+u*(4+20/(9-9*v))-3*v+ &
    (5-5/(9._dp*(omv*omv)))*(u*u))+ &
    G2(1,0)*(108*zeta2-18*H2(0,1)-18*H2(1,0)-36*H2(1,1)+ &
    H1(0)*(18-nf+7*u+5*(u*u))+H1(1)*(18-nf+7*u+5*(u*u)))+ &
    G3(2,3,0)*(4*nf+18*H1(0)-(164*H1(1))/9._dp+ &
    (-20*omv*u-3*(107+9*u+9*v)*(omv*omv)+5*(u*u))/ &
    (9._dp*(omv*omv)))+(H3(1,0,0)* &
    (96*nf+(-20*omv*u-3*(437+9*u+9*v)*(omv*omv)+5*(u*u))/ &
    (omv*omv)))/9._dp+G2(2,3)* &
    ((-332*zeta2)/9._dp-(80*H2(0,0))/9._dp+(164*H2(0,1))/9._dp+(164*H2(1,0))/9._dp+ &
    H1(0)*(107._dp/3._dp-4*nf+u*(3+20/(9-9*v))+3*v- &
    (5*(u*u))/(9._dp*(omv*omv))))+ &
    G2(3,2)*((622*zeta2)/9._dp-(100*H2(0,0))/9._dp-(70*H2(0,1))/9._dp- &
    (10*H2(1,0))/9._dp+(H1(0)*(-12-83*v-u*(83+80*v)-40*(u*u)- &
    &  40*(v*v)))/9._dp)+G3(2,1,2)* &
    (36*H1(0)+(-1431+108*nf+u*(-83+20/omv-80*v)-83*v+ &
    (-40-5/(omv*omv))*(u*u)-40*(v*v))/9._dp)+ &
    G3(0,1,2)*(18*H1(0)+36*H1(1)+ &
    (3+9*nf+u*(-173-20/omv-80*v)-110*v+ &
    (-85+5/(omv*omv))*(u*u)-40*(v*v))/9._dp)+ &
    G3(3,0,0)*((128*H1(0))/9._dp+(164*H1(1))/9._dp+ &
    (2*(-252+12*nf-55*v-5*u*(11+8*v)-20*(u*u)-20*(v*v)))/9._dp)+ &
    G2(0,0)*((128*H2(0,0))/9._dp+(164*H2(0,1))/9._dp+(244*H2(1,0))/9._dp+18*H2(1,1)+ &
    (480*(17+omw)-738*nf-4872*zeta2+15*(nf*nf))/108._dp+ &
    (H1(1)*(-417+27*nf-110*v-10*u*(11+8*v)-40*(u*u)-40*(v*v)))/ &
    &  9._dp+(2*H1(0)*(-252+12*nf-55*v-5*u*(11+8*v)-20*(u*u)- &
    &  20*(v*v)))/9._dp)+G3(3,3,0)* &
    ((-128*H1(0))/9._dp+(8*(-9+2*nf-7*v-u*(7+10*v)-5*(u*u)- &
    &  5*(v*v)))/9._dp)+G2(3,3)* &
    ((-256*zeta2)/9._dp+(128*H2(0,0))/9._dp- &
    (8*H1(0)*(-9+2*nf-7*v-u*(7+10*v)-5*(u*u)-5*(v*v)))/9._dp)+ &
    (H3(0,0,1)*(12+83*v+u*(83+80*v)+40*(u*u)+40*(v*v)))/9._dp+ &
    G3(2,0,0)*((-244*H1(0))/9._dp-18*H1(1)+ &
    (417-27*nf+110*v+10*u*(11+8*v)+40*(u*u)+40*(v*v))/9._dp)+ &
    G3(0,0,2)*(-20*H1(0)-18*H1(1)+ &
    (417-27*nf+110*v+10*u*(11+8*v)+40*(u*u)+40*(v*v))/9._dp)+ &
    G3(0,2,0)*(-10*H1(0)-18*H1(1)+ &
    (417-27*nf+110*v+10*u*(11+8*v)+40*(u*u)+40*(v*v))/9._dp)+ &
    G3(3,2,0)*((92*H1(0))/9._dp+(12+83*v+u*(83+80*v)+40*(u*u)+40*(v*v))/ &
    &  9._dp)+G3(3,0,2)*((34*H1(0))/3._dp+ &
    (12+83*v+u*(83+80*v)+40*(u*u)+40*(v*v))/9._dp)+ &
    G2(0,1)*(-48*zeta2-18*H2(0,1)-18*H2(1,0)-36*H2(1,1)+ &
    (H1(0)*(-3-9*nf+110*v+u*(173+20/omv+80*v)+ &
    (85-5/(omv*omv))*(u*u)+40*(v*v)))/9._dp+ &
    (H1(1)*(-3-9*nf+110*v+u*(173+20/omv+80*v)+ &
    (85-5/(omv*omv))*(u*u)+40*(v*v)))/9._dp)+ &
    G2(2,1)*(-36*zeta2+36*H2(0,1)-36*H2(1,0)+ &
    (H1(0)*(1431-108*nf+83*v+u*(83-20/omv+80*v)+ &
    &  5*(8+1/(omv*omv))*(u*u)+40*(v*v)))/9._dp+ &
    (H1(1)*(1431-108*nf+83*v+u*(83-20/omv+80*v)+ &
    &  5*(8+1/(omv*omv))*(u*u)+40*(v*v)))/9._dp)+ &
    H1(1)*((-4*H2(0,0))/3._dp-(5*H2(1,0))/3._dp+(40*(nf*nf))/27._dp+ &
    nf*(-5086._dp/81._dp+6*u+3*v+(2*v)/omu-12*u*v-(107*zeta2)/9._dp- &
    (3*(omv*omv))/u-6*(u*u)-6*(v*v)+(38*(v*v))/(3._dp*(omu*omu))) &
    +(49396-4116*v-(6480*v)/omu+u*(-6558+9768*v)-2628*zeta3+ &
    (2442*(omv*omv))/u+4884*(u*u)+4884*(v*v)- &
    (5184*(v*v))/(omu*omu)+ &
    &  6*zeta2*(1347+440*v+40*u*(11+1/omv+8*v)+ &
    (160-10/(omv*omv))*(u*u)+160*(v*v)))/108._dp)+ &
    H2(0,0)*((5*(nf*nf))/3._dp+nf*(-1253._dp/27._dp-v+(8-12*v)*(u*u)- &
    &  6*(u*u*u)+u*(-5+6*v-6*(v*v))-2*(v*v))+ &
    (-60*u*(1+u)*w+omv*(-144+24*v*(6+w*(86-157*u+407*(u*u)))+ &
    w*(34509+6624*u-21648*zeta2-6546*(u*u)+4884*(u*u*u))+ &
    &  6*(463+814*u)*w*(v*v)))/(108._dp*omv*w))+ &
    (48000*(nf*nf)+3*(6794370*zeta4+ &
    &  360*zeta2*(36497+v*(-3242-48/w)+48/w+ &
    (8728+80/omv-13024*v)*(u*u)-6512*(u*u*u)+ &
    u*(-9442+20/omv+5024*v-6512*(v*v))-3704*(v*v))+ &
    &  80*(158795-17559*zeta3+36*v*(-676+216/omu+83*zeta3)+ &
    &  36*u*(-676+(83+10/omv)*zeta3+80*v*zeta3)+ &
    &  90*zeta3*(16-1/(omv*omv))*(u*u)+1440*zeta3*(v*v)))+ &
    (40*nf*(-2*(59573-180*v+u*(-65909+6336*v-10602*zeta3)+ &
    &  10602*zeta3+6336*(u*u))+ &
    &  162*omu*zeta2*(-359+24*v+96*(-2+3*v)*(u*u)+144*(u*u*u)+ &
    &  48*(v*v)+24*u*(5-6*v+6*(v*v)))))/omu)/38880._dp+ &
    H2(1,0)*((10*(nf*nf))/9._dp+nf* &
    (-995._dp/54._dp-v+(8-12*v)*(u*u)-6*(u*u*u)+ &
    u*(-5+6*v-6*(v*v))-2*(v*v))+ &
    (-17724*zeta2+3*((-2182-20/omv+3256*v)*(u*u)+1628*(u*u*u)+ &
    &  4*u*(497-314*v+407*(v*v))+ &
    (48+247*w+8*v*(-6+41*w)+926*w*(v*v))/w))/108._dp)+ &
    G3(0,3,0)*((64*H1(0))/9._dp-(8*H1(1))/9._dp+ &
    (10*(omv*omv)*(-48+2*nf-11*v-4*(v*v))+ &
    u*(-27-2*(3+5*u)*v+(113+5*u)*(v*v)-80*(v*v*v)))/ &
    (9._dp*(omv*omv)))+(H3(1,0,1)* &
    (v*(83+40*v)*(omv*omv)+5*(u*u)*(9-16*v+8*(v*v))+ &
    u*(63-66*v-77*(v*v)+80*(v*v*v))))/(9._dp*(omv*omv))+ &
    G2(0,3)*((-176*zeta2)/9._dp-(128*H2(0,0))/9._dp+(8*H2(0,1))/9._dp+ &
    (8*H2(1,0))/9._dp+(H1(0)*(10*(omv*omv)*(48-2*nf+11*v+4*(v*v))+ &
    u*(27+2*(3+5*u)*v-(113+5*u)*(v*v)+80*(v*v*v))))/ &
    (9._dp*(omv*omv)))+((33-2*nf)/6._dp-(17*H1(0))/2._dp-(9*H1(1))/2._dp)* &
    (G1(0)*G1(0)*G1(0))+(9*(G1(0)*G1(0)*G1(0)*G1(0)))/8._dp+ &
    ((4*(-12+nf))/27._dp+(17*H1(0))/18._dp-H1(1)/18._dp)*(G1(2)*G1(2)*G1(2))+ &
    (G1(2)*G1(2)*G1(2)*G1(2))/72._dp+((68*(-12+nf))/27._dp-(17*H1(1))/18._dp)* &
    (H1(0)*H1(0)*H1(0))+(289*(H1(0)*H1(0)*H1(0)*H1(0)))/72._dp+ &
    G1(3)*((-622*zeta2*H1(1))/9._dp+(64*H3(0,0,0))/3._dp+20*H3(0,0,1)+ &
    (100*H3(0,1,0))/9._dp+(100*H3(1,0,0))/9._dp+(70*H3(1,0,1))/9._dp+ &
    (10*H3(1,1,0))/9._dp+(8*(60*zeta3+ &
    &  6*zeta2*(-108+8*nf-7*u-7*v-10*u*v-5*(u*u)-5*(v*v))))/ &
    &  27._dp+(8*H2(0,0)*(-90+4*nf+7*v+u*(7+10*v)+5*(u*u)+5*(v*v)))/ &
    &  9._dp+(H2(0,1)*(12+83*v+u*(83+80*v)+40*(u*u)+40*(v*v)))/9._dp+ &
    (H2(1,0)*(12+83*v+u*(83+80*v)+40*(u*u)+40*(v*v)))/9._dp+ &
    (76._dp/3._dp-(8*nf)/9._dp-(8*H1(1))/9._dp)*(H1(0)*H1(0))+ &
    (68*(H1(0)*H1(0)*H1(0)))/9._dp+ &
    H1(0)*((-4*H1(1))/3._dp+(-8232*zeta2- &
    &  2*nf*(-232+27*v+108*(-2+3*v)*(u*u)+162*(u*u*u)+ &
    &  54*(v*v)+27*u*(5-6*v+6*(v*v)))+ &
    (3*(-10*u*(1+u)*w+ &
    omv*(-24+4*v*(6+w*(86-157*u+407*(u*u)))+ &
    w*(-2416+1104*u-1091*(u*u)+814*(u*u*u))+ &
    (463+814*u)*w*(v*v))))/(omv*w))/54._dp- &
    (4*(H1(1)*H1(1)))/9._dp))+ &
    G2(2,0)*((-147._dp/2._dp+6*nf)*H1(1)+18*H2(0,0)-(2*H2(0,1))/9._dp- &
    (82*H2(1,0))/9._dp+H1(0)*(-H1(1)+ &
    (-969+81*nf-110*v-10*u*(11+8*v)-40*(u*u)-40*(v*v))/9._dp)+ &
    (-13542*zeta2+12*(nf*nf)- &
    &  27*nf*(-21+4*v+16*(-2+3*v)*(u*u)+24*(u*u*u)+8*(v*v)+ &
    &  4*u*(5-6*v+6*(v*v)))+ &
    &  6*(-2130+254*v+(-1091-10/omv+1628*v)*(u*u)+814*(u*u*u)+ &
    &  463*(v*v)+2*u*(502-314*v+407*(v*v))))/108._dp+ &
    (17*(H1(0)*H1(0)))/2._dp-(H1(1)*H1(1))/2._dp)+ &
    G2(0,2)*((-147._dp/2._dp+6*nf)*H1(1)+(332*H2(0,0))/9._dp+(232*H2(0,1))/9._dp+ &
    (242*H2(1,0))/9._dp+H1(0)*(-H1(1)+ &
    (-471+45*nf-110*v-10*u*(11+8*v)-40*(u*u)-40*(v*v))/9._dp)+ &
    (-12462*zeta2+12*(nf*nf)- &
    &  27*nf*(-21+4*v+16*(-2+3*v)*(u*u)+24*(u*u*u)+8*(v*v)+ &
    &  4*u*(5-6*v+6*(v*v)))+ &
    &  6*(-2130+254*v+(-1091-10/omv+1628*v)*(u*u)+814*(u*u*u)+ &
    &  463*(v*v)+2*u*(502-314*v+407*(v*v))))/108._dp+ &
    (17*(H1(0)*H1(0)))/2._dp-(H1(1)*H1(1))/2._dp)+ &
    H1(0)*H1(0)*((-47._dp/3._dp+(14*nf)/9._dp)*H1(1)+(68*H2(0,0))/9._dp+(85*H2(1,0))/9._dp+ &
    (4996*nf+3*(-19687+340*u+340*v-2058*zeta2)-84*(nf*nf))/216._dp- &
    (5*(H1(1)*H1(1)))/12._dp)+G1(0)*G1(0)* &
    (-98+(27*nf)/4._dp+(5*u)/2._dp+(5*v)/2._dp-(105*zeta2)/4._dp+ &
    (9*G2(0,2))/2._dp-9*G2(1,2)+(9*G2(2,0))/2._dp-4*G2(3,0)+4*G1(3)*H1(0)+ &
    G1(2)*(21._dp/2._dp-nf-(17*H1(0))/2._dp+H1(1)/2._dp)+((9+nf)*H1(1))/2._dp+ &
    H1(0)*(-41._dp/3._dp+(14*nf)/9._dp+(17*H1(1))/2._dp)+G1(1)*(9*H1(0)+9*H1(1))+ &
    &  4*H2(0,0)+5*H2(1,0)-(5*(nf*nf))/72._dp-(G1(2)*G1(2))/4._dp+ &
    (67*(H1(0)*H1(0)))/4._dp-(H1(1)*H1(1))/4._dp)+ &
    G1(2)*G1(2)*((-188*nf+15*(181-4*u-4*v+42*zeta2))/216._dp-G2(0,2)/2._dp+ &
    G2(1,2)-G2(2,0)/2._dp+(4*G2(3,0))/9._dp-(4*G1(3)*H1(0))/9._dp+ &
    H1(0)*(7._dp/3._dp-(4*nf)/9._dp-(11*H1(1))/6._dp)-(4*(-12+nf)*H1(1))/9._dp- &
    (4*H2(0,0))/9._dp-(5*H2(1,0))/9._dp-(29*(H1(0)*H1(0)))/12._dp+(H1(1)*H1(1))/12._dp) &
    +G2(3,0)*(H1(0)*(128._dp/3._dp-(20*nf)/9._dp+3*u+3*v+(8*H1(1))/9._dp)- &
    (64*H2(0,0))/3._dp-(172*H2(0,1))/9._dp-(92*H2(1,0))/9._dp+ &
    (H1(1)*(-(v*(83+40*v))-u*(83+80*v)-40*(u*u)))/9._dp+ &
    (30*u*(1+u)*w+omv*(72+ &
    &  6*v*(-12+w*(-172+314*u-814*(u*u)+ &
    &  9*nf*(1-6*u+12*(u*u))))+ &
    w*(7248-3312*u+8232*zeta2+3273*(u*u)-2442*(u*u*u)+ &
    nf*(-464+270*u-432*(u*u)+324*(u*u*u)))+ &
    &  3*(-463-814*u+36*nf*(1+3*u))*w*(v*v)))/(54._dp*omv*w)- &
    (68*(H1(0)*H1(0)))/9._dp+(4*(H1(1)*H1(1)))/9._dp)+ &
    ((-188*nf+15*(181-4*u-4*v+42*zeta2))/216._dp-(4*H2(0,0))/9._dp- &
    (5*H2(1,0))/9._dp)*(H1(1)*H1(1))+ &
    G2(1,2)*((147-12*nf)*H1(1)-36*H2(0,0)-72*H2(1,0)+ &
    H1(0)*(135-12*nf+7*u+2*H1(1)+5*(u*u))+ &
    (8*u*(-67+10*v)*(omv*omv)+(-407+54*nf)*(omv*omv*omv)+ &
    omv*(u*u)*(5726-484*v+2898*zeta2+2*(-407+54*nf)*(u*u*u)- &
    &  463*(v*v))+u*u*(-6*nf*(-1+v)* &
    (-58+3*v+12*(-2+3*v)*(u*u)+6*(v*v)+ &
    &  3*u*(5-6*v+6*(v*v)))+ &
    u*(u*(1101-2719*v+1628*(v*v))+ &
    &  2*(-537+846*v-721*(v*v)+407*(v*v*v)))))/(18._dp*omv*(u*u)) &
    -17*(H1(0)*H1(0))+H1(1)*H1(1))+ &
    G1(1)*(6*(-33+2*nf)*H2(0,0)+12*(-12+nf)*H2(1,1)+36*H3(0,0,0)+ &
    &  36*H3(0,0,1)+36*H3(0,1,0)+36*H3(1,0,0)+72*H3(1,1,0)+ &
    H2(1,0)*(-192+14*nf-7*u-5*(u*u))+ &
    H2(0,1)*(-150+10*nf+7*u+5*(u*u))+ &
    (-216*zeta3+6*zeta2*(-177+10*nf+7*u+5*(u*u)))/3._dp+ &
    (H1(1)*(8*u*(67-10*v)*(omv*omv)+(407-54*nf)*(omv*omv*omv)+ &
    omv*(u*u)*(-5726+484*v-2898*zeta2+(814-108*nf)*(u*u*u)+ &
    &  463*(v*v))+u*u*(6*nf*(-1+v)* &
    (-58+3*v+12*(-2+3*v)*(u*u)+6*(v*v)+ &
    &  3*u*(5-6*v+6*(v*v)))+ &
    u*(u*(-1101+2719*v-1628*(v*v))+ &
    &  2*(537-846*v+721*(v*v)-407*(v*v*v))))))/ &
    (18._dp*omv*(u*u))+(-H1(0)-H1(1))*(G1(2)*G1(2))+ &
    (57-2*nf+15*H1(1))*(H1(0)*H1(0))+17*(H1(0)*H1(0)*H1(0))+ &
    H1(0)*((54-2*nf)*H1(1)+(8*u*(67-10*v)*(omv*omv)+ &
    (407-54*nf)*(omv*omv*omv)+ &
    omv*(u*u)*(-5726+484*v-1602*zeta2+(814-108*nf)*(u*u*u)+ &
    &  463*(v*v))+u*u* &
    (6*nf*(-1+v)*(-58+3*v+12*(-2+3*v)*(u*u)+6*(v*v)+ &
    &  3*u*(5-6*v+6*(v*v)))+ &
    u*(u*(-1101+2719*v-1628*(v*v))+ &
    &  2*(537-846*v+721*(v*v)-407*(v*v*v)))))/ &
    (18._dp*omv*(u*u))-3*(H1(1)*H1(1)))-3*(H1(1)*H1(1))+ &
    G1(2)*(3*H1(1)+H1(0)*(3+4*H1(1))+2*(H1(0)*H1(0))+2*(H1(1)*H1(1)))- &
    H1(1)*H1(1)*H1(1))+G1(2)*(G2(1,2)*(-3-2*H1(0)-2*H1(1))+ &
    G2(3,0)*(-4._dp/3._dp-(8*H1(0))/9._dp-(8*H1(1))/9._dp)+ &
    G2(0,2)*(3._dp/2._dp+H1(0)+H1(1))+G2(2,0)*(3._dp/2._dp+H1(0)+H1(1))+ &
    (316._dp/3._dp-8*nf)*H2(1,0)-(404*H3(0,0,0))/9._dp-(160*H3(0,0,1))/9._dp- &
    (160*H3(0,1,0))/9._dp+(160*H3(1,0,1))/9._dp-(568*H3(1,1,0))/9._dp+ &
    H1(1)*((8*H2(0,0))/9._dp+(10*H2(1,0))/9._dp+ &
    (60*(-250+omw)+2384*nf-3534*zeta2-96*(nf*nf))/108._dp)- &
    (40*(nf*nf))/27._dp+H2(0,0)* &
    (235-16*nf+u*(3+20/(9-9*v))+3*v-(5*(u*u))/(9._dp*(omv*omv))) &
    +(H2(0,1)*(-135+u*(-83+20/omv-80*v)-83*v+ &
    (-40-5/(omv*omv))*(u*u)-40*(v*v)))/9._dp+ &
    nf*(5086._dp/81._dp-6*u-3*v-(2*v)/omu+12*u*v+(91*zeta2)/9._dp+ &
    (3*(omv*omv))/u+6*(u*u)+6*(v*v)-(38*(v*v))/(3._dp*(omu*omu))) &
    +(-49396+u*(6558-9768*v)+4116*v+(6480*v)/omu+3084*zeta3- &
    (2442*(omv*omv))/u-4884*(u*u)+ &
    &  6*zeta2*(-819-440*v-(40*u*(1+omv*(11+8*v)))/omv+ &
    &  10*(-16+1/(omv*omv))*(u*u)-160*(v*v))-4884*(v*v)+ &
    (5184*(v*v))/(omu*omu))/108._dp+ &
    G1(3)*(H1(0)*(4._dp/3._dp+(8*H1(1))/9._dp)+(8*(H1(0)*H1(0)))/9._dp)+ &
    ((4*(-93+nf))/9._dp+(17*H1(1))/6._dp)*(H1(0)*H1(0))- &
    (289*(H1(0)*H1(0)*H1(0)))/18._dp+ &
    H1(0)*(((-69+8*nf)*H1(1))/9._dp+(8*H2(0,0))/9._dp+(10*H2(1,0))/9._dp+ &
    (48*u*(-67+10*v)*(omv*omv)+6*(-407+54*nf)*(omv*omv*omv)+ &
    omv*(u*u)*(10830*zeta2-96*(nf*nf)+648*nf*(u*u*u)- &
    &  6*(-3226+394*v+814*(u*u*u)+463*(v*v)))- &
    &  2*(u*u)*(2*nf*(-1+v)* &
    (74+27*v+108*(-2+3*v)*(u*u)+54*(v*v)+ &
    &  27*u*(5-6*v+6*(v*v)))- &
    &  3*u*(u*(1101-2719*v+1628*(v*v))+ &
    &  2*(-532+846*v-721*(v*v)+407*(v*v*v)))))/ &
    (108._dp*omv*(u*u))+(5*(H1(1)*H1(1)))/6._dp)+ &
    (4*(-12+nf)*(H1(1)*H1(1)))/9._dp-(H1(1)*H1(1)*H1(1))/18._dp)+ &
    H1(0)*((76._dp/3._dp-(8*nf)/9._dp)*H2(0,0)-(5*(-57+2*nf)*H2(1,0))/9._dp+ &
    H1(1)*((-8*H2(0,0))/9._dp-(10*H2(1,0))/9._dp+ &
    (766*nf-15*(287+4*u+4*v)+6462*zeta2-24*(nf*nf))/108._dp)+ &
    (480*(nf*nf)+2*nf*(-10010+u*(360-972*v)+990*v+(324*v)/omu- &
    (144*v)/omw+(486*omv*v)/u-2502*zeta2-972*(v*v)+ &
    (2052*(v*v))/(omu*omu))+ &
    &  3*(38860-7314*v-(6480*v)/omu+(1776*v)/omw- &
    (2442*omv*v)/u+12*u*(-270+407*v)+9516*zeta3+4884*(v*v)- &
    (5184*(v*v))/(omu*omu)+ &
    &  6*zeta2*(4047+440*v+32*u*(-2-5/omv+10*v)+ &
    &  40*(-5+1/(omv*omv))*(u*u)+160*(v*v))))/324._dp- &
    (4*(-12+nf)*(H1(1)*H1(1)))/9._dp+(H1(1)*H1(1)*H1(1))/18._dp)+ &
    G1(0)*(G2(0,2)*((-30+nf)/2._dp-9*H1(0))+G2(2,0)*((-30+nf)/2._dp-9*H1(0))+ &
    G2(3,0)*((-4*(-30+nf))/9._dp+8*H1(0))+G2(1,2)*(30-nf+18*H1(0))+ &
    (72-6*nf)*H2(1,1)-(64*H3(0,0,0))/3._dp-(172*H3(0,0,1))/9._dp- &
    (92*H3(0,1,0))/9._dp-(332*H3(1,0,0))/9._dp-(20*H3(1,0,1))/3._dp+ &
    (10*H3(1,1,0))/9._dp+H2(0,1)* &
    (487._dp/6._dp-5*nf+(-7-20/(9._dp*omv))*u+ &
    (-5+5/(9._dp*(omv*omv)))*(u*u))+ &
    (H2(1,0)*(1857-162*nf+220*v+20*u*(11+8*v)+80*(u*u)+ &
    &  80*(v*v)))/18._dp+(H2(0,0)* &
    (omv*omv*(1152-64*nf+110*v+40*(v*v))+ &
    u*(27+2*(3+5*u)*v-(113+5*u)*(v*v)+80*(v*v*v))))/ &
    (9._dp*(omv*omv))+(-(omv*omv* &
    (3*(883-1628*v+8*nf*(-11+27*v))*(u*u)+ &
    &  6*(-407+54*nf)*(u*u*u)- &
    &  3*u*(1305-680*v-322*zeta3+814*(v*v))+ &
    v*(-3027-609*v+12*nf*(5+9*v)+966*zeta3+ &
    &  12*zeta2*(-123+220*v+80*(v*v)))))+ &
    &  6*zeta2*(u*u)*(-294+213*v+466*(v*v)- &
    &  5*u*(15-28*v+14*(v*v))-390*(v*v*v))+ &
    &  2*omv*u*(6*nf*(-9+22*v-40*(v*v)+27*(v*v*v))+ &
    &  6*zeta2*(123-490*v+137*(v*v)+240*(v*v*v))))/ &
    (54._dp*omw*(omv*omv))+ &
    (H1(1)*(omv*(13422*zeta2-4884*(u*u*u)+ &
    &  36*nf*(-29+3*v+12*(-2+3*v)*(u*u)+18*(u*u*u)+6*(v*v)+ &
    &  3*u*(5-6*v+6*(v*v))))+ &
    &  6*(2715-2969*v-209*(v*v)+u*u*(1101-2719*v+1628*(v*v))+ &
    &  463*(v*v*v)+2*u*(-502+816*v-721*(v*v)+407*(v*v*v)))))/ &
    (108._dp*omv)+((17*H1(0))/18._dp+H1(1)/2._dp)*(G1(2)*G1(2))+ &
    G1(1)*(H1(0)*(-30+nf-18*H1(1))+(-30+nf)*H1(1)-18*(H1(0)*H1(0)))+ &
    G1(3)*((4*(-30+nf)*H1(0))/9._dp-8*(H1(0)*H1(0)))+ &
    (73._dp/6._dp-(19*nf)/9._dp-(119*H1(1))/18._dp)*(H1(0)*H1(0))- &
    (289*(H1(0)*H1(0)*H1(0)))/18._dp+ &
    G1(2)*(-((-30+nf)*(-39+4*nf))/36._dp+ &
    H1(0)*(23._dp/3._dp+nf-(26*H1(1))/9._dp)-(3*H1(1))/2._dp+ &
    (145*(H1(0)*H1(0)))/9._dp-H1(1)*H1(1))+(3*(H1(1)*H1(1)))/2._dp+ &
    H1(0)*((-37._dp/6._dp-nf)*H1(1)-8*H2(0,0)-10*H2(1,0)+ &
    (-6*(omv*omv)*(-4528+514*v+463*(v*v))+ &
    omv*(16302*w*zeta2+4884*(u*u*u*u)- &
    &  4*nf*(412-439*v+54*(-7+9*v)*(u*u*u)+162*(u*u*u*u)- &
    &  27*(v*v)+27*(u*u)*(13-26*v+18*(v*v))+54*(v*v*v)+ &
    u*(-547+324*v-270*(v*v)+162*(v*v*v))))- &
    &  6*u*(5688-8014*v+3305*(v*v)+ &
    u*u*(1915-4347*v+2442*(v*v))-1793*(v*v*v)+ &
    u*(-2285+5632*v-5789*(v*v)+2442*(v*v*v))+814*(v*v*v*v)) &
    )/(108._dp*omv*w)+(35*(H1(1)*H1(1)))/18._dp)+(H1(1)*H1(1)*H1(1))/2._dp)- &
    (4*(-12+nf)*(H1(1)*H1(1)*H1(1)))/27._dp+(H1(1)*H1(1)*H1(1)*H1(1))/72._dp
    return
    end function Hgamma_4a2re
          
          
    function Hgamma_4a2im(u,v, H1,H2,H3,H4,G1,G2,G3,G4,nf) &
    result(fxn)
    use constants
    implicit  none
    include 'src/Inc/hpls.f'
    real(dp)::fxn
    real(dp)::u,v,w,omu,omv,omw,nf
    w=1-u-v
    omu=1-u
    omv=1-v
    omw=1-w
    fxn=(-16*pi*G3(0,0,1))/9._dp+18*pi*G3(0,0,2)+(16*pi*G3(0,0,3))/9._dp+ &
    &  26*pi*G3(0,1,0)-10*pi*G3(0,1,2)+18*pi*G3(0,2,0)-(92*pi*G3(0,2,1))/9._dp+ &
    &  18*pi*G3(0,2,2)-(152*pi*G3(0,2,3))/9._dp+(8*pi*G3(0,3,0))/9._dp+ &
    (8*pi*G3(0,3,2))/9._dp+(128*pi*G3(0,3,3))/9._dp-36*pi*G3(1,0,2)- &
    &  72*pi*G3(1,1,1)+72*pi*G3(1,1,2)-36*pi*G3(1,2,0)+72*pi*G3(1,2,1)- &
    &  36*pi*G3(1,2,2)+18*pi*G3(2,0,0)-18*pi*G3(2,0,1)+18*pi*G3(2,0,2)- &
    &  18*pi*G3(2,0,3)-18*pi*G3(2,1,0)+72*pi*G3(2,1,1)-18*pi*G3(2,1,2)+ &
    &  18*pi*G3(2,2,0)-36*pi*G3(2,2,1)+36*pi*G3(2,2,3)-(80*pi*G3(2,3,0))/9._dp+ &
    (164*pi*G3(2,3,2))/9._dp+(80*pi*G3(2,3,3))/9._dp-(164*pi*G3(3,0,0))/9._dp- &
    (70*pi*G3(3,0,1))/9._dp-(164*pi*G3(3,0,2))/9._dp+(16*pi*G3(3,0,3))/3._dp- &
    (164*pi*G3(3,2,0))/9._dp-(70*pi*G3(3,2,1))/9._dp-(80*pi*G3(3,2,3))/9._dp+ &
    (128*pi*G3(3,3,0))/9._dp-(128*pi*G3(3,3,3))/9._dp+ &
    G2(2,2)*(-6*(-12+nf)*pi+36*pi*H1(0))+ &
    G2(3,0)*(-(pi*(-16+4*nf+9*u+9*v))/3._dp+(76*pi*H1(0))/3._dp+ &
    (164*pi*H1(1))/9._dp)-6*(-12+nf)*pi*H2(1,1)-(124*pi*H3(0,0,0))/3._dp- &
    (172*pi*H3(0,0,1))/9._dp-18*pi*H3(0,1,0)-(104*pi*H3(1,0,0))/3._dp- &
    (20*pi*H3(1,0,1))/3._dp+(80*pi*H3(1,1,0))/3._dp+ &
    G2(1,0)*(18*pi*H1(0)+36*pi*H1(1)+pi*(-18+nf-7*u-5*(u*u)))+ &
    G2(1,2)*(36*pi*H1(1)+pi*(-27+4*nf-7*u-5*(u*u)))+ &
    (pi*H2(0,1)*(1461-90*nf-(2*(20+63*omv)*u)/omv+ &
    &  10*(-9+1/(omv*omv))*(u*u)))/18._dp+ &
    G2(1,1)*(-72*pi*H1(1)+2*pi*(21-2*nf+7*u+5*(u*u)))+ &
    G2(2,3)*((-2*pi*H1(0))/9._dp-(164*pi*H1(1))/9._dp+ &
    (pi*(-321+36*nf-27*u-(20*u)/omv-27*v+(5*(u*u))/(omv*omv)))/ &
    &  9._dp)+G2(0,1)*(-8*pi*H1(0)+10*pi*H1(1)+ &
    (pi*(3+9*nf+u*(-173-20/omv-80*v)-110*v+ &
    (-85+5/(omv*omv))*(u*u)-40*(v*v)))/9._dp)- &
    (pi*H2(1,0)*((-807+82*nf-54*v)*(omv*omv)+ &
    &  2*u*(-47+5*u+74*v-27*(v*v))))/(18._dp*(omv*omv))+ &
    G2(3,3)*((-128*pi*H1(0))/9._dp+(8*pi* &
    (-9+2*nf-7*v-u*(7+10*v)-5*(u*u)-5*(v*v)))/9._dp)+ &
    G2(0,0)*((-164*pi*H1(0))/9._dp-18*pi*H1(1)+ &
    (pi*(417-27*nf+110*v+10*u*(11+8*v)+40*(u*u)+40*(v*v)))/9._dp) &
    +G2(3,2)*(18*pi*H1(0)+(pi*(12+83*v+u*(83+80*v)+40*(u*u)+ &
    &  40*(v*v)))/9._dp)-(pi*H2(0,0)* &
    ((32*nf-9*(80+3*v))*(omv*omv)+ &
    u*(56-92*v+36*(v*v)+5*u*(8-18*v+9*(v*v)))))/(9._dp*(omv*omv)) &
    +G2(0,2)*((-241*pi*H1(0))/9._dp-18*pi*H1(1)+ &
    (pi*(-81*nf+5*(183+44*v+4*u*(11+8*v)+16*(u*u)+16*(v*v))))/ &
    &  18._dp)+G2(2,0)*((-161*pi*H1(0))/9._dp-18*pi*H1(1)+ &
    (pi*(-81*nf+5*(183+44*v+4*u*(11+8*v)+16*(u*u)+16*(v*v))))/ &
    &  18._dp)+G2(0,3)*((56*pi*H1(0))/9._dp-(8*pi*H1(1))/9._dp+ &
    (pi*(10*(omv*omv)*(-48+2*nf-11*v-4*(v*v))+ &
    u*(-27-2*(3+5*u)*v+(113+5*u)*(v*v)-80*(v*v*v))))/ &
    (9._dp*(omv*omv)))+G2(2,1)* &
    (-18*pi*H1(0)+18*pi*H1(1)+ &
    (pi*(omv*omv*(-1431+108*nf-83*v-40*(v*v))+ &
    u*(-63+66*v+77*(v*v)-5*u*(9-16*v+8*(v*v))-80*(v*v*v)))) &
    /(9._dp*(omv*omv)))+(pi*H1(1)* &
    (5430-5938*v-418*(v*v)+u*u*(2202-5438*v+3256*(v*v))+ &
    omv*(3450*zeta2-1628*(u*u*u)+ &
    &  12*nf*(-29+3*v+12*(-2+3*v)*(u*u)+18*(u*u*u)+6*(v*v)+ &
    &  3*u*(5-6*v+6*(v*v))))+926*(v*v*v)+ &
    &  4*u*(-502+816*v-721*(v*v)+407*(v*v*v))))/(36._dp*omv)- &
    (pi*(-3*(-407+54*nf)*v*(omv*omv*omv*omv)+ &
    u*(omv*omv)*(1221-2643*v-6156*v*zeta2+2334*v*zeta3+ &
    &  1818*(v*v)+996*zeta2*(v*v)+ &
    &  3*(u*u*u)*(-209-1647*v+778*zeta3+1628*(v*v))- &
    &  2592*(v*v*v)+480*zeta2*(v*v*v)+ &
    &  3*(u*u)*(-165+1865*v-1556*zeta3+778*v*zeta3-3198*(v*v)+ &
    &  814*(v*v*v))-3*u* &
    (363+558*v-778*zeta3+1556*v*zeta3-914*(v*v)+ &
    &  1221*(v*v*v))-2*nf* &
    (81-327*v-204*v*zeta2+6*(5+27*v)*(u*u*u*u)- &
    u*u*u*(87+369*v+204*zeta2-324*(v*v))+276*(v*v)- &
    &  342*(v*v*v)+u*u* &
    (-204*(-2+v)*zeta2+3*(55+47*v-214*(v*v)+54*(v*v*v))) &
    +u*(204*(-1+2*v)*zeta2- &
    &  3*(63-104*v+46*(v*v)+81*(v*v*v)))))+ &
    &  12*omv*zeta2*(u*u)*(-513+1622*v-1175*(v*v)-34*(v*v*v)+ &
    &  80*(v*v*v*v))+2*(u*u*u)* &
    (33*(15+37*v)*(u*u)*((-1+v)*(-1+v))+ &
    &  6*zeta2*(1026-2671*v-5*(-2+v)*v*(u*u*u)+2137*(v*v)+ &
    &  5*v*(u*u)*(8-24*v+15*(v*v))-290*(v*v*v)-237*(v*v*v*v)+ &
    u*(-513+999*v-324*(v*v)-287*(v*v*v)+120*(v*v*v*v))+ &
    &  40*(v*v*v*v*v)))))/(54._dp*omw*u*(omu*omu)*(omv*omv))+ &
    (((-30+nf)*pi)/2._dp-9*pi*G1(1)+(9*pi*G1(2))/2._dp-4*pi*G1(3)- &
    (17*pi*H1(0))/2._dp-(9*pi*H1(1))/2._dp)*(G1(0)*G1(0))+ &
    (9*pi*(G1(0)*G1(0)*G1(0)))/2._dp+ &
    ((3*pi)/2._dp+(4*pi*G1(3))/9._dp+(17*pi*H1(0))/18._dp+(3*pi*H1(1))/2._dp)* &
    (G1(2)*G1(2))-(pi*(G1(2)*G1(2)*G1(2)))/2._dp+ &
    (((-1131+70*nf)*pi)/18._dp-(155*pi*H1(1))/18._dp)*(H1(0)*H1(0))+ &
    (17*pi*(H1(0)*H1(0)*H1(0)))/18._dp+(3*pi*(H1(1)*H1(1)))/2._dp+ &
    G1(2)*(3*(-25+2*nf)*pi*H1(1)+ &
    G1(3)*((-4*pi)/3._dp-(8*pi*H1(0))/9._dp-(8*pi*H1(1))/9._dp)+ &
    (322*pi*H2(0,0))/9._dp-(2*pi*H2(0,1))/9._dp-(242*pi*H2(1,0))/9._dp+ &
    H1(0)*((-17*pi*H1(1))/9._dp-(pi* &
    (9*(-62+2*nf+3*v)*(omv*omv)+ &
    u*(47-5*u-74*v+27*(v*v))))/(9._dp*(omv*omv)))- &
    (pi*(5430-5938*v-418*(v*v)+u*u*(2202-5438*v+3256*(v*v))+ &
    omv*(3210*zeta2-1628*(u*u*u)+ &
    &  12*nf*(-29+3*v+12*(-2+3*v)*(u*u)+18*(u*u*u)+6*(v*v)+ &
    &  3*u*(5-6*v+6*(v*v))))+926*(v*v*v)+ &
    &  4*u*(-502+816*v-721*(v*v)+407*(v*v*v))))/(36._dp*omv)- &
    (169*pi*(H1(0)*H1(0)))/18._dp-(3*pi*(H1(1)*H1(1)))/2._dp)+ &
    G1(0)*(G1(3)*((-4*(-30+nf)*pi)/9._dp+8*pi*H1(0))+ &
    G1(1)*(-((-30+nf)*pi)+18*pi*H1(0))+ &
    G1(2)*(((-27+nf)*pi)/2._dp-8*pi*H1(0)+pi*H1(1))+(100*pi*H2(0,0))/3._dp+ &
    (164*pi*H2(0,1))/9._dp+(304*pi*H2(1,0))/9._dp+18*pi*H2(1,1)+ &
    (2*pi*H1(1)*(-168+18*nf-55*v-5*u*(11+8*v)-20*(u*u)- &
    &  20*(v*v)))/9._dp+H1(0)* &
    (-(pi*H1(1))+(pi*(2*(omv*omv)*(-159+10*nf-55*v-20*(v*v))+ &
    u*(-27-2*(3+5*u)*v+(113+5*u)*(v*v)-80*(v*v*v))))/ &
    (9._dp*(omv*omv)))-(pi* &
    (5070*omv*zeta2+36*omv*nf* &
    (-29+3*v+12*(-2+3*v)*(u*u)+18*(u*u*u)+6*(v*v)+ &
    &  3*u*(5-6*v+6*(v*v)))+ &
    &  6*(2795-3129*v-814*omv*(u*u*u)-129*(v*v)+ &
    u*u*(1101-2719*v+1628*(v*v))+463*(v*v*v)+ &
    &  2*u*(-542+856*v-721*(v*v)+407*(v*v*v)))))/(108._dp*omv)- &
    (pi*(G1(2)*G1(2)))/2._dp-(pi*(H1(0)*H1(0)))/2._dp-(pi*(H1(1)*H1(1)))/2._dp)+ &
    G1(3)*((-124*pi*H2(0,0))/3._dp-(172*pi*H2(0,1))/9._dp-18*pi*H2(1,0)- &
    (pi*H1(1)*(v*(83+40*v)+u*(83+80*v)+40*(u*u)))/9._dp+ &
    H1(0)*((8*pi*H1(1))/9._dp+(8*pi* &
    (-21+2*nf-7*v-u*(7+10*v)-5*(u*u)-5*(v*v)))/9._dp)+ &
    (pi*(30*u*(1+u)*w+omv* &
    (72+6*v*(-12+w*(-172+314*u-814*(u*u)+ &
    &  9*nf*(1-6*u+12*(u*u))))+ &
    w*(7248-3312*u+5280*zeta2+3273*(u*u)-2442*(u*u*u)+ &
    nf*(-464+270*u-432*(u*u)+324*(u*u*u)))+ &
    &  3*(-463-814*u+36*nf*(1+3*u))*w*(v*v))))/(54._dp*omv*w)+ &
    (4*pi*(H1(0)*H1(0)))/9._dp+(4*pi*(H1(1)*H1(1)))/9._dp)+ &
    H1(0)*((-289._dp/6._dp+3*nf)*pi*H1(1)+8*pi*H2(0,0)+10*pi*H2(1,0)- &
    (pi*(6*(omv*omv)*(-352+254*v+463*(v*v))+ &
    omv*(-4890*w*zeta2-4884*(u*u*u*u)+ &
    &  4*nf*(52-79*v+54*(-7+9*v)*(u*u*u)+162*(u*u*u*u)- &
    &  27*(v*v)+27*(u*u)*(13-26*v+18*(v*v))+54*(v*v*v)+ &
    u*(-187+324*v-270*(v*v)+162*(v*v*v))))+ &
    &  6*u*(1332-3228*v+2875*(v*v)+ &
    u*u*(1915-4347*v+2442*(v*v))-1793*(v*v*v)+ &
    u*(-2105+5462*v-5789*(v*v)+2442*(v*v*v))+814*(v*v*v*v))) &
    )/(108._dp*omv*w)+(17*pi*(H1(1)*H1(1)))/18._dp)+ &
    G1(1)*(G1(2)*(-3*pi-2*pi*H1(0)-2*pi*H1(1))+H1(0)*(-24*pi+20*pi*H1(1))- &
    &  36*pi*H2(0,0)-18*pi*H2(0,1)-18*pi*H2(1,0)-36*pi*H2(1,1)+ &
    pi*H1(1)*(30-4*nf+7*u+5*(u*u))+ &
    (pi*(8*u*(-67+10*v)*(omv*omv)+(-407+54*nf)*(omv*omv*omv)+ &
    omv*(u*u)*(5726-484*v+2250*zeta2+2*(-407+54*nf)*(u*u*u)- &
    &  463*(v*v))+u*u*(-6*nf*(-1+v)* &
    (-58+3*v+12*(-2+3*v)*(u*u)+6*(v*v)+ &
    &  3*u*(5-6*v+6*(v*v)))+ &
    u*(u*(1101-2719*v+1628*(v*v))+ &
    &  2*(-537+846*v-721*(v*v)+407*(v*v*v))))))/ &
    (18._dp*omv*(u*u))+pi*(G1(2)*G1(2))+pi*(H1(0)*H1(0))+pi*(H1(1)*H1(1)))+ &
    (pi*(H1(1)*H1(1)*H1(1)))/2._dp
    return
    end function Hgamma_4a2im

end module
