!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine coupling
          use topwidth
      implicit none
      include 'types.f'
c--- initialize electroweak couplings and calculate alpha-s; this
c--- must be called at the beginning of "chooser" to enable it to
c--- set up all the variables (e.g. twidth). Once nflav is set,
c--- alpha-s will be determined again (in coupling2).

      include 'cplx.h'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'couple.f'
      include 'zcouple_cms.f'
      include 'scale.f'
      include 'nlooprun.f'
      include 'ewinput.f'
      include 'nflav.f'
      include 'b0.f'
      include 'kpart.f'
      include 'mpicommon.f'
      include 'zerowidth.f'
      include 'nproc.f'
      include 'blha.f'
      integer:: i
      real(dp):: aemmz,alphas,cmass,bmass
      real(dp):: wwidth_cms,zwidth_cms
      complex(dp):: zwmass2,zzmass2
      character(len=3):: inlabel(10)
      character(len=5):: tworder
      common/qmass/cmass,bmass
      common/em/aemmz
      logical, save :: first = .true.

c--- blank out labels that indicate input parameters
      do i=1,10
        inlabel(i)='   '
      enddo
      inlabel(3)='(+)'
      inlabel(4)='(+)'
      inlabel(8)='(+)'

      if (ewscheme == -1) then
c--- This is the MCFM default, corresponding to an effective
c--- field theory approach valid for scales below the top-mass
c--- (see Georgi, Nucl. Phys. B 363 (1991) 301).
c--- There are 4 inputs here instead of the usual 3 ...
         Gf = Gf_inp
         aemmz  = aemmz_inp
         wmass  = wmass_inp
         zmass  = zmass_inp
         inlabel(5)='(+)'
         inlabel(6)='(+)'
         inlabel(1)='(+)'
         inlabel(2)='(+)'
         inlabel(8)='   '
c--- ... and as result, both xw and mtop are derived
c--- (using wmass=zmass*sqrt(rho)*cos(theta_w)
c---    and rho=1+3._dp*aemmz/16._dp/pi/xw*(mt/wmass)**2 )
         xw  = fourpi*aemmz/(8._dp*wmass**2*Gf/rt2)
         mt  = sqrt(16._dp*pisq/3._dp/rt2/Gf*(
     &          wmass**2/zmass**2/(1._dp-xw)-1._dp))

      elseif (ewscheme == 0) then
c------------------------------------------------------------
c     option=0 : MadEvent default (= AlpGen with iewopt=2)
c------------------------------------------------------------

c-- equal to the input values
         xw  = xw_inp
         aemmz  = aemmz_inp
         zmass  = zmass_inp
         inlabel(7)='(+)'
         inlabel(6)='(+)'
         inlabel(1)='(+)'
c-- derived
         wmass  = zmass * sqrt( One - xw )
         Gf = aemmz * Fourpi/xw/(8._dp*wmass**2/Rt2)

      elseif (ewscheme == 1) then
c-----------------------------------------------------
c     option=1 : LUSIFER and AlpGen (iewopt=3) default
c-----------------------------------------------------

c-- equal to the input values
         zmass  = zmass_inp
         wmass  = wmass_inp
         Gf = Gf_inp
         inlabel(1)='(+)'
         inlabel(2)='(+)'
         inlabel(5)='(+)'
c-- derived
         xw  = One-(wmass/zmass)**2
         aemmz  = Rt2*Gf*wmass**2*xw/pi

      elseif (ewscheme == 2) then
c-------------------------------------------------------------------
c     option=2 : W and Z mass are derived from couplings
c-------------------------------------------------------------------

c-- equal to the input values
         Gf = Gf_inp
         aemmz  = aemmz_inp
         xw  = xw_inp
         inlabel(5)='(+)'
         inlabel(6)='(+)'
         inlabel(7)='(+)'
c-- derived
         wmass  = sqrt(aemmz*pi/xw/Gf/Rt2)
         zmass  = wmass/sqrt(One-xw)

      elseif (ewscheme == 3) then
c-----------------------------------------------------------------
c     option=3 : USER choice : you should know what you're doing!!
c-----------------------------------------------------------------
         Gf = Gf_inp
         aemmz  = aemmz_inp
         xw  = xw_inp
         wmass  = wmass_inp
         zmass  = zmass_inp
         inlabel(5)='(+)'
         inlabel(6)='(+)'
         inlabel(7)='(+)'
         inlabel(1)='(+)'
         inlabel(2)='(+)'

      elseif (ewscheme == 4) then
c-----------------------------------------------------
c     option=4 : Complex mass scheme with GF,MW,MZ +widths as inputs
c-----------------------------------------------------

c-- equal to the input values
         zmass  = zmass_inp
         wmass  = wmass_inp
         Gf = Gf_inp
         inlabel(1)='(+)'
         inlabel(2)='(+)'
         inlabel(5)='(+)'
c-- derived
         if (zerowidth) then
                 wwidth_cms = 0d0
                 zwidth_cms = 0d0
         else
                 wwidth_cms = wwidth
                 zwidth_cms = zwidth
         endif

       zwmass2=cplx2(wmass**2,-wmass*wwidth_cms)
       zzmass2=cplx2(zmass**2,-zmass*zwidth_cms)
       zxw=cone-zwmass2/zzmass2
       zaemmz=cmplx(Rt2*Gf/pi,kind=dp)*wmass**2*(1._dp-(wmass/zmass)**2)

       zesq  = cplx1(fourpi)*zaemmz

       call couplz_cms(zxw)

c---Set up parameters for call to couplz and definition of esq,gwsq
      aemmz=real(zaemmz,kind=dp)   ! this is real anyway
      xw=abs(zxw)
      elseif (ewscheme == 5) then
c-----------------------------------------------------
c     option=5 : Complex mass scheme with GF,MW,MZ +widths as inputs
c-----------------------------------------------------

c-- equal to the input values
         zmass  = zmass_inp
         wmass  = wmass_inp
         Gf = Gf_inp
         inlabel(1)='(+)'
         inlabel(2)='(+)'
         inlabel(5)='(+)'
c-- derived
         zwmass2=cplx2(wmass**2,-wmass*wwidth)
         zzmass2=cplx2(zmass**2,-zmass*zwidth)
         zxw=cone-zwmass2/zzmass2
         zaemmz  = Rt2*Gf/pi*abs(zwmass2*(1-zwmass2/zzmass2))

         zesq  = fourpi*zaemmz

c---Set up parameters for call to couplz and definition of esq,gwsq
         aemmz=real(zaemmz,kind=dp)   ! this is real anyway
         xw=abs(zxw)
      else
         write(6,*) 'ewscheme=',ewscheme,' is not a valid input.'
         stop
      endif

c---- Output of all of all branches of if statement is aemmz and xw
c--- Now set up the other derived parameters
      gwsq=fourpi*aemmz/xw
      esq=fourpi*aemmz
      gw=sqrt(gwsq)
      call couplz(xw)

c set up complexified versions of real parameters if not in CMS
      if (ewscheme < 4) then
        zesq=cplx2(esq,0d0)
        zxw=cplx2(xw,0d0)
      endif
      call couplz_cms(zxw)

c Weak mixing angle in the complex mass scheme

c--- Calculate the appropriate Higgs vacuum expectation value.
c--- This vevsq is defined so that gwsq/(4*wmass**2)=Gf*rt2=1/vevsq
c--- (ie differs from definition in ESW)
      vevsq=1._dp/rt2/Gf

c--- set up the beta-function
      b0=(xn*11._dp-2._dp*nflav)/6._dp

c--- initialize the pdf set
      nlooprun=0
      if (useblha==0) then
         call pdfwrap
      endif

      cmass=sqrt(mcsq)
      bmass=sqrt(mbsq)
      musq=scale**2

c--- set the number of loops to use in the running of alpha_s
c--- if it hasn't been set by pdfwrap already
      if (nlooprun == 0) then
        if (kpart==klord) then
          nlooprun=1
        else
          nlooprun=2
        endif
      endif

c--- initialize alpha_s
      as=alphas(abs(scale),amz,nlooprun)

      ason2pi=as/twopi
      ason4pi=as/fourpi
      gsq=fourpi*as

c--- Set-up twidth, using LO formula everywhere
      if ((nproc == 1610) .or. (nproc == 1650)) then
          if (rank == 0) then
              write (*,*) "Top decay width fixed at LO with mb=0, wwidth=0"
              write (*,*) "Assembly routine will take into account higher order pieces."
          endif
          twidth=lotopdecaywidth(mt,0d0,wmass,0d0)
          tworder='(LO) '
          !twidth = 1.5415d0
          ! from 1210.2808
          ! LO 1.4806
          ! NLO 1.5109
          ! NNLO 1.5415
      else
          twidth=lotopdecaywidth(mt,mb,wmass,wwidth)
          tworder='(LO) '
          if (kpart /= klord) then
              twidth=twidth*(1d0 + nloratiotopdecay(mt,mb,wmass,wwidth,mt))
              tworder = '(NLO)'
          endif
      endif

c Corrections to top width
c      write(6,*) 'relative NLO', nloratiotopdecay(mt,mb,wmass,wwidth)
c      write(6,*) 'relative NNLO', nnlotopdecay(mt,wmass)
c      stop

c      if (verbose .and. rank == 0 .and. first ) then
      if (rank == 0 .and. first ) then
      write(6,*) '************** Electroweak parameters **************'
      write(6,*) '*                                                  *'
      write(6,75) 'zmass',inlabel(1),zmass,'wmass',inlabel(2),wmass
      write(6,75) 'zwidth',inlabel(3),zwidth,'wwidth',inlabel(4),wwidth
      write(6,76) 'Gf',inlabel(5),gf,'1/aemmz',inlabel(6),1._dp/aemmz
      write(6,75) 'xw',inlabel(7),xw,'mtop',inlabel(8),mt
      write(6,75) 'gwsq',inlabel(9),gwsq,'esq',inlabel(10),esq
      write(6,77) 'top width',twidth,'     at order  ',tworder
      write(6,78) 'mb',mb,'mc',mc
      write(6,*) '*                                                  *'
      write(6,*) '* Parameters marked (+) are input, others derived  *'
      write(6,*) '****************************************************'
      endif

   75 format(' * ',a6,a3,f13.7,3x,a7,a3,f12.7,'  *')
   76 format(' * ',a6,a3,d13.6,3x,a7,a3,f12.7,'  *')
   77 format(' * ',a9,f13.7,1x,a15,a5,'      *')
   78 format(' * ',a5,4x,f13.7,6x,a4,2x,f13.7,'  *')

      first = .false.

      return
      end


      block data wsalam1
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'ewcharge.f'
      data Q(-5)/+0.333333333333333_dp/
      data Q(-4)/-0.666666666666667_dp/
      data Q(-3)/+0.333333333333333_dp/
      data Q(-2)/-0.666666666666667_dp/
      data Q(-1)/+0.333333333333333_dp/
      data Q(0)/+0._dp/
      data Q(+1)/-0.333333333333333_dp/
      data Q(+2)/+0.666666666666667_dp/
      data Q(+3)/-0.333333333333333_dp/
      data Q(+4)/+0.666666666666667_dp/
      data Q(+5)/-0.333333333333333_dp/
      data tau/1._dp,-1._dp,1._dp,-1._dp,1._dp,0._dp,-1._dp,1._dp,-1._dp,1._dp,-1._dp/
      end



