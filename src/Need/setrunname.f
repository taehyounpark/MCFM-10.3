!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine setrunname(scalestart,fscalestart,werkdir)
        use singletop2_scale_m, only: use_DDIS
      implicit none
      include 'types.f'
      include 'flags.f'
      include 'masses.f'
      include 'kprocess.f'
      include 'jetcuts.f'
      include 'pdlabel.f'
      include 'kpart.f'
      include 'vdecayid.f'
      include 'runstring.f'
      include 'dynamicscale.f'
      include 'taucut.f'
      include 'ewcorr.f'
      real(dp):: scalestart,fscalestart
      integer:: nlength
      character(len=255):: outlabel1
      character(len=1024):: runname
      character(len=*), intent(in) :: werkdir
      character(len=3):: strmh,getstr,strpt
      character(len=6):: strtaucut
      character(len=9):: strscale
      character(len=20):: part,kpartstring
      character(len=6):: case,kcasestring
      common/runname/runname
      common/nlength/nlength

c      if (abs(scalestart-fscalestart) < 1._dp) then
c--- if the scales are the same, use this scale as the label
c        strscale=getstr(int(scalestart))
c       else
c--- .... otherwise, use the percentage of (muR/muF)
c        strscale=getstr(int(scalestart/fscalestart*100._dp))
c      endif
      strscale=getstr(int(scalestart))//'__'
     &       //getstr(int(fscalestart))//'_'

      if (dynamicscale) then
        write(strscale,'(F4.2,"_",F4.2)') scalestart,fscalestart
      endif

      if (kcase==ktopanom .and. use_DDIS) then
        write (strscale,'(A)') "DDIS"
      endif

c--- convert kpart and kcase to strings
      part=kpartstring(kpart)
      case=kcasestring(kcase)

      if     ( (kcase==kWHbbar)
     &    .or. (kcase==kZHbbar)
     &    .or. (kcase==kqq_tth)
     &    .or. (kcase==ktottth)
     &    .or. (kcase==kHWW_4l)
     &    .or. (kcase==kHWW_tb)
     &    .or. (kcase==kHWWint)
     &    .or. (kcase==kHWWHpi)
     &    .or. (kcase==kggWW4l)
     &    .or. (kcase==kggVV4l)
     &    .or. (kcase==kHZZ_4l)
     &    .or. (kcase==kHZZ_tb)
     &    .or. (kcase==kHZZint)
     &    .or. (kcase==kHZZHpi)
     &    .or. (kcase==kggZZ4l)
     &    .or. (kcase==kggfus0)
     &    .or. (kcase==kggfus1)
     &    .or. (kcase==khjetma)
     &    .or. (kcase==kh2jmas)
     &    .or. (kcase==kggfus2)
     &    .or. (kcase==kggfus3) ) then
        strmh=getstr(int(hmass))
        outlabel1=case//'_'//trim(part)//'_'//trim(pdlabel)//'_'//strscale//
     &   '_'//strmh
      elseif (  (kcase==kH_1jet) ) then
        strmh=getstr(int(hmass))
        strpt=getstr(int(ptjetmin))
        outlabel1=case//'_'//trim(part)//'_'//trim(pdlabel)//'_'//strscale//
     &   '_'//strmh//'_pt'//strpt(1:2)
      elseif ( (kcase==kW_2jet)
     &    .or. (kcase==kZ_2jet) ) then
        if     (Gflag .eqv. .false.) then
          outlabel1=case//'_'//trim(part)//'_'//trim(pdlabel)//'_'//strscale//'_qrk'
        elseif (Qflag .eqv. .false.) then
          outlabel1=case//'_'//trim(part)//'_'//trim(pdlabel)//'_'//strscale//'_glu'
        else
          outlabel1=case//'_'//trim(part)//'_'//trim(pdlabel)//'_'//strscale
        endif
      else
        if     (kewcorr == kexact) then
          outlabel1=case//'_'//trim(part)//'_exact_'//trim(pdlabel)//'_'//strscale
        elseif (kewcorr == ksudakov) then
          outlabel1=case//'_'//trim(part)//'_sudakov_'//trim(pdlabel)//'_'//strscale
        else
          outlabel1=case//'_'//trim(part)//'_'//trim(pdlabel)//'_'//strscale
        endif
      endif

      if (usescet) then
        write(strtaucut,'(ES6.1E1)') taucut
        runname=trim(outlabel1)//'_'//strtaucut//'_'//trim(runstring)
      else
        if (vdecayid) then
         runname=trim(outlabel1)//'_'//v34id//v56id//'_'//trim(runstring)
        else
         runname=trim(outlabel1)//'_'//trim(runstring)
        endif
      endif

c--- add working directory, if necessary
      if (werkdir  /=  '') then
        runname=trim(werkdir)//'/'//trim(runname)
      endif

      nlength=len(trim(runname))

      return
      end


      function getstr(no)
      character(len=3)::getstr
c returns a string of length 3 from an integer::
      integer:: no,i1,i2,i3,izero

      izero=ichar('0')

      i1=abs(no)/100
      i2=(abs(no)-i1*100)/10
      i3=abs(no)-i1*100-i2*10

      if    (i1==0.and.i2==0) then
        if (no < 0) then
        getstr='-'//char(i3+izero)//'_'
        else
        getstr=char(i3+izero)//'__'
        endif
      elseif(i1==0) then
        getstr=char(i2+izero)//char(i3+izero)//'_'
      else
        getstr=char(i1+izero)//char(i2+izero)//char(i3+izero)
      endif

      return
      end

