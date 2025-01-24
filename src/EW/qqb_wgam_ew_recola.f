!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_wgam_ew_recola(pmcfm,msqnlo)
c Each program, which uses RECOLA must have:
      use recola
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nwz.f'
      include 'masses.f'
      include 'scale.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'first.f'
      include 'epinv.f'
      include 'scheme.f'
      include 'zcouple_cms.f'
c Variables for this demo file
      real (dp)          :: p(0:3,1:5),MW,MZ,WW,WZ,GetDeltarr
      integer::colour(5),hel(5),j,k,ii
      complex(dp)::Anlo(2)
      real(dp) :: msqlo(-nf:nf,-nf:nf),msqnlo(-nf:nf,-nf:nf),
     & msqnlosp(-nf:nf,-nf:nf),msqnlodp(-nf:nf,-nf:nf)
      integer, save:: i0,ip,jj(60:70),kk(60:70)
      real(dp), parameter::
     & Qfsq=(xn*(3._dp*(-1._dp/3._dp)**2+2._dp*(2._dp/3._dp)**2)+3._dp)

      real(dp)::pmcfm(mxpart,4)

c      logical:: recola_use_gfermi
c      common /recola_use_gfermi/recola_use_gfermi

      call set_delta_uv_rcl(epinv)
      call set_delta_ir_rcl(epinv,epinv**2)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c Step 1
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c Set the inputs for the computation.
c The variables which can be set and the subroutines to set them are
c in the file "input.f90".
c All variables have default values, so this step is optional.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c      MZ = 9.1153480619183000D+01;WZ = 2.4942663787727999D+00
c      MZ=sqrt(2d0)*MZ
c      WZ=sqrt(2d0)*WZ
c      call set_pole_mass_w_rcl(MW,WW)
c      call set_pole_mass_z_rcl(MZ,WZ)

      scheme = 'tH-V'

      if (first) then
        first=.false.
c The standard output is selected
        call set_output_file_rcl('wgam_ew')
        call set_dynamic_settings_rcl(1)
c        call use_dim_reg_soft_rcl
c Map MCFM parameters onto those used in Recola
        call set_pole_mass_w_rcl(wmass,wwidth)
        call set_pole_mass_z_rcl(zmass,zwidth)
        call use_gfermi_scheme_rcl(a=real(zaemmz,kind=dp))
c        call use_alphaZ_scheme_rcl(a=real(zaemmz,kind=dp))
c        if (recola_use_gfermi) then
c          call use_gfermi_scheme_rcl(a=real(zaemmz,kind=dp))
c        else
c          call use_alpha0_scheme_rcl(a=real(zaemmz,kind=dp))
c        endif

c Let's print the squared amplitude
        call set_print_level_parameters_rcl(0)
        call set_print_level_amplitude_rcl(0)
        call set_print_level_squared_amplitude_rcl(0)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c Step 2
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c Define the processes to be computed and select the power of g_s
c (i.e. the strong coupling: g_s^2 = 4*pi*alpha_s), by calling
c the subroutines of "process_definition.f90".
c The processes are defined calling subroutine "define_process_rcl"
c successively with different process number argument. At least one
c call of define_process_rcl must be present.
c All power of g_s are selected by default, the call of the
c subroutines for their selection are optional.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        if (nwz==+1) then
        ip=60; i0=ip
        call define_process_rcl(ip,'u d~ -> nu_e e+ A','NLO')
        call unselect_all_gs_powers_BornAmpl_rcl(ip)
        call unselect_all_gs_powers_LoopAmpl_rcl(ip)
        call select_gs_power_BornAmpl_rcl(ip,0)
        call select_gs_power_LoopAmpl_rcl(ip,0)
        jj(ip) = 2; kk(ip) = -1

        ip=ip+1
        call define_process_rcl(ip,'d~ u -> nu_e e+ A','NLO')
        call unselect_all_gs_powers_BornAmpl_rcl(ip)
        call unselect_all_gs_powers_LoopAmpl_rcl(ip)
        call select_gs_power_BornAmpl_rcl(ip,0)
        call select_gs_power_LoopAmpl_rcl(ip,0)
        jj(ip) = -1; kk(ip) = 2

        elseif (nwz==-1) then
        ip=65; i0=ip
        call define_process_rcl(ip,'d u~ -> e- nu_e~ A','NLO')
        call unselect_all_gs_powers_BornAmpl_rcl(ip)
        call unselect_all_gs_powers_LoopAmpl_rcl(ip)
        call select_gs_power_BornAmpl_rcl(ip,0)
        call select_gs_power_LoopAmpl_rcl(ip,0)
        jj(ip) = 1; kk(ip) = -2

        ip=ip+1
        call define_process_rcl(ip,'u~ d -> e- nu_e~ A','NLO')
        call unselect_all_gs_powers_BornAmpl_rcl(ip)
        call unselect_all_gs_powers_LoopAmpl_rcl(ip)
        call select_gs_power_BornAmpl_rcl(ip,0)
        call select_gs_power_LoopAmpl_rcl(ip,0)
        jj(ip) = -2; kk(ip) = 1
        endif

        call generate_processes_rcl

      endif

      call set_mu_uv_rcl(scale)
      call set_mu_ir_rcl(scale)
      call set_alphas_rcl(as,scale,5)

      p(:,1)=-[pmcfm(1,4),pmcfm(1,1),pmcfm(1,2),pmcfm(1,3)]
      p(:,2)=-[pmcfm(2,4),pmcfm(2,1),pmcfm(2,2),pmcfm(2,3)]
      p(:,3)=+[pmcfm(3,4),pmcfm(3,1),pmcfm(3,2),pmcfm(3,3)]
      p(:,4)=+[pmcfm(4,4),pmcfm(4,1),pmcfm(4,2),pmcfm(4,3)]
      p(:,5)=+[pmcfm(5,4),pmcfm(5,1),pmcfm(5,2),pmcfm(5,3)]


c     qqb
c      call compute_process_rcl(i0,p,'NLO')
c      colour(1:5)=[0,1,0,0,0]
c      hel(1:5)=[-1,+1,-1,+1,-1]
c      call get_Amplitude_rcl(i0,0,'NLO',colour,hel,Anlo(1))
c      write(6,*) 'i0,Anlo',i0,Anlo(1)

      msqlo=0._dp
      msqnlo=0._dp
      msqnlosp=0._dp
      msqnlodp=0._dp

      do ii = i0,ip

c      call set_delta_uv_rcl(0d0)
c      call set_delta_ir_rcl(0d0,0d0)
c      call compute_process_rcl(ii,p,'NLO')
c      call get_squared_amplitude_rcl(ii, 0, 'LO', msqlo(jj(ii),kk(ii)))
c      call get_squared_amplitude_rcl(ii, 0, 'NLO', msqnlo(jj(ii),kk(ii)))

c      call set_delta_uv_rcl(1d0)
c      call set_delta_ir_rcl(1d0,0d0)
c      call compute_process_rcl(ii,p,'NLO')
c      call get_squared_amplitude_rcl(ii, 0, 'NLO', msqnlosp(jj(ii),kk(ii)))

c      call set_delta_uv_rcl(0d0)
c      call set_delta_ir_rcl(0d0,1d0)
c      call compute_process_rcl(ii,p,'NLO')
c      call get_squared_amplitude_rcl(ii, 0, 'NLO', msqnlodp(jj(ii),kk(ii)))

c      call set_delta_uv_rcl(1d0)
c      call set_delta_ir_rcl(1d0,0d0)
      call compute_process_rcl(ii,p,'NLO')
      call get_squared_amplitude_rcl(ii, 0, 'LO', msqlo(jj(ii),kk(ii)))
      call get_squared_amplitude_rcl(ii, 0, 'NLO', msqnlo(jj(ii),kk(ii)))

      enddo

c      msqnlosp=msqnlosp-msqnlo
c      msqnlodp=msqnlodp-msqnlo

c      call reset_recola_rcl

c Repeat from 1st generation to 2nd
      do j=3,4
      do k=3,4
      msqlo(j,-k)=msqlo(j-2,-k+2)
      msqlo(-j,k)=msqlo(-j+2,k-2)
      msqnlo(j,-k)=msqnlo(j-2,-k+2)
      msqnlo(-j,k)=msqnlo(-j+2,k-2)
      enddo
      enddo

c Translate from alpha(MZ) scheme to a hybrid scheme in which one power
c of alpha is treated as alpha(0) due to having one external photon
c      msqnlo(:,:)=msqnlo(:,:)
c     &-msqlo(:,:)*(abs(zesq)/fourpi)/twopi*Qfsq
c     & *1d0/3d0*(2d0*epinv+2d0*log(musq/zmass**2)+16d0/3d0)

c Translate from Gmu scheme to a hybrid scheme in which one power
c of alpha is treated as alpha(0) due to having one external photon
      msqnlo(:,:)=msqnlo(:,:)
     &-msqlo(:,:)*((abs(zesq)/fourpi)/twopi*Qfsq
     & *1d0/3d0*(2d0*epinv+2d0*log(musq/zmass**2))-GetDeltarr())

c Translation to overall prefactor 1/Gamma(1-e) instead of Gamma(1+e)
c using the fact that double pole = (-14/9)*(alpha/2/pi)*LO
      msqnlo(:,:)=msqnlo(:,:)
     &+msqlo(:,:)*(abs(zesq)/fourpi)/twopi*(-14d0/9d0)*pi**2/6._dp

c      if (recola_use_gfermi .eqv. .false.) then
c! this appears to convert poles from alpha0 scheme to expected result ... ?
c        msqnlo(:,:)=msqnlo(:,:)
c     &  +msqlo(:,:)*2d0*abs(zesq)/fourpi**2*(4d0/3d0*(3d0*(2d0*4d0/9d0+3d0*1d0/9d0)+3d0))
c     &   *(epinv+log(musq/wmass**2))
c      else
c        msqnlo(:,:)=msqnlo(:,:)
c     &  -msqlo(:,:)*abs(zesq)/fourpi**2*(4d0/3d0*(3d0*(2d0*4d0/9d0+3d0*1d0/9d0)+3d0))
c     &   *(epinv+log(musq/wmass**2))
c      endif

      return
      end
