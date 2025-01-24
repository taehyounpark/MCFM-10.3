!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine lumxmsq_wgamma(p,xx,z1,z2,QB,order,xmsq,central)
          use SCET
          use LHAPDF
          implicit none
          include 'types.f'
          include 'nf.f'
          include 'mxpart.f'
          include 'taucut.f'
          include 'facscale.f'
          include 'beamtype.f'

          real(dp), intent(in) :: p(mxpart,4), xx(2), z1, z2, QB(2)
          integer, intent(in) :: order
          real(dp), intent(out) :: xmsq
          logical, intent(in) :: central

          real(dp) :: soft1(-1:1), soft2(-1:3)
          real(dp) :: beama0(-5:5), beamb0(-5:5),
     &                beama1(-5:5, -1:1), beamb1(-5:5, -1:1),
     &                beama2(-5:5, -1:3), beamb2(-5:5, -1:3)
          real(dp) :: assemble_wgamma

          real(dp) :: msq(-nf:nf, -nf:nf, 0:2)
          real(dp) :: origtaucut

          integer:: m

          call softqqbis(order,soft1,soft2)

          if (order >= 0) then
              call fdist(ih1,xx(1),facscale,beama0,1)
              call fdist(ih2,xx(2),facscale,beamb0,2)
          endif
          if (order >= 1) then
              call xbeam1bis(ih1,z1,xx(1),QB(1),beama1,1)
              call xbeam1bis(ih1,z2,xx(2),QB(2),beamb1,2)
          endif
          if (order >= 2) then
              call xbeam2bis(ih1,z1,xx(1),QB(1),beama2,1)
              call xbeam2bis(ih2,z2,xx(2),QB(2),beamb2,2)
          endif

          call wgam_mat(p,msq)

          xmsq = assemble_wgamma(p,xx,order,soft1,soft2,
     &              beama0,beamb0,beama1,beamb1,beama2,beamb2,msq)

          ! only fill scetreweight for the central value binning
          ! not for pdf uncertainties and scale variation
          if (central .and. doMultitaucut) then
              scetreweight(:) = 0._dp
              ! is there any reasonable way for the other taucut values to give a non-zero result
              ! when xmsq is zero? I don't think so..
              if (xmsq /= 0._dp) then
                  origtaucut = taucut
                  do m=1,size(tcutarray)
                    taucut = tcutarray(m)
                    scetreweight(m) = assemble_wgamma(p,xx,order,soft1,soft2,
     &                      beama0,beamb0,beama1,beamb1,beama2,beamb2,msq)
                  enddo
                  taucut = origtaucut
                  scetreweight(:) = scetreweight(:) / xmsq
              endif
          endif

      end subroutine


      function assemble_wgamma(p,xx,order,soft1,soft2,beama0,beamb0,beama1,beamb1,beama2,beamb2,msq)
          use SCET
          implicit none
          include 'types.f'
          include 'nf.f'
          include 'constants.f'
          include 'mxpart.f'
          include 'qcdcouple.f'
          include 'tiny.f'
          include 'taucut.f'

          real(dp) :: assemble_wgamma
          real(dp), intent(in) :: p(mxpart,4), xx(2)
          integer, intent(in) :: order
          real(dp), intent(in) :: soft1(-1:1), soft2(-1:3)
          real(dp), intent(in) :: beama0(-5:5), beamb0(-5:5),
     &                beama1(-5:5, -1:1), beamb1(-5:5, -1:1),
     &                beama2(-5:5, -1:3), beamb2(-5:5, -1:3)
          real(dp), intent(in) :: msq(-nf:nf, -nf:nf, 0:2)

          real(dp) :: assemble, dot
          real(dp) :: Q, msqpow(-5:5,-5:5)
          real(dp) :: bit, hard(2), tauc

          integer :: j,k

          if (dynamictau) then
            tauc=getdynamictau(p)
          else
            tauc=taucut
          endif

c compute power corrections if required
          if ((incpowcorr) .or. (onlypowcorr)) then
              Q=sqrt(two*dot(p,1,2))
              call powcorr_qa(order,tauc,xx(1),xx(2),Q,beama0,beamb0,msqpow)
          endif

          assemble_wgamma = 0._dp

          do j=-nf,nf
            do k=-nf,nf
              if (msq(j,k,0) < tiny) cycle

              hard(1) = msq(j,k,1)/msq(j,k,0)/ason2pi
              hard(2) = msq(j,k,2)/msq(j,k,0)/ason2pi**2

              bit = assemble(order,tauc,beama0(j),beamb0(k),
     &                         beama1(j,:),beamb1(k,:),
     &                         beama2(j,:),beamb2(k,:),soft1,soft2,hard)

              if (incpowcorr) then
                bit=bit+msqpow(j,k)+msqpow(j,0)+msqpow(0,k)
              endif
              if (onlypowcorr) then
                bit=msqpow(j,k)+msqpow(j,0)+msqpow(0,k)
              endif

              bit = bit * msq(j,k,0)

              assemble_wgamma = assemble_wgamma + bit

            enddo
          enddo

      end function
