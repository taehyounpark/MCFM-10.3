!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

c--- converts integer kpart to corresponding 4-character string
      function kpartstring(k)
      implicit none
      include 'kpart.f'
      character(len=20):: kpartstring
      integer k

      if     (k == klord) then
        kpartstring='lo'
      elseif (k == kvirt) then
        kpartstring='virt'
      elseif (k == kreal) then
        kpartstring='real'
      elseif (k == ktota) then
        kpartstring='nlo'
      elseif (k == kfrag) then
        kpartstring='frag'
      elseif (k == ktodk) then
        kpartstring='todk'
      elseif (k == ksnlo) then
        kpartstring='snlo'
      elseif (k == knnlo) then
        kpartstring='nnlo'
      elseif (k == kn3lo) then
        kpartstring='n3lo'
      elseif (k == kresummed .and. kresorder == 2) then
          if (krespart == kresexp) then
              kpartstring = 'resexpLO'
          elseif (krespart == kresonly) then
              kpartstring = 'resonlyLO'
          elseif (krespart == kresabove) then
              kpartstring = 'resaboveLO'
          elseif (krespart == kresmatchcorr) then
              kpartstring = 'resmatchcorrLO'
          else
              kpartstring='resLO'
          endif
      elseif (k == kresummed .and. kresorder == 3) then
          if (krespart == kresexp) then
              kpartstring = 'resexpLOp'
          elseif (krespart == kresonly) then
              kpartstring = 'resonlyLOp'
          elseif (krespart == kresabove) then
              kpartstring = 'resaboveLOp'
          elseif (krespart == kresmatchcorr) then
              kpartstring = 'resmatchcorrLOp'
          else
              kpartstring='resLOp'
          endif
      elseif (k == kresummed .and. kresorder == 4) then
          if (krespart == kresexp) then
              kpartstring = 'resexpNLO'
          elseif (krespart == kresonly) then
              kpartstring = 'resonlyNLO'
          elseif (krespart == kresabove) then
              kpartstring = 'resaboveNLO'
          elseif (krespart == kresmatchcorr) then
              kpartstring = 'resmatchcorrNLO'
          else
              kpartstring='resNLO'
          endif
      elseif (k == kresummed .and. kresorder == 5) then
          if (krespart == kresexp) then
              kpartstring = 'resexpNLOp'
          elseif (krespart == kresonly) then
              kpartstring = 'resonlyNLOp'
          elseif (krespart == kresabove) then
              kpartstring = 'resaboveNLOp'
          elseif (krespart == kresmatchcorr) then
              kpartstring = 'resmatchcorrNLOp'
          else
              kpartstring='resNLOp'
          endif
      elseif (k == kresummed .and. kresorder == 6) then
          if (krespart == kresexp) then
              kpartstring = 'resexpNNLO'
          elseif (krespart == kresonly) then
              kpartstring = 'resonlyNNLO'
          elseif (krespart == kresabove) then
              kpartstring = 'resaboveNNLO'
          elseif (krespart == kresmatchcorr) then
              kpartstring = 'resmatchcorrNNLO'
          else
              kpartstring='resNNLO'
          endif
      elseif (k == kresummed .and. kresorder == 7) then
          if (krespart == kresexp) then
              kpartstring = 'resexpNNLOp'
          elseif (krespart == kresonly) then
              kpartstring = 'resonlyNNLOp'
          elseif (krespart == kresabove) then
              kpartstring = 'resaboveNNLOp'
          elseif (krespart == kresmatchcorr) then
              kpartstring = 'resmatchcorrNNLOp'
          else
              kpartstring='resNNLOp'
          endif

      elseif (k == kresummed .and. kresorder == 8) then
          if (krespart == kresexp) then
              kpartstring = 'resexpN3LO'
          elseif (krespart == kresonly) then
              kpartstring = 'resonlyN3LO'
          elseif (krespart == kresabove) then
              kpartstring = 'resaboveN3LO'
          elseif (krespart == kresmatchcorr) then
              kpartstring = 'resmatchcorrN3LO'
          else
              kpartstring='resN3LO'
          endif
      else
        write(6,*) 'Unexpected kpart in kpartstring: ',k
        stop
      endif

      if (coeffonly) then
        kpartstring=trim(kpartstring)//'coeff'
      endif

      return
      end


