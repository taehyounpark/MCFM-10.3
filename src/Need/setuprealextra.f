!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine setuprealextra(nprocextra)
      implicit none
c Routines to setup process number corresponding to additional real radiation
c that must also run when using current value of nproc
      include 'nproc.f'
      integer nprocextra

      if     (nproc == 92) then
        nprocextra = 920
      elseif (nproc == 97) then
        nprocextra = 970
      elseif (nproc == 101) then
        nprocextra = 1010
      elseif (nproc == 114) then
        nprocextra = 115
      elseif (nproc == 141) then
        nprocextra = 142
      elseif (nproc == 146) then
        nprocextra = 147
      elseif (nproc == 161) then
        nprocextra = 162
      elseif (nproc == 166) then
        nprocextra = 167
      elseif (nproc == 171) then
        nprocextra = 172
      elseif (nproc == 176) then
        nprocextra = 177
      elseif (nproc == 181) then
        nprocextra = 182
      elseif (nproc == 186) then
        nprocextra = 187
      elseif (nproc == 233) then
        nprocextra = 234
      elseif (nproc == 238) then
        nprocextra = 239
      elseif (nproc == 501) then
        nprocextra = 502
      elseif (nproc == 511) then
        nprocextra = 512
      else
        write(6,*) 'Unexpected process in setuprealextra: ',nproc
        stop
      endif

      return
      end
