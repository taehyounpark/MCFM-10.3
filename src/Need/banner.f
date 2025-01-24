!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine banner
          use m_config, only: cfg_string_len, cfg_get
      implicit none
      include 'types.f'

      logical, parameter :: prelim = .false.

c--- warning message, if necessary
      if (prelim) then
        write(6,*) '*                                                *'
        write(6,*) '*           PRELIMINARY VERSION                  *'
        write(6,*) '*                                                *'
        write(6,*) '*  NOTE: This is a private release of the MCFM   *'
        write(6,*) '*  code that has not yet been made public on the *'
        write(6,*) '*  usual website. As such:                       *'
        write(6,*) '*                                                *'
        write(6,*) '*   + Please do not redistribute without the     *'
        write(6,*) '*     knowledge of the authors;                  *'
        write(6,*) '*                                                *'
        write(6,*) '*   + Please notify the authors of any bugs      *'
        write(6,*) '*     or problems so that they can be corrected  *'
        write(6,*) '*     before the next official release.          *'
        write(6,*) '*                                                *'
      endif

      write (6,*) '**************** MCFM - version 10.3 ***************'
      write (6,*) '*                                                  *'
      write (6,*) '*  MCFM, v10.3                     January 2023    *'
      write (6,*) '*  CuTe-MCFM, v1.2                                 *'
      write (6,*) '*                                                  *'
      write (6,*) '*  On the web: https://mcfm.fnal.gov/              *'
      write (6,*) '*                                                  *'
      write (6,*) '*  MCFM Authors:                                   *'
      write (6,*) '*                                                  *'
      write (6,*) '*   John Campbell <johnmc@fnal.gov>                *'
      write (6,*) '*   Keith Ellis <ellis@fnal.gov>                   *'
      write (6,*) '*   Tobias Neumann <tneumann@fnal.gov>             *'
      write (6,*) '*   Ciaran Williams <ciaranwi@buffalo.edu>         *'
      write (6,*) '*                                                  *'
      write (6,*) '*  CuTe-MCFM Authors:                              *'
      write (6,*) '*                                                  *'
      write (6,*) '*   Thomas Becher <becher@itp.unibe.ch>            *'
      write (6,*) '*   Tobias Neumann <tneumann@fnal.gov>             *'
      write (6,*) '*                                                  *'
      write (6,*) '*   See https://mcfm.fnal.gov/                     *'
      write (6,*) '*     for a full list of contributors.             *'
      write (6,*) '*                                                  *'
      write (6,*) '****************************************************'

      call writereference

      write (6,*) '****************************************************'
      write (6,*) '*   MCFM uses the libraries                        *'
      write (6,*) '*                                                  *'
      write (6,*) '*    AMOS (Netlib)                                 *'
      write (6,*) '*    Chaplin 1.2 (Buehler, Duhr)                   *'
      write (6,*) '*    HandyG 0.1.4 (Naterop, Signer, Ulrich)        *'
      write (6,*) '*    hplog6 1.6 (Gehrmann, Remiddi)                *'
      write (6,*) '*    LHAPDF 6.5.1 (Buckley, et al.)                *'
      write (6,*) '*    QCDLoop 2.0.9 (Carazza, Ellis, Zanderighi)    *'
      write (6,*) '*    OneLOop (van Hameren)                         *'
      write (6,*) '*    Quad Double 2.3.22 (Hida, Li, Bailey)         *'
      write (6,*) '*                                                  *'
      write (6,*) '****************************************************'

      return
      end








