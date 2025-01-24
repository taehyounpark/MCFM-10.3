!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine parser(pre,post,helname,KCD,KCC,KCB,rattot,tot)
      implicit none
      include 'types.f'
      include 'WWjetlabels.f'
      include 'KCdef.f'
      character(len=4):: helname
      character*(*) pre,post
      character(len=40):: name
      character(len=99):: line
      integer k,ii
      real(dp):: xr, xi

      KCD(:)=0._dp
      KCC(:)=0._dp
      KCB(:)=0._dp

      name=trim(pre)//helname//trim(post)//'.f'

c      write(6,*) name

      open(unit=31,file=name,status='unknown')

      read(31,'(14x,D25.15,x,D25.15,x)') xr,xi
      rattot=cmplx(xr,xi, kind=dp)
      read(31,'(14x,D25.15,x,D25.15,x)') xr,xi
      tot=cmplx(xr,xi, kind=dp)

c      write(6,*) 'rat',rattot
c      write(6,*) 'tot',tot

      do k=1,5
      read(31,*)
      enddo

   91 continue

      read(31,'(a99)') line
c      write(6,*) line

      if (index(line,'KCD') > 0) then
        ii=0
        if     (index(line,'KCD(d7x2x56    )') > 0) then
          ii=d7x2x56
        elseif (index(line,'KCD(d7x2x34    )') > 0) then
          ii=d7x2x34
        elseif (index(line,'KCD(d2x17x56   )') > 0) then
          ii=d2x17x56
        elseif (index(line,'KCD(d17x2x56   )') > 0) then
          ii=d17x2x56
        elseif (index(line,'KCD(d7x1x56    )') > 0) then
          ii=d7x1x56
        elseif (index(line,'KCD(d7x1x34    )') > 0) then
          ii=d7x1x34
        elseif (index(line,'KCD(d27x1x34   )') > 0) then
          ii=d27x1x34
        elseif (index(line,'KCD(d1x27x34   )') > 0) then
          ii=d1x27x34
        elseif (index(line,'KCD(d1x7x2     )') > 0) then
          ii=d1x7x2
        elseif (index(line,'KCD(d2x56x7sl  )') > 0) then
          ii=d2x56x7sl
        elseif (index(line,'KCD(d2x7x56sl  )') > 0) then
          ii=d2x7x56sl
        elseif (index(line,'KCD(d2x34x7sl  )') > 0) then
          ii=d2x34x7sl
        elseif (index(line,'KCD(d2x7x34sl  )') > 0) then
          ii=d2x7x34sl
        elseif (index(line,'KCD(d17x2x34sl )') > 0) then
          ii=d17x2x34sl
        elseif (index(line,'KCD(d17x2x56sl )') > 0) then
          ii=d17x2x56sl
        elseif (index(line,'KCD(d1x56x7sl  )') > 0) then
          ii=d1x56x7sl
        elseif (index(line,'KCD(d1x7x56sl  )') > 0) then
          ii=d1x7x56sl
        elseif (index(line,'KCD(d1x34x7sl  )') > 0) then
          ii=d1x34x7sl
        elseif (index(line,'KCD(d1x7x34sl  )') > 0) then
          ii=d1x7x34sl
        elseif (index(line,'KCD(d27x1x34sl )') > 0) then
          ii=d27x1x34sl
        elseif (index(line,'KCD(d1x27x34sl )') > 0) then
          ii=d1x27x34sl
        elseif (index(line,'KCD(d2x1x34sl  )') > 0) then
          ii=d2x1x34sl
        elseif (index(line,'KCD(d7x34x12sl )') > 0) then
          ii=d7x34x12sl
        elseif (index(line,'KCD(d1x2x56sl  )') > 0) then
          ii=d1x2x56sl
        elseif (index(line,'KCD(d7x12x34sl )') > 0) then
          ii=d7x12x34sl
        elseif (index(line,'KCD(d1x2x7sl   )') > 0) then
          ii=d1x2x7sl
        elseif (index(line,'KCD(d2x1x56sl  )') > 0) then
          ii=d2x1x56sl
        elseif (index(line,'KCD(d1x2x34sl  )') > 0) then
          ii=d1x2x34sl
        elseif (index(line,'KCD(d12x7x34sl )') > 0) then
          ii=d12x7x34sl
        elseif (index(line,'KCD(d2x1x7sl   )') > 0) then
          ii=d2x1x7sl
        elseif (index(line,'KCD(d23x1x45q  )') > 0) then
          ii=d23x1x45q
        elseif (index(line,'KCD(d23x1x67q  )') > 0) then
          ii=d23x1x67q
        elseif (index(line,'KCD(d45x1x67q  )') > 0) then
          ii=d45x1x67q
        endif
        read(line,'(10x,11x,9x,D25.15,1x,D25.15,1x)') xr,xi
        if (ii > 0) KCD(ii) = cmplx(xr,xi, kind=dp)
c        write(6,*) k,ii,KCD(ii)
      endif

      if (index(line,'KCD') > 0) goto 91

   92 continue

      read(31,'(a99)') line
c      write(6,*) line

      if (index(line,'KCC') > 0) then
        ii=0
        if     (index(line,'KCC(c7x134   )') > 0) then
          ii=c7x134
        elseif (index(line,'KCC(c27x56   )') > 0) then
          ii=c27x56
        elseif (index(line,'KCC(c7x156   )') > 0) then
          ii=c7x156
        elseif (index(line,'KCC(c34x27   )') > 0) then
          ii=c34x27
        elseif (index(line,'KCC(c56x17   )') > 0) then
          ii=c56x17
        elseif (index(line,'KCC(c17x34   )') > 0) then
          ii=c17x34
        elseif (index(line,'KCC(c2x7         )') > 0) then
          ii=c2x7
        elseif (index(line,'KCC(c2x56         )') > 0) then
          ii=c2x56
        elseif (index(line,'KCC(c2x34         )') > 0) then
          ii=c2x34
        elseif (index(line,'KCC(c56x34   )') > 0) then
          ii=c56x34
        elseif (index(line,'KCC(c1x7         )') > 0) then
          ii=c1x7
        elseif (index(line,'KCC(c1x56         )') > 0) then
          ii=c1x56
        elseif (index(line,'KCC(c1x34         )') > 0) then
          ii=c1x34
        elseif (index(line,'KCC(c1x27         )') > 0) then
          ii=c1x27
        elseif (index(line,'KCC(c27x56sl )') > 0) then
          ii=c27x56sl
        elseif (index(line,'KCC(c34x27sl )') > 0) then
          ii=c34x27sl
        elseif (index(line,'KCC(c17x56sl )') > 0) then
          ii=c17x56sl
        elseif (index(line,'KCC(c17x34sl )') > 0) then
          ii=c17x34sl
        elseif (index(line,'KCC(c34x56sl )') > 0) then
          ii=c34x56sl
        elseif (index(line,'KCC(c12x34sl )') > 0) then
          ii=c12x34sl
        elseif (index(line,'KCC(c12x56sl )') > 0) then
          ii=c12x56sl
        elseif (index(line,'KCC(c45x67q  )') > 0) then
          ii=c45x67q
        elseif (index(line,'KCC(c23x45q  )') > 0) then
          ii=c23x45q
        elseif (index(line,'KCC(c23x67q  )') > 0) then
          ii=c23x67q
        elseif (index(line,'KCC(c1x23q   )') > 0) then
          ii=c1x23q
        elseif (index(line,'KCC(c1x45q   )') > 0) then
          ii=c1x45q
        elseif (index(line,'KCC(c1x67q   )') > 0) then
          ii=c1x67q
        endif
        read(line,'(10x,9x,9x,D25.15,1x,D25.15,1x)') xr,xi
        if (ii > 0) KCC(ii) = cmplx(xr,xi, kind=dp)
c        write(6,*) k,ii,KCC(ii)
      endif

      if (index(line,'KCC') > 0) goto 92
      backspace(31)

   93 continue

      read(31,'(a99)') line

      if (index(line,'KCB') > 0) then
        ii=0
        if     (index(line,'KCB(b256  )') > 0) then
          ii=b256
        elseif (index(line,'KCB(b134  )') > 0) then
          ii=b134
        elseif (index(line,'KCB(b56   )') > 0) then
          ii=b56
        elseif (index(line,'KCB(b234  )') > 0) then
          ii=b234
        elseif (index(line,'KCB(b156  )') > 0) then
          ii=b156
        elseif (index(line,'KCB(b34   )') > 0) then
          ii=b34
        elseif (index(line,'KCB(b17   )') > 0) then
          ii=b17
        elseif (index(line,'KCB(b127  )') > 0) then
          ii=b127
        elseif (index(line,'KCB(b27   )') > 0) then
          ii=b27
        elseif (index(line,'KCB(b256sl)') > 0) then
          ii=b256sl
        elseif (index(line,'KCB(b134sl)') > 0) then
          ii=b134sl
        elseif (index(line,'KCB(b56sl )') > 0) then
          ii=b56sl
        elseif (index(line,'KCB(b27sl )') > 0) then
          ii=b27sl
        elseif (index(line,'KCB(b234sl)') > 0) then
          ii=b234sl
        elseif (index(line,'KCB(b156sl)') > 0) then
          ii=b156sl
        elseif (index(line,'KCB(b34sl )') > 0) then
          ii=b34sl
        elseif (index(line,'KCB(b567sl)') > 0) then
          ii=b567sl
        elseif (index(line,'KCB(b347sl)') > 0) then
          ii=b347sl
        elseif (index(line,'KCB(b127sl)') > 0) then
          ii=b127sl
        elseif (index(line,'KCB(b17sl )') > 0) then
          ii=b17sl
        elseif (index(line,'KCB(b45q  )') > 0) then
          ii=b45q
        elseif (index(line,'KCB(b145q )') > 0) then
          ii=b145q
        elseif (index(line,'KCB(b67q  )') > 0) then
          ii=b67q
        elseif (index(line,'KCB(b123q )') > 0) then
          ii=b123q
        elseif (index(line,'KCB(b23q  )') > 0) then
          ii=b23q
        elseif (index(line,'KCB(b167q )') > 0) then
          ii=b167q
        endif
        read(line,'(10x,6x,9x,D25.15,1x,D25.15,1x)') xr,xi
        if (ii > 0) KCB(ii) = cmplx(xr,xi, kind=dp)
c        write(6,*) k,ii,KCB(ii)
      endif

      if (index(line,'KCB') > 0) goto 93

      close(unit=31)
c      pause

      return
      end

