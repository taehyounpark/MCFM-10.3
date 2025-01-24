
      integer :: lbjscheme
      common/lbjscheme_com/lbjscheme
      logical:: jetsopphem
      common /cjetsopphem/jetsopphem
      real(dp) :: leptptmin,leptptmax,leptrapmin,leptrapmax
      real(dp) :: misspt,Rjlmin,Rllmin,delyjjmin
      real(dp) :: leptpt2min,leptpt2max,leptrap2min,leptrap2max
      real(dp) :: leptpt3min,leptpt3max,leptrap3min,leptrap3max
      real(dp) :: gammptmin,gammptmax,gammrapmin,gammrapmax
      real(dp) :: Rgalmin,mtrans34cut
      real(dp) :: gammpt2,Rgagamin,gammpt3,Rgajetmin
      real(dp) :: leptveto1min,leptveto1max,leptveto2min,leptveto2max
      real(dp) :: gammvetomin,gammvetomax
      real(dp) :: missrelpt,mllmin,mllmax,ptllmin,gammptprod,m3lmin
      common/leptcuts0/leptptmin,leptptmax,leptrapmin,leptrapmax
      common/leptcuts1/misspt,Rjlmin,Rllmin,delyjjmin
      common/leptcuts2a/leptpt2min,leptpt2max,leptrap2min,leptrap2max
      common/leptcuts2b/leptpt3min,leptpt3max,leptrap3min,leptrap3max
      common/leptcuts3/gammptmin,gammptmax,gammrapmin,gammrapmax
      common/leptcuts4/Rgalmin,mtrans34cut
      common/leptcuts5/gammpt2,Rgagamin,gammpt3,Rgajetmin
      common/leptcuts6/leptveto1min,leptveto1max,leptveto2min,leptveto2max
      common/leptcuts7/gammvetomin,gammvetomax
      common/leptcuts8/missrelpt,mllmin,mllmax,ptllmin,gammptprod,m3lmin
      real(dp) :: elptmin,muptmin,Relelmin,Relmumin,Rmumumin,elrapmax,murapmax
      real(dp) :: elvetomin,elvetomax,muvetomin,muvetomax
      common/elmucuts0/elptmin,muptmin,Relelmin,Relmumin,Rmumumin,elrapmax,murapmax
      common/elmucuts1/elvetomin,elvetomax,muvetomin,muvetomax
      real(dp) :: Rlepiso,fraclepiso
      common/lepiso/Rlepiso,fraclepiso

      real(dp) :: y34min, y34max
      common/leptcuts7/y34min,y34max
