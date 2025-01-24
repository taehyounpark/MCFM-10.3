      function Li3(x)
      implicit none
      include 'types.f'
      real(dp):: Li3
      
C     returns Li_3(x) for real x, minf < x < 1
      real(dp):: x
      complex(dp):: wgplg
      
      if (x > 1._dp) then
        write(6,*) 'x>1 in Li3 function, src/Lib/Li3.f', x
        Li3 = 1d0
        return
      endif
      
      Li3 = real(wgplg(2,1,x))
      end

      function ReLi3(x)
      implicit none
      include 'types.f'
      real(dp):: ReLi3
      
      real(dp):: x
      complex(dp):: wgplg
      
      ReLi3 = real(wgplg(2,1,x))
      end
