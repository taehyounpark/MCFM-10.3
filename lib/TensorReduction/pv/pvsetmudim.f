      subroutine pvsetmudim(mu)
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/TRscale.f'
      real(dp):: mu
      scale=mu
      musq=scale**2
      return
      end
