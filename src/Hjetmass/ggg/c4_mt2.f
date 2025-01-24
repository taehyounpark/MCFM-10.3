c  Copyright (C) 2019-2022, respective authors of MCFM.

c  This program is free software: you can redistribute it and/or modify it under
c  the terms of the GNU General Public License as published by the Free Software
c  Foundation, either version 3 of the License, or (at your option) any later
c  version.

c  This program is distributed in the hope that it will be useful, but WITHOUT ANY
c  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
c  PARTICULAR PURPOSE. See the GNU General Public License for more details.

c  You should have received a copy of the GNU General Public License along with
c  this program. If not, see <http://www.gnu.org/licenses/>

      t1 = 0.1D1 / sman
      t2 = cdlogwrap(real(tman * t1,dp))
      t3 = cdlogwrap(real(uman * t1,dp))
      t4 = 0.1D1 / uman
      t5 = cdlogwrap(real(tman * t4,dp))
      t6 = mH ** 2
      t7 = sman + tman + uman
      t8 = uman ** 2
      t9 = cA * t8
      t10 = t7 * LogMumu
      t11 = t7 * LogMumt
      t2 = pi ** 2 - t2 ** 2 - t3 ** 2 - t5 ** 2
      t1 = dilogc(real(-t1 * t6 + 1,dp)) + dilogc(real(-t6 / tman + 1,dp)) + dilogc(real(-t4 *
     &t6 + 1,dp))
      t3 = sman + tman
      t4 = tman * sman
      c4mt2 = (0.7D1 / 0.540D3) * cA * sman * tman * t7 + (0.7D1 / 0.60D
     &2) * Kg * t8 * t7 - (0.7D1 / 0.18D2) * tr * t8 * t7 + (0.7D1 / 0.1
     &80D3) * LogMuMtop * t8 * t7 + (0.289D3 / 0.1080D4) * t9 * t7 + (0.
     &77D2 / 0.1080D4) * t9 * (LogMums * t7 + t10 + t11) - (0.7D1 / 0.36
     &0D3) * t9 * (sman * t2 + t2 * tman + t2 * uman) + (0.7D1 / 0.90D2)
     & * t8 * (cA * (sman * t1 + tman * t1 + uman * t1) + (-sman - tman
     &- uman) * LogMuMtop * tr) - (0.7D1 / 0.54D2) * tr * ((t3 * t8 + um
     &an * t8) * LogMums + t10 * t8 + t11 * t8 + t4 * t3) + (0.61D2 / 0.
     &135D3) * cF * t8 * t7 - (0.14D2 / 0.45D2) * t4 * tr * uman

