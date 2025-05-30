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
      t2 = mH ** 2
      t3 = 0.1D1 / uman
      t4 = cdlogwrap(real(tman * t1,dp))
      t5 = cdlogwrap(real(uman * t1,dp))
      t6 = cdlogwrap(real(tman * t3,dp))
      t1 = dilogc(real(-t1 * t2 + 1,dp)) + dilogc(real(-t2 / tman + 1,dp)) + dilogc(real(-t2 *
     &t3 + 1,dp))
      t2 = sman ** 2
      t3 = -pi ** 2 + t4 ** 2 + t5 ** 2 + t6 ** 2
      t4 = cA * t2
      t5 = sman + tman + uman
      t6 = tman + uman
      t7 = t5 * LogMumu
      t8 = t5 * LogMumt
      t9 = uman * tman
      c2mt2 = (0.7D1 / 0.90D2) * t2 * (cA * (sman * t1 + tman * t1 + uma
     &n * t1) + (-sman - tman - uman) * LogMuMtop * tr) + (0.7D1 / 0.360
     &D3) * t4 * (sman * t3 + tman * t3 + t3 * uman) + (0.7D1 / 0.540D3)
     & * cA * tman * uman * t5 - (0.7D1 / 0.54D2) * tr * ((sman * t2 + t
     &2 * t6) * LogMums + t7 * t2 + t8 * t2 + t9 * t6) + (0.77D2 / 0.108
     &0D4) * t4 * (LogMums * t5 + t7 + t8) + (0.289D3 / 0.1080D4) * t4 *
     & t5 + (0.7D1 / 0.60D2) * Kg * t2 * t5 + (0.7D1 / 0.180D3) * LogMuM
     &top * t2 * t5 - (0.7D1 / 0.18D2) * t2 * tr * t5 + (0.61D2 / 0.135D
     &3) * cF * t2 * t5 - (0.14D2 / 0.45D2) * t9 * sman * tr

