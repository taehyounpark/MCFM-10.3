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
      t4 = cdlogwrap(real(uman * t1,dp))
      t5 = cdlogwrap(real(tman * t3,dp))
      t6 = cdlogwrap(real(tman * t1,dp))
      t4 = -pi ** 2 + t4 ** 2 + t5 ** 2 + t6 ** 2
      t5 = tman ** 2
      t6 = cA * t5
      t7 = sman + tman + uman
      t1 = dilogc(real(-t1 * t2 + 1,dp)) + dilogc(real(-t2 / tman + 1,dp)) + dilogc(real(-t2 *
     &t3 + 1,dp))
      t2 = t7 * LogMumu
      t3 = t7 * LogMumt
      t8 = sman + uman
      t9 = sman * uman
      c3mt2 = (0.7D1 / 0.360D3) * t6 * (sman * t4 + t4 * tman + t4 * uma
     &n) + (0.289D3 / 0.1080D4) * t6 * t7 + (0.7D1 / 0.90D2) * t5 * (cA
     &* (sman * t1 + tman * t1 + uman * t1) + (-sman - tman - uman) * Lo
     &gMuMtop * tr) + (0.77D2 / 0.1080D4) * t6 * (LogMums * t7 + t2 + t3
     &) - (0.7D1 / 0.54D2) * tr * ((t5 * t8 + tman * t5) * LogMums + t2
     &* t5 + t3 * t5 + t9 * t8) + (0.7D1 / 0.540D3) * t9 * cA * t7 - (0.
     &7D1 / 0.18D2) * t5 * tr * t7 + (0.7D1 / 0.60D2) * Kg * t5 * t7 + (
     &0.7D1 / 0.180D3) * LogMuMtop * t5 * t7 + (0.61D2 / 0.135D3) * cF *
     & t5 * t7 - (0.14D2 / 0.45D2) * t9 * tman * tr

