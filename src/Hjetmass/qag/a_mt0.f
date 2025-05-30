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

      t1 = tman ** 2
      t2 = uman ** 2
      t3 = t2 + t1
      t4 = 0.1D1 / sman
      t5 = mH ** 2
      t6 = 0.1D1 / uman
      t7 = cdlog(t4 * uman)
      t8 = cdlog(t4 * tman)
      t9 = cdlog(tman * t6)
      t10 = 0.1D1 / t3
      t11 = t10 * t4
      t3 = t11 * t3
      t12 = LogMumt + LogMumu
      t12 = t1 * t12 + t12 * t2
      t13 = t3 * LogMums
      t1 = t2 + t1
      amt0 = -(0.844D3 / 0.9D1) * t3 + (0.40D2 / 0.9D1) * t10 * (uman +
     &tman) - (0.40D2 / 0.3D1) * t11 * t12 - (0.52D2 / 0.3D1) * t13 + (0
     &.68D2 / 0.27D2) * t3 * pi ** 2 + (0.104D3 / 0.27D2) * t3 * nl - 16
     & * t3 * cli2(-t4 * t5 + 1) - (0.64D2 / 0.9D1) * t11 * (t1 * cli2
     &(-t5 / tman + 1) + t1 * cli2(-t5 * t6 + 1)) - 4 * t11 * (t1 *
     &t7 ** 2 + t1 * t8 ** 2) + (0.4D1 / 0.9D1) * t11 * (t1 * t9 ** 2 +
     &nl * t12) - (0.16D2 / 0.3D1) * t3 * Kquark - (0.8D1 / 0.3D1) * t3
     &* Kgluon + (0.16D2 / 0.9D1) * t13 * nl

