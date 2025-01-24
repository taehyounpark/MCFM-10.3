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
      t5 = cdlogwrap(real(tman * t3,dp))
      t6 = cdlogwrap(real(uman * t1,dp))
      t7 = sman ** 2
      t8 = (LogMums + LogMumt + LogMumu) * t7
      t9 = tman * uman
      c2mt0 = t7 * ((0.2D1 / 0.3D1) * LogMuMtop + 7 * cA - (0.20D2 / 0.3
     &D1) * tr) + (0.4D1 / 0.3D1) * t7 * ((dilogc(real(-t1 * t2 + 1,dp)) + dilogc
     &(real(-t2 / tman + 1,dp)) + dilogc(real(-t2 * t3 + 1,dp))) * cA - LogMuMtop * tr) +
     &(0.11D2 / 0.9D1) * t8 * cA - cA * t7 * (pi ** 2 - t4 ** 2 - t5 **
     &2 - t6 ** 2) / 3 - (0.20D2 / 0.9D1) * tr * (t9 + t8) - 2 * t7 * (-
     &Kg + cF) + (0.2D1 / 0.9D1) * t9 * cA

