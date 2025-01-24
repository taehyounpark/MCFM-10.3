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
      t3 = dilogc(real(-t1 * t2 + 1,dp))
      t4 = dilogc(real(-t2 / tman + 1,dp))
      t5 = 0.1D1 / uman
      t2 = dilogc(real(-t2 * t5 + 1,dp))
      t6 = cdlogwrap(real(uman * t1,dp))
      t5 = cdlogwrap(real(tman * t5,dp))
      t1 = cdlogwrap(real(tman * t1,dp))
      t7 = (sman + tman) * uman + sman * tman
      t8 = tman + uman
      t9 = sman ** 2
      t10 = sman * t9
      t11 = sman * t8
      t12 = t11 * LogMuMtop
      t13 = tman ** 2
      t14 = uman ** 2
      t15 = t9 + t13 + t14
      t16 = pi ** 2
      t6 = t6 ** 2
      t5 = t5 ** 2
      t1 = t1 ** 2
      t17 = t15 * LogMuMtop
      t18 = tman + uman
      t19 = t3 + t4 + t2
      t20 = uman * tman
      t21 = t20 * cA
      t22 = t14 + t13
      t23 = t15 * LogMumu
      t24 = t15 * LogMumt
      t25 = t9 * cA
      t26 = tr * t9
      t8 = t10 * t8
      t27 = -t1 + t16 - t6 - t5
      t27 = t25 * (t13 * t27 + t14 * t27 + t27 * t9)
      t28 = cA * t10
      t29 = LogMums + LogMumt + LogMumu
      t30 = t28 * (t29 * tman + t29 * uman)
      t31 = t21 * t9 * (LogMumt + LogMumu)
      t1 = (0.11D2 / 0.1512D4) * t25 * (LogMums * t15 + t23 + t24) - t27
     & / 504 + (0.167D3 / 0.3780D4) * t28 * t18 + (0.11D2 / 0.756D3) * t
     &30 + (0.43D2 / 0.15120D5) * t21 * sman * t18 - (0.1751D4 / 0.37800
     &D5) * t20 * sman * tr * t18 + (0.211D3 / 0.15120D5) * t31 - (0.23D
     &2 / 0.756D3) * t26 * t20 * t29 + (0.2216D4 / 0.14175D5) * cF * t10
     & * t18 + t9 * (cA * (t2 * t7 + t3 * t7 + t4 * t7) - t12 * tr) / 63
     & + t9 * ((t1 * t7 - t11 * t16 + t5 * t7 + t6 * t7) * cA + t17) / 2
     &52 + (0.1108D4 / 0.14175D5) * t9 * cF * t15 + Kg * t10 * t18 / 42
     &+ t9 * ((t13 * t19 + t14 * t19 + t19 * t9 + t20 * LogMuMtop) * cA
     &+ t12 - t17 * tr) / 126 + t21 * (t14 + t13) / 756 - (0.5D1 / 0.63D
     &2) * t10 * tr * t18 - (0.5D1 / 0.378D3) * tr * ((t22 * t9 + t9 **
     &2) * LogMums + t20 * t22 + t23 * t9 + t24 * t9) + (0.167D3 / 0.756
     &0D4) * t25 * t15 + t9 * Kg * t15 / 84 - (0.5D1 / 0.126D3) * t26 *
     &t15
      c2mt4 = -(0.5D1 / 0.189D3) * tr * (t8 * LogMums + t8 * LogMumt + t
     &8 * LogMumu + t13 * t14) + t1 + t20 * ((t9 * ((0.31D2 / 0.2160D4)
     &* LogMums - (0.13D2 / 0.2520D4) * t16 + (0.28429D5 / 0.635040D6))
     &+ t20 / 378) * cA + t9 * ((0.23D2 / 0.840D3) * Kg + (0.3721D4 / 0.
     &22680D5) * cF - (0.1411D4 / 0.9450D4) * tr) + t9 * (-(0.23D2 / 0.1
     &260D4) * tr + (0.23D2 / 0.2520D4)) * LogMuMtop)

