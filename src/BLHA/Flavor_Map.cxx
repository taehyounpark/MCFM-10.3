// Copyright (C) 2021 John Campbell, Stefan Hoeche, Christian T Preuss.
//
//  This program is free software: you can redistribute it and/or modify it under
//  the terms of the GNU General Public License as published by the Free Software
//  Foundation, either version 3 of the License, or (at your option) any later
//  version.
//
//  This program is distributed in the hope that it will be useful, but WITHOUT ANY
//  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
//  PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along with
//  this program. If not, see <http://www.gnu.org/licenses/>
 

#include "MCFM/Flavor_Map.h"

#include <iostream>

MCFM::Flavor_Table MCFM::s_flavors;

using namespace MCFM;

Flavor_Table::Flavor_Table() {
  for (int i(1);i<=6;++i) {
    m_spin[-i]=m_spin[i]=1;
    m_strong[-i]=-(m_strong[i]=3);
    m_anti[m_anti[i]=-i]=i;
  }
  for (int i(11);i<=16;++i) {
    m_spin[-i]=m_spin[i]=1;
    m_strong[-i]=m_strong[i]=0;
    m_anti[m_anti[i]=-i]=i;
  }
  m_spin[-24]=m_spin[24]=2;
  m_strong[-24]=m_strong[24]=0;
  m_anti[m_anti[24]=-24]=24;
  m_spin[23]=m_spin[22]=2;
  m_anti[23]=23;
  m_anti[22]=22;
  m_strong[23]=m_strong[22]=0;
  m_spin[21]=2;
  m_strong[21]=8;
  m_anti[21]=21;
  m_spin[25]=0;
  m_strong[25]=0;
  m_anti[25]=25;
}

std::ostream &MCFM::operator<<(std::ostream &s,const Leg &l) {
  return s<<"(id="<<l.m_id<<",is="<<l.m_is<<",fl="<<l.m_fl<<")";
}
