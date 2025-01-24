// Copyright (C) 2021 John Campbell, Stefan Hoeche, Christian T Preuss.

#ifndef MCFM__Flavor_Map_H
#define MCFM__Flavor_Map_H

// Standard includes. 
#include <map>
#include <vector>
#include <iostream>

namespace MCFM {

  struct Flavor_Table {
    std::map<int,int> m_spin, m_strong, m_anti;
    std::map<int,double> m_mass, m_width;
    Flavor_Table();
  };// end of struct Flavor_Table

  extern Flavor_Table s_flavors;

  struct Leg {
    int m_fl, m_id, m_is;
    inline Leg(int fl=0,int id=0,int is=0):
      m_fl(fl), m_id(id), m_is(is) {}
    inline Leg Cross() const
    { return Leg(s_flavors.m_anti[m_fl],m_id,!m_is); }
    inline int Kfcode() const { return abs(m_fl); }
    inline int IsAnti() const { return m_fl<0; }
    inline int Strong() const { return s_flavors.m_strong[m_fl]; }
    inline double Mass() const { return s_flavors.m_mass[m_fl]; }
    inline double Width() const { return s_flavors.m_width[m_fl]; }
    inline bool IsScalar() const  { return s_flavors.m_spin[m_fl]==0; }
    inline bool IsFermion() const { return s_flavors.m_spin[m_fl]==1; }
    inline bool IsVector() const  { return s_flavors.m_spin[m_fl]==2; }
  };// end of struct Leg

  std::ostream &operator<<(std::ostream &s,const Leg &l); 

  typedef std::map<int,int> FlavorMulti_Map;

  class Order_Flavor {
    FlavorMulti_Map *p_fmm;
    int Order_Multi(const Leg &a,const Leg &b)
    {
      if ((*p_fmm)[int(a.Kfcode())]==0 ||
	  (*p_fmm)[int(b.Kfcode())]==0) return 0;
      if ((*p_fmm)[int(a.Kfcode())]>
	  (*p_fmm)[int(b.Kfcode())]) return 1;
      return 0;
    }
    int Order_SVFT(const Leg &a,const Leg &b) {
      if (a.IsScalar() && !b.IsScalar()) return 1;
      if (a.IsVector() && !b.IsScalar() && !b.IsVector()) return 1;
      if (a.IsFermion() && !b.IsFermion() && 
	  !b.IsScalar() && !b.IsVector()) return 1;
      return 0;
    }
  public:
    Order_Flavor(FlavorMulti_Map *const fmm): p_fmm(fmm) {}
    int operator()(const Leg &a,const Leg &b) {
      if (!a.Strong()&&b.Strong()) return 1;
      if (a.Strong()&&!b.Strong()) return 0;
      if (a.Mass()>b.Mass()) return 1;
      if (a.Mass()<b.Mass()) return 0;
      if (Order_SVFT(a,b)) return 1;
      if (Order_SVFT(b,a)) return 0;
      if (!a.IsAnti()&&b.IsAnti()) return 1;
      if (a.IsAnti()&&!b.IsAnti()) return 0;
      if (p_fmm) {
	if (Order_Multi(a,b)) return 1;
	if (Order_Multi(b,a)) return 0;
      }
      return a.Kfcode()<b.Kfcode();
    }
  };// end of class Order_Flavor

}// end of namespace MCFM

#endif
