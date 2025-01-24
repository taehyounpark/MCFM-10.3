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
 

// MCFM routine for loop MEs.
extern "C" { void qqb_gam2j_v_(double *p,double *msqv); }

namespace MCFM {

  class qqb_gam2j: public Process {
  private:

    int m_type, m_id[4];
    
  public:

    static bool InitializeProcess(CXX_Interface *const interface,
				  const Process_Info &pi,
				  const std::vector<Leg> &legs) {
      if (pi.m_oqcd!=3 || pi.m_oew!=1) return false;
      if (legs.size()!=5) return false;
      if (pi.m_decids.size()) return false;
      for (size_t i(0);i<legs.size();++i)
	if (legs[i].Mass()) return false;
      if (legs[0].m_fl!=22) return false;
      if (legs[1].m_fl==21 && legs[2].m_fl==21) {
	if (legs[3].m_fl>0 && legs[3].m_fl<6 && legs[4].m_fl==-legs[3].m_fl)// ggqa
	  return interface->AddProcess(pi,new qqb_gam2j(legs,1))>-1;
      }
      if (legs[1].m_fl>0 && legs[1].m_fl<6 && legs[3].m_fl==-legs[1].m_fl) {
       	if (legs[2].m_fl==legs[1].m_fl && legs[4].m_fl==-legs[2].m_fl)// aaaa
	  return interface->AddProcess(pi,new qqb_gam2j(legs,2))>-1;
	if (legs[2].m_fl>0 && legs[2].m_fl<6 && legs[4].m_fl==-legs[2].m_fl)// qrqr
	  return interface->AddProcess(pi,new qqb_gam2j(legs,3))>-1;
      }
      return false;
    }

    inline void SwapIDs(int i1,int i2) {
      std::swap<int>(m_id[i1],m_id[i2]);
      m_cfac=ISSymmetryFactor(m_legs,0);
    }

    qqb_gam2j(const std::vector<Leg> &legs,int type):
      Process(legs,3,4), m_type(type), m_id{3,4,1,2} {
      static int first(true);
      if (first) {
	first=false;
        nproc_.nproc=280;
	blha_.useblha=1;
	chooser_();
      }
      m_res.resize(4);
      if (m_type==1) {
	if (m_legs[1].m_is && m_legs[3].m_is) SwapIDs(1,2);
	if (m_legs[2].m_is && m_legs[3].m_is) SwapIDs(1,3);
	if (m_legs[1].m_is && m_legs[4].m_is) SwapIDs(0,2);
	if (m_legs[2].m_is && m_legs[4].m_is) SwapIDs(0,3);
      }
    }

    void Calc(const std::vector<FourVec> &p,int oqcd) {
      // Convert momenta.
      SetMom(p_p,p,m_legs[m_id[0]],0);// IS
      SetMom(p_p,p,m_legs[m_id[1]],1);// IS
      SetMom(p_p,p,m_legs[0],2);// FS y
      SetMom(p_p,p,m_legs[m_id[2]],3);// FS
      SetMom(p_p,p,m_legs[m_id[3]],4);// FS
      // Calculate result.
      blha_.blhatype=m_type;
      if (m_type==1) {
	flags_.qflag=false;
	flags_.gflag=true;
	blha_.blhafl[0]=m_legs[m_id[0]].m_fl;
	blha_.blhafl[1]=m_legs[m_id[1]].m_fl;
      }
      else {
	flags_.gflag=false;
	flags_.qflag=true;
      }
      epinv_.epinv=epinv2_.epinv2=0.0;
      qqb_gam2j_v_(p_p,p_msqv);
      double res = p_msqv[MSQId(m_legs[m_id[0]],m_legs[m_id[1]])];
      m_res[0] = res*m_cfac;
      if (!m_polecheck) return;
      epinv_.epinv=1.0;
      qqb_gam2j_v_(p_p,p_msqv);
      double res1 = p_msqv[MSQId(m_legs[m_id[0]],m_legs[m_id[1]])];
      epinv2_.epinv2=1.0;
      qqb_gam2j_v_(p_p,p_msqv);
      double res2 = p_msqv[MSQId(m_legs[m_id[0]],m_legs[m_id[1]])];
      m_res[1]=(res1-res)*m_cfac;
      m_res[2]=(res2-res1)*m_cfac;
      m_res[3]=m_res[2]/(-qcdcouple_.ason2pi*m_bfac);
    }

    int GetScheme() const { return 1; }

  };// end of class qqb_gam2j

}// end of namespace MCFM
