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
extern "C" {
  void qqb_ww_v_(double *p,double *msqv);
  void qqb_ww_(double *p,double *msqv);}

namespace MCFM {

  class qqb_ww: public Process {
  private:

    int m_type;

  public:

    static int CheckDecays(const Process_Info &pi,
			   const std::vector<Leg> &legs,
			   int d0,int d1,int d2,int d3) {
      if (pi.m_decids.size()==0) {
	if (legs[d0].m_fl!=-legs[d2].m_fl &&
	    legs[d0].m_fl!=-legs[d3].m_fl &&
	    legs[d1].m_fl!=-legs[d2].m_fl &&
	    legs[d1].m_fl!=-legs[d3].m_fl) return 1;
	return 0;
      }
      if (pi.m_decids.size()!=2) return 0;
      int decids[2]{((1<<legs[d0].m_id)|(1<<legs[d1].m_id)),
	  ((1<<legs[d2].m_id)|(1<<legs[d3].m_id))};
      int decfls[2]{legs[(legs[d0].m_fl%2)?d0:d1].m_fl>0?-24:24,
	  legs[(legs[d2].m_fl%2)?d2:d3].m_fl>0?-24:24};
      for (int i(0);i<2;++i)
	if (pi.m_decids[0]==decids[i] && pi.m_decids[1]==decids[1-i] &&
	    pi.m_decfls[0]==decfls[i] && pi.m_decfls[1]==decfls[1-i]) return 1;
      return 0;
    }

    static bool InitializeProcess(CXX_Interface *const interface,
				  const Process_Info &pi,
				  const std::vector<Leg> &legs) {
      if (pi.m_oqcd!=1 || pi.m_oew!=4) return false;
      if (legs.size()!=6) return false;
      for (size_t i(0);i<legs.size();++i)
	if (legs[i].Mass()) return false;
      if (legs[0].m_fl>10 && legs[0].m_fl<17 &&
	  legs[2].m_fl==-legs[0].m_fl+1 &&
          legs[1].m_fl>10 && legs[1].m_fl<17 &&
	  legs[3].m_fl==-legs[1].m_fl-1 &&
	  legs[4].m_fl>0 && legs[4].m_fl<6 && legs[5].m_fl==-legs[4].m_fl &&
	  CheckDecays(pi,legs,0,2,1,3))
	return interface->AddProcess(pi,new qqb_ww(legs,1))>-1;
      if (legs[0].m_fl>10 && legs[0].m_fl<17 &&
	  legs[2].m_fl==-legs[0].m_fl-1 &&
          legs[1].m_fl>10 && legs[1].m_fl<17 &&
	  legs[3].m_fl==-legs[1].m_fl+1 &&
	  legs[4].m_fl>0 && legs[4].m_fl<6 && legs[5].m_fl==-legs[4].m_fl &&
	  CheckDecays(pi,legs,0,2,1,3))
	return interface->AddProcess(pi,new qqb_ww(legs,2))>-1;
      if (legs[0].m_fl>10 && legs[0].m_fl<17 &&
	  legs[3].m_fl==-legs[0].m_fl-1 &&
          legs[1].m_fl>10 && legs[1].m_fl<17 &&
	  legs[2].m_fl==-legs[1].m_fl+1 &&
	  legs[4].m_fl>0 && legs[4].m_fl<6 && legs[5].m_fl==-legs[4].m_fl &&
	  CheckDecays(pi,legs,0,3,1,2))
	return interface->AddProcess(pi,new qqb_ww(legs,2))>-1;
      return false;
    }

    qqb_ww(const std::vector<Leg> &legs,int type):
      Process(legs,5,4), m_type(type) {
      static int first(true);
      if (first) {
	first=false;
	nproc_.nproc=61;
	blha_.useblha=1;
	chooser_();
      }
      m_res.resize(4);
    }

    void Calc(const std::vector<FourVec> &p,int oqcd) {
      // Convert momenta.
      SetMom(p_p,p,m_legs[5],0);// IS q
      SetMom(p_p,p,m_legs[4],1);// IS qb'
      if (m_type==1) {
	SetMom(p_p,p,m_legs[0],2);// FS ve
	SetMom(p_p,p,m_legs[2],3);// FS mu-
	SetMom(p_p,p,m_legs[1],4);// FS e+
	SetMom(p_p,p,m_legs[3],5);// FS vmub
      }
      else {
	SetMom(p_p,p,m_legs[0],4);// FS e-
	SetMom(p_p,p,m_legs[2],5);// FS veb
	SetMom(p_p,p,m_legs[1],2);// FS vmu
	SetMom(p_p,p,m_legs[3],3);// FS mu+
      }
      // Calculate result.
      epinv_.epinv=epinv2_.epinv2=0.0;
      qqb_ww_v_(p_p,p_msqv);
      double res = p_msqv[MSQId(m_legs[5],m_legs[4])];
      m_res[0] = res*m_cfac;
      if (!m_polecheck) return;
      epinv_.epinv=1.0;
      qqb_ww_v_(p_p,p_msqv);
      double res1 = p_msqv[MSQId(m_legs[5],m_legs[4])];
      epinv2_.epinv2=1.0;
      qqb_ww_v_(p_p,p_msqv);
      double res2 = p_msqv[MSQId(m_legs[5],m_legs[4])];
      m_res[1] = (res1-res)*m_cfac;
      m_res[2] = (res2-res1)*m_cfac;
      m_res[3] = m_res[2]/(-qcdcouple_.ason2pi*m_bfac);
    }

    int GetScheme() const { return 1; }

  };// end of class qqb_ww

}// end of namespace MCFM
