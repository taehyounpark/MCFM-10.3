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
  void qqb_zzg_v_(double *p,double *msqv);
  void qqb_zzg_(double *p,double *msqv);}

namespace MCFM {

  class qqb_zzg: public Process {
  private:

    int m_type, m_dec;
    double m_ccfac;

  public:

    static int CheckDecays(const Process_Info &pi,
			   const std::vector<Leg> &legs,
			   int d0,int d1,int d2,int d3) {
      return 1; // JC hack!
      if (pi.m_decids.size()==0) {
	if (legs[d0].m_fl==-legs[d3].m_fl-1 &&
	    legs[d2].m_fl==-legs[d1].m_fl+1) return 0;
	return 1;
      }
      if (pi.m_decids.size()!=2) return 0;
      int decids[2]{((1<<legs[d0].m_id)|(1<<legs[d1].m_id)),
	  ((1<<legs[d2].m_id)|(1<<legs[d3].m_id))};
      for (int i(0);i<2;++i)
	if (pi.m_decids[0]==decids[i] && pi.m_decids[1]==decids[1-i] &&
	    pi.m_decfls[0]==23 && pi.m_decfls[1]==23) return 1;
      return 0;
    }

    static bool InitializeProcess(CXX_Interface *const interface,
				  const Process_Info &pi,
				  const std::vector<Leg> &legs) {
      if (pi.m_oqcd!=2 || pi.m_oew!=4) return false;
      if (legs.size()!=7) return false;
      for (size_t i(0);i<legs.size();++i)
	if (legs[i].Mass()) return false;
      if ((legs[0].m_fl%2)==0 && (legs[1].m_fl%2)==0) return false;
      if (legs[0].m_fl>10 && legs[0].m_fl<17 &&
	  legs[2].m_fl==-legs[0].m_fl &&
          legs[1].m_fl>10 && legs[1].m_fl<17 &&
	  legs[3].m_fl==-legs[1].m_fl &&
	  legs[5].m_fl>0 && legs[5].m_fl<6 && legs[6].m_fl==-legs[5].m_fl &&
	  legs[4].m_fl==21  &&
          CheckDecays(pi,legs,0,2,1,3))
	return interface->AddProcess
	  (pi,new qqb_zzg(legs,(legs[0].m_fl%2)?0:1,pi.m_decids.size()?1:0))>-1;
      return false;
    }

    qqb_zzg(const std::vector<Leg> &legs,int type,int dec):
      Process(legs,6,5), m_type(type), m_dec(dec),
      m_ccfac(1.) {
      static int first(true);
      if (first) {
	first=false;
	nproc_.nproc=(legs[1-m_type].m_fl%2)?481:482;
	blha_.useblha=1;
	chooser_();
        if (legs[0].m_fl == legs[1].m_fl && !m_dec) {
            interference_.interference=true;
            vsymfact_.vsymfact=0.25;
        }
	if ((legs[0].m_fl==11 && legs[1].m_fl==12
          && legs[2].m_fl==-11 && legs[3].m_fl==-12) ||
	    (legs[0].m_fl==13 && legs[1].m_fl==14
          && legs[2].m_fl==-13 && legs[3].m_fl==-14)) {
      // Set blhatype = 1 for process with ZZ/WW interference (e- e+ nue nue~)
          blha_.blhatype = 1;
        }
        else {
          blha_.blhatype = 0;
        }
      }
      m_res.resize(4);
    }

    void Calc(const std::vector<FourVec> &p,int oqcd) {
      // Convert momenta.
      SetMom(p_p,p,m_legs[6],0);// IS q
      SetMom(p_p,p,m_legs[5],1);// IS qb'
      SetMom(p_p,p,m_legs[4],6);// FS g
      if (m_type) {
	SetMom(p_p,p,m_legs[1],2);// FS e-
	SetMom(p_p,p,m_legs[3],3);// FS e+
	SetMom(p_p,p,m_legs[0],4);// FS mu-/v
	SetMom(p_p,p,m_legs[2],5);// FS mu+/v~
      }
      else {
	SetMom(p_p,p,m_legs[0],2);// FS e-
	SetMom(p_p,p,m_legs[2],3);// FS e+
	SetMom(p_p,p,m_legs[1],4);// FS mu-/v
	SetMom(p_p,p,m_legs[3],5);// FS mu+/v~
      }
      // Set blhatype = 1 for process with ZZ/WW interference (e- e+ nue nue~)
      // blha_.blhatype = 1;
      // Calculate result.
      if (blha_.blhatype > 0) {
         blha_.blhatype=(m_legs[6].m_fl%2)?1:2;
      }
      epinv_.epinv=epinv2_.epinv2=0.0;
      qqb_zzg_v_(p_p,p_msqv);
      double res = p_msqv[MSQId(m_legs[6],m_legs[5])];
      m_res[0] = res*m_cfac/m_ccfac;
      if (!m_polecheck) return;
      epinv_.epinv=1.0;
      qqb_zzg_v_(p_p,p_msqv);
      double res1 = p_msqv[MSQId(m_legs[6],m_legs[5])];
      epinv2_.epinv2=1.0;
      qqb_zzg_v_(p_p,p_msqv);
      double res2 = p_msqv[MSQId(m_legs[6],m_legs[5])];
      m_res[1] = (res1-res)*m_cfac/m_ccfac;
      m_res[2] = (res2-res1)*m_cfac/m_ccfac;
      m_res[3] = m_res[2]/(-qcdcouple_.ason2pi*m_bfac);
    }

    int GetScheme() const { return 0; }

  };// end of class qqb_zzg

}// end of namespace MCFM
