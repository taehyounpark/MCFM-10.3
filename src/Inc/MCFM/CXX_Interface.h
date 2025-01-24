// Copyright (C) 2021 John Campbell, Stefan Hoeche, Christian T Preuss.

#ifndef MCFM__CXX_Interface_H
#define MCFM__CXX_Interface_H

// MCFM includes.
#include "MCFM/CXX_Wrapper.h"
#include "MCFM/Flavor_Map.h"

// Standard includes. 
#include <cmath>
#include <map>
#include <sstream>

// MCFM interface and individual process interfaces.
namespace MCFM {

  template <typename __Tp>
  std::ostream &operator<<(std::ostream &str,const std::vector<__Tp> &v);

  std::vector<int> ID(size_t id);
  size_t ID(const std::vector<int> &id);

  // Process parameters.
  struct Process_Info {
    std::vector<int> m_ids, m_decids, m_decfls;
    int m_nin, m_oqcd, m_oew;
    std::string m_model;
    inline Process_Info(const std::vector<int> &ids=std::vector<int>(),
			const int nin=0,const int oqcd=0,const int oew=0):
      m_ids(ids), m_nin(nin), m_oqcd(oqcd), m_oew(oew), m_model("sm") {}
    inline int UID() const { return m_oqcd+100*m_oew; }
  };// end of struct Process_Info

  // Forward declarations.
  class CXX_Interface;

  // Base class for interfacing MCFM processes.
  class Process {
  protected:

    double *p_p, *p_msqv;

    std::vector<Leg> m_legs;
    std::vector<double> m_res;

    double m_sf, m_cfac, m_bfac;
    int m_polecheck, m_is[2];

    inline void SetMom(double *m,const std::vector<FourVec> &p,
		       const Leg &l,const int n) const
    { GetMom(m,n,l.m_is?-p[l.m_id]:p[l.m_id]); }
    inline int MSQId(const Leg &l1,const Leg &l2) const
    { return mr(-MCFMId(l1.m_fl),-MCFMId(l2.m_fl)); }

    double ISSymmetryFactor(const std::vector<Leg> &legs,
			    const int mode=0) const;
    double FSSymmetryFactor(const std::vector<Leg> &legs,
			    const int mode=0) const;

  public:

    Process(const std::vector<Leg> &legs,int is1,int is2);

    ~Process();

    double GetSchemeFactor(int scheme,int mode=1);

    virtual void Calc(const std::vector<FourVec> &pn,int oqcd) = 0;

    // Regularization scheme: 0 = tHV/CDR, 1 = DR
    virtual int GetScheme() const = 0;

    // Process specific citations.
    virtual std::string GetReferences() const;

    // Symmetry factor, including colour and helicity average.
    inline double GetSymmetryFactor() { return m_sf; }

    inline const std::vector<double> &GetResult() const { return m_res; }

    inline void SetPoleCheck(int check) { m_polecheck=check; }

  };// end of class Process

  // The MCFM interface.
  class CXX_Interface {
  private:

    // List of processes.
    std::vector<Process*> m_processes;
    std::map<int,std::map<std::vector<int>,int> > m_procmap;
    // List of references.
    std::vector<std::string> m_ref_list, m_cite_list;
    // MCFM version number.
    std::string m_vno;
    int m_banner;
    // Debug level.
    static int s_verbose;

    // Get parameter from list.
    inline std::string ReadParam(std::string key,
      std::map<std::string,std::string>& params);

    // Convert strings.
    template <class _T> _T ToType
    (const std::string &v,const size_t prec=16) {
      std::stringstream c; _T r;
      c.precision(prec); c<<v; c>>r;
      return r;
    }

  public:

    CXX_Interface(const int banner=1);
    ~CXX_Interface();

    std::string GetReferences() const;
    std::string GetStartupMessage() const;
    std::string GetFinishMessage(const int aggregate=0) const;

    bool Initialize(std::map<std::string,std::string>& parameters);
    bool CheckInput(std::map<std::string,std::string>& parameters);
    void PrintSettings();

    // Convert list of PDG IDs to MCFM process ID.
    int InitializeProcess(const Process_Info &pi);
    inline int AddProcess(const Process_Info &pi,Process *const proc) {
      m_processes.push_back(proc);
      return m_procmap[pi.UID()][pi.m_ids]=m_processes.size()-1; }
    inline int GetProcessID(const Process_Info &pi) const
    { std::map<int,std::map<std::vector<int>,int> >::
	const_iterator iit(m_procmap.find(pi.UID()));
      if (iit==m_procmap.end()) return -1;
      std::map<std::vector<int>,int>::
	const_iterator pit(iit->second.find(pi.m_ids));
      return pit!=iit->second.end()?pit->second:-1; }
    inline Process *GetProcess(const size_t pid,const int lo=1) const
    { return pid<m_processes.size()?m_processes[pid]:NULL; }

    // Set pole check mode.
    void SetPoleCheck(int check)
    { for (int i(0);i<m_processes.size();++i)
        if(m_processes[i]) m_processes[i]->SetPoleCheck(check); }

    // Set scale and alphaS value.
    void SetMuR2(const double &mur2)
    { mcfmscale_.scale=sqrt(mcfmscale_.musq=mur2); }
    void SetAlphaS(const double &as)
    { qcdcouple_.as=as; qcdcouple_.ason2pi=as/(2.*M_PI);
      qcdcouple_.ason4pi=qcdcouple_.ason2pi/2.;
      qcdcouple_.gsq=4.*M_PI*as; }

    // Get alphaS value.
    double GetAlphaS() { return qcdcouple_.as; }

    // Call MCFM matrix element routines.
    inline void Calc(int procID,const std::vector<FourVec> &p,int oqcd)
    { m_processes[procID]->Calc(p,oqcd); }
    inline void Calc(const Process_Info &pi,
		     const std::vector<FourVec> &p,int oqcd)
    { return Calc(GetProcessID(pi),p,oqcd); }

    // Access results.
    inline const std::vector<double> &GetResult(int procID) const
    { return m_processes[procID]->GetResult(); }
    inline const std::vector<double> &GetResult
    (const Process_Info &pi) const
    { return GetResult(GetProcessID(pi)); }

    // Get pointer to process interface.
    inline Process *GetProcess(int procID) const
    { return m_processes[procID]; }
    inline Process *GetProcess(const Process_Info &pi) const
    { return GetProcess(GetProcessID(pi)); }

    inline static void SetVerbose(int v)
    { verbose_.verbose=s_verbose=v; }

    inline static int Verbose() { return s_verbose; }

  };// end of class CXX_Interface

} // end of namespace MCFM

#endif
