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
 
// Program to test MCFM interface against MadLoop.

// Standard includes.
#include <ctime>
#include <fstream>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

// MCFM includes
#include "MCFM/CXX_Interface.h"
#include "params.cxx"
#include "Rambo.cxx"

#include <sys/stat.h>
#include <dlfcn.h>
#include <unistd.h>

void set_sm_parameters_ml(std::map<std::string,std::string> &params,
			  std::string cn, std::string pn);

using namespace MCFM;

inline int mpml(const int id,const int i)
{ return i+id*4; }

inline void GetMomML(double *m,const int n,const FourVec &p)
{ m[mpml(n,0)]=p[0]; for (int i(1);i<4;++i) m[mpml(n,i)]=p[i]; }

extern "C" {

  inline void MakeFortranString
  (char *output,std::string input,unsigned int length)
  {
    for (unsigned int i=0;i<length;++i) 
      output[i]=(char)32;
    for (size_t j=0;j<input.length();++j) 
      output[j]=(char)input[j];
  }

  void setpara2_(char *name);
  void setmadlooppath_(char *name);
  void setparamlog_(int *log);

}

typedef void (*ME_Function)(double *,double *,double *,double *,int *);
typedef void (*DD_Function)(double *,double *);
typedef void (*C_Function)(char *);
typedef void (*I_Function)(int *);
typedef void (*Void_Function)(void);

template <class _T> std::string ToString
(const _T &v,const size_t prec=16) {
  std::stringstream c; std::string r;
  c.precision(prec); c<<v; c>>r;
  return r;
}

class Library_Loader {
private:
  std::vector<std::string>    m_paths;
  std::map<std::string,void*> m_libs;
public:
  Library_Loader(): m_paths{"./"} {}
  inline void AddPath(const std::string &path) { m_paths.push_back(path); }
  void *LoadLibrary(const std::string &name) {
    std::map<std::string,void*>::iterator lit(m_libs.find(name));
    if (lit!=m_libs.end()) return lit->second;
    // std::cout<<"loading library 'lib"<<name<<LIB_SUFFIX<<"' {"<<std::endl;
    for (size_t i(0);i<m_paths.size();++i) {
      std::string libname(m_paths[i]+"/lib"+name+LIB_SUFFIX);
      struct stat buffer;
      // std::cout<<"checking for '"<<libname<<"' ... "<<std::flush;
      if (stat(libname.c_str(),&buffer)) {
	// std::cout<<" not found"<<std::endl;
	continue;
      }
      // std::cout<<" found"<<std::endl;
      void *module(dlopen(libname.c_str(),RTLD_LAZY|RTLD_GLOBAL));
      if (module!=nullptr) {
	// std::cout<<"} found in '"<<m_paths[i]<<"'"<<std::endl;
	m_libs[name]=module;
	return module;
      }
      char *err(dlerror());
      if (err!=nullptr) {
	std::cout<<err<<std::endl;
	break;
      }
    }
    // std::cout<<"} failed"<<std::endl;
    std::cout<<"Failed to load 'lib"<<name<<LIB_SUFFIX<<"'."<<std::endl;
    return nullptr;
  }
  void *GetLibraryFunction(const std::string &libname,
			   const std::string &funcname,
			   void *&module) {
    // std::cout<<"executing library function '"<<funcname
    // 		   <<"' from 'lib"<<libname<<LIB_SUFFIX<<"' ... "<<std::flush;
    if (module==nullptr) module=LoadLibrary(libname);
    if (module==nullptr) return nullptr;
    void *func(dlsym(module,funcname.c_str()));
    char *error(dlerror());
    if (error!=nullptr) {
      // std::cout<<"failed"<<std::endl;
      std::cout<<error<<std::endl;
      std::cout<<"Failed to load function '"<<funcname<<"'."<<std::endl;
      return nullptr;
    }
    // std::cout<<"done"<<std::endl;
    return func;
  }
};//end of class Library_Loader

int convert_pid(const char* in) {
  if (strcmp("d~",in)==0) return -1;
  if (strcmp("d",in)==0) return 1;
  if (strcmp("u~",in)==0) return -2;
  if (strcmp("u",in)==0) return 2;
  if (strcmp("s~",in)==0) return -3;
  if (strcmp("s",in)==0) return 3;
  if (strcmp("c~",in)==0) return -4;
  if (strcmp("c",in)==0) return 4;
  if (strcmp("b~",in)==0) return -5;
  if (strcmp("b",in)==0) return 5;
  if (strcmp("t~",in)==0) return -6;
  if (strcmp("t",in)==0) return 6;
  if (strcmp("e+",in)==0) return -11;
  if (strcmp("e-",in)==0) return 11;
  if (strcmp("mu+",in)==0) return -13;
  if (strcmp("mu-",in)==0) return 13;
  if (strcmp("tau+",in)==0) return -15;
  if (strcmp("tau-",in)==0) return 15;
  if (strcmp("ve~",in)==0) return -12;
  if (strcmp("ve",in)==0) return 12;
  if (strcmp("vmu~",in)==0) return -14;
  if (strcmp("vmu",in)==0) return 14;
  if (strcmp("vtau~",in)==0) return -16;
  if (strcmp("vtau",in)==0) return 16;
  if (strcmp("h",in)==0) return 25;
  if (strcmp("g",in)==0) return 21;
  if (strcmp("a",in)==0) return 22;
  if (strcmp("Z",in)==0) return 23;
  if (strcmp("W+",in)==0) return 24;
  if (strcmp("W-",in)==0) return -24;
  return 0;
}

std::string convert_pid_to_ml(int pid,const std::string &model) {
  switch (pid) {
  case -1:  return "d~";
  case 1:   return "d";
  case -2:  return "u~";
  case 2:   return "u";
  case -3:  return "s~";
  case 3:   return "s";
  case -4:  return "c~";
  case 4:   return "c";
  case -5:  return "b~";
  case 5:   return "b";
  case -6:  return "t~";
  case 6:   return "t";
  case -11: return "e+";
  case 11:  return "e-";
  case -13: return "mu+";
  case 13:  return "mu-";
  case -15: return "ta+";
  case 15:  return "ta-";
  case -12: return "ve~";
  case 12:  return "ve";
  case -14: return "vm~";
  case 14:  return "vm";
  case -16: return "vt~";
  case 16:  return "vt";
  case 25:  return model=="heft"?"x0":"h";
  case 21:  return "g";
  case 22:  return "a";
  case 23:  return "z";
  case 24:  return "w+";
  case -24: return "w-";
  default: break;
  }
  return "X";
}

// Main program.
int main(int argc, char** argv) {

  double Ecm = 350., mcfm_nptsfac=1.;
  int foEW = -1, foQCD = -1, diagonalCKM = 1;
  int nPts = 10, checkPoles = 0;
  int printMoms = 0, ampType = 11, fampType = -1;
  std::string filename;
  std::vector<int> decids, decfls;

  // Collect settings.
  std::map<std::string,std::string> params;
  read_parameters("params.lh",params);

  // Optional arguments.
  int opt, verbosity = 0;
  while ((opt = getopt(argc, argv, "pvmA:n:e:o:S:W:CD:P:h")) != -1) {
    switch (opt) {
    case 'p': checkPoles = 1; break;
    case 'v': verbosity += 1; break;
    case 'm': printMoms = 1; break;
    case 'n': nPts = atoi(optarg); break;
    case 'e': mcfm_nptsfac = atoi(optarg); break;
    case 'o': filename = optarg; break;
    case 'A': fampType = atoi(optarg); break;
    case 'W': foEW = atoi(optarg); break;
    case 'S': foQCD = atoi(optarg); break;
    case 'C': diagonalCKM = 0; break;
    case 'P': {
      size_t pos(0);
      std::string arg(optarg);
      for (;pos<arg.length();++pos) if (arg[pos]=='=') break;
      params[arg.substr(0,pos)]=arg.substr(pos+1,std::string::npos);
      break;
    }
    case 'D': {
      size_t pos(0);
      std::string arg(optarg);
      for (;pos<arg.length();++pos)
	if (arg[pos]==',') break;
      decids.push_back(atoi(arg.substr(0,pos).c_str()));
      decfls.push_back(convert_pid(arg.substr(pos+1,std::string::npos).c_str()));
      break;
    }
    case 'h':
      std::cout<<argv[0]<<": Command line options\n"
	       <<"  -h        -- print this help message\n"
	       <<"  -v        -- increase verbosity\n"
	       <<"  -p        -- print ME values and enable pole check\n"
	       <<"  -m        -- print momentum configuration\n"
	       <<"  -C        -- enable non-diagonal CKM\n"
	       <<"  -n points -- set number of phase space points to 'points'\n"
	       <<"  -e factor -- set factor for MCFM extra points to 'points'\n"
	       <<"  -S oQCD   -- set QCD order to 'oQCD'\n"
	       <<"  -W oEW    -- set EW order to 'oEW'\n"
	       <<"  -A type   -- set OpenLoops amplitude type to 'type'\n"
	       <<"  -P p=val  -- set parameter 'p' to 'val'\n"
	       <<"  -D id,fl  -- add decay id and flavor 'id' and 'fl'\n";
      exit(0);							
    default:
      std::cerr<<"Usage: "<<argv[0]
	       <<" [-vpmC] [-n #pts] [-S oQCD] [-W oEW] [-A type] [-P p,val] [-D id,fl] [process]"
	       <<std::endl;
      exit(EXIT_FAILURE);
    }
  }
  if (optind >= argc) {
    std::cerr<<"You did not specify a process"<<std::endl;
    exit(EXIT_FAILURE);
  }

  //---------------------------------------------------------------------------
  // Process and initialization.

  // Process.
  std::vector<int> ids;
  int oEW = 0, oQCD = 1, nh = 0; // Needed for OpenLoops.
  for (int i(optind);i<argc;++i) {
    ids.push_back(convert_pid(argv[i]));
    if (ids.back()==21 || (abs(ids.back()))>0 && abs(ids.back())<7) ++oQCD;
    else {
      if (ids.back()==25) ++nh;
      ++oEW;
    }
  }
  // Note: in general, this does not capture all cases.
  oQCD-=2;
  oQCD = std::max(oQCD,0);

  int nIn = 2, nOut = ids.size()-nIn;

  // Fix coupling orders and amplitude types for special cases.
  if (nh && oEW == nh) {
    params["H_width"] = "0";
    if (params["model"] != "heft") {
      ampType = 12;
      oQCD -= 1;
    }
    oQCD += 2;
  }
  if (foEW>=0) oEW=foEW;
  if (foQCD>=0) oQCD=foQCD;
  if (fampType>0) ampType=fampType;

  //---------------------------------------------------------------------------
  // MCFM.

  // Create instance of MCFM interface.
  CXX_Interface mcfm;

  // Initialize MCFM.
  mcfm.Initialize(params);
  mcfm.SetVerbose(verbosity);

  // Print settings.
  if (verbosity) mcfm.PrintSettings();

  // Initialize process.
  Process_Info pi(ids,2,oQCD,oEW);
  pi.m_decids=decids;
  pi.m_decfls=decfls;
  pi.m_model=params["model"];
  int pid_mcfm(mcfm.InitializeProcess(pi));
  if (pid_mcfm < 0) {
    std::cout << " Process not available in MCFM." << std::endl;
    exit(EXIT_FAILURE);
  }

  // Enable calculation of pole tems.
  mcfm.SetPoleCheck(checkPoles);

  //---------------------------------------------------------------------------
  // MadLoop.
  std::string procstring, filestring;
  for (size_t i(optind);i<optind+nIn;++i) {
    procstring+=convert_pid_to_ml(convert_pid(argv[i]),params["model"])+" ";
    filestring+=ToString(convert_pid(argv[i]))+"_";
  }
  procstring+="> ";
  for (size_t i(optind+nIn);i<argc;++i) {
    procstring+=convert_pid_to_ml(convert_pid(argv[i]),params["model"])+" ";
    filestring+=ToString(convert_pid(argv[i]))+"_";
  }
  filestring.erase(filestring.length()-1,1);
  if (params["model"]!="sm") {
    filestring+="_"+params["model"]; 
    procstring+="/ t ";
  }

  // Set MadLoop parameters.
  std::string mlPath("./ML5/"+filestring), cardName("OLP.mg5");
  std::string paramName(mlPath+"/SubProcesses/MadLoop5_resources/param_card.dat");
  set_sm_parameters_ml(params,cardName,paramName);

  // Initialize process.
  int pid_ml(0);
  std::ofstream card(cardName,std::ios::app);
  card<<"add process "<<procstring
      <<"QCD="<<oQCD-1<<" QED="<<oEW
      <<" [virt=QCD] QCD="<<oQCD<<" QED="<<oEW
      <<" @"<<pid_ml<<std::endl;
  card<<"output "<<mlPath<<std::endl;

  Library_Loader loader;
  loader.AddPath(mlPath+"/lib");
  loader.AddPath(mlPath+"/lib/collier_lib");
  loader.AddPath(mlPath+"/lib/ninja_lib");
  loader.LoadLibrary("collier");
  loader.LoadLibrary("ninja");
  loader.LoadLibrary("cts");
  loader.LoadLibrary("iregi");
  loader.LoadLibrary("ML5_DHELAS");
  void *gmodule(loader.LoadLibrary("ML5_MODEL"));
  if (gmodule==nullptr) exit(1);
  std::string bp(mlPath+"/SubProcesses/MadLoop5_resources");
  std::string libname("ML5_P"+ToString(pid_ml));
  void *module(loader.LoadLibrary(libname));
  if (module==nullptr) exit(1);

  char name[512];
  int log(0), onoff(1), info;
  MakeFortranString(name,bp+"/param_card.dat",512);
  ((I_Function)loader.GetLibraryFunction("","setparamlog_",gmodule))(&log);
  ((C_Function)loader.GetLibraryFunction("","setpara2_",gmodule))(name);
  MakeFortranString(name,bp,512);
  ((C_Function)loader.GetLibraryFunction("","setmadlooppath_",module))(name);
  ((I_Function)loader.GetLibraryFunction
   ("","ml5_"+ToString(pid_ml)+"_force_stability_check_",module))(&onoff);
  double mur(std::stod(params["scale"])), as(std::stod(params["alpha_S"]));
  ((DD_Function)loader.GetLibraryFunction("","update_as_param2_",gmodule))(&mur,&as);
  // ((V_Function)loader.GetLibraryFunction("","printout_",gmodule))();
  ((I_Function)loader.GetLibraryFunction
   ("","ml5_"+ToString(pid_ml)+"_get_nsqso_loop_",module))(&info);
  ME_Function me = (ME_Function)loader.GetLibraryFunction
    ("","ml5_"+ToString(pid_ml)+"_sloopmatrix_thres_",module);

  double *res_ml = new double[4*(1+info)];
  double *prec_ml = new double[1+info];

  // Compare CPU time.
  clock_t t;
  double  cput_mcfm;
  double  cput_ml;
  double  cput_ps;

  std::cout << "\n Generating " << nPts << " flat phase space points\n\n";

  double *pn_ml = new double[4*(nIn+nOut)];
  if (checkPoles) {
    for (int iPt(0);iPt<nPts;++iPt) {
      // Generate phase space point with OpenLoops.
      std::vector<FourVec> pn;
      GenPoint(Ecm, pn, ids, nIn, params);
      for (size_t n(0);n<pn.size();++n) GetMomML(pn_ml,n,pn[n]);
      
      // Calculate result with MadLoop.
      t = clock();
      int retcode;
      double prec(-1.0);
      me(pn_ml,res_ml,&prec,prec_ml,&retcode);
      if (retcode/100==4) {
	std::cout<<"Unstable point {\n";
	std::cout.precision(16);
	for (size_t i(0);i<pn.size();++i)
	  std::cout<<"  pn["<<i<<"]=FourVec"<<pn[i]<<";\n";
	std::cout<<"}"<<std::endl;
      }
      t = clock() - t;
      cput_ml += double(t)/CLOCKS_PER_SEC;
    
      // Calculate result with MCFM.
      mcfm.SetAlphaS(std::stod(params["alpha_S"]));
      t = clock();
      mcfm.Calc(pid_mcfm, pn, 1);
      t = clock() - t;
      cput_mcfm += double(t)/CLOCKS_PER_SEC;
      std::vector<double> res_mcfm = mcfm.GetResult(pid_mcfm);
      // Translate MCFM result to same normalization and scheme conventions as OpenLoops.
      for (int i(0);i<4;++i)
	res_mcfm[i] /= mcfm.GetProcess(pid_mcfm)->GetSymmetryFactor();
      res_mcfm[0] -= res_mcfm[3]*mcfm.GetProcess(pid_mcfm)->GetSchemeFactor(0);

      //  Print result and optionally phase space point.
      if (ampType == 11) {
        std::cout << " Finite: MCFM = " << std::setprecision(17) << res_mcfm[0] << " MadLoop = " << res_ml[1]
                  << " Ratio = " << res_mcfm[0]/res_ml[1] << std::endl
		  << " IR:     MCFM = " << res_mcfm[1] << " MadLoop = " << res_ml[2]
		  << " Ratio = " << res_mcfm[1]/res_ml[2] << std::endl
		  << " IR2:    MCFM = " << res_mcfm[2] << " MadLoop = " << res_ml[3]
		  << " Ratio = " << res_mcfm[2]/res_ml[3] << std::endl;
      }
      double born_ml = ampType==11?res_ml[0]:res_ml[1];
      std::cout << " Born:   MCFM = " << std::setprecision(17) << res_mcfm[3] << " MadLoop = " << born_ml
                << " Ratio = " << res_mcfm[3]/born_ml
                << std::endl << std::endl;
      if (printMoms) {
	for (size_t i(0);i<pn.size();++i)
	  std::cout << " p[" << i << "] = " << std::setprecision(17)
		    << pn[i] << " ( p^2 = " << m2(pn[i]) << " )\n";
	std::cout << std::endl;
      }
    }
  } else {
    // Phase space.
    t = clock();
    for (int iPt(0);iPt<nPts;++iPt) {
      // Generate phase space point with Rambo.
      std::vector<FourVec> pn;
      GenPoint(Ecm, pn, ids, nIn, params);
      for (size_t n(0);n<pn.size();++n) GetMomML(pn_ml,n,pn[n]);
    }
    t = clock() - t;
    cput_ps = double(t)/CLOCKS_PER_SEC;
    
    // MadLoop.
    t = clock();
    for (int iPt(0);iPt<nPts;++iPt) {
      // Generate phase space point with Rambo.
      std::vector<FourVec> pn;
      GenPoint(Ecm, pn, ids, nIn, params);
      for (size_t n(0);n<pn.size();++n) GetMomML(pn_ml,n,pn[n]);
      
      // Calculate result with MadLoop.
      int retcode;
      double prec(-1.0);
      me(pn_ml,res_ml,&prec,prec_ml,&retcode);
      if (retcode/100==4) {
	std::cout<<"Unstable point {\n";
	std::cout.precision(16);
	for (size_t i(0);i<pn.size();++i)
	  std::cout<<"  pn["<<i<<"]=FourVec"<<pn[i]<<";\n";
	std::cout<<"}"<<std::endl;
      }
    }
    cput_ml=(double(clock()-t)-cput_ps)/CLOCKS_PER_SEC;
    
    // MCFM.
    t = clock();
    for (int iPt(0);iPt<nPts;++iPt) {
      // Generate phase space point with Rambo.
      std::vector<FourVec> pn;
      double pn_rcl[nIn+nOut][4];
      GenPoint(Ecm, pn, ids, nIn, params);
      for (size_t n(0);n<pn.size();++n) GetMomML(pn_ml,n,pn[n]);

      // Calculate result with MCFM.
      for (int i(0);i<mcfm_nptsfac;++i) {
	mcfm.SetAlphaS(std::stod(params["alpha_S"]));
	mcfm.Calc(pid_mcfm, pn, 1);
      }
      std::vector<double> res_mcfm = mcfm.GetResult(pid_mcfm);
      // Translate MCFM result to same normalization and scheme conventions as OpenLoops.
      for (int i(0);i<4;++i)
	res_mcfm[i] /= mcfm.GetProcess(pid_mcfm)->GetSymmetryFactor();
      res_mcfm[0] -= res_mcfm[3]*mcfm.GetProcess(pid_mcfm)->GetSchemeFactor(0);
    }
    cput_mcfm=(double(clock()-t)-cput_ps)/CLOCKS_PER_SEC/mcfm_nptsfac;
  }

  // Print average CPU time.
  std::cout << " ----------------------------------------------------\n"
	    << "  Average CPU time per event:\n"
	    << "    MCFM:      " << std::setprecision(17) << std::scientific << cput_mcfm/double(nPts) << "s\n"
	    << "    MadLoop:   " << std::setprecision(17) << std::scientific << cput_ml/double(nPts) << "s\n"
	    << "    Ratio:     " << cput_ml/cput_mcfm << "\n"
	    << " ----------------------------------------------------\n";

  // Done.
  return 0;
}

void set_sm_parameters_ml(std::map<std::string,std::string> &params,
			  std::string cn, std::string pn) {
  std::string model = std::stod(params["bottom_mass"])>0.?"loop_sm":"loop_sm-no_b_mass";
  if (params["model"]=="heft") model="HC_NLO_X0_UFO-heft";
  std::string cardname(cn);
  int res(std::system(("rm -f "+cardname).c_str()));
  int nofprocs(0);
  std::ofstream card(cardname,std::ios::app);
  if (model.substr(0,model.find('-'))=="loop_sm")
    card<<"set complex_mass_scheme True\n";
  card<<"import model "<<model<<"\n";
  //SubProcesses/MadLoop5_resources
  std::string fname(pn);
  res=std::system(("rm -f "+fname).c_str());
  std::ofstream param(fname,std::ios::out);
  param.precision(12);
  param<<"Block LOOP\n";
  param<<"    1 "<<params["scale"]<<"\n";
  param<<"Block SMINPUTS\n";
  param<<"    1 "<<1.0/std::stod(params["alpha_EM"])<<"\n";
  param<<"    3 "<<params["alpha_S"]<<"\n";
  param<<"Block MASS\n";
  param<<"    1 "<<params["down_mass"]<<"\n";
  param<<"    2 "<<params["up_mass"]<<"\n";
  param<<"    3 "<<params["strange_mass"]<<"\n";
  param<<"    4 "<<params["charm_mass"]<<"\n";
  param<<"    5 "<<params["bottom_mass"]<<"\n";
  param<<"    6 "<<params["top_mass"]<<"\n";
  param<<"    11 "<<params["electron_mass"]<<"\n";
  param<<"    12 "<<0.0<<"\n";
  param<<"    13 "<<params["muon_mass"]<<"\n";
  param<<"    14 "<<0.0<<"\n";
  param<<"    15 "<<params["tau_mass"]<<"\n";
  param<<"    16 "<<0.0<<"\n";
  param<<"    21 "<<0.0<<"\n";
  param<<"    22 "<<0.0<<"\n";
  param<<"    23 "<<params["Z_mass"]<<"\n";
  param<<"    24 "<<params["W_mass"]<<"\n";
  param<<"    25 "<<params["H_mass"]<<"\n";
  param<<"Block YUKAWA\n";
  param<<"    4 "<<params["charm_yukawa"]<<"\n";
  param<<"    5 "<<params["bottom_yukawa"]<<"\n";
  param<<"    6 "<<params["top_yukawa"]<<"\n";
  param<<"    15 "<<params["tau_yukawa"]<<"\n";
  param<<"DECAY  6 "<<params["top_width"]<<"\n";
  param<<"DECAY 23 "<<params["Z_width"]<<"\n";
  param<<"DECAY 24 "<<params["W_width"]<<"\n";
  param<<"DECAY 25 "<<params["H_width"]<<"\n";
  param<<"Block QNUMBERS 82\n";
  param<<"    1 0\n";
  param<<"    2 1\n";
  param<<"    3 8\n";
  param<<"    4 1\n";
}
