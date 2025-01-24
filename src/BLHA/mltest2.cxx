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
 
// Program to test MCFM interface against OpenLoops.

// Standard includes.
#include <ctime>
#include <fstream>
#include <algorithm>
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


using namespace MCFM;

// Main program.
int main(int argc, char** argv) {

  // Create instance of MCFM interface.
  CXX_Interface mcfm;

  // Read settings.
  std::map<std::string,std::string> params;
  read_parameters("params.lh",params);

  // Command-line options.
  int foEW = -1, foQCD = -1, diagonalCKM = 1, nIn = 2;
  int nPts = 1000000, checkPoles = 1;
  int printMom = 0, mode = 0, ampType = 11, fampType = -1;
  std::vector<int> decids, decfls;
  int opt, verbosity = 0;
  while ((opt = getopt(argc, argv, "vM:m:A:n:i:o:S:W:C:P:h")) != -1) {
    switch (opt) {
    case 'v': verbosity += 1; break;
    case 'M': printMom = atoi(optarg); break;
    case 'm': mode = atoi(optarg); break;
    case 'n': nPts = atoi(optarg); break;
    case 'i': nIn = atoi(optarg); break;
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
    case 'h':
      std::cout<<argv[0]<<": Command line options\n"
	       <<"  -h        -- print this help message\n"
	       <<"  -v        -- increase verbosity\n"
	       <<"  -C        -- enable non-diagonal CKM\n"
	       <<"  -m mode   -- set stability check mode to 'mode'\n"
	       <<"  -n nin    -- set number of incoming particles to 'nin'\n"
	       <<"  -n points -- set number of phase space points to 'points'\n"
	       <<"  -S oQCD   -- set QCD order to 'oQCD'\n"
	       <<"  -W oEW    -- set EW order to 'oEW'\n"
	       <<"  -M print  -- set print mode to 'print'\n"
	       <<"  -P p=val  -- set parameter 'p' to 'val'\n"<<std::endl;
      exit(0);							
    default:
      std::cerr<<"Usage: "<<argv[0]<<" [-vmC] [-n #pts] [-S oQCD] [-W oEW] [-P p,val]"<<std::endl;
      exit(EXIT_FAILURE);
    }
  }
  if (optind >= argc) {
    std::cerr<<" Please provide an input file."<<std::endl;
    exit(EXIT_FAILURE);
  }

  //---------------------------------------------------------------------------
  // Read phase-space points.

  std::string infilename(argv[optind]);
  std::cout << "Input file: " << infilename << std::endl;
  std::ifstream infile(infilename);
  if (!infile.good())
    std::cout << "Input file " << infilename << " not found.\n";
  int nPt(0), nans(0), pid_mcfm(0), pid_rcl(1);
  std::vector<int> ids;
  std::vector<FourVec> pn;
  std::string line;
  std::vector<double> accs;
  double *res_ml, *prec_ml;
  ME_Function me;
  while (std::getline(infile,line)) {
    line = reduce(line);
    line = chop(line,"!#");
    // Read configuration.
    if (line != "" && line != "eof") {
      int id(0);
      double e(0.0), px(0.0), py(0.0), pz(0.0);
      std::stringstream linestream(line);
      linestream.precision(16);
      linestream >> id >> e >> px >> py >> pz;
      ids.push_back(id);
      pn.push_back(FourVec(e,px,py,pz));
    }

    // At end of phase-space point do calculations.
    if (line == "" || line == "eof") {
      if (ids.size() > 0 && pn.size() > 0 && ids.size() == pn.size()) {
	// Count number of points and print this point.
	++nPt;
	if (printMom&2) {
	  std::cout << "Phase-space point " << nPt << std::endl;
	  FourVec psum(0.,0.,0.,0.);
	  for (size_t iPtcl(0); iPtcl<ids.size(); ++iPtcl) {
	    std::cout << ids[iPtcl] << " " << pn[iPtcl] << " "
		      << (pn[iPtcl][0]*pn[iPtcl][0]-pn[iPtcl][1]*pn[iPtcl][1]
			  -pn[iPtcl][2]*pn[iPtcl][2]-pn[iPtcl][3]*pn[iPtcl][3]) << std::endl;
	    psum[0]+=iPtcl<nIn?-pn[iPtcl][0]:pn[iPtcl][0];
	    psum[1]+=iPtcl<nIn?-pn[iPtcl][1]:pn[iPtcl][1];
	    psum[2]+=iPtcl<nIn?-pn[iPtcl][2]:pn[iPtcl][2];
	    psum[3]+=iPtcl<nIn?-pn[iPtcl][3]:pn[iPtcl][3];
	  }
	  std::cout<<"Sum "<<psum<<"\n";
	}
	int nOut = ids.size()-nIn;
	if (nPt==1) {
	  // Set up process parameters.
	  int oEW = 0, oQCD = 1, nh = 0;
	  for (size_t iPtcl(0); iPtcl<ids.size(); ++iPtcl) {
	    if (ids[iPtcl]==21 || (abs(ids[iPtcl])>0 && abs(ids[iPtcl])<7)) ++oQCD;
	    else {
	      if (ids[iPtcl]==25) ++nh;
	      ++oEW;
	    }
	  }
	  oQCD-=2;
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
	  
	  // Initialize MCFM.
	  mcfm.Initialize(params);
	  Process_Info pi(ids,nIn,oQCD,oEW);
	  pi.m_decids=decids;
	  pi.m_decfls=decfls;
	  pi.m_model=params["model"];
	  pid_mcfm=mcfm.InitializeProcess(pi);
	  if (pid_mcfm < 0) {
	    std::cout << " Process not available in MCFM." << std::endl;
	    exit(EXIT_FAILURE);
	  }
	  mcfm.SetPoleCheck(checkPoles);
	
	  //---------------------------------------------------------------------------
	  // MadLoop.
	  std::string procstring, filestring;
	  for (int i(0);i<nIn;++i) {
	    procstring += convert_pid_to_ml(ids[i],params["model"])+" ";
	    filestring+=ToString(ids[i])+"_";
	  }
	  procstring += "> ";
	  for (int i(nIn);i<nIn+nOut;++i) {
	    procstring += " "+convert_pid_to_ml(ids[i],params["model"])+" ";
	    filestring+=ToString(ids[i])+"_";
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
	  me = (ME_Function)loader.GetLibraryFunction
	    ("","ml5_"+ToString(pid_ml)+"_sloopmatrix_thres_",module);
	  res_ml = new double[4*(1+info)];
	  prec_ml = new double[1+info];
	}
	std::vector<double> res_ol, res_mcfm;
	if (mode&1) {
	  int retcode;
	  double prec(-1.0);
	  double pn_ml[4*(nIn+nOut)];
	  for (size_t n(0);n<pn.size();++n) GetMomML(pn_ml,n,pn[n]);
	  me(pn_ml,res_ml,&prec,prec_ml,&retcode);
	  res_ol={res_ml[1],res_ml[2],res_ml[3]};
	}
	else {
	  mcfm.SetAlphaS(std::stod(params["alpha_S"]));
	  mcfm.Calc(pid_mcfm, pn, 1);
	  res_ol = mcfm.GetResult(pid_mcfm);
	  for (int i(0);i<4;++i)
	    res_ol[i] /= mcfm.GetProcess(pid_mcfm)->GetSymmetryFactor();
	  res_ol[0] -= res_ol[3]*mcfm.GetProcess(pid_mcfm)->GetSchemeFactor(0);
	}
	if (mode==0 || mode==3) {
	  for (size_t i(0);i<pn.size();++i)
	    pn[i] = FourVec(pn[i][0],pn[i][2],pn[i][3],pn[i][1]);
	}
	if (mode&2) {
	  int retcode;
	  double prec(-1.0);
	  double pn_ml[4*(nIn+nOut)];
	  for (size_t n(0);n<pn.size();++n) GetMomML(pn_ml,n,pn[n]);
	  me(pn_ml,res_ml,&prec,prec_ml,&retcode);
	  res_mcfm={res_ml[1],res_ml[2],res_ml[3]};
	}
	else {
	  mcfm.SetAlphaS(std::stod(params["alpha_S"]));
	  mcfm.Calc(pid_mcfm, pn, 1);
	  res_mcfm = mcfm.GetResult(pid_mcfm);
	  for (int i(0);i<4;++i)
	    res_mcfm[i] /= mcfm.GetProcess(pid_mcfm)->GetSymmetryFactor();
	  res_mcfm[0] -= res_mcfm[3]*mcfm.GetProcess(pid_mcfm)->GetSchemeFactor(0);
	}
	// Calculate accuracy.
	double acc = -log10(fabs((res_ol[0]-res_mcfm[0])/res_ol[0]));
	if (printMom&1) {
	  std::cout << std::setprecision(17) 
		    << ((mode&2)?"ML5:  ":"MCFM: ") << res_mcfm[0] << std::endl
		    << ((mode&1)?"ML5:  ":"MCFM: ") << res_ol[0] << std::endl
		    << "Acc:  " << acc << std::endl << std::endl;
	}
	if (std::isnan(acc) || std::isinf(acc)) ++nans;
	else accs.push_back(acc);
      }
      ids.clear();
      pn.clear();
      if (accs.size()>=nPts) break;
    }
  }
  double sum(0.), sum2(0.), n(accs.size());
  for (size_t i(0);i<n;++i) {
    sum+=accs[i];
    sum2+=accs[i]*accs[i];
  }
  std::cout << "Checked "
	    << ((mode&2)?"ML5":"MCFM") << " against "
	    << ((mode&1)?"ML5":"MCFM") << std::endl;
  std::cout<<"Average Acc: "<<sum/n<<" +- "
	   <<sqrt((sum2/n-pow(sum/n,2))/(n-1.0))<<"\n";
  std::sort(accs.begin(),accs.end());
  std::cout<<"Median  Acc: "<<accs[n/2]
	   <<" - "<<accs[n/2]-accs[n/4]
	   <<" + "<<accs[3*n/4]-accs[n/2]<<"\n";
  std::cout<<"Worst   Acc: "<<accs[0]<<" ("<<nans<<" nans)\n";
  std::cout<<"Best    Acc: "<<accs[n-1]<<"\n\n";
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
