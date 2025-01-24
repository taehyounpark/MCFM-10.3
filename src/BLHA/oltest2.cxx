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

// OpenLoops includes.
#include "openloops.h"
void set_sm_parameters_ol(std::map<std::string,std::string> &params,
			  const int diagonalCKM=1);

using namespace MCFM;

// Convert FourVec struct to array used by OpenLoops.
void GetMomOL(const std::vector<FourVec>& pn, double* pn_ol);

// Main program.
int main(int argc, char** argv) {

  // Create instance of MCFM interface.
  CXX_Interface mcfm;

  // Read settings.
  std::map<std::string,std::string> params;
  read_parameters("params.lh",params);

  // Command-line options.
  int foEW = -1, foQCD = -1, diagonalCKM = 1, nIn = 2;
  int nPts = 1000000, checkPoles = 1, loopacc = 8;
  int printMom = 0, mode = 0, ampType = 11, fampType = -1;
  std::vector<int> decids, decfls;
  int opt, verbosity = 0;
  while ((opt = getopt(argc, argv, "vM:m:a:A:n:i:o:S:W:C:P:h")) != -1) {
    switch (opt) {
    case 'v': verbosity += 1; break;
    case 'M': printMom = atoi(optarg); break;
    case 'm': mode = atoi(optarg); break;
    case 'n': nPts = atoi(optarg); break;
    case 'i': nIn = atoi(optarg); break;
    case 'a': loopacc = atoi(optarg); break;
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
	       <<"  -A type   -- set OpenLoops amplitude type to 'type'\n"
	       <<"  -M print  -- set print mode to 'print'\n"
	       <<"  -a acc    -- set OpenLoops loop accuracy to 'acc'\n"
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
  int nPt(0), nans(0), pid_mcfm(0), pid_ol(0);
  std::vector<int> ids;
  std::vector<FourVec> pn;
  std::string line;
  std::vector<double> accs;
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
	    ol_setparameter_string("model", params["model"].c_str());
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
	
	  // Initialize OpenLoops.
	  set_sm_parameters_ol(params,diagonalCKM);
	  std::string process;
	  for (int i(0);i<nIn;++i) process += std::to_string(ids[i])+" ";
	  process += "->";
	  for (int i(nIn);i<nIn+nOut;++i) process += " "+std::to_string(ids[i]);
	  ol_setparameter_int("order_ew", oEW);
	  pid_ol=ol_register_process(process.c_str(), ampType);
	  if (pid_ol < 0) {
	    std::cout << " Process not available in OpenLoops." << std::endl;
	    exit(EXIT_FAILURE);
	  }
	  ol_start();
	}
	std::vector<double> res_ol, res_mcfm;
	if (mode&1) {
	  double pn_ol[5*(nIn+nOut)];
	  GetMomOL(pn, pn_ol);
	  double res_ol_tree, acc_ol, res_ol_loop[3] = {0., 0., 0.};
	  ol_setparameter_int("psp_cleaning", 1);
	  ol_setparameter_int("hp_mode", 2);
	  ol_setparameter_int("hp_loopacc", loopacc);
	  if (ampType == 12) ol_evaluate_loop2(pid_ol, pn_ol, &res_ol_tree, &acc_ol);
	  else ol_evaluate_loop(pid_ol, pn_ol, &res_ol_tree, res_ol_loop, &acc_ol);
	  res_ol = {res_ol_loop[0], res_ol_loop[1], res_ol_loop[2], res_ol_tree};
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
	  double pn_ol[5*(nIn+nOut)];
	  GetMomOL(pn, pn_ol);
	  double res_ol_tree, acc_ol, res_ol_loop[3] = {0., 0., 0.};
	  ol_setparameter_int("psp_cleaning", 1);
	  ol_setparameter_int("hp_mode", 2);
	  ol_setparameter_int("hp_loopacc", loopacc);
	  if (ampType == 12) ol_evaluate_loop2(pid_ol, pn_ol, &res_ol_tree, &acc_ol);
	  else ol_evaluate_loop(pid_ol, pn_ol, &res_ol_tree, res_ol_loop, &acc_ol);
	  res_mcfm = {res_ol_loop[0], res_ol_loop[1], res_ol_loop[2], res_ol_tree};
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
		    << ((mode&2)?"OL2:  ":"MCFM: ") << res_mcfm[0] << std::endl
		    << ((mode&1)?"OL2:  ":"MCFM: ") << res_ol[0] << std::endl
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
	    << ((mode&2)?"OL2":"MCFM") << " against "
	    << ((mode&1)?"OL2":"MCFM") << std::endl;
  std::cout<<"Average Acc: "<<sum/n<<" +- "
	   <<sqrt((sum2/n-pow(sum/n,2))/(n-1.0))<<"\n";
  std::sort(accs.begin(),accs.end());
  std::cout<<"Median  Acc: "<<accs[n/2]
	   <<" - "<<accs[n/2]-accs[n/4]
	   <<" + "<<accs[3*n/4]-accs[n/2]<<"\n";
  std::cout<<"Worst   Acc: "<<accs[0]<<" ("<<nans<<" nans)\n";
  std::cout<<"Best    Acc: "<<accs[n-1]<<"\n\n";
  // Done.
  ol_finish();
  return 0;
}

// Get momenta as needed for OpenLoops.
void GetMomOL(const std::vector<FourVec>& pn, double* pn_ol) {
  for (size_t i(0); i<pn.size(); ++i) {
    pn_ol[5*i]   = pn[i].e();
    pn_ol[5*i+1] = pn[i].px();
    pn_ol[5*i+2] = pn[i].py();
    pn_ol[5*i+3] = pn[i].pz();
    pn_ol[5*i+4] = sqrt(m2(pn[i]));
  }
}

void set_sm_parameters_ol(std::map<std::string,std::string> &params,
			  const int diagonalCKM) {
  ol_setparameter_string("model", params["model"].c_str());
  ol_setparameter_double("mass(1)", std::stod(params["down_mass"]));
  ol_setparameter_double("mass(2)", std::stod(params["up_mass"]));
  ol_setparameter_double("mass(3)", std::stod(params["strange_mass"]));
  ol_setparameter_double("mass(4)", std::stod(params["charm_mass"]));
  ol_setparameter_double("yuk(4)", std::stod(params["charm_yukawa"]));
  ol_setparameter_double("mass(5)", std::stod(params["bottom_mass"]));
  ol_setparameter_double("yuk(5)", std::stod(params["bottom_yukawa"]));
  ol_setparameter_double("mass(6)", std::stod(params["top_mass"]));
  ol_setparameter_double("yuk(6)", std::stod(params["top_yukawa"]));
  ol_setparameter_double("width(6)", std::stod(params["top_width"]));
  ol_setparameter_double("mass(11)", std::stod(params["electron_mass"]));
  ol_setparameter_double("mass(13)", std::stod(params["muon_mass"]));
  ol_setparameter_double("mass(15)", std::stod(params["tau_mass"]));
  ol_setparameter_double("yuk(15)", std::stod(params["tau_yukawa"]));
  ol_setparameter_double("width(15)", std::stod(params["tau_width"]));
  ol_setparameter_double("mass(23)", std::stod(params["Z_mass"]));
  ol_setparameter_double("width(23)", std::stod(params["Z_width"]));
  ol_setparameter_double("mass(24)", std::stod(params["W_mass"]));
  ol_setparameter_double("width(24)", std::stod(params["W_width"]));
  ol_setparameter_double("mass(25)", std::stod(params["H_mass"]));
  ol_setparameter_double("width(25)", std::stod(params["H_width"]));
  // EW parameters.
  ol_setparameter_int("ew_scheme",2);
  ol_setparameter_double("alpha_qed_mz", std::stod(params["alpha_EM"]));
  // CKM elements.
  if (diagonalCKM) ol_setparameter_int("ckmorder", 0);
  else {
    ol_setparameter_double("VCKMdu", std::stod(params["CKM_u_d"]));
    ol_setparameter_double("VCKMsu", std::stod(params["CKM_u_s"]));
    ol_setparameter_double("VCKMbu", std::stod(params["CKM_u_b"]));
    ol_setparameter_double("VCKMdc", std::stod(params["CKM_c_d"]));
    ol_setparameter_double("VCKMsc", std::stod(params["CKM_c_s"]));
    ol_setparameter_double("VCKMbc", std::stod(params["CKM_c_b"]));
  }
  // Scales and couplings.
  ol_setparameter_double("alpha_s", std::stod(params["alpha_S"]));
  ol_setparameter_double("mu", std::stod(params["scale"]));
  ol_setparameter_double("muren", std::stod(params["scale"]));
}
