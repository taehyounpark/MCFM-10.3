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
 
// Parameter setting for test programs.

// Standard includes.
#include <fstream>
#include <sstream>

std::string trim(const std::string& str,
  const std::string& whitespace = " \t") {
  const size_t strBegin = str.find_first_not_of(whitespace);
  if (strBegin == std::string::npos) return "";
  const size_t strEnd   = str.find_last_not_of(whitespace);
  const size_t strRange = strEnd - strBegin + 1;
  return str.substr(strBegin, strRange);
}

std::string reduce(const std::string& str,
  const std::string& fill = " ",
  const std::string& whitespace = " \t") {
  std::string res   = trim(str, whitespace);
  size_t beginSpace = res.find_first_of(whitespace);
  while (beginSpace != std::string::npos) {
    const size_t endSpace = res.find_first_not_of(whitespace, beginSpace);
    res.replace(beginSpace, endSpace-beginSpace, fill);
    const size_t newStart = beginSpace + fill.length();
    beginSpace = res.find_first_of(whitespace, newStart);
  }
  return res;
}

std::string chop(const std::string& str, const std::string com) {
  std::string res = str;
  const size_t beginComment = res.find_first_of(com);
  if (beginComment == std::string::npos) return res;
  return res.substr(0,beginComment);
}

void read_parameters(const std::string fname,
  std::map<std::string,std::string>& params) {
  params.clear();
  std::ifstream paramfile(fname);
  if (!paramfile.good()) return;
  std::string line;
  while (std::getline(paramfile,line)) {
    line = reduce(line);
    line = chop(line,"!#");
    if (line != "") {
      const size_t sep = line.find_first_of(" ");
      std::string  key = line.substr(0,sep);
      std::string  arg = line.substr(sep+1,line.size()-sep);
      params[key] = arg;
    }
  }
}

void get_sm_parameters(std::map<std::string,std::string> &params,
		       const int diagonalCKM=1) {
  // Masses and widths.
  params["n_flav"]             = "5";
  params["down_mass"]          = "0";
  params["up_mass"]            = "0";
  params["strange_mass"]       = "0";
  params["charm_mass"]         = "0";
  params["charm_yukawa"]       = "2.0164";
  params["bottom_mass"]        = "0";
  params["bottom_yukawa"]      = "23.04";
  params["top_mass"]           = "173.21";
  params["top_width"]          = "0";
  params["electron_mass"]      = "0";
  params["muon_mass"]          = "0";
  params["tau_mass"]           = "0";
  params["tau_yukawa"]         = "3.157729";
  params["tau_width"]          = "2.26735e-12";
  params["H_mass"]             = "125.";
  params["H_width"]            = "0.00407";
  params["Z_mass"]             = "91.1876";
  params["Z_width"]            = "2.4952";
  params["W_mass"]             = "80.385";
  params["W_width"]            = "2.085";
  // EW parameters.
  params["ew_scheme"]          = "5";
  params["alpha_EM"]           = "7.7676170824406213e-3";
  params["Gf"]                 = "1.197450526084678e-5";
  params["sin2_thetaW"]        = "0.2228972225239183";
  // Scales and couplings.
  params["order_alpha_S"]      = "1";
  params["alpha_S"]            = "0.118";
  params["scale"]              = params["Z_mass"];
  // CKM elements.
  if (diagonalCKM) {
    params["CKM_u_d"]          = "1";
    params["CKM_u_s"]          = "0";
    params["CKM_u_b"]          = "0";
    params["CKM_c_d"]          = "0";
    params["CKM_c_s"]          = "1";
    params["CKM_c_b"]          = "0";
  } else {
    params["CKM_u_d"]          = "0.97383";
    params["CKM_u_s"]          = "0.2272";
    params["CKM_u_b"]          = "0.00396";
    params["CKM_c_d"]          = "0.2271";
    params["CKM_c_s"]          = "0.97296";
    params["CKM_c_b"]          = "0.04221";
  }
}
