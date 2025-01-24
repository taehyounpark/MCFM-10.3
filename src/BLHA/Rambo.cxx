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
 

#include <functional>
#include "MCFM/CXX_Interface.h"

namespace MCFM {

// Solve f(x) = target for x in the specified range.
// Adapted from Pythia8304/src/MathTools.
bool brent(double& solutionOut, std::function<double(double)> f,  double target,
  double xLo, double xHi) {

  double tol = 1.e-10;
  int maxIter = 10000;

  // Range checks.
  if (xLo > xHi) return false;

  // Evaluate function - targetValue at lower boundary.
  double f1 = f(xLo) - target;
  if (abs(f1) < tol) {
    solutionOut = xLo;
    return true;
  }
  // Evaluate function - targetValue at upper boundary.
  double f2 = f(xHi) - target;
  if (abs(f2) < tol) {
    solutionOut = xHi;
    return true;
  }

  // Check if root is bracketed.
  if ( f1 * f2 > 0.0) return false;

  // Start searching for root.
  double x1 = xLo;
  double x2 = xHi;
  double x3 = 0.5 * (xLo + xHi);

  int iter=0;
  while(++iter < maxIter) {
    // Now check at x = x3.
    double f3 = f(x3) - target;
    // Check if tolerance on f has been reached.
    if (abs(f3) < tol) {
      solutionOut = x3;
      return true;
    }
    // Is root bracketed in lower or upper half?
    if (f1 * f3 < 0.0) xHi = x3;
    else xLo = x3;
    // Check if tolerance on x has been reached.
    if ((xHi - xLo) < tol * (abs(xHi) < 1.0 ? xHi : 1.0)) {
      solutionOut = 0.5 * (xLo + xHi);
      return true;
    }

    // Work out next step to take in x.
    double den = (f2 - f1) * (f3 - f1) * (f2 - f3);
    double num = x3 * (f1 - f2) * (f2 - f3 + f1) + f2 * x1 * (f2 - f3)
               + f1 * x2 * (f3 - f1);
    double dx = xHi - xLo;
    if (den != 0.0) dx = f3 * num / den;

    // First attempt, using gradient
    double x = x3 + dx;

    // If this was too far, just step to the middle
    if ((xHi - x) * (x - xLo) < 0.0) {
      dx = 0.5 * (xHi - xLo);
      x = xLo + dx;
    }
    if (x < x3) {
        x2 = x3;
        f2 = f3;
    }
    else {
        x1 = x3;
        f1 = f3;
    }
    x3 = x;
  }

  // Maximum number of iterations exceeded.
  return false;
}

// Function to generate massless final-state momenta using RAMBO.
// Adapted from Pythia8.302/src/Vincia.cc.
double Rambo(double Ecm, int n, std::vector<FourVec>& pn) {
  // Set size of output vector
  pn.resize(n);
  // Special case of one outgoing momentum.
  if (n == 1) {
    pn[0] = FourVec(Ecm, 0., 0., 0.);
    return 1.0;
  }
  // Create momentum-sum four-vector
  FourVec R;
  // Generate nParticles independent massless 4-momenta with isotropic angles
  for (int i = 0; i < n; ++i) {
    // Cos(theta), sin(theta), and phi
    double c   = 2.0*(double(rand())/RAND_MAX) - 1.0;
    double s   = sqrt(1.0-c*c);
    double phi = 2.0*M_PI*(double(rand())/RAND_MAX);
    // Norm
    double r12 = 0.0;
    while (r12 == 0.0) {
      double r1 = double(rand())/RAND_MAX;
      double r2 = double(rand())/RAND_MAX;
      r12 = r1*r2;
    }
    double En = -log(r12);
    pn[i].e(En);
    pn[i].pz(En*c);
    pn[i].py(En*s*cos(phi));
    pn[i].px(En*s*sin(phi));
    // Add to vector and add to sum
    for (int j(0);j<4;++j) R[j] += pn[i][j];
  }
  // Compute ECM and normalise to unity (with sign flip)
  double Rmass = sqrt(R.e()*R.e()-R.px()*R.px()-R.py()*R.py()-R.pz()*R.pz());
  for (int j(0);j<4;++j) R[j] /= -Rmass;
  // Transform momenta so add up to (eCM, 0, 0, 0)
  double a = 1.0/(1.0-R.e());
  double x = Ecm/Rmass;
  for (int i = 0; i < n; ++i) {
    double bq = R.px()*pn[i].px()+R.py()*pn[i].py()+R.pz()*pn[i].pz();
    pn[i].px( x * (pn[i].px()+R.px()*(pn[i].e()+a*bq)) );
    pn[i].py( x * (pn[i].py()+R.py()*(pn[i].e()+a*bq)) );
    pn[i].pz( x * (pn[i].pz()+R.pz()*(pn[i].e()+a*bq)) );
    pn[i].e(  x * (-R.e()*pn[i].e()+bq) );
  }
  // The weight is always unity for the massless algorithm.
  return 1.0;

}

// Function to generate massive final-state momenta using RAMBO.
// Adapted from Pythia8.304/src/PhaseSpace.cc.
double Rambo(double Ecm, std::vector<double> mOut, std::vector<FourVec>& pOut) {

  // Call the massless genPoint, initializing weight.
  int nOut = mOut.size();
  if (nOut < 1 || Ecm <= 0.) return 0.;
  if (nOut == 1) {
    pOut.resize(1);
    pOut[0] = FourVec(Ecm, 0., 0., sqrt(Ecm*Ecm-mOut[0]*mOut[0]));
    return 1.0;
  }
  double weight = Rambo(Ecm, nOut, pOut);
  bool massesnonzero = false;

  // Set up the function determining the rescaling parameter xi.
  std::vector<double> energies;
  for (int i(0); i<nOut;++i) {
    energies.push_back(pOut[i].e());
    if ((mOut[i]/Ecm)*(mOut[i]/Ecm) > 1e-9) massesnonzero = true;
  }

  // If none of the reduced masses is > 1e-9, return.
  if (!massesnonzero) return weight;

  // Set up the mass and energy vectors.
  std::vector<double> mXi(0), energiesXi(0);
  if (mOut.size() == energies.size()) {mXi = mOut; energiesXi = energies;}

  // Define the Xi function.
  std::function<double(double)> rhs = [&mXi, &energiesXi](double xi) -> double{
    double retval = 0.;
    for (size_t i(0); i<mXi.size();++i)
      retval += sqrt(mXi[i]*mXi[i] + xi*xi*energiesXi[i]*energiesXi[i]);
    return retval;
  };
  
  // Rescale all the momenta.
  double xi(0);
  brent(xi, rhs, Ecm, 0., 1.);
  for (int i(0);i<nOut;++i) {
    pOut[i].m_p[1] *= xi;
    pOut[i].m_p[2] *= xi;
    pOut[i].m_p[3] *= xi;
    pOut[i].e(sqrt(mOut[i]*mOut[i] + xi*xi*pOut[i].e()*pOut[i].e()) );
  }

  // Determine the quantities needed for the calculation of the weight.
  double sumP(0.), prodPdivE(1.), sumP2divE(0.);
  for (int i(0); i<nOut; ++i) {
    double pAbs2 = pOut[i].px()*pOut[i].px()+pOut[i].py()*pOut[i].py()
      +pOut[i].pz()*pOut[i].pz();
    double pAbs  = sqrt(pAbs2);
    sumP      += pAbs;
    prodPdivE *= pAbs/pOut[i].e();
    sumP2divE += pAbs2/pOut[i].e();
  }

  // There's a typo in eq. 4.11 of the Rambo paper by Kleiss, Stirling
  // and Ellis, the Ecm below is not present there.
  weight *= pow(sumP/Ecm,2*nOut-3)*prodPdivE*Ecm/sumP2divE;
  return weight;
}
  
// Function to ensure precision in massive PS points.
// Taken from Sherpa 2.3.0/PHASIC++/Channels/Rambo.C
void MassivePoint(std::vector<FourVec>& p, double ET, std::vector<double> ms, int nin, int nout) {
  int itmax = 6;
  double accu = ET * 1.e-14;
  double xmt = 0.;
  double x;
  double xm2[nin+nout+1];
  double p2[nin+nout+1];
  double E[nin+nout+1];
  for (short int i=nin;i<nin+nout;i++) {
    xmt   += ms[i];
    xm2[i] = ms[i]*ms[i];
    p2[i]  = p[i].e()*p[i].e();
  }
  x = sqrt(1.-(xmt/ET)*(xmt/ET));

  // Massive particles: rescale their momenta by a common factor x.
  // Loop to calculate x.
  double f0,g0,x2;
  short int iter = 0;
  for (;;) {
    f0 = -ET;
    g0 = 0.;
    x2 = x*x;
    for (short int i=nin;i<nin+nout;i++) {
      E[i] = sqrt(xm2[i]+x2*p2[i]);
      f0  += E[i];
      g0  += p2[i]/E[i];
    }
    if (fabs(f0)<accu) break;
    iter++;
    if (iter>itmax) break;
    x -= f0/(x*g0);
  }

  // Construct momenta.
  for (short int i=nin;i<nin+nout;i++)
    p[i] = FourVec(E[i],x*p[i].px(),x*p[i].py(),x*p[i].pz());
}

// Function to sample a phase space point flatly.
double GenPoint(double Ecm, std::vector<FourVec>& pn,
  int nIn, const std::vector<double>& masses) {
  double sqrts;
  std::vector<double> mOut;
  for (size_t i(nIn);i<masses.size();++i) mOut.push_back(masses[i]);
  if (nIn == 2) {
    double x = 1.;//double(rand())/RAND_MAX;
    sqrts = x*Ecm;
    double a(0.5), b(0.5);
    // Special treatment if one outgoing particle.
    if (mOut.size() == 1) {
      sqrts = mOut[0]+(sqrts-mOut[0])*double(rand())/RAND_MAX;
      a = sqrt(sqrts*sqrts-mOut[0]*mOut[0])/(2.*sqrts)+1./2.;
      b = 1.-a;
    }
    pn.push_back(FourVec(a*sqrts, 0., 0., a*sqrts));
    pn.push_back(FourVec(b*sqrts, 0., 0., -b*sqrts));
  } else if (nIn == 1) {
    pn.push_back(FourVec(Ecm, 0., 0., 0.));
    sqrts = Ecm;
  } else return 0.;
  std::vector<FourVec> pOut;
  double wt = Rambo(sqrts, mOut, pOut);
  for (size_t i(0);i<masses.size()-nIn;++i) pn.push_back(pOut[i]);
  //  MassivePoint(pn, sqrts, masses, nIn, mOut.size());
  return wt;
}

// Function to sample a phase space point flatly.
double GenPoint(double Ecm, std::vector<FourVec>& pn, std::vector<int>& ids, int nIn,
		std::map<std::string,std::string>& params) {
  // Special handling of massive Higgs and top.
  const int nOut = ids.size()-nIn;
  std::vector<double> masses(nIn+nOut, 0.);
  for (int i(nIn);i<nIn+nOut;++i) {
    if (ids[i] == 25) masses[i] = std::stod(params["H_mass"]);
    if (abs(ids[i]) == 6) masses[i] = std::stod(params["top_mass"]);
  }
  return GenPoint(Ecm, pn, nIn, masses);
}

}//end of namespace MCFM
