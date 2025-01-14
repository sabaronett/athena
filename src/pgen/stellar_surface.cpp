//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 *
 * This program is free software: you can redistribute and/or modify it under the terms
 * of the GNU General Public License (GPL) as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of GNU GPL in the file LICENSE included in the code
 * distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

// C headers

// C++ headers
#include <algorithm>  // min
#include <cmath>      // sqrt
#include <fstream>
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../mesh/mesh.hpp"
#include "../nr_radiation/radiation.hpp"
#include "../parameter_input.hpp"

// File scope variables
static Real theta;
static int radial;

int RefinementCondition(MeshBlock *pmb);

//======================================================================================
//! \file stellar_surface.cpp
//  \brief Radiative transfer test of a stellar surface
//======================================================================================

void TwoBeams(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
              const AthenaArray<Real> &w, FaceField &b,
              AthenaArray<Real> &ir, Real time, Real dt,
              int is, int ie, int js, int je, int ks, int ke, int ngh);

void TwoBeamHydro(
    MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (adaptive)
    EnrollUserRefinementCondition(RefinementCondition);

  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, TwoBeamHydro);
  if (NR_RADIATION_ENABLED || IM_RADIATION_ENABLED)
    EnrollUserRadBoundaryFunction(BoundaryFace::inner_x1, TwoBeams);
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  theta = pin->GetOrAddReal("problem", "theta", 0.0);
  radial = pin->GetOrAddInteger("problem", "radial", 0);
  return;
}


//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief stellar surface test
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real gamma = peos->GetGamma();

  // Initialize hydro variable
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        phydro->u(IDN,k,j,i) = 1.0;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = 1.0/(gamma-1.0);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
        }
      }
    }
  }

  // Now initialize opacity and specific intensity
  if (NR_RADIATION_ENABLED || IM_RADIATION_ENABLED) {
    int nfreq = pnrrad->nfreq;
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          for (int ifr=0; ifr < nfreq; ++ifr) {
            pnrrad->sigma_s(k,j,i,ifr) = 0.0;
            pnrrad->sigma_a(k,j,i,ifr) = 0.0;
            pnrrad->sigma_pe(k,j,i,ifr) = 0.0;
            pnrrad->sigma_p(k,j,i,ifr) = 0.0;
          }
          for (int n=0; n<pnrrad->n_fre_ang; ++n) {
              pnrrad->ir(k,j,i,n) = 0.0;
          }
        }
      }
    }
  }
  return;
}

void TwoBeams(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
              const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir, Real time,
              Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {
  int nang = prad->nang;       // total n-hat angles N
  Real max_mu_x = 0;           // max(\mu_x)
  int nact = 0;                // active angles
  
  if (radial == 1) {
    for (int n=0; n<nang; ++n) { // find independent of Rad_angles.txt order
      if (prad->mu(0,0,0,0,n) > max_mu_x) max_mu_x = prad->mu(0,0,0,0,n);
    }
    for (int n=0; n<nang; ++n) { // count active angles
      if (prad->mu(0,0,0,0,n) == max_mu_x) ++nact;
    }
  }

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
        for (int n=0; n<nang; ++n) {
          ir(k,j,is-i,n) = 0.0;

          if (radial == 1 && prad->mu(0,k,j,is-i,n) == max_mu_x) { // radial emission only
            Real ir_adj = 1/(prad->wmu(n)*nact);
            
            if (theta > 0.0) {
              Real const &x2 = pco->x2v(j);
              Real dis = std::abs(x2 - theta);
              
              if (dis < pco->dx2v(j)) ir(k,j,is-i,n) = ir_adj;
            }
            else {
              ir(k,j,is-i,n) = ir_adj;
            }
          }
          else {                                             // isotropic surface emission
            Real ir_adj = 1/(prad->wmu(n)*nang/2);

            if (theta > 0.0) {
              Real const &x2 = pco->x2v(j);
              Real dis = std::abs(x2 - theta);
              
              if (dis < pco->dx2v(j) && prad->mu(0,k,j,is-i,n) > 0) {
                ir(k,j,is-i,n) = ir_adj;
              }
            }
            else {
              if (prad->mu(0,k,j,is-i,n) > 0) ir(k,j,is-i,n) = ir_adj;
            }
          }
        }
      }
    }
  }
  return;
}

void TwoBeamHydro(
    MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {
  for (int n=0; n<NHYDRO; ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          prim(n,k,j,is-i) = prim(n,k,j,is);
        }
      }
    }
  }
}

// refinement condition: density curvature
int RefinementCondition(MeshBlock *pmb) {
  Coordinates *pco = pmb->pcoord;
  Real xmin = pco->x1f(pmb->is);
  Real xmax = pco->x1f(pmb->ie+1);
  Real ymin = pco->x2f(pmb->js);
  Real ymax = pco->x2f(pmb->je+1);
  Real zmin = pco->x3f(pmb->ks);
  Real zmax = pco->x3f(pmb->ke+1);
  // refine : delta rho > 0.9*amp
  if (ymin > -1 && ymax < 1) return 1;
  return 0;
}
