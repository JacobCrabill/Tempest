#include "solver.hpp"

void solver::applyBoundaryConditions(void)
{
  for (int bnd=0; bnd<Geo->nBounds; bnd++) {

#pragma omp parallel for
    for (int i=0; i<Geo->nBndPts[bnd]; i++) {
      int iv = bndPts(bnd,i);

      switch (bcList[bnd]) {
      case SUP_IN:
        U(iv,0) = params->rhoBound;
        U(iv,1) = params->rhoBound*params->uBound;
        U(iv,2) = params->rhoBound*params->vBound;
        if (params->nDims == 3)
          U(iv,3) = params->rhoBound*params->wBound;
        U(iv,params->nDims+1) = params->pBound/(params->gamma-1) + 0.5*U(iv,0)*params->Uinf*params->Uinf;

      case SUP_OUT:
        // Extrapolated, so do nothing
        break;

      case SLIP_WALL: {
        Vec3 norm = bndNorm(bnd,i);
        double rho = U(iv,0);
        double u = U(iv,1)/rho;
        double v = U(iv,2)/rho;
        double w = 0;
        if (nDims == 3)
          w = U(iv,3)/rho;
        double vn = u*norm[0]+v*norm[1]+w*norm[2];
        // Reflect velocity at wall (Strong enforcement, not weak)
        u -= vn*norm[0];
        v -= vn*norm[1];
        w -= vn*norm[2];
        U(iv,1) = rho*u;
        U(iv,2) = rho*v;
        if (nDims == 3)
          U(iv,3) = rho*w;
        break;
      }
      case CHAR: {
        double gamma = params->gamma;
        Vec3 norm = bndNorm(bnd,i);
        double rho = U(iv,0);
        double u = U(iv,1)/rho;
        double v = U(iv,2)/rho;
        double w = 0;
        if (nDims == 3)
          w = U(iv,3)/rho;
        double p = (gamma-1)*(U(iv,nDims+1) - 0.5*rho*(u*u+v*v+w*w));
        double vn = u*norm[0]+v*norm[1]+w*norm[2];

        double vnBound = params->uBound*norm[0]+params->vBound*norm[1]+params->wBound*norm[2];

        double r_plus  = vn + 2./(gamma-1.)*sqrt(gamma*p/rho);
        double r_minus = vnBound - 2./(gamma-1.)*sqrt(gamma*params->pBound/params->rhoBound);

        double cStar = 0.25*(gamma-1.)*(r_plus-r_minus);
        double vn_star = 0.5*(r_plus+r_minus);

        // Inflow
        if (vn<0) {
          // HACK
          double one_over_s = pow(params->rhoBound,gamma)/params->pBound;

          // freestream total enthalpy
          double vSq = params->uBound*params->uBound+params->vBound*params->vBound+params->wBound*params->wBound;
          double h_free_stream = gamma/(gamma-1.)*params->pBound/params->rhoBound + 0.5*vSq;

          U(iv,0) = pow(1./gamma*(one_over_s*cStar*cStar),1./(gamma-1.));

          // Compute velocity on the right side
          U(iv,1) = U(iv,0)*(vn_star*norm[0] + params->uBound - vnBound*norm[0]);
          U(iv,2) = U(iv,0)*(vn_star*norm[1] + params->vBound - vnBound*norm[1]);
          U(iv,3) = U(iv,0)*(vn_star*norm[2] + params->wBound - vnBound*norm[2]);

          double pR = U(iv,0)/gamma*cStar*cStar;
          U(iv,4) = U(iv,0)*h_free_stream - pR;
        }
        // Outflow
        else {
          double one_over_s = pow(rho,gamma)/p;

          // freestream total enthalpy
          U(iv,0) = pow(1./gamma*(one_over_s*cStar*cStar), 1./(gamma-1.));

          // Compute velocity on the right side
          double uR = vn_star*norm[0] + (u - vn*norm[0]);
          double vR = vn_star*norm[1] + (v - vn*norm[1]);
          double wR = vn_star*norm[2] + (w - vn*norm[2]);

          U(iv,1) = U(iv,0)*uR;
          U(iv,2) = U(iv,0)*vR;
          U(iv,3) = U(iv,0)*wR;

          double pR = U(iv,0)/gamma*cStar*cStar;
          double vSq = uR*uR+vR*vR+wR*wR;
          U(iv,4) = (pR/(gamma-1.0)) + 0.5*U(iv,0)*vSq;
        }
        break;
      }
      }
    }
  }
}
