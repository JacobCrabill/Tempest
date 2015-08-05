/*!
 * \file solver.cpp
 * \brief Class to store all solution data & apply FR operators
 *
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * \version 0.0.1
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2014 Jacob Crabill.
 *
 */

#include "../include/solver.hpp"

#include <sstream>
#include <omp.h>

#include "input.hpp"
#include "flux.hpp"
#include "geo.hpp"

solver::solver()
{
  //tg = NULL;
}

solver::~solver()
{
//  if (tg != NULL)
//    delete tg;
}

void solver::setup(input *params, geo *Geo)
{
  this->params = params;
  this->Geo = Geo;
  //this->tg = Geo->tg; // Geo will have initialized this already if needed

  params->time = 0.;

  if (params->equation == NAVIER_STOKES) {
    if (params->nDims == 2)
      params->nFields = 4;
    else if (params->nDims == 3)
      params->nFields = 5;
  }

  nDims = params->nDims;
  nFields = params->nFields;

  gridID = Geo->gridID;
  gridRank = Geo->gridRank;
  nprocPerGrid = Geo->nprocPerGrid;

  /* Setup the final dual-mesh data structures */
  setupDualMesh();

  if (params->restart)
    readRestartFile();
  else
    initializeSolution();

  if (params->meshType == OVERSET_MESH)
    setupOverset();

  /* Additional Setup */

  // Time advancement setup
  switch (params->timeType) {
    case 0:
      nRKSteps = 1;
      RKb = {1};
      break;
    case 4:
      nRKSteps = 4;
      RKa = {.5, .5, 1.};
      RKb = {1./6., 1./3., 1./3., 1./6.};
      break;
    default:
      FatalError("Time-Stepping type not supported.");
  }

  F.setup(nEdges,nDims,nFields);

  tempFL.setup(nDims,nFields);
  tempFR.setup(nDims,nFields);

  waveSp.resize(nEdges);

  divF.setup(nRKSteps,nVerts,nFields);

  if (params->viscous) {
    gradU.setup(nVerts,nDims,nFields);
  }

//#ifndef _NO_MPI
//  finishMpiSetup();
//#endif
}

void solver::update(void)
{
  // For RK time-stepping, store the starting solution values
  if (nRKSteps>1)
    copyU_U0();

  /* Intermediate residuals for Runge-Kutta time integration */

  for (int step=0; step<nRKSteps-1; step++) {

    if (step == 0)
      params->rkTime = params->time;
    else
      params->rkTime = params->time + RKa[step-1]*params->dt;

    moveMesh(step);

    calcResidual(step);

    /* If in first stage, compute CFL-based timestep */
    if (step == 0 && params->dtType == 1) calcDt();  // -- NOT CONSISTENT WITH MOVING MESH SEQUENCE --

    timeStepA(step);

  }

  /* Final Runge-Kutta time advancement step */

  if (nRKSteps == 1) {
    params->rkTime = params->time;
    /* Calculate CFL-based timestep */
    if (params->dtType == 1) calcDt();
  }
  else {
    params->rkTime = params->time + params->dt;
  }

  moveMesh(nRKSteps-1);

  calcResidual(nRKSteps-1);

  // Reset solution to initial-stage values
  if (nRKSteps>1)
    copyU0_U();

  for (int step=0; step<nRKSteps; step++)
    timeStepB(step);

  params->time += params->dt;
}


void solver::calcResidual(int step)
{

#ifndef _NO_MPI
  doCommunication();
#endif

  applyBoundaryConditions();

  if (params->viscous)
    calcGradU();

  calcInviscidFlux();

#ifndef _NO_MPI
  calcInviscidFlux_mpi();
#endif

  if (params->viscous) {

    calcViscousFlux();

#ifndef _NO_MPI
    calcViscousFlux_mpi();
#endif

  }

  calcFluxDivergence(step);

}

void solver::setupDualMesh(void) {

  if (params->rank==0) cout << "Solver: Setting up dual mesh arrays" << endl;

  nVerts = Geo->nVerts;
  nEdges = Geo->nEdges;

  xv = Geo->xv;
  v2ne = Geo->v2nv;
  v2e = Geo->v2e;
  v2v = Geo->v2v;
  Ae = Geo->e2A;
  vol = Geo->v2vol;
  bcList = Geo->bcList;
  bndNorm = Geo->bndNorm;
  bndArea = Geo->bndArea;
  bndPts = Geo->bndPts;
  nBndPts = Geo->nBndPts;

  // Solution at vertices
  U.setup(nVerts,nFields);

  // Flux at faces
  F.setup(nEdges,nFields);
  Fn.setup(nEdges,nFields);

  // Dual-mesh face areas
  A.setup(nVerts,getMax(v2ne));
  normDir.setup(nVerts,getMax(v2ne));
  normDir.initializeToValue(1);
  for (int i=0; i<nVerts; i++) {
    for (int j=0; j<v2ne[i]; j++) {
      A(i,j) = Ae[v2e(i,j)];
      if (i == Geo->e2v(v2e(i,j),1)) {
        // This node is on "right" of edge, so "flip the normal"
        normDir(i,j) = -1;
      }
    }
  }
}

void solver::finishMpiSetup(void)
{
//  if (params->rank==0) cout << "Solver: Setting up MPI face communicataions" << endl;
}


void solver::calcInviscidFlux(void)
{
#pragma omp parallel for
  for (int i=0; i<nEdges; i++) {
    int ivL = Geo->e2v(i,0);
    int ivR = Geo->e2v(i,1);

    inviscidFlux(U[ivL],tempFL,params);
    inviscidFlux(U[ivR],tempFR,params);

    double tempFnL[5] = {0,0,0,0,0};
    double tempFnR[5] = {0,0,0,0,0};

    // Get primitive variables
    double rhoL = U(ivL,0);     double rhoR = U(ivR,0);
    double uL = U(ivL,1)/rhoL;  double uR = U(ivR,1)/rhoR;
    double vL = U(ivL,2)/rhoL;  double vR = U(ivR,2)/rhoR;

    double wL, pL, vnL=0.;
    double wR, pR, vnR=0.;

    // Calculate pressure
    if (params->nDims==2) {
      pL = (params->gamma-1.0)*(U(ivL,3)-rhoL*(uL*uL+vL*vL));
      pR = (params->gamma-1.0)*(U(ivR,3)-rhoR*(uR*uR+vR*vR));
    }
    else {
      wL = U(ivL,3)/rhoL;   wR = U(ivR,3)/rhoR;
      pL = (params->gamma-1.0)*(U(ivL,4)-rhoL*(uL*uL+vL*vL+wL*wL));
      pR = (params->gamma-1.0)*(U(ivR,4)-rhoR*(uR*uR+vR*vR+wR*wR));
    }

    // Get normal fluxes, normal velocities
    Vec3 norm = point(xv[ivR]) - point(xv[ivL]);
    norm /= norm.norm();

    for (int dim=0; dim<params->nDims; dim++) {
      vnL += norm[dim]*U(ivL,dim+1)/rhoL;
      vnR += norm[dim]*U(ivR,dim+1)/rhoR;
      for (int i=0; i<params->nFields; i++) {
        tempFnL[i] += norm[dim]*tempFL(dim,i);
        tempFnR[i] += norm[dim]*tempFR(dim,i);
      }
    }

    // Get maximum eigenvalue for diffusion coefficient
    double csqL = max(params->gamma*pL/rhoL,0.0);
    double csqR = max(params->gamma*pR/rhoR,0.0);
    double eigL = std::fabs(vnL) + sqrt(csqL);
    double eigR = std::fabs(vnR) + sqrt(csqR);
    waveSp[i] = max(eigL,eigR);

    // Calculate Rusanov flux
    for (int k=0; k<params->nFields; k++) {
      Fn(i,k) = 0.5*(tempFnL[k]+tempFnR[k] - waveSp[i]*(U(ivR,k)-U(ivL,k)));
    }
  }
}

void solver::doCommunication()
{

}

void solver::calcInviscidFlux_mpi()
{

}

void solver::calcViscousFlux(void)
{

}

void solver::calcViscousFlux_mpi()
{

}


void solver::calcFluxDivergence(int step)
{
#pragma omp parallel for
  for (int i=0; i<nVerts; i++) {
    //if (Geo->v2b[i]) continue;
    for (int k=0; k<nFields; k++) divF(step,i,k) = 0;

    for (int j=0; j<v2ne[i]; j++) {
      for (int k=0; k<nFields; k++) {
        divF(step,i,k) += Fn(v2e(i,j),k)*A(i,j)*normDir(i,j);
      }
    }
  }

  matrix<double> tmpF(nDims,nFields);
#pragma omp parallel for collapse(2)
  for (int ib=0; ib<nBounds; ib++) {
    for (int i=0; i<nBndPts[ib]; i++) {
      int iv = bndPts(ib,i);
      inviscidFlux(U[iv],tmpF,params);
      for (int dim=0; dim<nDims; dim++) {
        for (int k=0; k<nFields; k++) {
          divF(step,iv,k) += tmpF(dim,k)*bndNorm(ib,i)[dim]*bndArea(ib,i);
        }
      }
    }
  }
}

void solver::calcGradU(void)
{

}

//void solver::calcEntropyErr(void)
//{

//}

void solver::moveMesh(int step)
{

}

void solver::calcDt(void)
{
  double dt = INFINITY;

#pragma omp parallel for reduction(min:dt)
  for (int i=0; i<nVerts; i++) {
    //if (Geo->v2b[i]) continue;

    double rho = U(i,0);
    double rhovMagSq = U(i,1)*U(i,1) + U(i,2)*U(i,2);
    if (nDims == 3)
      rhovMagSq += U(i,3)*U(i,3);
    double p = (params->gamma-1)*U(i,nDims+1) - 0.5*rhovMagSq/rho;
    double a = sqrt(max(params->gamma*p/rho,0.));

    double minA = INFINITY;
    for (int j=0; j<v2ne[i]; j++)
      minA = min(minA,A(i,j));
    double dx = vol[i]/minA;

    dt = std::min(dt, params->CFL*dx/a);
  }

#ifndef _NO_MPI
double dtTmp = dt;
MPI_Allreduce(&dtTmp, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif

  params->dt = dt;
}


void solver::timeStepA(int step)
{
#pragma omp parallel for
  for (int i=0; i<nVerts; i++) {
    //if (Geo->v2b[i]) continue;
    for (int j=0; j<nFields; j++) {
      U(i,j) = U0(i,j) - RKa[step]*params->dt*divF(step,i,j)/vol[i];
    }
  }
}


void solver::timeStepB(int step)
{
#pragma omp parallel for
  for (int i=0; i<nVerts; i++) {
    //if (Geo->v2b[i]) continue;
    for (int j=0; j<nFields; j++) {
      U(i,j) -= RKb[step]*params->dt*divF(step,i,j)/vol[i];
    }
  }
}

void solver::copyU_U0(void)
{
  U0 = U;
}

void solver::copyU0_U(void)
{
  U = U0;
}

void solver::readRestartFile(void) {

  ifstream dataFile;
  dataFile.precision(15);

  // Get the file name & open the file
  char fileNameC[50];
  string fileName = params->dataFileName;
#ifndef _NO_MPI
  /* --- All processors write their solution to their own .vtu file --- */
  sprintf(fileNameC,"%s_%.09d/%s_%.09d_%d.vtu",&fileName[0],params->restartIter,&fileName[0],params->restartIter,params->rank);
#else
  sprintf(fileNameC,"%s_%.09d.vtu",&fileName[0],params->restartIter);
#endif

  if (params->rank==0) cout << "Solver: Restarting from " << fileNameC << endl;

  dataFile.open(fileNameC);

  if (!dataFile.is_open())
    FatalError("Cannont open restart file.");

  // Find the start of the UnstructuredData region
  bool found = false;
  string str;
  while (getline(dataFile,str)) {
    stringstream ss;
    ss.str(str);
    ss >> str;
    if (str.compare("<UnstructuredGrid>")==0) {
      found = true;
      break;
    }
  }

  if (!found)
    FatalError("Cannot fine UnstructuredData tag in restart file.");

  dataFile.close();

  if (params->rank==0) cout << "Solver: Done reading restart file." << endl;

  // Finish setting up MPI faces
  if (params->rank==0 && params->nproc>1) cout << "MPIFace: Setting up MPI face communications." << endl;

}

void solver::initializeSolution()
{
  if (params->rank==0) cout << "Solver: Initializing Solution... " << flush;

  for (uint i=0; i<nVerts; i++) {
    switch (params->icType) {
    case 0: // Uniform Flow
      U(i,0) = params->rhoBound;
      U(i,1) = params->rhoBound*params->uBound;
      U(i,2) = params->rhoBound*params->vBound;
      if (nDims == 3)
        U(i,3) = params->rhoBound*params->wBound;
      U(i,nDims+1) = params->pBound/(params->gamma-1) + 0.5*params->rhoBound*params->Uinf*params->Uinf;
      break;

    case 1: // Isentropic Vortex test case
      // Grab from HiFiLES or Flurry
      break;

    case 2: // Fake shock tube
      if (xv(i,2) > 0) {
        U(i,0) = 10;
        U(i,1) = 0;
        U(i,2) = 0;
        U(i,3) = 0;
        U(i,4) = 1000;
      }
      else {
        U(i,0) = .1;
        U(i,1) = 0;
        U(i,2) = 0;
        U(i,3) = 0;
        U(i,4) = 10;
      }
      break;

    }
  }

  /* If running a moving-mesh case and using CFL-based time-stepping,
   * calc initial dt for grid velocity calculation */
  /* If running a moving-mesh case and using CFL-based time-stepping,
   * calc initial dt for grid velocity calculation */
  //if ( (params->motion!=0 || params->slipPenalty==1) && params->dtType == 1 ) {
  if (params->dtType == 1)
    calcDt();

  if (params->rank == 0) cout << "done." << endl;
}

/* ---- Overset-Grid Functions ---- */

void solver::setupOverset(void)
{
  //setupOversetData();

  // Give TIOGA a pointer to this solver for access to callback functions
  //tg->setcallback(this);

  //Geo->updateOversetConnectivity();
}

void solver::callDataUpdateTIOGA(void)
{
  //tg->dataUpdate(params->nFields,U.getData(),0);
}
