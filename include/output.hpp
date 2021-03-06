/*!
 * \file output.hpp
 * \brief Header file for restart & visualization data output functions
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
#pragma once

#include "global.hpp"
#include "solver.hpp"
#include "geo.hpp"

/*! Write solution to file (of type params->plotType) */
void writeData(solver *Solver, input *params);

/*! Write solution data to a CSV file. */
void writeCSV(solver *Solver, input *params);

/*! Write solution data to a Paraview .vtu file. */
void writeParaview(solver *Solver, geo* Geo, input *params);

/*! Compute the residual and print to the screen. */
void writeResidual(solver *Solver, input *params);

/*! Write a Tecplot mesh file compatible with TIOGA's testTioga FORTRAN interface */
void writeMeshTecplot(solver* Solver, input* params);
