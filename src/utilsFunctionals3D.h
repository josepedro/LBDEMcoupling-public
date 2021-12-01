/*
 * This file is part of the LBDEMcoupling software.
 *
 * LBDEMcoupling is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Copyright 2014 Johannes Kepler University Linz
 *
 * Author: Philippe Seil (philippe.seil@jku.at)
 */

#ifndef UTILS_FUNCTIONALS_3D
#define UTILS_FUNCTIONALS_3D

namespace plb {
  
  template<typename T>
  class CoutteProfile {
  public:
    CoutteProfile(T uMax, plint nx_, plint ny_, plint nz_, 
      plint dimension_, bool positiveDirection, plint verticalDirection);
    void operator() (plint iX, plint iY, plint iZ, T& density, Array<T,3>& velocity) const;
  private:
    T uMax;
    plint nx,ny,nz,dimension,verticalDirection;
    bool positiveDirection;
  };

}; /* namespace plb */

#include "utilsFunctionals3D.hh"

#endif /* UTILS_FUNCTIONALS_3D */
