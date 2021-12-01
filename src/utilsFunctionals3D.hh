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

namespace plb {

  template<typename T>
  CoutteProfile<T>::CoutteProfile(T uMax_, 
    plint nx_, plint ny_, plint nz_,
    plint dimension_, bool positiveDirection_,
    plint verticalDirection_)
    : uMax(uMax_), nx(nx_), ny(ny_), nz(nz_), dimension(dimension_), 
    positiveDirection(positiveDirection_), verticalDirection(verticalDirection_)
  { }
  template<typename T>
  void CoutteProfile<T>::operator() 
    (plint iX, plint iY, plint iZ, T& density, Array<T,3>& velocity) const
  {
    velocity.resetToZero();
    density = 1.0;
    T h = 0;
    T hMax = 0;
    switch(verticalDirection){
    case 0:
      h = iX;
      hMax = nx - 1;
      break;
    case 1:
      h = iY;
      hMax = ny - 1;
      break;
    case 2:
      h = iZ;
      hMax = nz - 1;
      break;
    }

    switch(dimension){
    case 0:
      if (positiveDirection) {
        velocity[0] = uMax*h/hMax;
        velocity[1] = 0.;
        velocity[2] = 0.;
      } else {
        velocity[0] = uMax*(hMax - h)/hMax;
        velocity[1] = 0.;
        velocity[2] = 0.;
      }
      break;
    case 1:
      if (positiveDirection) {
        velocity[0] = 0.;
        velocity[1] = uMax*h/hMax;
        velocity[2] = 0.;
      } else {
        velocity[0] = 0.;
        velocity[1] = uMax*(hMax - h)/hMax;
        velocity[2] = 0.;
      }
      break;
    case 2:
      if (positiveDirection) {
        velocity[0] = 0.;
        velocity[1] = 0.;
        velocity[2] = uMax*h/hMax;
      } else {
        velocity[0] = 0.;
        velocity[1] = 0.;
        velocity[2] = uMax*(hMax - h)/hMax;
      }
      break;
    }
  }
};
