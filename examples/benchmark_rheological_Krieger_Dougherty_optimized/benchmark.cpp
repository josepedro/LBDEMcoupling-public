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

#include "palabos3D.h"
#include "palabos3D.hh"

#include "plb_ib.h"

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>

// necessary LAMMPS/LIGGGHTS includes
#include "lammps.h"
#include "input.h"
#include "library.h"
#include "library_cfd_coupling.h"

#include "periodicPressureFunctionals3D.h"
#include "liggghtsCouplingWrapper.h"
#include "latticeDecomposition.h"

using namespace plb;
using namespace std;

typedef double T;

#define DESCRIPTOR descriptors::D3Q27Descriptor
#define DYNAMICS IBcompositeDynamics<T,DESCRIPTOR>(new CompleteRegularizedBGKdynamics<T,DESCRIPTOR>(parameters.getOmega()))
//#define DESCRIPTOR descriptors::D3Q19Descriptor
//#define DYNAMICS IBcompositeDynamics<T,DESCRIPTOR>(new RRdynamics<T,DESCRIPTOR>(parameters.getOmega()))

void writeVTK(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
              IncomprFlowParam<T> const& parameters,
              PhysUnits3D<T> const& units, plint iter)
{
  
  T p_fact = units.getPhysForce(1)/pow(units.getPhysLength(1),2)/3.;
  
  std::string fname(createFileName("vtk", iter, 6));
  
  VtkImageOutput3D<T> vtkOut(fname, units.getPhysLength(1));
  vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", units.getPhysVel(1));  
  vtkOut.writeData<float>(*computeDensity(lattice), "density",units.getPhysDensity(1)); 
  
  MultiScalarField3D<T> p(*computeDensity(lattice));
  subtractInPlace(p,1.);
  vtkOut.writeData<float>(p,"pressure",p_fact ); 
  
  IBscalarQuantity sf = SolidFraction;
  applyProcessingFunctional(new GetScalarQuantityFromDynamicsFunctional<T,DESCRIPTOR,T>(sf),
                            lattice.getBoundingBox(),lattice,p);

  vtkOut.writeData<float>(p,"solidfraction",1. ); 
  
 
  pcout << "wrote " << fname << std::endl;
}

int main(int argc, char* argv[]) {

    plbInit(&argc, &argv);
    
    const T lx = 2., ly = 2., lz = 2.;
    const T strain_rate = 0.00005;
    const T velocity_imposed = strain_rate*ly;
    const T d_part = 0.1; // 1 - particle diameter
    const plint N = 10; // 2 - number of grid points per particle diameter
    // const plint N = 5;
    const T v_frac = 0.5; // 3 - the solid fraction in the insertion region
    const T nu_f = 1e-4; // 4 - kinematic viscosity (m^2/s)
    // const T v_inf = 0.15;
    // const T v_inf = 0.15; // 5 - velocity of the particle (v_inf is an estimated settling velocity)
    // const T v_inf = velocity_imposed; // try_2
    const T v_inf = velocity_imposed/2.0; // try_3
    const T optimization_factor = 2.0/1000.0;
    
    const T uMax = 0.02 * optimization_factor; /* 6 - maybe this is the ULB 
    - uMax is the maximum velocity in LB units, proportional to the Mach number -
    Characteristic velocity in lattice units (proportional to Mach number). */
    const std::string outDir = "outDir/"; /* 7 - directory where your output shall be stored. It must
    contain the subdirectories outDir/post and outDir/tmp. */

    const T rho_f = 1000;
    const T r_ = d_part/2.;
    const T rho_s = 1100.;
    const T m = r_*r_*r_*4./3.*3.14*rho_s;
    const T Rep = d_part*velocity_imposed/nu_f;

    std::string lbOutDir(outDir), demOutDir(outDir);
    lbOutDir.append("tmp/"); demOutDir.append("post/");
    global::directories().setOutputDir(lbOutDir);

    LiggghtsCouplingWrapper wrapper(argv,global::mpi().getGlobalCommunicator());

    // particle size and volume fraction are handed over to LIGGGHTS 
    // as variables (see LIGGGHTS docu for details)
    wrapper.setVariable("r_part",d_part/2);
    wrapper.setVariable("v_frac",v_frac);
    wrapper.execFile("in.lbdem");
    
    PhysUnits3D<T> units(2.*r_,v_inf,nu_f,lx,ly,lz,N,uMax,rho_f);

    IncomprFlowParam<T> parameters(units.getLbParam());
    
    ///// Initialize the diagonal relaxation matrix from the .xml file.
    // Q19 and Q27 do not have the same number of relaxation parameters!
    Array<T, DESCRIPTOR<T>::numRelaxationTimes> allOmega;
    allOmega[0] = parameters.getOmega();     // relaxation of M200 and cyclic permutations
    allOmega[1] = parameters.getOmega();     // relaxation of M110 and cyclic permutations
    allOmega[2] = parameters.getOmega();    // relaxation of M210 and cyclic permutations
    allOmega[3] = parameters.getOmega();    // relaxation of M220 and cyclic permutations
    allOmega[4] = 1.0; // relaxation of bulk moment (M200 + M020 + M002)
    RRdynamics<T,DESCRIPTOR>::allOmega = allOmega;

    plint nx = parameters.getNx(), ny = parameters.getNy(), nz = parameters.getNz();

    // get lattice decomposition from LIGGGHTS and create lattice according to parallelization
    // given in the LIGGGHTS input script
    LatticeDecomposition lDec(parameters.getNx(),parameters.getNy(),parameters.getNz(),
                              wrapper.lmp);
    SparseBlockStructure3D blockStructure = lDec.getBlockDistribution();
    ExplicitThreadAttribution* threadAttribution = lDec.getThreadAttribution();
    plint envelopeWidth = 1;

    MultiBlockLattice3D<T, DESCRIPTOR> 
      lattice (MultiBlockManagement3D (blockStructure, threadAttribution, envelopeWidth ),
               defaultMultiBlockPolicy3D().getBlockCommunicator(),
               defaultMultiBlockPolicy3D().getCombinedStatistics(),
               defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),
               new DYNAMICS );

    defineDynamics(lattice,lattice.getBoundingBox(),new DYNAMICS);
    
    const T maxT = ceil(3.*lz/v_inf);
    // const T vtkT = 0.1;
    // const T vtkT = 0.000000001;
    // const T vtkT = 1.4;
    const T vtkT = 14.0;
    const T logT = 0.0000001;

    // const plint maxSteps = units.getLbSteps(maxT);
    // const plint maxSteps = 10000000000000;
    // const plint maxSteps = 10;
    const plint maxSteps = 100;
    const plint vtkSteps = max<plint>(units.getLbSteps(vtkT),1);
    const plint logSteps = max<plint>(units.getLbSteps(logT),1);

    writeLogFile(parameters, "sedimenting spheres benchmark");

    lattice.periodicity().toggle(0,true); // periodic boundary condition on axis X
    lattice.periodicity().toggle(2,true); // periodic boundary condition on axis Z
    Box3D inlet(0,0,0,ny-1,0,nz-1), outlet(nx-1,nx-1,0,ny-1,0,nz-1);

    // set strain rate
    T vel = units.getLbVel(velocity_imposed);
    //T vel = units.getLbVel(v_inf);
    Box3D lid(0, nx - 1, ny - 1, ny - 1, 0, nz - 1);
    Box3D bottom(0, nx - 1, 0, 0, 0, nz - 1);
    // Yann's approach
    OnLatticeBoundaryCondition3D<T, DESCRIPTOR>* boundaryCondition = 
      createLocalBoundaryCondition3D<T, DESCRIPTOR>();
    boundaryCondition->addVelocityBoundary1P(bottom, lattice);
    boundaryCondition->addVelocityBoundary1N(lid, lattice);
    setBoundaryVelocity(lattice, lid, Array<T, 3>(vel, 0., 0.));
    setBoundaryVelocity(lattice, bottom, Array<T, 3>(0., 0., 0.));

    // initializeAtEquilibrium( lattice, lattice.getBoundingBox(), 
    //                          PressureGradient<T>(rhoHi,rhoLo,nz,0) );
    /*
    initializeAtEquilibrium( lattice, lattice.getBoundingBox(), 
                             PoiseuilleProfileAndPressureGradient<T>(rhoHi,rhoLo,uMax,nx,ny,nz,0) );
    */
    lattice.initialize();
    T dt_phys = units.getPhysTime(1); // dt_dem = t_step = 0.000266667
    plint demSubsteps = 10; // dt_phys = dt_dem*demSubsteps; 0.000266667*10 = 0.00266667 = dt_phys
    T dt_dem = dt_phys/(T)demSubsteps;


    pcout << "------------------------------\n"
          << "omega: " << parameters.getOmega() << "\n"
          << "dt_phys: " << dt_phys << "\n"
          << "dx_phys: " << units.getPhysLength(1) << "\n"
          << "maxT: " << maxT << " | maxSteps: " << maxSteps << "\n"
          << "v_inf: " << v_inf << "\n"
          << "Rep : " << Rep << "\n"
          << "vtkT: " << vtkT << " | vtkSteps: " << vtkSteps << "\n"
          << "grid size: " << nx << " " << ny << " " << nz << "\n"
          << "strain_rate: " << strain_rate << "\n"
          << "velocity_imposed: " << velocity_imposed << "\n"
          << "uMax (ULB) in physical units: " << units.getPhysVel(uMax) << "\n"
          << "uMax (ULB) in lattice units: " << uMax << "\n"
          << "------------------------------" << std::endl;
    // set timestep and output directory
    wrapper.setVariable("t_step",dt_dem);
    wrapper.setVariable("dmp_stp",vtkSteps*demSubsteps);
    wrapper.setVariable("dmp_dir",demOutDir);


    wrapper.execFile("in2.lbdem");
    wrapper.runUpto(demSubsteps-1);

    clock_t start = clock();
    clock_t loop = clock();
    clock_t end = clock();

    // Outputing information
    // Energy
    std::string fname_energy(global::directories().getOutputDir() + "lattice_average_energy.csv");
    plb_ofstream ofile_energy(fname_energy.c_str());
    ofile_energy << "iT" << "," << "average_energy" << std::endl;

    // Physical Velocity
    std::string fname_velocity_x(global::directories().getOutputDir() + "lattice_average_velocity_x.csv");
    plb_ofstream ofile_velocity_x(fname_velocity_x.c_str());
    ofile_velocity_x << "iT";
    for (plint iY = 0; iY < ny; ++iY) {
      ofile_velocity_x << "," << iY;
    }
    ofile_velocity_x << std::endl;

    // Physical Velocity Gradient
    std::string fname_velocity_gradient_x(
      global::directories().getOutputDir() + "lattice_average_velocity_gradient_x.csv");
    plb_ofstream ofile_velocity_gradient_x(fname_velocity_gradient_x.c_str());
    ofile_velocity_gradient_x << "iT";
    for (plint iY = 0; iY < ny - 1; ++iY) {
      ofile_velocity_gradient_x << "," << iY;
    }
    ofile_velocity_gradient_x << std::endl;
    // --------------------------------------

    // Relative Apparent Viscosity
    std::string fname_relative_apparent_viscosity(
      global::directories().getOutputDir() + "relative_apparent_viscosity.csv");
    plb_ofstream ofile_relative_apparent_viscosity(fname_relative_apparent_viscosity.c_str());
    ofile_relative_apparent_viscosity << "iT";
    for (plint iY = 0; iY < ny - 1; ++iY) {
      ofile_relative_apparent_viscosity << "," << iY;
    }
    ofile_relative_apparent_viscosity << std::endl;
    // --------------------------------------

    // Relative Apparent Viscosity Averaged
    T measurements_relative_apparent_viscosity_averaged = 0;
    T averaged_velocities[ny]; averaged_velocities[ny - 1] = 0;
    std::string fname_relative_apparent_viscosity_averaged(
      global::directories().getOutputDir() + "relative_apparent_viscosity_averaged.csv");
    plb_ofstream ofile_relative_apparent_viscosity_averaged(
      fname_relative_apparent_viscosity_averaged.c_str());
    ofile_relative_apparent_viscosity_averaged << "iT";
    for (plint iY = 0; iY < ny - 1; ++iY) {
      ofile_relative_apparent_viscosity_averaged << "," << iY;
      averaged_velocities[iY] = 0;
    }
    ofile_relative_apparent_viscosity_averaged << std::endl;
    // --------------------------------------

    // Loop over main time iteration.
    for (plint iT=0; iT<=maxSteps; ++iT) {

      bool initWithVel = false;
      setSpheresOnLattice(lattice,wrapper,units,initWithVel);
      

      if(iT%vtkSteps == 0 /*&& iT > 0*/) { // LIGGGHTS does not write at timestep 0
        // writeVTK(lattice,parameters,units,iT);
        // writing files here
        // Energy
        ofile_energy << iT << "," << setprecision(10) << getStoredAverageEnergy<T>(lattice) << std::endl;

        // Physical Velocity
        ofile_velocity_x << iT;
        for (plint iY = 0; iY < ny; ++iY) {
          ofile_velocity_x << "," << setprecision(10) << units.getPhysVel(computeAverage(*computeVelocityComponent(lattice,
                                                         Box3D(0, nx - 1, iY, iY, 0, nz - 1),
                                                         0)));
        }
        ofile_velocity_x << std::endl;

        // Physical Velocity Gradient
        ofile_velocity_gradient_x << iT;
        for (plint iY = 0; iY < ny - 1; ++iY) {
          T velocity_gradient = ( units.getPhysVel(computeAverage(*computeVelocityComponent(lattice,
                                                    Box3D(0, nx - 1, iY + 1, iY + 1, 0, nz - 1),
                                                    0))) - 
                                  units.getPhysVel(computeAverage(*computeVelocityComponent(lattice,
                                                    Box3D(0, nx - 1, iY, iY, 0, nz - 1),
                                                    0))) )/units.getPhysLength(1);
          ofile_velocity_gradient_x << "," << setprecision(10) << velocity_gradient;
        }
        ofile_velocity_gradient_x << std::endl;

        // Relative Apparent Viscosity
        ofile_relative_apparent_viscosity << iT;
        for (plint iY = 0; iY < ny - 1; ++iY) {
          T velocity_gradient = ( units.getPhysVel(computeAverage(*computeVelocityComponent(lattice,
                                                    Box3D(0, nx - 1, iY + 1, iY + 1, 0, nz - 1),
                                                    0))) - 
                                  units.getPhysVel(computeAverage(*computeVelocityComponent(lattice,
                                                    Box3D(0, nx - 1, iY, iY, 0, nz - 1),
                                                    0))) )/units.getPhysLength(1);
          T relative_apparent_viscosity = velocity_gradient/strain_rate;
          ofile_relative_apparent_viscosity << "," << setprecision(10) << relative_apparent_viscosity;
        }
        ofile_relative_apparent_viscosity << std::endl;

        // Relative Apparent Viscosity Averaged
        for (plint iY = 0; iY < ny; ++iY) {
          averaged_velocities[iY] += units.getPhysVel(computeAverage(*computeVelocityComponent(lattice,
                                                    Box3D(0, nx - 1, iY, iY, 0, nz - 1),
                                                    0)));          
        }
        measurements_relative_apparent_viscosity_averaged++;
        ofile_relative_apparent_viscosity_averaged << iT;
        for (plint iY = 0; iY < ny - 1; ++iY) {
          T velocity_gradient = averaged_velocities[iY + 1] - averaged_velocities[0];
          velocity_gradient /= measurements_relative_apparent_viscosity_averaged; // averaging in fact
          T distance = (iY + 1)*units.getPhysLength(1); //delta position in physical units ((iY + 1) is the number of cells (LB units))
          T meanStrainRate = velocity_gradient / distance; //measured strain rate in physical units
          T relative_apparent_viscosity = meanStrainRate/strain_rate;
          ofile_relative_apparent_viscosity_averaged << "," << setprecision(10) << relative_apparent_viscosity;
        }
        ofile_relative_apparent_viscosity_averaged << std::endl;
      }

      T rhoAvgIn = computeAverageDensity(lattice,inlet);
      T rhoAvgOut = computeAverageDensity(lattice,outlet);

      lattice.collideAndStream();

      getForcesFromLattice(lattice,wrapper,units);

      wrapper.run(demSubsteps);


      if(iT%logSteps == 0){
        end = clock();
        T time = difftime(end,loop)/((T)CLOCKS_PER_SEC);
        T totaltime = difftime(end,start)/((T)CLOCKS_PER_SEC);
        T mlups = ((T) (lattice.getNx()*lattice.getNy()*lattice.getNz()*logSteps))/time/1e6;
        pcout << "time: " << time << " " ;
        pcout << "calculating at " << mlups << " MLU/s"
              << " | total time running: " << totaltime << std::endl;
        loop = clock();
      }
    }
    T totaltime = difftime(end,start)/((T)CLOCKS_PER_SEC);
    T totalmlups = ((T) (lattice.getNx()*lattice.getNy()*lattice.getNz()*(maxSteps+1)))/totaltime/1e6;
    T totalPhysicalSimulationTime = maxSteps*dt_phys;
    pcout << " ********************** \n"
          << "total time: " << totaltime << "\n"
	  << "totalPhysicalSimulationTime: " << totalPhysicalSimulationTime << " seconds" << "\n"
          << "calculating at " << totalmlups << " MLU/s" << std::endl;

    ofile_energy.close();
    ofile_velocity_x.close();
    ofile_velocity_gradient_x.close();
    ofile_relative_apparent_viscosity.close();
    ofile_relative_apparent_viscosity_averaged.close();
}
