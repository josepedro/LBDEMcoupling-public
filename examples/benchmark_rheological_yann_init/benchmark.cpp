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

/*
This case was used as a benchmark. It consists of a tank with fluid 
where the upper third is filled with particles that then settle under
gravity. The command line parameter v_inf is an estimated settling velocity.
If chosen too large, the timestep will be very small. If chosen too small,
the simulation will crash.

A few working parameter sets:

./benchmark 0.1 5 0.1 1e-4 0.15 0.02 your/out/dir
results in ~70 particles on a 51x51x101 grid. Feasible on a desktop computer.
Particle Reynolds number ~150

./benchmark 0.06 5 0.1 1e-4 0.10 0.02 your/out/dir
~340 particles on a 84x84x164 grid - might need 1-2h or so on a desktop computer.
Particle Reynolds number ~90


./benchmark 0.02 5 0.1 1e-5 0.15 0.02 your/out/dir 
gives you ~10k particles on a 251x251x501 grid - computationally demanding!
Particle Reynolds number ~300

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

#include "utilsFunctionals3D.h"
#include "liggghtsCouplingWrapper.h"
#include "latticeDecomposition.h"

using namespace plb;
using namespace std;

typedef double T;

const T omega = 1.0;

#define DESCRIPTOR descriptors::D3Q27Descriptor
#define DYNAMICS IBcompositeDynamics<T,DESCRIPTOR>(new CompleteRegularizedBGKdynamics<T,DESCRIPTOR>(omega))
//#define DESCRIPTOR descriptors::D3Q19Descriptor
//#define DYNAMICS IBcompositeDynamics<T,DESCRIPTOR>(new RRdynamics<T,DESCRIPTOR>(omega))

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

    const T d_part = 2*0.0001;
    const T dx = 2*0.00001;
    //const T dx = 4*0.00001;
    const T dt = 6*0.0000001;
    const T lx = 0.01, ly = 0.00102, lz = 0.002;
    const T rho_f = 1000;
    const T r_ = d_part/2.;
    const T dynamic_viscosity = 0.1; // (Pa.s or N.s/m^2 or kg/(m.s)) m0
    const T nu_f =  dynamic_viscosity/rho_f; // kinematic viscosity (m^2/s)
    const T rho_s = 2*1000;
    //const T strain_rate = 100; // (s^-1) this comes from Yann thesis
    //const T strain_rate = 100.0/10.0; // (s^-1) this comes from an attempt
    const T strain_rate = 5000.0; // (s^-1) this comes from Yann thesis highest value
    const plint N = d_part/dx; // N is the number of grid points per particle diameter
    const T v_inf = ly*strain_rate;
    const T uMax = v_inf*(dt/dx);
    const std::string outDir = "outDir/";
    const T Rep = rho_f*r_*r_*strain_rate/dynamic_viscosity; // Particle Reynolds number
    const T volume_fraction = 0.36945129606215965; // non-newtonian -> 0.25 <= volume_fraction < volume_fraction_critic (0.5)

    std::string lbOutDir(outDir), demOutDir(outDir);
    lbOutDir.append("tmp/"); demOutDir.append("post/");
    global::directories().setOutputDir(lbOutDir);

    LiggghtsCouplingWrapper wrapper(argv,global::mpi().getGlobalCommunicator());

    // particle size and volume fraction are handed over to LIGGGHTS 
    // as variables (see LIGGGHTS docu for details)
    wrapper.setVariable("r_part",d_part/2.0);
    
    wrapper.execFile("in.lbdem");
    
    PhysUnits3D<T> units(2.*r_,v_inf,nu_f,lx,ly,lz,N,uMax,rho_f);

    IncomprFlowParam<T> parameters(units.getLbParam());
    
    ///// Initialize the diagonal relaxation matrix from the .xml file.
    // Q19 and Q27 do not have the same number of relaxation parameters!
    Array<T, DESCRIPTOR<T>::numRelaxationTimes> allOmega;
    allOmega[0] = omega;     // relaxation of M200 and cyclic permutations
    allOmega[1] = omega;     // relaxation of M110 and cyclic permutations
    allOmega[2] = omega;    // relaxation of M210 and cyclic permutations
    allOmega[3] = omega;    // relaxation of M220 and cyclic permutations
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

    defineDynamics(lattice, lattice.getBoundingBox(), new DYNAMICS);    
    
    const T maxT = 30.0/strain_rate; // Yann's thesis = 20.0/strain_rate
    //const T vtkT = 0.006/8.0;
    const T vtkT = 0.001;
    const T logT = 0.0000001;

    const plint maxSteps = units.getLbSteps(maxT);
    const plint vtkSteps = max<plint>(units.getLbSteps(vtkT),1);
    const plint logSteps = max<plint>(units.getLbSteps(logT),1);

    writeLogFile(parameters, "Rheological Study from Yann's Thesis");

    lattice.periodicity().toggle(0,true); // periodic boundary condition on axis X
    lattice.periodicity().toggle(2,true); // periodic boundary condition on axis Z

    // set strain rate
    //T vel = units.getLbVel(0.0);
    T vel = units.getLbVel(v_inf);
    Box3D lid(0, nx - 1, ny - 1, ny - 1, 0, nz - 1);
    Box3D bottom(0, nx - 1, 0, 0, 0, nz - 1);
    // Yann's approach
    /*
    OnLatticeBoundaryCondition3D<T, DESCRIPTOR>* boundaryCondition = 
      createLocalBoundaryCondition3D<T, DESCRIPTOR>();
    boundaryCondition->addVelocityBoundary1P(bottom, lattice);
    boundaryCondition->addVelocityBoundary1N(lid, lattice);
    setBoundaryVelocity(lattice, lid, Array<T, 3>(vel, 0., 0.));
    */
    // VelocityBounceBack's approach
    VelocityBounceBack<T,DESCRIPTOR> vbbDynamics_bottom = VelocityBounceBack<T,DESCRIPTOR>(1.0, Array<T, 3>(0., 0., 0.));
    defineDynamics(lattice, bottom, vbbDynamics_bottom.clone());
    VelocityBounceBack<T,DESCRIPTOR> vbbDynamics_lid = VelocityBounceBack<T,DESCRIPTOR>(1.0, Array<T, 3>(vel, 0., 0.));
    defineDynamics(lattice, lid, vbbDynamics_lid.clone());



/*
    initializeAtEquilibrium( lattice, lattice.getBoundingBox(), 
                             CoutteProfile<T>(vel,nx,ny,nz,0,true,2) );
*/  

    // -----------------------------------------------------------
    lattice.initialize();
    T dt_phys = units.getPhysTime(1);
    plint demSubsteps = 10;
    T dt_dem = dt_phys/(T)demSubsteps;

    pcout << "------------------------------\n"
          << "omega: " << omega << "\n" 
          << "dt_phys: " << dt_phys << "\n"
          << "maxT: " << maxT << " | maxSteps: " << maxSteps << "\n"
          << "v_inf: " << v_inf << "\n"
          << "Re : " << parameters.getRe() << "\n"
          << "Particle Reynolds number: " << Rep << "\n"
          << "vtkT: " << vtkT << " | vtkSteps: " << vtkSteps << "\n"
          << "grid size: " << nx << " " << ny << " " << nz << "\n"
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
    
    // ---------------------------------------
    // Loop over main time iteration.
    for (plint iT=0; iT<=maxSteps; ++iT) {
      bool initWithVel = false;
      setSpheresOnLattice(lattice,wrapper,units,initWithVel);

      if(iT%vtkSteps == 0 /*&& iT > 0*/) { // LIGGGHTS does not write at timestep 0
        writeVTK(lattice,parameters,units,iT);
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
                                                    0))) )/units.getPhysLength(1);;
          ofile_velocity_gradient_x << "," << setprecision(10) << velocity_gradient;
        }
        ofile_velocity_gradient_x << std::endl;
      }

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
    pcout << " ********************** \n"
          << "total time: " << totaltime
          << " calculating at " << totalmlups << " MLU/s" << std::endl;

    ofile_energy.close();
    ofile_velocity_x.close();
    ofile_velocity_gradient_x.close();
}
