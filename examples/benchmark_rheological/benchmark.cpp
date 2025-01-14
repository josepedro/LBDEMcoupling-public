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

    T uMax;

    plint N;
    
    T nu_f,d_part,v_frac, v_inf;

    std::string outDir;
    
    try {
        global::argv(1).read(d_part);
        global::argv(2).read(N);
        global::argv(3).read(v_frac);
        global::argv(4).read(nu_f);
        global::argv(5).read(v_inf);
        global::argv(6).read(uMax);
        global::argv(7).read(outDir);
    } catch(PlbIOException& exception) {
        pcout << exception.what() << endl;
        pcout << "Command line arguments:\n";
        pcout << "1 : d_part\n";
        pcout << "2 : N per particle diameter\n";
        pcout << "3 : particle volume fraction\n";
        pcout << "4 : nu_fluid\n";
        pcout << "5 : estimated v_inf\n";
        pcout << "6 : uMax\n";
        pcout << "7 : outDir\n";
        exit(1);
    }

    std::string lbOutDir(outDir), demOutDir(outDir);
    lbOutDir.append("tmp/"); demOutDir.append("post/");
    global::directories().setOutputDir(lbOutDir);

    const T rho_f = 1000;

    LiggghtsCouplingWrapper wrapper(argv,global::mpi().getGlobalCommunicator());

    // particle size and volume fraction are handed over to LIGGGHTS 
    // as variables (see LIGGGHTS docu for details)
    wrapper.setVariable("r_part",d_part/2);
    wrapper.setVariable("v_frac",v_frac);
    
    wrapper.execFile("in.lbdem");


    T g = 9.81;

    const T lx = 2., ly = 2., lz = 2.;


    T r_ = d_part/2.;
    T rho_s = 1100.;
    T m = r_*r_*r_*4./3.*3.14*rho_s;
    
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
    //const T vtkT = 0.1;
    const T vtkT = 1.4;
    const T logT = 0.0000001;

    const plint maxSteps = units.getLbSteps(maxT);
    const plint vtkSteps = max<plint>(units.getLbSteps(vtkT),1);
    const plint logSteps = max<plint>(units.getLbSteps(logT),1);

    writeLogFile(parameters, "sedimenting spheres benchmark");

    // periodic boundary condition on axis X
    lattice.periodicity().toggle(0,true);
    Box3D inlet(0,0,1,ny-2,1,nz-2), outlet(nx-1,nx-1,1,ny-2,1,nz-2);
    
    T deltaRho = 0.000001;
    // T rhoHi = 1.+0.5*deltaRho, rhoLo = 1.-0.5*deltaRho;
    T rhoHi = 1., rhoLo = 1.-deltaRho;

    // initializeAtEquilibrium( lattice, lattice.getBoundingBox(), 
    //                          PressureGradient<T>(rhoHi,rhoLo,nz,0) );
    initializeAtEquilibrium( lattice, lattice.getBoundingBox(), 
                             PoiseuilleProfileAndPressureGradient<T>(rhoHi,rhoLo,uMax,nx,ny,nz,0) );
    
    lattice.initialize();
    T dt_phys = units.getPhysTime(1);
    plint demSubsteps = 10;
    T dt_dem = dt_phys/(T)demSubsteps;


    pcout << "------------------------------\n"
          << "omega: " << parameters.getOmega() << "\n" 
          << "dt_phys: " << dt_phys << "\n"
          << "maxT: " << maxT << " | maxSteps: " << maxSteps << "\n"
          << "v_inf: " << v_inf << "\n"
          << "Re : " << parameters.getRe() << "\n"
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
    // Loop over main time iteration.
    for (plint iT=0; iT<=maxSteps; ++iT) {

      bool initWithVel = false;
      setSpheresOnLattice(lattice,wrapper,units,initWithVel);
      

      if(iT%vtkSteps == 0 && iT > 0) // LIGGGHTS does not write at timestep 0
        writeVTK(lattice,parameters,units,iT);

      T rhoAvgIn = computeAverageDensity(lattice,inlet);
      T rhoAvgOut = computeAverageDensity(lattice,outlet);

      lattice.collideAndStream();

      // applying a pressure gradient across the periodic boundary 
      applyProcessingFunctional
      (new ZhangPeriodicPressureFunctional3D<T,DESCRIPTOR>(rhoHi, 
                       rhoAvgOut,0,1), inlet, lattice);
      applyProcessingFunctional
      (new ZhangPeriodicPressureFunctional3D<T,DESCRIPTOR>(rhoLo, 
                       rhoAvgIn,0,-1), outlet, lattice);

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

}
