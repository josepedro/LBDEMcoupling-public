# LBDEMcoupling

* [About](#about)
* [Requirements](#requirements)
* [Setting Up a Simulation](#setting_up)
* [Implicit Assumptions, Known Issues](#assumptions)
* [Gallery](#gallery)
* [Documented Compiler Switches](#compilerswitches)
* [Citing LBDEMcoupling](#citing)
* [References](#references)
* [License and Copyright](#license)

<a name="about"></a>
## About

LBDEMcoupling is a coupling between the Lattice Boltzmann (LB) library
Palabos (http://www.palabos.org) and the Discrete Element Method code
LIGGGHTS® (http://www.ligggghts.com). It implements the model of Noble
and Torczinsky [[1]](#ref1) for resolved coupling between particles
and a fluid phase.

<a name="requirements"></a>
## Requirements

### LIGGGHTS installation

1. First of all, you need to move the files from LBDEMcoupling to LIGGGHTS source code directory:
         fix_lb_coupling_onetoone.cpp
         fix_lb_coupling_onetoone.h

2. Clone repository:
```console
git clone git@github.com:CFDEMproject/LIGGGHTS-PUBLIC.git && cd LIGGGHTS-PUBLIC && git pull
```

3. Compile the code in order to generate the executable lmp_auto:
```console
cd src && sed -i 's/#AUTOINSTALL_VTK = "OFF"/AUTOINSTALL_VTK = "ON"/g' MAKE/Makefile.user_default && make auto
```

4. Add VTK library in the path LD_LIBRARY_PATH:
```console
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<full_path>/LIGGGHTS-PUBLIC/lib/vtk/install/lib
```

5. Compile the shared library:
```console
make makeshlib && make -f Makefile.shlib auto
```

6. Add LIGGGHTS libraries in the path LD_LIBRARY_PATH:
```console
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<full_path>/LIGGGHTS-PUBLIC/src
```

<a name="setting_up"></a>
## Setting Up a Simulation

### CMakeLists.txt settings

In order to run the examples you just need to modify the following fields with the correct paths:
```console
set(PALABOS_ROOT "/home/pedro/singularity/singularity-ce-3.8.1/workspace/palabos")
set(LIGGGHTS_ROOT "/home/pedro/singularity/singularity-ce-3.8.1/workspace/LIGGGHTS-PUBLIC")
set(LBDEM_ROOT "/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public")
```

And then you just need to go inside 'build' folder and proceed with:
```console
cmake .. && make -j
```

And finally you have the binary created.

### Coupling API

As for now, the coupling API is rather limited, and several tasks are coded explicitly in each case. This might change in the future. For now, the most important API elements are

* the header file `plb_ib.h`: an aggregate header that provides all functionality implemented in LBDEMcoupling
* `class PhysUnits3D`: a conversion utility to convert between Palabos simulation units and physical units. See [this excellent document by Jonas Latt](http://wiki.palabos.org/_media/howtos:lbunits.pdf) for how unit conversion is performed, and the file `${LBDEM_ROOT}/src/physunits.h` for the detailed interface.
* `class LIGGGHTScouplingWrapper`: a wrapper for an instance of LIGGGHTS. Important methods
  * `LIGGGHTScouplingWrapper(char **argv, MPI_COMM communicator)`: The constructor needs the arguments array, and a valid MPI communicator. The latter is usually supplied via the Palabos call `global::mpi().getGlobalCommunicator()`. Currently, no command line argument passing from Palabos to LIGGGHTS is implemented.
  * `void execFile(char* const fname)` executes the commands given in a file
  * `void execCommand(char* const command)` executes a single LIGGGHTS command. There is also a version of this function accepting a `std::stringstream` for convenience.
  * `void run(plint nStep)` and `void runUpto(plint nStep)`: equivalent to the LIGGGHTS commands `run` and `run upto`
  * `int getNumParticles()` returns the total number of particles in the domain
  * `void setVariable(char const *name, double value)` and `void setVariable(char const *name, std::string &value)` define a LIGGGHTS variable just like the `variable` command. The latter creates a `variable string` (see LIGGGHTS docu).
* wrapper functions for the actual coupling
  * `void setSpheresOnLattice(MultiBlockLattice3D &lattice,LIGGGHTScouplingWrapper &wrapper,PhysUnits3D &units, bool initWithVel)` writes particle information to the lattice. If `initWithVel` is `true`, the velocity at all cells covered by a particle is set to the particle velocity. This can be useful for initialization purposes.
  * `void getForcesFromLattice(MultiBlockLattice3D &lattice,LIGGGHTScouplingWrapper &wrapper,PhysUnits3D &units)` collects the hydrodynamic forces on the particles.


<a name="gallery"></a>
## Gallery

### Settling spheres

<img src="doc/img/settling.png" alt="10000 settling spheres">

### Square Channel with Particles

<img src="doc/img/showcaseRectChannel.png">


<a name="assumptions"></a>
## Implicit Assumptions, Known Issues

### Sphere size

The code implicitly assumes that all your particles are larger than
four grid spacings and smaller than half the smallest extent of any
partition.

### Multiple Definition Errors

If you come across a multiple definition error during compilation,
please consult ![this
thread: http://www.palabos.org/forum/read.php?11,6746,7581](http://www.palabos.org/forum/read.php?11,6746,7581) 
in the Palabos forum (copy the URL if the link does not work).

<a name="compilerswitches"></a>
## Documented Compiler Switches

LBDEM_USE_MULTISPHERE switches on support for the (non-public)
multisphere model of LIGGGHTS

<a name="citing"></a>
## Citing LBDEMcoupling

If you found LBDEMcoupling useful (which we hope you do), you might want to cite it. If so, please cite the following conference paper:

Seil, P., & Pirker, S. (2017). LBDEMcoupling: Open-Source Power for Fluid-Particle Systems. In *Proceedings of the 7th International Conference on Discrete Element Methods* (pp. 679-686). Springer Singapore.

<a name="references"></a>
## References

<a name="ref1">[1]</a> Noble, D. R., & Torczynski, J. R. (1998). A
lattice-Boltzmann method for partially saturated computational
cells. *International Journal of Modern Physics C*, 9(08), 1189-1201.

<a name="license"></a>
## License and Copyright

(c) Johannes Kepler University Linz, Austria

released under the GPLv3

main author: Philippe Seil (philippe.seil@jku.at)

LIGGGHTS® is a registered trade mark of DCS Computing GmbH, the
producer of the LIGGGHTS® software.

Palabos is a registered trademark of FlowKit Ltd., the developer of the Palabos software.
