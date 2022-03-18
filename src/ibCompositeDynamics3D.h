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
 * composite dynamics that implement the method of Noble and Torczinsky
 */

#ifndef IB_COMPOSITE_DYNAMICS_H_LBDEM
#define IB_COMPOSITE_DYNAMICS_H_LBDEM

#include "ibDynamicsParticleData.h"

namespace plb {

  template<typename T, template<typename U> class Descriptor>
  class IBcompositeDynamics : public CompositeDynamics<T,Descriptor>
  {
  public:

    IBcompositeDynamics(Dynamics<T,Descriptor>* baseDynamics_,
                        bool automaticPrepareCollision_ = true);
    IBcompositeDynamics(const IBcompositeDynamics &orig);
    IBcompositeDynamics(HierarchicUnserializer &unserializer);
    ~IBcompositeDynamics();

    virtual IBcompositeDynamics<T,Descriptor>* clone() const;    
    
    virtual int getId() const;

    virtual void serialize(HierarchicSerializer& serializer) const;
    virtual void unserialize(HierarchicUnserializer& unserializer);

    virtual void prepareCollision(Cell<T,Descriptor>& cell);
    virtual void collide(Cell<T,Descriptor>& cell,  BlockStatistics& statistics);

    virtual void defineVelocity(Cell<T,Descriptor>& cell, 
                                Array<T,Descriptor<T>::d> const& u);

    IBdynamicsParticleData<T,Descriptor> particleData;
  private:
    static int id;  

    // performance tweaks
    // // have static variables for some temporary quantities
    // // to save memory
    static Array<T,Descriptor<T>::q> fEqSolid;
    static Array<T,Descriptor<T>::q> fEq;
    static Array<T,Descriptor<T>::q> fPre;
  };

  template<typename T, template<typename U> class Descriptor>
  class PartialBBTRTdynamics : public BaseTRTdynamics<T,Descriptor>{
  public:
    // inherit constructors
    using BaseTRTdynamics<T, Descriptor>::BaseTRTdynamics;

    //PartialBBTRTdynamics();
    PartialBBTRTdynamics(const PartialBBTRTdynamics &orig);
    PartialBBTRTdynamics(HierarchicUnserializer &unserializer);
    //~PartialBBTRTdynamics();

    /// Clone the object on its dynamic type.
    PartialBBTRTdynamics<T,Descriptor>* clone() const override;

    /// Return a unique ID for this class.
    int getId() const override;

    virtual void serialize(HierarchicSerializer& serializer) const;

    virtual void unserialize(HierarchicUnserializer& unserializer);

    /* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics& statistics_) override;

    /// Implementation of the collision step, with imposed macroscopic variables
    void collideExternal(Cell<T,Descriptor>& cell, T rhoBar,
                                 Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& stat) override;

    /// Compute equilibrium distribution function
    T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const override;

    /// if desired, initialize interior of sphere with sphere velocity
    void defineVelocity(Cell<T,Descriptor>& cell, 
                        Array<T,Descriptor<T>::d> const& u);

    IBdynamicsParticleData<T,Descriptor> particleData;

  private:
    static int id;

    // performance tweaks
    // // have static variables for some temporary quantities
    // // to save memory
    static Array<T,Descriptor<T>::q> fEqSolid;
    static Array<T,Descriptor<T>::q> fEq;
    static Array<T,Descriptor<T>::q> fPre;
  };

}; /* namespace plb */

#include "ibCompositeDynamics3D.hh"

#endif /* IB_COMPOSITE_DYNAMICS3D_H_LBDEM */
