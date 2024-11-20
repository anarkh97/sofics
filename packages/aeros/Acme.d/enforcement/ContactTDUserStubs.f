c----------------------------------------------------------------------
c user_initialize_model is called once to allow for set-up and 
c parameter calculation through the userquery_x functions
c----------------------------------------------------------------------
      subroutine user_initialize_model(enf, id, rdata, idata, istat)
      integer idata(*), id, istat
      real*8  rdata(*)
      istat = 0
      return
      end
      
c----------------------------------------------------------------------
c user_initialize_time_step is called at the beginning of each timestep
c to allow for resetting state data, etc through the userset_x functions.
c----------------------------------------------------------------------
      subroutine user_initialize_time_step(enf, id, rdata, idata, istat)
      integer idata(*), id, istat
      real*8  rdata(*)
      istat = 0
      return
      end
      
c----------------------------------------------------------------------
c user_initialize_node_state_data is called to initialize state
c data passed in "sdata".
c----------------------------------------------------------------------
      subroutine user_init_node_state_data(enf, id, rdata, idata,
     *                                           sdata, istat)
      integer idata(*), id, istat
      real*8  rdata(*), sdata(*)
      istat = 0
      return
      end
      
c----------------------------------------------------------------------
c user_interaction_type returns the classification that the interaction
c currently falls in. The types are enumerated in ContactTDEnfModel.h 
c and are repeated here :
c                 TDEM_FRICTIONLESS=0,    // return value 0
c                 TDEM_FRICTIONAL,        // 1 
c                 TDEM_TIED,              // 2
c                 TDEM_ADHESIVE,          // 3
c                 TDEM_SPRINGY,           // 4
c                 TDEM_ADHESIVE_FRICTION, // 5
c The control what constraints the enforcement routine tries to impose
c and what predictors will be used to generate a trial force for the
c user-defined material model. A FRICTIONLESS interaction will only
c have impenetrability enforced. A FRICTIONAL interaction will have
c impenetrability enforced and, in addition, a trial force due to 
c trying to prevent tangential motion i.e. stick. A TIED interaction
c will enforce that the node and surface are coincident at the same 
c (convected) point for all time. An ADHESIVE interaction is like a 
c FRICTIONLESS interaction with the addition of a trial force for 
c seperation that attempts to bring the node and surface together. 
c A SPRINGY interaction is similar to a TIED interaction in that the
c node and surface are always bound together at the same "contact
c point" but in this case with a penalty-like trial force instead
c of a Lagrange multiplier. An ADHESIVE_FRICTION interaction is 
c the combination of the FRICTIONAL and the ADHESIVE types so that 
c the trial force is composed of a normal force (in the case of 
c separation) and a tangential force that attempts to stick the 
c node and surface together.
c----------------------------------------------------------------------
      subroutine user_interaction_type(enf, ni, id, rdata, itype, 
     *                                idata, index, istat)
      integer istat, id, itype, index, idata(*)
      real*8  rdata(*)
      istat = 0
      return
      end

c----------------------------------------------------------------------
c user_is_active simply takes the current value of the gap ("gap" in
c the arguement list and the material parameters in rdata and idata to
c determine whether this interaction is "active" i.e. will generate a 
c force based on proximity of the two surfaces.
c----------------------------------------------------------------------
      subroutine user_is_active(enf, ni, id, rdata, idata,
     *                          itype, index, gap, istat)
      integer istat, id, itype, index, idata(*)
      real*8  gap, rdata(*)
      istat = 0
      return
      end
  
c----------------------------------------------------------------------
c user_limit_force takes the value of the trial force passed in by 
c "force" and returns the constitutively determined value using the 
c the arguements:
c gap : the normal closest-point  distance between the node and surface, 
c rel_disp : the absolute vectorial distance between node & contact point 
c slip : the incremental vectorial distance between node & contact point
c normal : outward unit normal to the surface
c dt : timestep
c area : tributary area of the node
c The last arguement is 0 if the force was left as is, and nonzero if
c it was modified.
c----------------------------------------------------------------------
      subroutine user_limit_force(enf, ni, id, rdata, idata, itype, 
     *                            index, gap, rel_disp, slip, 
     *                            normal, dt, area, force, 
     *                            istat)
      integer istat, id, itype, index, idata(*)
      real*8 gap, dt, area
      real*8 rel_disp(*), slip(*), normal(*), force(*), rdata(*)
      force(1) = 0.0
      force(2) = 0.0
      force(3) = 0.0
      istat = 1
      return
      end
