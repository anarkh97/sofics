Piston is a simple spring–mass system simulator designed for studying basic
fluid–structure interaction (FSI) problems, where a spring–mass system
interacts with a compressible fluid flow.

To capture the two-way coupled FSI phenomenon, Piston communicates with the
fluid dynamics solver M2C. The coupling is implemented using a staggered
partitioned scheme, in which the structural and fluid solvers operate on
time grids offset by half a time step. At each time step, Piston provides
the displacements and velocities of the fluid–structure interface, and in
return receives the corresponding loads computed by M2C.
