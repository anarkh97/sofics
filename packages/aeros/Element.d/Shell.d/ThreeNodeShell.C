#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <Corotational.d/GeomState.h>
#include <Corotational.d/Shell3Corotator.h>
#include <Corotational.d/utilities.h>
#include <Driver.d/PolygonSet.h>
#include <Element.d/Shell.d/ThreeNodeShell.h>
#include <Element.d/State.h>
#include <Hetero.d/InterpPoint.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/Vector.h>
#include <Math.d/matrix.h>
#include <Utils.d/Conwep.d/BlastLoading.h>
#include <Utils.d/dbg_alloca.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/pstress.h>

// tria3d   - three node shell stiffness routine
// mass8    - three node shell mass routine
// sands8   - three node shell stress routine
// trithmfr - three node shell thermal load routine 

extern int verboseFlag;
extern "C"      {
void _FORTRAN(tria3d)(int&, double*, double*, double*, double&, double&,
                      double*, double*);
void _FORTRAN(mass8)(double*, double*, double*, double*, double&,
                     double*, const int&, double*, double*, const int&,
                     double&, const int&);
void _FORTRAN(sands8)(double*, double*, double*, double&, double&, double*,
                      double*, double*, const int&, const int&, const int&, 
                      const int&, const int&, const int&, double&);
void _FORTRAN(trithmfr)(double*, double*, double*, const double&, const double&, const double&,
                        const double&, const double&, const double&, double*, const int&);
}

ThreeNodeShell::ThreeNodeShell(int* nodenums, double _w)
{
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
	nn[2] = nodenums[2];
	w = _w;
        pbc = NULL;
}

Element *
ThreeNodeShell::clone()
{
	return new ThreeNodeShell(*this);
}

void
ThreeNodeShell::renum(const int *table)
{
	nn[0] = table[nn[0]];
	nn[1] = table[nn[1]];
	nn[2] = table[nn[2]];
}

void
ThreeNodeShell::renum(EleRenumMap& table)
{
	nn[0] = table[nn[0]];
	nn[1] = table[nn[1]];
	nn[2] = table[nn[2]];
}

void
ThreeNodeShell::getVonMises(Vector& stress, Vector& weight, CoordSet &cs,
		       	    Vector& elDisp, int strInd, int surface,
                            double *ndTemps, double ylayer, double zlayer, int avgnum)
{
	weight = 1.0;
  if(strInd == -1) return;

        auto &nd1 = cs.getNode(nn[0]);
        auto &nd2 = cs.getNode(nn[1]);
        auto &nd3 = cs.getNode(nn[2]);

        double x[3], y[3], z[3], h[3];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

        // Thickness
        h[0] = h[1] = h[2] = prop->eh;
	
	// determine average nodal temperature difference relative to ambient
        double thermalStrain = 0.0;
        if(ndTemps) {
	  double dt = (ndTemps[0]+ndTemps[1]+ndTemps[2])/3 - prop->Ta;
	  thermalStrain = (prop->W)*dt;
        }

        int maxsze = 1;
        int maxstr = 7;
        int maxgus = 3;
        int elm    = 1;

        double elStress[3][7];

	// SET STRAIN FLAG IF USER WANTS STRAIN OUTPUT
	int strainFlg = 0;
	if(strInd > 6) strainFlg = 1;

       _FORTRAN(sands8)(x,y,z,prop->E,prop->nu,h,elDisp.data(),
                        (double*)elStress,
                        strainFlg,maxsze,maxstr,maxgus,elm,surface,thermalStrain);

        if(strInd < 7) {
          stress[0] = elStress[0][strInd];
          stress[1] = elStress[1][strInd];
          stress[2] = elStress[2][strInd];
        } else {
          stress[0] = elStress[0][strInd-7];
          stress[1] = elStress[1][strInd-7];
          stress[2] = elStress[2][strInd-7];
        }
}

void
ThreeNodeShell::getAllStress(FullM& stress, Vector& weight, CoordSet &cs,
                             Vector& elDisp, int strInd, int surface,
                             double *ndTemps)
{
        weight = 1.0;

        auto &nd1 = cs.getNode(nn[0]);
        auto &nd2 = cs.getNode(nn[1]);
        auto &nd3 = cs.getNode(nn[2]);

        double x[3], y[3], z[3], h[3];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

        // Thickness
        h[0] = h[1] = h[2] = prop->eh;
	
	// determine average nodal temperature difference relative to ambient
	double dt = 0.0;
	if(ndTemps)
          dt = (ndTemps[0]+ndTemps[1]+ndTemps[2])/3 - prop->Ta;
        double thermalStrain = (prop->W)*dt;

        int maxsze = 1;
        int maxstr = 7;
        int maxgus = 3;
        int elm    = 1;

        double elStress[3][7];

       _FORTRAN(sands8)(x,y,z,prop->E,prop->nu,h,elDisp.data(),
                        (double*)elStress,
                        strInd,maxsze,maxstr,maxgus,elm,surface,thermalStrain);

// Store all Stress or all Strain as defined by strInd
        int i,j;
        for (i=0; i<3; ++i) {
          for (j=0; j<6; ++j) {
            stress[i][j] = elStress[i][j];
          }
        }

// Get Element Principals
        double svec[6], pvec[3];
        for (j=0; j<6; ++j) {
          svec[j] = stress[0][j];
        }
// Convert Engineering to Tensor Strains
        if(strInd != 0) {
          svec[3] /= 2;
          svec[4] /= 2;
          svec[5] /= 2;
        }
        pstress(svec,pvec);
        for (i=0; i<3; ++i) {
          for (j=0; j<3; ++j) {
            stress[i][j+6] = pvec[j];
          }
        }
}

double
ThreeNodeShell::getMass(const CoordSet& cs) const
{
        if (prop == NULL) return 0.0;

        auto &nd1 = cs.getNode(nn[0]);
        auto &nd2 = cs.getNode(nn[1]);
        auto &nd3 = cs.getNode(nn[2]);

        double r1[3], r2[3], r3[3], v1[3], v2[3], v3[3];

        r1[0] = nd1.x; r1[1] = nd1.y; r1[2] = nd1.z;
        r2[0] = nd2.x; r2[1] = nd2.y; r2[2] = nd2.z;
        r3[0] = nd3.x; r3[1] = nd3.y; r3[2] = nd3.z;

        v1[0] = r3[0] - r1[0];
        v1[1] = r3[1] - r1[1];
        v1[2] = r3[2] - r1[2];

        v2[0] = r2[0] - r1[0];
        v2[1] = r2[1] - r1[1];
        v2[2] = r2[2] - r1[2];

        crossprod(v1, v2, v3);

        double area = 0.5*sqrt(v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2]);
        double mass = area*prop->rho*prop->eh;

        return mass;

}

void
ThreeNodeShell::getGravityForce(CoordSet& cs, double *gravityAcceleration, 
                                Vector& gravityForce, int gravflg, GeomState*)
{
        double mass = getMass(cs);
        double massPerNode = mass/3.0;

        double fx = massPerNode*gravityAcceleration[0];
        double fy = massPerNode*gravityAcceleration[1];
        double fz = massPerNode*gravityAcceleration[2];
        double mx[3],my[3],mz[3];
        int i;   
        for(i = 0; i < 3; ++i) {
          mx[i] = 0.0;
          my[i] = 0.0;
          mz[i] = 0.0;
        }

        // Lumped with no fixed-end moments
        if (gravflg == 0) {

        }
        // Consistent or lumped with fixed end moments.  Compute treating shell as 3 beams.
        else {

          auto &nd1 = cs.getNode(nn[0]);
          auto &nd2 = cs.getNode(nn[1]);
          auto &nd3 = cs.getNode(nn[2]);
          double x[3], y[3], z[3];
          x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
          x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
          x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

          double T1[3],T2[3],T3[3];
          // Vector 1 from Node 1->2
          T1[0] = x[1] - x[0];
          T1[1] = y[1] - y[0];
          T1[2] = z[1] - z[0];
          normalize( T1 );
          // Vector 2 from Node 1->3
          T2[0] = x[2] - x[0];
          T2[1] = y[2] - y[0];
          T2[2] = z[2] - z[0];
          normalize( T2 );
          // Local Z-axis as cross between V1 and V2
          crossprod( T1, T2, T3 );
          normalize( T3);

          int beam, beamnode[3][2];
          beamnode[0][0] = 0;
          beamnode[0][1] = 1;
          beamnode[1][0] = 0;
          beamnode[1][1] = 2;
          beamnode[2][0] = 1;
          beamnode[2][1] = 2;

          for(beam=0; beam<3; ++beam) {
            double length, dx, dy, dz, localg[3];
            int n1, n2;
            n1 = beamnode[beam][0];
            n2 = beamnode[beam][1];
            dx = x[n2] - x[n1];
            dy = y[n2] - y[n1];
            dz = z[n2] - z[n1];
            length = sqrt(dx*dx + dy*dy + dz*dz);
            // Local X-axis from Node 1->2
            T1[0] = x[n2] - x[n1];
            T1[1] = y[n2] - y[n1];
            T1[2] = z[n2] - z[n1];
            normalize( T1 );
            // Local Y-axis as cross between Z and X
            crossprod( T3, T1, T2 );
            normalize( T2);

            for(i = 0; i < 3; ++i)
              localg[i] = 0.0;
            for(i = 0; i < 3; ++i) {
              localg[0] += T1[i]*gravityAcceleration[i];
              localg[1] += T2[i]*gravityAcceleration[i];
              localg[2] += T3[i]*gravityAcceleration[i];
            }
            double lmy,lmz;
            if (gravflg == 2) { // consistent
              lmy = -massPerNode*localg[2]*length/12.0;
              lmz = massPerNode*localg[1]*length/12.0;
            }
            else { // lumped with fixed-end moments
              lmy = -massPerNode*localg[2]*length/16.0;
              lmz = massPerNode*localg[1]*length/16.0;
            }
            mx[n1] += ((T2[0]*lmy) + (T3[0]*lmz));
            my[n1] += ((T2[1]*lmy) + (T3[1]*lmz));
            mz[n1] += ((T2[2]*lmy) + (T3[2]*lmz));
            mx[n2] -= ((T2[0]*lmy) + (T3[0]*lmz));
            my[n2] -= ((T2[1]*lmy) + (T3[1]*lmz));
            mz[n2] -= ((T2[2]*lmy) + (T3[2]*lmz));
          }
        }
        
        gravityForce[0]  = fx;
        gravityForce[1]  = fy;
        gravityForce[2]  = fz;
        gravityForce[3]  = mx[0];
        gravityForce[4]  = my[0];
        gravityForce[5]  = mz[0];
        gravityForce[6]  = fx;
        gravityForce[7]  = fy;
        gravityForce[8]  = fz;
        gravityForce[9]  = mx[1];
        gravityForce[10] = my[1];
        gravityForce[11] = mz[1];
        gravityForce[12] = fx;
        gravityForce[13] = fy;
        gravityForce[14] = fz;
        gravityForce[15] = mx[2];
        gravityForce[16] = my[2];
        gravityForce[17] = mz[2];
}

FullSquareMatrix
ThreeNodeShell::massMatrix(const CoordSet &cs, double *mel, int cmflg) const
{
        // Check for phantom element, which has no mass
        if(prop == NULL) {
           FullSquareMatrix ret(18,mel);
	   ret.zero();
           return ret;
        }

        auto &nd1 = cs.getNode(nn[0]);
        auto &nd2 = cs.getNode(nn[1]);
        auto &nd3 = cs.getNode(nn[2]);

        double x[3], y[3], z[3], h[3];
        double *gravityAcceleration = NULL, *grvfor = NULL, totmas = 0.0;

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

        h[0] = h[1] = h[2] = prop->eh;

        int grvflg = 0, masflg = 0;

	const int numdof = 18;

       _FORTRAN(mass8)(x,y,z,h,prop->rho,(double *)mel,numdof,
               gravityAcceleration,grvfor,grvflg,totmas,masflg);

        FullSquareMatrix ret(18,mel);

        return ret;
}

FullSquareMatrix
ThreeNodeShell::stiffness(const CoordSet &cs, double *d, int flg) const
{
        // Check for phantom element, which has no stiffness
        if (prop == NULL) {
           FullSquareMatrix ret(18,d);
           ret.zero();
           return ret;
        }

        auto &nd1 = cs.getNode(nn[0]);
        auto &nd2 = cs.getNode(nn[1]);
        auto &nd3 = cs.getNode(nn[2]);

        double x[3], y[3], z[3], h[3];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

        h[0] = h[1] = h[2] = prop->eh;

        // Check for a zero thickness
        if(h[0] <= 0.0) {
          fprintf(stderr," *** ERROR: Shell element # %d has zero or negative thickness. Exiting...\n", getGlNum()+1);
          exit(-1);
        }

        _FORTRAN(tria3d)(flg, x, y, z, prop->E, prop->nu, h, (double *)d);

        FullSquareMatrix ret(18,d);
       
        return ret;
}

int
ThreeNodeShell::numNodes() const
{
 	return 3;
}

int*
ThreeNodeShell::nodes(int *p) const
{
 	if(p == 0) p = new int[3];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
	return p;
}

int
ThreeNodeShell::numDofs() const
{
 	return 18;
}

int*
ThreeNodeShell::dofs(DofSetArray &dsa, int *p) const
{
 	if(p == 0) p = new int[18];

        dsa.number(nn[0], DofSet::XYZdisp | DofSet::XYZrot, p  );
        dsa.number(nn[1], DofSet::XYZdisp | DofSet::XYZrot, p+6);
        dsa.number(nn[2], DofSet::XYZdisp | DofSet::XYZrot, p+12);

	return p;
}

void
ThreeNodeShell::markDofs(DofSetArray &dsa) const
{
        dsa.mark(nn, 3, DofSet::XYZdisp | DofSet::XYZrot);
}

Corotator *
ThreeNodeShell::getCorotator(CoordSet &cs, double *kel, int fitAlgShell, int )
{
 int flag = 0; // signals stiffness routine to keep local matrix 
 FullSquareMatrix myStiff = stiffness(cs, kel, flag);
 corot = new Shell3Corotator(nn[0], nn[1], nn[2], myStiff, fitAlgShell);
 return corot;
}

void
ThreeNodeShell::computeDisp(CoordSet &cs, State &state, const InterpPoint &ip,
                            double *res, GeomState *geomState)
{
 const double *gp = ip.xy;
 double gap = std::sqrt(ip.gap[0]*ip.gap[0]+ip.gap[1]*ip.gap[1]+ip.gap[2]*ip.gap[2]);
 if(gap < 1e-12 || geomState == NULL) {

   double xyz[3][6];
   state.getDV(nn[0], xyz[0], xyz[0]+3);
   state.getDV(nn[1], xyz[1], xyz[1]+3);
   state.getDV(nn[2], xyz[2], xyz[2]+3);

   int j;
   for(j=0; j<6; ++j)
      res[j] = (1.0-gp[0]-gp[1]) * xyz[0][j] + gp[0]*xyz[1][j] + gp[1]*xyz[2][j];
 
 }
 else {

   double xyz[3][12];
   double w[3];
   state.getDVRot(nn[0], xyz[0], xyz[0]+6);
   state.getDVRot(nn[1], xyz[1], xyz[1]+6);
   state.getDVRot(nn[2], xyz[2], xyz[2]+6);

   // Compute the translation and linear velocity at the interpolation point
   int j;
   for(j=0; j<3; ++j)
      res[j] = (1.0-gp[0]-gp[1]) * xyz[0][j] + gp[0]*xyz[1][j] + gp[1]*xyz[2][j];
   for(j=0; j<3; ++j)
      res[j+3] = (1.0-gp[0]-gp[1]) * xyz[0][j+6] + gp[0]*xyz[1][j+6] + gp[1]*xyz[2][j+6];

   // Compute the angular velocity at the interpolation point
   for(j=0; j<3; ++j)
      w[j] = (1.0-gp[0]-gp[1]) * xyz[0][j+9] + gp[0]*xyz[1][j+9] + gp[1]*xyz[2][j+9];

   // Get Nodes original coordinates (C0 configuration)
   Node &node1 = cs.getNode(nn[0]);
   Node &node2 = cs.getNode(nn[1]);
   Node &node3 = cs.getNode(nn[2]);

   // Get Nodes current coordinates (C0n configuration)
   NodeState &ns1 = (*geomState)[nn[0]];
   NodeState &ns2 = (*geomState)[nn[1]];
   NodeState &ns3 = (*geomState)[nn[2]];   

   double xl0[3][3], xln[3][3], t0[3][3], t0n[3][3], vld[18];

   // C0    = initial configuration
   // C0n   = nth configuration
   // xl0   = C0 local coordinates
   // xln   = C0n local coordinates
   // t0    = transformation matrix to C0
   // t0n   = transformation matrix to C0n
   // vld   = local deformation vector

   // Compute transformation matrices and 
   // extract deformational displacement from C0 to C0n configurations

   corot->extractDefDisp(node1,node2,node3, ns1,ns2,ns3, xl0,xln, t0,t0n, vld );

   double locGap[3], gn[3], dv[3];
   // initial configuration gap vector in C0 local coordinates
   for(int i=0; i<3; ++i) {
     locGap[i] = 0;
     for(j=0; j<3; ++j)
       locGap[i] += t0[i][j]*ip.gap[j];
   }
   // nth configuration gap vector in global coordinates
   for(int i=0; i<3; ++i) {
     gn[i] = 0;
     for(j=0; j<3; ++j)
       gn[i] += t0n[j][i]*locGap[j];
   }

   // correction to displacement of fluid node due to gap = t0n^T*t0*ip.gap - ip.gap
   for(int i=0; i<3; ++i) {
     res[i] += (gn[i] - ip.gap[i]);
   }

   // correction to velocity of fluid node due to gap = w x gn
   crossprod(w, gn, dv);
   for(int i=0; i<3; ++i) { 
     res[i+3] += dv[i];
   }

   // TODO: 
   // 1. correction to displacement due to deformational rotation (small)
   // 2. correction to fluid load due to gap (this will involve adding moments in ThreeNodeShell::getFlLoad)
 }
}

void
ThreeNodeShell::printInfo(CoordSet&, State &state, double gp[2])
{
 double xyz[3][6];
 state.getDV(nn[0], xyz[0], xyz[0]+3);
 state.getDV(nn[1], xyz[1], xyz[1]+3);
 state.getDV(nn[2], xyz[2], xyz[2]+3);

 fprintf(stderr, "INfo %e %e %e %e %e %e %e %e %e\n",
     xyz[0][0], xyz[0][1], xyz[0][2],
     xyz[1][0], xyz[1][1], xyz[1][2],
     xyz[2][0], xyz[2][1], xyz[2][2]);
}

void
ThreeNodeShell::getFlLoad(CoordSet &, const InterpPoint &ip, double *flF, 
                          double *resF, GeomState *gs)
{
 const double *gp = ip.xy;
 int i;
 for(i = 0; i < 3; ++i) {
   resF[i+3]  = resF[i+9] = resF[i+15] = 0.0;
   resF[i]    = (1.0-gp[0]-gp[1]) * flF[i];
   resF[6+i]  = gp[0] * flF[i];
   resF[12+i] = gp[1] * flF[i];
  }
}

int
ThreeNodeShell::getTopNumber() const
{
  return 108;
}

void
ThreeNodeShell::computePressureForce(CoordSet& cs, Vector& elPressureForce,
                                     GeomState *geomState, int cflg, double time) {
     double pressure = pbc->val;
     // Check if Conwep is being used. If so, add the pressure from the blast loading function.
     if (pbc->conwep && pbc->conwepswitch) {
       double* CurrentElementNodePositions = (double*) dbg_alloca(sizeof(double)*3*4);
       int Offset;
       for(int i = 0; i < 4; ++i) {
         Offset = i*3;
         if (i==3) {
           CurrentElementNodePositions[Offset+0] = cs[nn[2]]->x;
           CurrentElementNodePositions[Offset+1] = cs[nn[2]]->y;
           CurrentElementNodePositions[Offset+2] = cs[nn[2]]->z;
         }
         else {
           CurrentElementNodePositions[Offset+0] = cs[nn[i]]->x;
           CurrentElementNodePositions[Offset+1] = cs[nn[i]]->y;
           CurrentElementNodePositions[Offset+2] = cs[nn[i]]->z;
         }
       }
       pressure += BlastLoading::ComputeShellPressureLoad(CurrentElementNodePositions, time, *(pbc->conwep));
     }
     double px = 0.0;
     double py = 0.0;
     double pz = 0.0;

     double mx[3],my[3],mz[3];
     int i;
     for(i=0; i<3; ++i) {
       mx[i]=0.0;
       my[i]=0.0;
       mz[i]=0.0;
     }

     // Compute area of shell
     auto &nd1 = cs.getNode(nn[0]);
     auto &nd2 = cs.getNode(nn[1]);
     auto &nd3 = cs.getNode(nn[2]);

     double r1[3], r2[3], r3[3], v1[3], v2[3], normal[3];

     r1[0] = nd1.x; r1[1] = nd1.y; r1[2] = nd1.z;
     r2[0] = nd2.x; r2[1] = nd2.y; r2[2] = nd2.z;
     r3[0] = nd3.x; r3[1] = nd3.y; r3[2] = nd3.z;

     v1[0] = r3[0] - r1[0];
     v1[1] = r3[1] - r1[1];
     v1[2] = r3[2] - r1[2];

     v2[0] = r2[0] - r1[0];
     v2[1] = r2[1] - r1[1];
     v2[2] = r2[2] - r1[2];

     // Compute normal to shell using vector cross product rule
     crossprod(v2, v1, normal);

     double magnitude = sqrt(normal[0]*normal[0] + normal[1]*normal[1] 
                                                  + normal[2]*normal[2]);
     double area = 0.5*magnitude;

     // compute pressure force per node
     double pressureForce = pressure * area / 3.0;

     double xn[3][3];
       for(i=0; i<3; ++i) {
         xn[0][i] = r1[i];
         xn[1][i] = r2[i];
         xn[2][i] = r3[i];
       }

     double mcoef = 8; // 12 ???

     if(!geomState) {

       // compute unit normal to shell surface

       normal[0] /= magnitude;
       normal[1] /= magnitude;
       normal[2] /= magnitude;

       px = pressureForce*normal[0];
       py = pressureForce*normal[1];
       pz = pressureForce*normal[2];

       if (cflg == 1) {

         int beam, beamnode[3][2];
         beamnode[0][0] = 0;
         beamnode[0][1] = 1;
         beamnode[1][0] = 0;
         beamnode[1][1] = 2;
         beamnode[2][0] = 1;
         beamnode[2][1] = 2;

         for(beam=0; beam<3; ++beam) {
           double length, dx, dy, dz;
           int n1, n2;
           n1 = beamnode[beam][0];
           n2 = beamnode[beam][1];
           dx = xn[n2][0] - xn[n1][0];
           dy = xn[n2][1] - xn[n1][1];
           dz = xn[n2][2] - xn[n1][2];
           length = sqrt(dx*dx + dy*dy + dz*dz);
           // Local X-axis from Node 1->2
           for(i=0; i<3; i++ ) v1[i] = xn[n2][i] - xn[n1][i];
           normalize( v1 );
           // Local Y-axis as cross between Z and X
           crossprod( normal, v1, v2 );
           normalize( v2 );

           double lmy = -pressureForce*length/mcoef; // mcoef = 8 seems to work better for cylinder shell model

           mx[n1] += (v2[0]*lmy);
           my[n1] += (v2[1]*lmy);
           mz[n1] += (v2[2]*lmy);
           mx[n2] -= (v2[0]*lmy);
           my[n2] -= (v2[1]*lmy);
           mz[n2] -= (v2[2]*lmy);
         }
       }

       elPressureForce[0]  = px;
       elPressureForce[1]  = py;
       elPressureForce[2]  = pz;
       elPressureForce[3]  = mx[0];
       elPressureForce[4]  = my[0];
       elPressureForce[5]  = mz[0];

       elPressureForce[6]  = px;
       elPressureForce[7]  = py;
       elPressureForce[8]  = pz;
       elPressureForce[9]  = mx[1];
       elPressureForce[10] = my[1];
       elPressureForce[11] = mz[1];

       elPressureForce[12] = px;
       elPressureForce[13] = py;
       elPressureForce[14] = pz;
       elPressureForce[15] = mx[2];
       elPressureForce[16] = my[2];
       elPressureForce[17] = mz[2];
     }

     //if local coordinates are needed for nonlinear analysis
     if (geomState) {

        //compute centroid
        double xc0[3];
        double t0[3][3],xl0[3][3];
        int nod;
        for (i=0;i<3;i++)
        xc0[i] = ( xn[0][i] + xn[1][i] + xn[2][i] )/3.0;

        // Compute t0 transformation matrix with x axis along side 1-2 
        for( i=0; i<3; i++ ) t0[0][i] = xn[1][i] - xn[0][i];
        normalize( t0[0] );

        // local y axis
        for( i=0; i<3; i++ ) t0[1][i] = xn[2][i] - xn[0][i];
        crossprod( t0[0], t0[1], t0[2] );
        normalize( t0[2] );

        // local z axis
        crossprod( t0[2], t0[0], t0[1] );
        normalize( t0[1] );

        // Compute local coordinates of undeformed element
        for( nod=0; nod<3; nod++) {
             for( i=0; i<3; i++ ) {
                 xl0[nod][i] = t0[i][0]*(xn[nod][0] - xc0[0])
                              +t0[i][1]*(xn[nod][1] - xc0[1])
                              +t0[i][2]*(xn[nod][2] - xc0[2]);
                }
             }

        elPressureForce[0]  = 0;
        elPressureForce[1]  = 0;
        elPressureForce[2]  = pressureForce;
        elPressureForce[3]  = static_cast<double>(cflg)*pressureForce*( xl0[2][1] - xl0[0][1]
                                                  +xl0[1][1] - xl0[0][1])/mcoef;
        elPressureForce[4]  = static_cast<double>(cflg)*pressureForce*( xl0[0][0] - xl0[2][0]
                                                  +xl0[0][0] - xl0[1][0])/mcoef;
        elPressureForce[5]  = 0;

        elPressureForce[6]  = 0;
        elPressureForce[7]  = 0;
        elPressureForce[8]  = pressureForce;
        elPressureForce[9]  = static_cast<double>(cflg)*pressureForce*( xl0[0][1] - xl0[1][1]
                                                  +xl0[2][1] - xl0[1][1])/mcoef;
        elPressureForce[10] = static_cast<double>(cflg)*pressureForce*( xl0[1][0] - xl0[0][0]
                                                  +xl0[1][0] - xl0[2][0])/mcoef;
        elPressureForce[11] = 0;

        elPressureForce[12] = 0;
        elPressureForce[13] = 0;
        elPressureForce[14] = pressureForce;
        elPressureForce[15] = static_cast<double>(cflg)*pressureForce*( xl0[1][1] - xl0[2][1]
                                                  +xl0[0][1] - xl0[2][1])/mcoef;
        elPressureForce[16] = static_cast<double>(cflg)*pressureForce*( xl0[2][0] - xl0[1][0]
                                                  +xl0[2][0] - xl0[0][0])/mcoef;
        elPressureForce[17] = 0;

     }

}

/* dec
int
ThreeNodeShell::facelist(PolygonSet &pgs, int *flist)
{
 if(flist != 0) {
    flist[0] = pgs.addLine2( nn[0], nn[1] );
    flist[1] = pgs.addLine2( nn[1], nn[2] );
    flist[2] = pgs.addLine2( nn[2], nn[0] );
 }
 return 3 ;
}

// end dec */

//------------------------------------------------------------------------

void
ThreeNodeShell::getThermalForce(CoordSet& cs, Vector& ndTemps,
				Vector &elThermalForce, int glflag, 
				GeomState *)
{
   if(prop == NULL) {
     elThermalForce.zero();
     return;
   }

   double Tref = prop->Ta;
   
   double x[3], y[3], z[3];
   
   // get nodal temperatures and relate to reference temperature
   double meant = 0.0; // determine the average nodal temperature
   for(int i=0; i<3; i++)
     meant += ndTemps[i]/3;

   meant -= Tref;

   // Compute Node's coordinates
   auto &nd1 = cs.getNode(nn[0]);
   auto &nd2 = cs.getNode(nn[1]);
   auto &nd3 = cs.getNode(nn[2]);
      
   x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
   x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
   x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;      
   
   double alpha = 1.5;
     
   //Call FORTRAN routine to determine elemental thermal force vector from
   //membrane effects -- returned in global coordinates if glflag = 0
   _FORTRAN(trithmfr)(x, y, z, meant, prop->W, prop->E, prop->nu,
		      prop->eh, alpha, (double *)elThermalForce.data(), 
                      glflag);
}

#ifdef USE_EIGEN3
#include <Element.d/FelippaShell.d/ShellElementTemplate.hpp>
#include <Element.d/FelippaShell.d/EffMembraneTriangle.hpp>
#include <Element.d/FelippaShell.d/AndesBendingTriangle.hpp>

typedef ShellElementTemplate<double,EffMembraneTriangle,AndesBendingTriangle> Impl;

double
ThreeNodeShell::getMassThicknessSensitivity(CoordSet& cs)
{
  if(prop == NULL) return 0.0;
  else return getMass(cs)/prop->eh;
}

void
ThreeNodeShell::getWeightNodalCoordinateSensitivity(Vector &dwdx, CoordSet& cs, double *gravityAcceleration)
{
  if(prop == NULL || gravityAcceleration == NULL) {
    dwdx.zeroAll();
    return;
  }

  double x[3] = { cs[nn[0]]->x, cs[nn[1]]->x, cs[nn[2]]->x };
  double y[3] = { cs[nn[0]]->y, cs[nn[1]]->y, cs[nn[2]]->y };
  double z[3] = { cs[nn[0]]->z, cs[nn[1]]->z, cs[nn[2]]->z };

  double rhoh = prop->rho*prop->eh;

  Impl::andesmsWRTcoord(glNum+1, x, y, z, dwdx.data(), rhoh);

  double gravAccNorm = sqrt(gravityAcceleration[0]*gravityAcceleration[0] +
                            gravityAcceleration[1]*gravityAcceleration[1] +
                            gravityAcceleration[2]*gravityAcceleration[2]);

  dwdx *= gravAccNorm;

}

void
ThreeNodeShell::getGravityForceThicknessSensitivity(CoordSet& cs, double *gravityAcceleration,
                                                       Vector& dGfdthick, int gravflg, GeomState *geomState)
{
  if(prop == NULL) {
    dGfdthick.zero();
    return;
  }

  getGravityForce(cs, gravityAcceleration, dGfdthick, gravflg, geomState);

  dGfdthick /= prop->eh;
}

void
ThreeNodeShell::getGravityForceNodalCoordinateSensitivity(CoordSet& cs, double *gravityAcceleration,
                                                          GenFullM<double> &dGfdx, int gravflg, GeomState*)
{
  if(prop == NULL) {
    dGfdx.zero();
    return;
  }

  double x[3] = { cs[nn[0]]->x, cs[nn[1]]->x, cs[nn[2]]->x };
  double y[3] = { cs[nn[0]]->y, cs[nn[1]]->y, cs[nn[2]]->y };
  double z[3] = { cs[nn[0]]->z, cs[nn[1]]->z, cs[nn[2]]->z };

  double rhoh = prop->rho*prop->eh;

  Impl::andesgfWRTcoord(glNum+1, x, y, z, dGfdx.data(),
                        gravityAcceleration, gravflg, rhoh);
}

void
ThreeNodeShell::getStiffnessThicknessSensitivity(CoordSet &cs, FullSquareMatrix &dStiffdThick, int flg)
{
  if(prop == NULL) {
    dStiffdThick.zero();
    return;
  }

  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);

  double x[3], y[3], z[3];

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

  ShellMaterial<double> *mat = new ShellMaterialType0<double>(prop->E, prop->eh, prop->nu,
                                                              prop->rho, prop->Ta, prop->W); 

  double disp[18]; for(int i=0; i<18; ++i) disp[i] = 0;
  Impl::andesstfWRTthic(glNum+1, dStiffdThick.data(), (double*)NULL, prop->nu,
                        x, y, z, disp, 0, mat, flg);
  delete mat;
}

void
ThreeNodeShell::getStiffnessNodalCoordinateSensitivity(FullSquareMatrix *&dStiffdx, CoordSet &cs)
{
  if(prop == NULL) {
    for(int i=0; i<9; ++i) dStiffdx[i].zero();
    return;
  }

  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);

  double x[3], y[3], z[3];

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

  double *data[9];
  for(int i=0; i<9; ++i) data[i] = dStiffdx[i].data();

  int flg = 1;

  Impl::andesstfWRTcoord(glNum+1, data, prop->E, prop->nu,
                         prop->rho, prop->eh, prop->Ta, prop->W,
                         (double*)NULL, x, y, z, 0, (double*)NULL,
                         flg);
}

void
ThreeNodeShell::getVonMisesThicknessSensitivity(Vector &dStdThick, Vector &weight, CoordSet &cs,
                                                Vector &elDisp, int, int surface, double *ndTemps,
                                                int avgnum, double, double)
{
  weight = 1.0;

  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);

  double x[3], y[3], z[3];

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

  ShellMaterial<double> *mat = new ShellMaterialType0<double>(prop->E, prop->eh, prop->nu, 
                                                              prop->rho, prop->Ta, prop->W);
  int sflg = 0;

  Impl::andesvmsWRTthic(glNum+1, prop->nu, x, y, z, elDisp.data(),
                        dStdThick.getData(), 0, mat, surface, sflg,
                        ndTemps);
  delete mat;
}

void
ThreeNodeShell::getVonMisesNodalCoordinateSensitivity(GenFullM<double> &dStdx, Vector &weight, CoordSet &cs,
                                                      Vector &elDisp, int, int surface, double *ndTemps,
                                                      int avgnum, double, double)
{
  weight = 1.0;

  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);

  double x[3], y[3], z[3];

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

  int sflg = 0;

  Impl::andesvmsWRTcoord(glNum+1, prop->E, prop->nu, prop->rho,
                         prop->eh, prop->Ta, prop->W, (double*)NULL,
                         x, y, z, elDisp.data(), dStdx.getData(),
                         0, (double*)NULL, surface, sflg, ndTemps);
}

void
ThreeNodeShell::getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, GenFullM<double> *dDispDisp, CoordSet &cs,
                                                   Vector &elDisp, int, int surface, double *ndTemps,
                                                   int avgnum, double, double)
{
  weight = 1.0;

  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);

  double x[3], y[3], z[3];

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

  double* disp = elDisp.data();

  ShellMaterial<double> *mat = new ShellMaterialType0<double>(prop->E, prop->eh, prop->nu, 
                                                              prop->rho, prop->Ta, prop->W);
  int sflg = 0;

  Impl::andesvmsWRTdisp(glNum+1, prop->nu, x, y, z, elDisp.data(),
                        dStdDisp.getData(), 0, mat, surface, sflg,
                        ndTemps);
  if(dDispDisp) dStdDisp ^= (*dDispDisp);
  delete mat;
}

int ThreeNodeShell::getMassType() const { return 0; }

#endif
