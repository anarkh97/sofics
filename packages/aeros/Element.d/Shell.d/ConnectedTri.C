#include <Math.d/FullSquareMatrix.h>
#include <Element.d/Shell.d/ConnectedTri.h>
#include <Utils.d/dofset.h>
#include <Hetero.d/FlExchange.h>
#include <Element.d/State.h>

ConnectedTri::ConnectedTri(int *n)
{
  for(int i = 0; i < 4; ++i)
    nn[i] = n[i];
}

int
ConnectedTri::numNodes() const
{ 
  return 4;
}

void
ConnectedTri::renum(const int *table)
{
	nn[0] = table[nn[0]];
	nn[1] = table[nn[1]];
	nn[2] = table[nn[2]];
	nn[3] = table[nn[3]];
}

void
ConnectedTri::renum(EleRenumMap& table)
{
	nn[0] = table[nn[0]];
	nn[1] = table[nn[1]];
	nn[2] = table[nn[2]];
	nn[3] = table[nn[3]];
}

int *
ConnectedTri::nodes(int *p) const
{ 
  if(p == 0) p = new int[4];
  for(int i = 0; i < 4; ++i)
    p[i] = nn[i];
  return p;
}

void
ConnectedTri::markDofs(DofSetArray &dsa) const
{
        dsa.mark(nn+3, 1,  DofSet::XYZdisp | DofSet::XYZrot);
}

int *
ConnectedTri::dofs(DofSetArray &dsa, int *p) const
{
 	if(p == 0) p = new int[6];

        dsa.number(nn[3],DofSet::XYZdisp | DofSet::XYZrot, p  );

	return p;
}

int
ConnectedTri::numDofs() const
{ 
  return 6;
}

FullSquareMatrix
ConnectedTri::massMatrix(const CoordSet &cs,double *mel,int cmflg) const
{
           FullSquareMatrix ret(6,mel);

	   ret.zero();

           return ret;
}

FullSquareMatrix
ConnectedTri::stiffness(const CoordSet &cs,double *d,int cmflg) const
{
           FullSquareMatrix ret(6,d);

	   ret.zero();

           return ret;
}

void
ConnectedTri::getFlLoad(CoordSet &cs, const InterpPoint &ip, double *flF, 
                          double *resF, GeomState *gs)
{
 const double *gp = ip.xy;
 int i;
 Node *n1 = cs[nn[0]];
 Node *n2 = cs[nn[1]];
 Node *n3 = cs[nn[2]];
 Node *n4 = cs[nn[3]];
 double xyzRef[3] = { n4->x, n4->y, n4->z };
 for(i = 0; i < 3; ++i)
   resF[i] = flF[i];
 double xyz[3] = {
   (1.0-gp[0]-gp[1]) * n1->x + gp[0]*n2->x+gp[1]*n3->x,
   (1.0-gp[0]-gp[1]) * n1->y + gp[0]*n2->y+gp[1]*n3->y,
   (1.0-gp[0]-gp[1]) * n1->z + gp[0]*n2->z+gp[1]*n3->z };
 resF[3] = (xyz[1]-xyzRef[1])*flF[2] -  (xyz[2]-xyzRef[2])*flF[1];
 resF[4] = (xyz[2]-xyzRef[2])*flF[0] -  (xyz[0]-xyzRef[0])*flF[2];
 resF[5] = (xyz[0]-xyzRef[0])*flF[1] -  (xyz[1]-xyzRef[1])*flF[0];
}

void
ConnectedTri::computeDisp(CoordSet&cs, State &state, const InterpPoint &ip,
                            double *res, GeomState *gs)
{
 const double *gp = ip.xy;
 double d[6], dv[6]; 
 state.getDVRot(nn[3], d, dv);
 Node *n1 = cs[nn[0]];
 Node *n2 = cs[nn[1]];
 Node *n3 = cs[nn[2]];
 Node *n4 = cs[nn[3]];
 double xyzRef[3] = { n4->x, n4->y, n4->z };
 double xyz[3] = {
   (1.0-gp[0]-gp[1]) * n1->x + gp[0]*n2->x+gp[1]*n3->x,
   (1.0-gp[0]-gp[1]) * n1->y + gp[0]*n2->y+gp[1]*n3->y,
   (1.0-gp[0]-gp[1]) * n1->z + gp[0]*n2->z+gp[1]*n3->z };

 res[0] = d[0] + (xyz[2]-xyzRef[2])*d[4] - (xyz[1]-xyzRef[1])*d[5];
 res[1] = d[1] + (xyz[0]-xyzRef[0])*d[5] - (xyz[2]-xyzRef[2])*d[3];
 res[2] = d[2] + (xyz[1]-xyzRef[1])*d[3] - (xyz[0]-xyzRef[0])*d[4];
 res[3] = dv[0] + (xyz[2]-xyzRef[2])*dv[4] - (xyz[1]-xyzRef[1])*dv[5];
 res[4] = dv[1] + (xyz[0]-xyzRef[0])*dv[5] - (xyz[2]-xyzRef[2])*dv[3];
 res[5] = dv[2] + (xyz[1]-xyzRef[1])*dv[3] - (xyz[0]-xyzRef[0])*dv[4];
}

int
ConnectedTri::getTopNumber() const
{
 return 4;
}
