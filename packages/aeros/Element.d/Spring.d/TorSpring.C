#include <Element.d/Spring.d/TorSpring.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>
#include <Corotational.d/SpringCorotator.h>
#include <Math.d/FullSquareMatrix.h>

extern "C"      {
void    _FORTRAN(mstf11)(double*,double&,double&,double&,double&,double&,
                         double&,double&,double&,double&,double&,double&,
                         double&);
}

TorSpring::TorSpring(int* nodenums)
{
        nn[0] = nodenums[0];
}

Element *
TorSpring::clone()
{
 return new TorSpring(*this);
}

void
TorSpring::renum(const int *table)
{
  nn[0] = table[nn[0]];
}

void
TorSpring::renum(EleRenumMap& table)
{
  nn[0] = table[nn[0]];
}

FullSquareMatrix
TorSpring::massMatrix(const CoordSet &,double *mel,int cmflg) const
{
        FullSquareMatrix elementMassMatrix(3,mel);

	// zero the element mass matrix
	elementMassMatrix.zero();

        return elementMassMatrix;
}

FullSquareMatrix
TorSpring::stiffness(const CoordSet &, double *d, int flg) const
{
        double sc1 = prop->A;
        double sc2 = prop->E;
        double sc3 = prop->nu;

        double u1  = prop->rho;
        double u2  = prop->c;
        double u3  = prop->k;

        double v1  = prop->eh;
        double v2  = prop->P;
        double v3  = prop->Ta;

        double w1  = prop->Q;
        double w2  = prop->W;
        double w3  = prop->Ixx;

	// KHP Sometime modify this so that there is no call to FORTRAN
	// the stiffness matrix is easy to do right here in C++

	_FORTRAN(mstf11)((double*)d,sc1,sc2,sc3,u1,u2,u3,v1,v2,v3,w1,w2,w3);

        FullSquareMatrix ret(3,d);

        return ret;
}

int
TorSpring::numNodes() const
{
        return 1;
}

int *
TorSpring::nodes(int *p) const
{
        if(p == 0) p = new int[1];
        p[0] = nn[0];
        return p;
}

int
TorSpring::numDofs() const
{
        return 3;
}

int *
TorSpring::dofs(DofSetArray &dsa, int *p) const
{
        if(p == 0) p = new int[3];
      
	dsa.number(nn[0],DofSet::XYZrot,p);

        return p;
}

void
TorSpring::markDofs(DofSetArray& dsa) const
{
        dsa.mark( nn[0], DofSet::XYZrot );
}

Corotator *
TorSpring::getCorotator(CoordSet &cs, double *kel, int , int )
{
 int flag = 0;
 int numnodes = 1;
 FullSquareMatrix myStiff = stiffness(cs, kel, flag);

 return new SpringCorotator( nn, cs, numnodes, prop->A, prop->E, prop->nu,
                                myStiff);

}

