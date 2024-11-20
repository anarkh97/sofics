#include <cmath>

#include <Element.d/FluidQuad.d/BarSloshFS.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>

extern "C"      {
void    _FORTRAN(barsloshfs)(double*, double*, double&,
                             double*, const int *, const int *);
}

BarSloshFS::BarSloshFS(int* nodenums)
{
        nn[0] = nodenums[0];
        nn[1] = nodenums[1];
}

Element *
BarSloshFS::clone()
{
	return new BarSloshFS(*this);
}

void
BarSloshFS::renum(const int *table)
{
	nn[0] = table[nn[0]];
	nn[1] = table[nn[1]];
}

void
BarSloshFS::renum(EleRenumMap& table)
{
	nn[0] = table[nn[0]];
	nn[1] = table[nn[1]];
}

FullSquareMatrix
BarSloshFS::massMatrix(const CoordSet &cs, double *Ms, int cmflg) const
{

        auto &nd1 = cs.getNode(nn[0]);
        auto &nd2 = cs.getNode(nn[1]);

        int i;
        double x[2], y[2], MassMat[4];

        x[0] = nd1.x; y[0] = nd1.y;
        x[1] = nd2.x; y[1] = nd2.y;

        // Calculate mass matrix

        const int numgauss = 2;
        const int numdof   = 2;

        double h = prop ->eh;
        if(h == 0) { // PJSA 12/2/2014
          std::cerr << " *** ERROR: BarSloshFS element (type 302) has zero thickness.\n";
        }

        _FORTRAN(barsloshfs)(x, y, h, MassMat, &numgauss, &numdof);

        for (i=0;i<4;i++)
         Ms[i] = MassMat[i];

        FullSquareMatrix ret(2, Ms);

        return ret;

}

FullSquareMatrix
BarSloshFS::stiffness(const CoordSet &cs, double *kel, int flg) const
{

        FullSquareMatrix k_e(2,kel);

        k_e.zero();
                            
        return k_e;
}

int
BarSloshFS::numNodes() const
{
        return 2;
}

int *
BarSloshFS::nodes(int *p) const
{
        if(p == 0) p = new int[2];
        p[0] = nn[0];
        p[1] = nn[1];
        return p;
}

int
BarSloshFS::numDofs() const
{
        return 2;
}

int *
BarSloshFS::dofs(DofSetArray &dsa, int *p) const
{
        if(p == 0) p = new int[2];

        p[0] = dsa.locate(nn[0],DofSet::Potential);
        p[1] = dsa.locate(nn[1],DofSet::Potential);

        return p;
}

void
BarSloshFS::markDofs(DofSetArray& dsa) const
{
        dsa.mark( nn, 2, DofSet::Potential);
}

int
BarSloshFS::getTopNumber() const
{
  return 147;//1;
}
