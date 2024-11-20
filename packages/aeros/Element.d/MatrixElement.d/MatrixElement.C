#include <Element.d/MatrixElement.d/MatrixElement.h>
#include <Utils.d/dofset.h>
#include <Math.d/matrix.h>
#include <Math.d/FullSquareMatrix.h>

MatrixElement::MatrixElement(int _nnodes, int *_nn)
      : nnodes(_nnodes), nn(_nn), ndofs(0), alldofs(0), k_real(0), k_complex(0)
{ 
  //cerr << "here in MatrixElement::MatrixElement, nn = "; for(int i=0; i<nnodes; ++i) cerr << nn[i] << " "; cerr << endl;
  prop = new StructProp(); 
}

MatrixElement::~MatrixElement() 
{ 
}

int
MatrixElement::getTopNumber() const {
  return -1;
}

Element *
MatrixElement::clone()
{
  return(new MatrixElement(*this));
}

void
MatrixElement::setDofs(DofSet *d) 
{ 
  alldofs = d;
  ndofs = 0;
  for(int i=0; i<nnodes; ++i) ndofs += alldofs[i].count();
  alldofs = d; 
}

void 
MatrixElement::setStiffness(GenAssembledFullM<double> *k)
{ 
  k_real = k; 
}

void 
MatrixElement::setStiffness(GenAssembledFullM<complex<double> > *k)
{ 
  k_complex = k;
}

FullSquareMatrix
MatrixElement::stiffness(const CoordSet &cs, double *d, int flg) const
{
  FullSquareMatrix K(ndofs,d);
  for(int i=0; i<ndofs; ++i)
    for(int j=0; j<ndofs; ++j)
      K[i][j] = (k_real) ? (*k_real)[i][j] : ScalarTypes::Real((*k_complex)[i][j]);

  return(K);
}

FullSquareMatrix
MatrixElement::imagStiffness(CoordSet &cs, double *d, int flg)
{
  FullSquareMatrix K(ndofs,d);
  for(int i=0; i<ndofs; ++i)
    for(int j=0; j<ndofs; ++j)
      K[i][j] = (k_real) ? 0.0 : ScalarTypes::Imag((*k_complex)[i][j]);

  return(K);
}

void 
MatrixElement::renum(const int *table)
{
  // renumber the nodes
  for(int i=0; i<nnodes; ++i) nn[i] = table[nn[i]];
}

void
MatrixElement::renum(EleRenumMap& m)
{
  // renumber the nodes
  for(int i=0; i<nnodes; ++i) nn[i] = m[nn[i]];
}

void 
MatrixElement::markDofs(DofSetArray &dsa) const
{
  for(int i=0; i<nnodes;  ++i)
    dsa.mark(nn[i], alldofs[i].list());
}

int* MatrixElement::dofs(DofSetArray &dsa, int *p) const
{
  if(!p) p = new int[ndofs];
  int offset = 0;
  for(int i=0; i<nnodes; ++i) {
    dsa.number(nn[i], alldofs[i], p+offset);
    offset += alldofs[i].count();
  }
  //cerr << "here in MatrixElement::dofs, p = "; for(int i=0; i<ndofs; ++i) cerr << p[i] << " "; cerr << endl;

  return(p);
}

int* MatrixElement::nodes(int *p) const
{
  if(!p) p = new int[nnodes];
  for(int i=0; i<nnodes; ++i) p[i] = nn[i];
  return(p);
}
