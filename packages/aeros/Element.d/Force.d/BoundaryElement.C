#include <Element.d/Force.d/BoundaryElement.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Corotational.d/utilities.h>
#include <Corotational.d/GeomState.h>


BoundaryElement::BoundaryElement(int _nNodes, DofSet nodalDofs, int* _nn)
 : nNodes(_nNodes)
{
  // this constructor is for a force involving the same DofSet on each node, used for both inputs and outputs
  nn = new int[nNodes];
  for(int i = 0; i < nNodes; ++i)
    nn[i] = _nn[i];
  addTerms(nodalDofs);
}

void
BoundaryElement::addTerms(DofSet nodalDofs)
{
  terms.clear();
  inputs.clear();
  outputs.clear();
  nterms = 0;
  for(int i = 0; i < nNodes; ++i) {
    int j = 0;
    for(int k = 0; ; ++k) 
      if(nodalDofs.contains(1 << k)) {
        BCond t; 
        t.setData(nn[i], k, 0.0);
        terms.push_back(t);
        inputs.push_back(nterms);
        outputs.push_back(nterms);
        nterms++;
        j++;
        if(j == nodalDofs.count()) break;
      }
  }
}

BoundaryElement::BoundaryElement(int _nNodes, DofSet *nodalDofs, int* _nn)
 : nNodes(_nNodes)
{
  // this constructor is for a force involving a different DofSet on each nodes, used for both inputs and outputs
  nn = new int[nNodes];
  for(int i = 0; i < nNodes; ++i)
    nn[i] = _nn[i];
  addTerms(nodalDofs);
}

void
BoundaryElement::addTerms(DofSet *nodalDofs)
{
  terms.clear();
  inputs.clear();
  outputs.clear();
  nterms = 0;
  for(int i = 0; i < nNodes; ++i) {
    int j = 0;
    for(int k = 0; ; ++k)
      if(nodalDofs[i].contains(1 << k)) {
        BCond t; 
        t.setData(nn[i], k, 0.0);
        terms.push_back(t);
        inputs.push_back(nterms);
        outputs.push_back(nterms);
        nterms++;
        j++;
        if(j == nodalDofs[i].count()) break;
      }
  }
}

BoundaryElement::BoundaryElement(int _nNodes, DofSet *nodalInputDofs, DofSet *nodalOutputDofs, int* _nn)
 : nNodes(_nNodes)
{
  // this constructor is for a force involving a different input and output DofSet on each node
  nn = new int[nNodes];
  for(int i = 0; i < nNodes; ++i)
    nn[i] = _nn[i];
  addTerms(nodalInputDofs, nodalOutputDofs);
}

void
BoundaryElement::addTerms(DofSet *nodalInputDofs, DofSet *nodalOutputDofs)
{
  terms.clear();
  nterms = 0;
  for(int i = 0; i < nNodes; ++i) {
    int j = 0;
    DofSet nodalDofs_i = nodalInputDofs[i] | nodalOutputDofs[i];
    for(int k = 0; ; ++k)
      if(nodalDofs_i.contains(1 << k)) {
        BCond t; 
        t.setData(nn[i], k, 0.0);
        terms.push_back(t);
        if(nodalInputDofs[i].contains(1 << k)) inputs.push_back(nterms);
        if(nodalOutputDofs[i].contains(1 << k)) outputs.push_back(nterms);
        nterms++;
        j++;
        if(j == nodalDofs_i.count()) break;
      }
  }
}

BoundaryElement::~BoundaryElement()
{
  delete [] nn;
}

int
BoundaryElement::numNodes() const
{
  return nNodes;
}

void
BoundaryElement::renum(const int *table)
{
  for(int i = 0; i < numNodes(); ++i)
    if(nn[i] > -1)
      nn[i] = table[nn[i]];
  for(int i = 0; i < nterms; ++i)
    terms[i].nnum = table[terms[i].nnum];
}

void
BoundaryElement::renum(EleRenumMap& table)
{
  for(int i = 0; i < numNodes(); ++i)
    if(nn[i] > -1)
      nn[i] = table[nn[i]];
  for(int i = 0; i < nterms; ++i)
    terms[i].nnum = table[terms[i].nnum];
}

int*
BoundaryElement::nodes(int* p) const
{
  if(p == 0) p = new int[numNodes()];
  for(int i = 0; i < numNodes(); ++i) p[i] = nn[i];
  return p;
}

int
BoundaryElement::numDofs() const
{
  return nterms;
}

int *
BoundaryElement::dofs(DofSetArray &dsa, int *p) const
{
  if(p == 0) p = new int[numDofs()];
  for(int i = 0; i < nterms; i++)
    dsa.number(terms[i].nnum, 1 << terms[i].dofnum, p+i);
  return p;
}

void
BoundaryElement::markDofs(DofSetArray &dsa) const
{
  for(int i = 0; i < nterms; i++)
    dsa.mark(terms[i].nnum, 1 << terms[i].dofnum);
}

Corotator*
BoundaryElement::getCorotator(CoordSet&, double*, int, int)
{
  return this;
}

void
BoundaryElement::getNLVonMises(Vector& stress, Vector& weight,
                               GeomState&, CoordSet&, int)
{
  stress.zero();
  weight.zero();
}

void
BoundaryElement::getGravityForce(CoordSet&, double*, Vector& force, int, GeomState*)
{
  force.zero();
}
