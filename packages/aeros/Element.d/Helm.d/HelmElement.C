#include <cstdlib>
#include <cstdio>
#include <Element.d/Helm.d/HelmElement.h>

void
HelmElement::addFaces(PolygonSet *) {
  fprintf (stderr, "HelmElement::addFaces not implemented for this element. Exiting.\n");
  exit (1);
}

void HelmElement::edgeShapeFunctions(int n1, int n2, int *ng,
                                       double **gw, double **N) {
  fprintf (stderr, "HelmElement::edgeShapeFunctions not implemented for this element. Exiting.\n");
  exit (1);
}


FullSquareMatrix HelmElement::acousticm(CoordSet&, double *kel) {
  fprintf (stderr, "HelmElement::acousticm not implemented for this element. Exiting.\n");
  exit (1);
  return FullSquareMatrix();
}




void
HelmElement::getHelmForce(CoordSet&, ComplexVector &vc, ComplexVector &force)
{
 fprintf(stderr,"HelmElement::getHelmForce not implemented\n");
 force.zero();
}

void HelmElement::wErrors(CoordSet&, double *l2e, double *h1e,
                          double *l2, double *h1, ComplexD *u,
                          double kappa, double *waveDir) {
 fprintf(stderr,"HelmElement::wErrors not implemented.\n");
}


void 
HelmElement::getNormalDeriv(CoordSet&,ComplexD *uel, int ns, int *s, ComplexD *) {
 fprintf(stderr,"HelmElement::getNormalDeriv not implemented.\n");
}


void 
HelmElement::getNormalDeriv(CoordSet&,ComplexD *uel, int ns, int *s, ComplexD *,
                            double kappa, double *wavedir) {
 fprintf(stderr,"HelmElement::getNormalDeriv not implemented.\n");
}
