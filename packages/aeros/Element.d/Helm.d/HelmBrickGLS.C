#include	<Element.d/Helm.d/HelmBrickGLS.h>
#include        <Math.d/matrix.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Utils.d/dofset.h>
#include        <Utils.d/linkfc.h>

extern "C"      {
void  _FORTRAN(helmbrik8v)(double*, double*, double*, const int&,
                       double*, const int&, const int& );
void  _FORTRAN(helmbr8mas)(const int&,double*,const int&,double*,double*, double*, double &);
}

HelmBrickGLS::HelmBrickGLS(int* nodenums)
{
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
	nn[2] = nodenums[2];
	nn[3] = nodenums[3];
	nn[4] = nodenums[4];
	nn[5] = nodenums[5];
	nn[6] = nodenums[6];
	nn[7] = nodenums[7];
}

Element *
HelmBrickGLS::clone()
{
 return new HelmBrickGLS(*this);
}

void
HelmBrickGLS::renum(const int *table)
{
        nn[0] = table[nn[0]];
        nn[1] = table[nn[1]];
        nn[2] = table[nn[2]];
        nn[3] = table[nn[3]];
        nn[4] = table[nn[4]];
        nn[5] = table[nn[5]];
        nn[6] = table[nn[6]];
        nn[7] = table[nn[7]];
}

void
HelmBrickGLS::renum(EleRenumMap& table)
{
        nn[0] = table[nn[0]];
        nn[1] = table[nn[1]];
        nn[2] = table[nn[2]];
        nn[3] = table[nn[3]];
        nn[4] = table[nn[4]];
        nn[5] = table[nn[5]];
        nn[6] = table[nn[6]];
        nn[7] = table[nn[7]];
}

double
HelmBrickGLS::getMass(const CoordSet& cs) const
{

 return 0.0;

}

FullSquareMatrix
HelmBrickGLS::massMatrix(const CoordSet &cs,double *d,int cmflg) const
{
     Node nd1 = cs.getNode(nn[0]);
     Node nd2 = cs.getNode(nn[1]);
     Node nd3 = cs.getNode(nn[2]);
     Node nd4 = cs.getNode(nn[3]);
     Node nd5 = cs.getNode(nn[4]);
     Node nd6 = cs.getNode(nn[5]);
     Node nd7 = cs.getNode(nn[6]);
     Node nd8 = cs.getNode(nn[7]);

     double x[8], y[8], z[8], mm[64];

     x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
     x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
     x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
     x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;
     x[4] = nd5.x; y[4] = nd5.y; z[4] = nd5.z;
     x[5] = nd6.x; y[5] = nd6.y; z[5] = nd6.z;
     x[6] = nd7.x; y[6] = nd7.y; z[6] = nd7.z;
     x[7] = nd8.x; y[7] = nd8.y; z[7] = nd8.z;

     int i;

     double totmas = 0.0;

     const int numgauss = 2;
     const int numdof   = 8;

     _FORTRAN(helmbr8mas)(numgauss,mm,numdof,x,y,z,totmas);

      for (i=0;i<64;i++){
       d[i] = mm[i];
      } 

      FullSquareMatrix ret(8,d);

      double kappa = prop ->kappaHelm;
      double h = 0.0;
      for(i=0;i<64;i++) h += d[i];
      h = pow(h,1.0/3.0);
      double tauksq = 1.0 - 6.0/(kappa*h*kappa*h)*(1-cos(kappa*h))/(2+cos(kappa*h));
      coef = (-kappa*kappa+kappa*kappa*tauksq);

      ret /= getProperty()->rho;
      return ret;
}


FullSquareMatrix
HelmBrickGLS::stiffness(const CoordSet &cs, double *Ks, int flg ) const
{
        Node nd1 = cs.getNode(nn[0]);
        Node nd2 = cs.getNode(nn[1]);
        Node nd3 = cs.getNode(nn[2]);
        Node nd4 = cs.getNode(nn[3]);
        Node nd5 = cs.getNode(nn[4]);
        Node nd6 = cs.getNode(nn[5]);
        Node nd7 = cs.getNode(nn[6]);
        Node nd8 = cs.getNode(nn[7]);

        double x[8], y[8], z[8];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
        x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;
        x[4] = nd5.x; y[4] = nd5.y; z[4] = nd5.z;
        x[5] = nd6.x; y[5] = nd6.y; z[5] = nd6.z;
        x[6] = nd7.x; y[6] = nd7.y; z[6] = nd7.z;
        x[7] = nd8.x; y[7] = nd8.y; z[7] = nd8.z;

        const int numgauss = 2;
        const int numdof   = 8;
        const int outerr   = 6;

        _FORTRAN(helmbrik8v)(x, y, z, numgauss, Ks, numdof, outerr);

        FullSquareMatrix ret(8,Ks);

        ret /= getProperty()->rho;
        return ret;
}


FullSquareMatrix
HelmBrickGLS::acousticm(CoordSet &cs, double *Ks)
{
	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);
	Node nd3 = cs.getNode(nn[2]);
	Node nd4 = cs.getNode(nn[3]);
	Node nd5 = cs.getNode(nn[4]);
	Node nd6 = cs.getNode(nn[5]);
	Node nd7 = cs.getNode(nn[6]);
	Node nd8 = cs.getNode(nn[7]);

	double x[8], y[8], z[8], Kstiff[64];

	x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z; 
	x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
	x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
	x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;
	x[4] = nd5.x; y[4] = nd5.y; z[4] = nd5.z;
	x[5] = nd6.x; y[5] = nd6.y; z[5] = nd6.z;
	x[6] = nd7.x; y[6] = nd7.y; z[6] = nd7.z;
	x[7] = nd8.x; y[7] = nd8.y; z[7] = nd8.z;

  	const int numgauss = 2;
  	const int numdof   = 8;
	const int outerr   = 6;

        // First get the stiffness
        _FORTRAN(helmbrik8v)(x, y, z, numgauss, Ks, numdof, outerr);

        // Put the stiffness temporarely in Kstiff
        int i;
        for (i=0;i<64;i++)
          Kstiff[i] = Ks[i];

        // Then get the mass
        double kappa = prop ->kappaHelm;
        double totmas;
        _FORTRAN(helmbr8mas)(numgauss,Ks,numdof,x,y,z,totmas);

/*
//RT: The following implementation is due to Antonini

	double kk = kappa*kappa;

        double x_centroid = (x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7])/8.0;
        double y_centroid = (y[0]+y[1]+y[2]+y[3]+y[4]+y[5]+y[6]+y[7])/8.0;

        double d[8];
        double dsum=0.0;
        for (i=0; i<8; i++) {
          d[i] = sqrt((x[i]-x_centroid)*(x[i]-x_centroid)
                     +(y[i]-y_centroid)*(y[i]-y_centroid));
          dsum += d[i];
        }

        double h_ele = 2.0*(dsum)/8.0;


        // Thompson and Pinsky tau definition
        double alpha = 0.3926990875; // pi/4 = 22.5 degrees
        double h = h_ele;
        double k = kappa;
        double eps1 = k*h*cos(alpha);
        double eps2 = k*h*sin(alpha);
        double a0 = 1/(k*k);
        double a1 = 6/(k*k*h*h);
        double a2 = 4-cos(eps1)-cos(eps2)-2*cos(eps1)*cos(eps2);
        double a3 = (2+cos(eps1))*(2+cos(eps2));
        double tau = (1-a1*a2/a3);

        // Now add together [K] + kappa^2 [M]
        for (i=0;i<64;i++)
          Ks[i] = Kstiff[i] - (1.0-tau)*kk*Ks[i];
*/

   double h = 0.0;
   for(i=0;i<64;i++) h += Ks[i];
   h = pow(h,1.0/3.0);

   double tauksq = 1.0 - 6.0/(kappa*h*kappa*h)*(1-cos(kappa*h))/(2+cos(kappa*h));

   coef = (-kappa*kappa+kappa*kappa*tauksq);
   for (i=0;i<16;i++) Ks[i] = Kstiff[i] + coef*Ks[i];

        FullSquareMatrix ret(8,Ks);

	//ret.print();

        ret /= getProperty()->rho;
        return ret;
}

int
HelmBrickGLS::numNodes() const
{
 	return 8;
}

int*
HelmBrickGLS::nodes(int *p) const
{
 	if(p == 0) p = new int[8];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
 	p[3] = nn[3];
 	p[4] = nn[4];
 	p[5] = nn[5];
 	p[6] = nn[6];
 	p[7] = nn[7];
	return p;
}

int
HelmBrickGLS::numDofs() const
{
 	return 8;
}

int*
HelmBrickGLS::dofs(DofSetArray &dsa, int *p) const
{
 	if(p == 0) p = new int[8];

        p[0] = dsa.locate(nn[0],DofSet::Helm);
        p[1] = dsa.locate(nn[1],DofSet::Helm);
        p[2] = dsa.locate(nn[2],DofSet::Helm);
        p[3] = dsa.locate(nn[3],DofSet::Helm);
        p[4] = dsa.locate(nn[4],DofSet::Helm);
        p[5] = dsa.locate(nn[5],DofSet::Helm);
        p[6] = dsa.locate(nn[6],DofSet::Helm);
        p[7] = dsa.locate(nn[7],DofSet::Helm);

	return p;
}

void
HelmBrickGLS::markDofs(DofSetArray &dsa) const
{
 	dsa.mark(nn[0],DofSet::Helm);
 	dsa.mark(nn[1],DofSet::Helm);
 	dsa.mark(nn[2],DofSet::Helm);
 	dsa.mark(nn[3],DofSet::Helm);
 	dsa.mark(nn[4],DofSet::Helm);
 	dsa.mark(nn[5],DofSet::Helm);
 	dsa.mark(nn[6],DofSet::Helm);
 	dsa.mark(nn[7],DofSet::Helm);
}

int
HelmBrickGLS::getTopNumber() const
{
  return 144;//3;
}



void
HelmBrickGLS::addFaces(PolygonSet *pset)
{

        pset->addQuad(this,nn[0], nn[1], nn[2], nn[3]);
        pset->addQuad(this,nn[4], nn[5], nn[6], nn[7]);
        pset->addQuad(this,nn[3], nn[0], nn[4], nn[7]);
        pset->addQuad(this,nn[0], nn[1], nn[5], nn[4]);
        pset->addQuad(this,nn[2], nn[1], nn[5], nn[6]);
        pset->addQuad(this,nn[3], nn[2], nn[6], nn[7]);

}

int HelmBrickGLS::getDecFace(int iFace, int *fn) {
 switch(iFace) {
  case 0: fn[0] = nn[0];  fn[1] = nn[1]; fn[2] = nn[2]; fn[3] = nn[3]; break;
  case 1: fn[0] = nn[4];  fn[1] = nn[5]; fn[2] = nn[6]; fn[3] = nn[7]; break;
  case 2: fn[0] = nn[3];  fn[1] = nn[0]; fn[2] = nn[4]; fn[3] = nn[7]; break;
  case 3: fn[0] = nn[0];  fn[1] = nn[1]; fn[2] = nn[5]; fn[3] = nn[4]; break;
  case 4: fn[0] = nn[2];  fn[1] = nn[1]; fn[2] = nn[5]; fn[3] = nn[6]; break;
  default:
  case 5: fn[0] = nn[3];  fn[1] = nn[2]; fn[2] = nn[6]; fn[3] = nn[7]; break;
 }
 return 4;
}

