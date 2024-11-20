#include        <cstdio>
#include        <cmath>
#include        <Utils.d/dbg_alloca.h>

#include        <Element.d/Helm.d/TetraHelmGal.h>
#include        <Math.d/matrix.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Driver.d/Domain.h>
#include        <Utils.d/dofset.h>
#include        <Utils.d/linkfc.h>


TetraHelmGal::TetraHelmGal(int* nodenums) {
    nn[0] = nodenums[0];
    nn[1] = nodenums[1];
    nn[2] = nodenums[2];
    nn[3] = nodenums[3];
}


int
TetraHelmGal::getTopNumber() const
{
    return 141;//5;
}

Element * TetraHelmGal::clone() {
    return new TetraHelmGal(*this);
}


void TetraHelmGal::renum(const int *table) {
    nn[0] = table[nn[0]];
    nn[1] = table[nn[1]];
    nn[2] = table[nn[2]];
    nn[3] = table[nn[3]];
}


void TetraHelmGal::renum(EleRenumMap& table) {
    nn[0] = table[nn[0]];
    nn[1] = table[nn[1]];
    nn[2] = table[nn[2]];
    nn[3] = table[nn[3]];
}


double TetraHelmGal::getMass(const CoordSet& cs) const {
    fprintf(stderr,"TetraHelmGal::getMass not implemented and should not be called.\n");
    return 0.0;
}


FullSquareMatrix TetraHelmGal::massMatrix(const CoordSet &cs, double *mel, int cmflg) const
{
    FullSquareMatrix mass(4,mel);

    double dxdxi[4][3][3], det[4];
    computedxdxi(cs,4,tetra4_derivatives,dxdxi,det);

    int i, j, k, nint = 4;
    for(j=0;j<4;j++) for(k=j;k<4;k++) {
            mass[j][k] = 0.0;
            for(i=0;i<nint;i++) mass[j][k] += tetra4_weights[i]*
                                              tetra4_values[i][j]*
                                              tetra4_values[i][k]*det[i];
        }

    for(j=0;j<4;j++) for(k=0;k<j;k++) mass[j][k] = mass[k][j];

    mass /= getProperty()->rho;
    return mass;
}


FullSquareMatrix TetraHelmGal::stiffness(const CoordSet &cs, double *Ks, int flg ) const
{
    FullSquareMatrix stif(4,Ks);

    double dxdxi[4][3][3], det[4];
    computedxdxi(cs,4,tetra4_derivatives,dxdxi,det);

    int i, j, k, nint = 4;
    for(j=0;j<4;j++) for(k=0;k<4;k++) stif[j][k] = 0.0;

    for(i=0;i<nint;i++) {
        Matrix33 &x = dxdxi[i];
        double inv[3][3];
        inv[0][0] = -x[1][2]*x[2][1] + x[1][1]*x[2][2];
        inv[0][1] = x[0][2]*x[2][1] - x[0][1]*x[2][2];
        inv[0][2] = -x[0][2]*x[1][1] + x[0][1]*x[1][2];
        inv[1][0] = x[1][2]*x[2][0] - x[1][0]*x[2][2];
        inv[1][1] = -x[0][2]*x[2][0] + x[0][0]*x[2][2];
        inv[1][2] = x[0][2]*x[1][0] - x[0][0]*x[1][2];
        inv[2][0] = -x[1][1]*x[2][0] + x[1][0]*x[2][1];
        inv[2][1] = x[0][1]*x[2][0] - x[0][0]*x[2][1];
        inv[2][2] = -x[0][1]*x[1][0] + x[0][0]*x[1][1];

        double grad[4][3];
        int l,m,n;
        for(l=0;l<4;l++) for(m=0;m<3;m++) {
                grad[l][m] = 0.0;
                for(n=0;n<3;n++) grad[l][m] += inv[n][m] * tetra4_derivatives[i][l][n];
            }

        for(l=0;l<4;l++) for(m=0;m<4;m++) for(n=0;n<3;n++)
                    stif[l][m] += tetra4_weights[i] * grad[l][n] * grad[m][n] / det[i];
    }
    stif /= getProperty()->rho;
    return stif;
}


FullSquareMatrix TetraHelmGal::acousticm(CoordSet &cs, double *d)
{
    FullSquareMatrix sm(4,d);

    double stif[4][4], mass[4][4];

    double dxdxi[4][3][3], det[4];
    computedxdxi(cs,4,tetra4_derivatives,dxdxi,det);

    int i,j,k,nint = 4;
    for(j=0;j<4;j++) for(k=j;k<4;k++) {
            mass[j][k] = 0.0;
            for(i=0;i<nint;i++) mass[j][k] += tetra4_weights[i]*
                                              tetra4_values[i][j]*
                                              tetra4_values[i][k]*det[i];
        }

    for(j=0;j<4;j++) for(k=0;k<j;k++) mass[j][k] = mass[k][j];

    for(j=0;j<4;j++) for(k=0;k<4;k++) stif[j][k] = 0.0;

    for(i=0;i<nint;i++) {
        Matrix33 &x = dxdxi[i];
        double inv[3][3];
        inv[0][0] = -x[1][2]*x[2][1] + x[1][1]*x[2][2];
        inv[0][1] = x[0][2]*x[2][1] - x[0][1]*x[2][2];
        inv[0][2] = -x[0][2]*x[1][1] + x[0][1]*x[1][2];
        inv[1][0] = x[1][2]*x[2][0] - x[1][0]*x[2][2];
        inv[1][1] = -x[0][2]*x[2][0] + x[0][0]*x[2][2];
        inv[1][2] = x[0][2]*x[1][0] - x[0][0]*x[1][2];
        inv[2][0] = -x[1][1]*x[2][0] + x[1][0]*x[2][1];
        inv[2][1] = x[0][1]*x[2][0] - x[0][0]*x[2][1];
        inv[2][2] = -x[0][1]*x[1][0] + x[0][0]*x[1][1];

        double grad[4][3];
        int l,m,n;
        for(l=0;l<4;l++) for(m=0;m<3;m++) {
                grad[l][m] = 0.0;
                for(n=0;n<3;n++) grad[l][m] += inv[n][m] * tetra4_derivatives[i][l][n];
            }

        for(l=0;l<4;l++) for(m=0;m<4;m++) for(n=0;n<3;n++)
                    stif[l][m] += tetra4_weights[i] * grad[l][n] * grad[m][n] / det[i];
    }

    double kappa = prop -> kappaHelm;
    double kk = kappa * kappa;

    for (i=0;i<4;i++)
        for (j=0;j<4;j++) sm[i][j] = stif[i][j] - kk*mass[i][j];

    sm /= getProperty()->rho;
    return sm;
}


int TetraHelmGal::numNodes() const {
    return 4;
}


int* TetraHelmGal::nodes(int *p) const {
    if(p == 0) p = new int[4];
    p[0] = nn[0];
    p[1] = nn[1];
    p[2] = nn[2];
    p[3] = nn[3];
    return p;
}


int TetraHelmGal::numDofs() const {
    return 4;
}


int* TetraHelmGal::dofs(DofSetArray &dsa, int *p) const  {
    if(p == 0) p = new int[4];

    p[0] = dsa.locate(nn[0],DofSet::Helm);
    p[1] = dsa.locate(nn[1],DofSet::Helm);
    p[2] = dsa.locate(nn[2],DofSet::Helm);
    p[3] = dsa.locate(nn[3],DofSet::Helm);

    return p;
}


void TetraHelmGal::markDofs(DofSetArray &dsa) const {

    dsa.mark(nn[0],DofSet::Helm);
    dsa.mark(nn[1],DofSet::Helm);
    dsa.mark(nn[2],DofSet::Helm);
    dsa.mark(nn[3],DofSet::Helm);
}


void TetraHelmGal::addFaces(PolygonSet *pset) {
    fprintf(stderr,"TetraHelmGal::addFaces not implemented.\n");
/*
 pset->addTri(this,nn[0], nn[2], nn[1]);
 pset->addTri(this,nn[0], nn[1], nn[3]);
 pset->addTri(this,nn[0], nn[3], nn[2]);
 pset->addTri(this,nn[2], nn[3], nn[1]);
*/
}

int TetraHelmGal::getDecFace(int iFace, int *fn) {
    switch(iFace) {
        case 0: fn[0] = nn[0];  fn[1] = nn[2]; fn[2] = nn[1]; break;
        case 1: fn[0] = nn[0];  fn[1] = nn[1]; fn[2] = nn[3]; break;
        case 2: fn[0] = nn[0];  fn[1] = nn[3]; fn[2] = nn[2]; break;
        default:
        case 3: fn[0] = nn[2];  fn[1] = nn[3]; fn[2] = nn[1]; break;
    }
    return 3;
}



void TetraHelmGal::computedxdxi(const CoordSet &cs, int nint,
                                double (*derivatives)[4][3], Matrix33 *dxdxi, double *det) const{

    int i,j,k;
    Node *nd = (Node*)dbg_alloca(sizeof(Node)*numNodes());
    for(i=0;i<numNodes();i++) nd[i] = cs.getNode(nn[i]);

    double (*coord)[3] = (double(*)[3])dbg_alloca(numNodes()*3*sizeof(double));
    for(i=0;i<numNodes();i++) {
        coord[i][0] = nd[i].x;
        coord[i][1] = nd[i].y;
        coord[i][2] = nd[i].z;
    }

    for(j=0;j<3;j++) for(k=0;k<3;k++) {
            for(i=0;i<nint;i++) {
                dxdxi[i][j][k] = 0.0;
                int m;
                for(m=0;m<numNodes();m++)
                    dxdxi[i][j][k] += derivatives[i][m][k] * coord[m][j];
            }
        }

    for(i=0;i<nint;i++) {
        Matrix33 &x = dxdxi[i];
        det[i] =
                -x[0][2]*x[1][1]*x[2][0] + x[0][1]*x[1][2]*x[2][0] +
                x[0][2]*x[1][0]*x[2][1] - x[0][0]*x[1][2]*x[2][1] -
                x[0][1]*x[1][0]*x[2][2] + x[0][0]*x[1][1]*x[2][2];
    }

}



void TetraHelmGal::getNormalDeriv(CoordSet&cs,ComplexD *uel, int nsc,
                                  int *sc, ComplexD *grad,
                                  double kappa, double *waveDir) {

// Gradient is constant. Just compute it anywhere
    double dN[1][4][3];
    dN[0][0][0] = -1.0;
    dN[0][0][1] = -1.0;
    dN[0][0][2] = -1.0;
    dN[0][1][0] = 1.0;
    dN[0][1][1] = 0.0;
    dN[0][1][2] = 0.0;
    dN[0][2][0] = 0.0;
    dN[0][2][1] = 1.0;
    dN[0][2][2] = 0.0;
    dN[0][3][0] = 0.0;
    dN[0][3][1] = 0.0;
    dN[0][3][2] = 1.0;
    double dxdxi[1][3][3],det[1];
    computedxdxi(cs,1,dN,dxdxi,det);
    Matrix33 &x = dxdxi[0];
    double inv[3][3];
    inv[0][0] = (-x[1][2]*x[2][1] + x[1][1]*x[2][2])/det[0];
    inv[0][1] = (x[0][2]*x[2][1] - x[0][1]*x[2][2])/det[0];
    inv[0][2] = (-x[0][2]*x[1][1] + x[0][1]*x[1][2])/det[0];
    inv[1][0] = (x[1][2]*x[2][0] - x[1][0]*x[2][2])/det[0];
    inv[1][1] = (-x[0][2]*x[2][0] + x[0][0]*x[2][2])/det[0];
    inv[1][2] = (x[0][2]*x[1][0] - x[0][0]*x[1][2])/det[0];
    inv[2][0] = (-x[1][1]*x[2][0] + x[1][0]*x[2][1])/det[0];
    inv[2][1] = (x[0][1]*x[2][0] - x[0][0]*x[2][1])/det[0];
    inv[2][2] = (-x[0][1]*x[1][0] + x[0][0]*x[1][1])/det[0];

    grad[0] = grad[1] = grad[2] = ComplexD(0.0,0.0);
    int l,m,n;
    for(l=0;l<4;l++) {
        for(m=0;m<3;m++) {
            for(n=0;n<3;n++)
                grad[m] += inv[n][m] * dN[0][l][n] * uel[l];
        }
    }

// Exponential is not constant, so find the midpoint of the face
    int iVertex[3];
    int i;

    for(i=0;i<4;i++) {
        if (sc[0]==nn[i]) iVertex[0] = i;
        if (sc[1]==nn[i]) iVertex[1] = i;
        if (sc[2]==nn[i]) iVertex[2] = i;
    }

    Node nd[4];
    for(i=0;i<4;i++) nd[i] = cs.getNode(nn[i]);

    double coord[4][3];
    for(i=0;i<4;i++) {
        coord[i][0] = nd[i].x;
        coord[i][1] = nd[i].y;
        coord[i][2] = nd[i].z;
    }

    double point[3] = {0.0,0.0,0.0};
    for(m=0;m<3;m++) {
        for(l=0;l<3;l++) {
            point[m] += coord[iVertex[l]][m];
        }
        point[m] /= 3.0;
    }

    ComplexD tmp = exp(ComplexD(0.0,kappa)*
                       (waveDir[0]*point[0]+waveDir[1]*point[1]+waveDir[2]*point[2]));
    for(m=0;m<3;m++)
        grad[m] -= -ComplexD(0.0,1.0)*kappa*waveDir[m]*tmp;

}

// Integration rule with 4 points
// r[0] = 0.58541020; s[0] = 0.13819660; t[0] = 0.13819660; and symmetric
// n[0] := 1-r-s-t
// n[1] := r
// n[2] := s
// n[3] := t

double TetraHelmGal::tetra4_weights[4] = {
        .04166666666666666666,.04166666666666666666,
        .04166666666666666666,.04166666666666666666
};

double TetraHelmGal::tetra4_values [4][4] = {
        { 0.13819660, 0.58541020, 0.13819660, 0.13819660 },
        { 0.13819660, 0.13819660, 0.58541020, 0.13819660 },
        { 0.13819660, 0.13819660, 0.13819660, 0.58541020 },
        { 0.58541020, 0.13819660, 0.13819660, 0.13819660 }
};

double TetraHelmGal::tetra4_derivatives[4][4][3] = {
        { { -1.0,-1.0,-1.0 }, { 1.0,0.0,0.0 }, { 0.0,1.0,0.0 }, { 0.0,0.0,1.0 } },
        { { -1.0,-1.0,-1.0 }, { 1.0,0.0,0.0 }, { 0.0,1.0,0.0 }, { 0.0,0.0,1.0 } },
        { { -1.0,-1.0,-1.0 }, { 1.0,0.0,0.0 }, { 0.0,1.0,0.0 }, { 0.0,0.0,1.0 } },
        { { -1.0,-1.0,-1.0 }, { 1.0,0.0,0.0 }, { 0.0,1.0,0.0 }, { 0.0,0.0,1.0 } }
};

//HB
/*
template<class Scalar>
double TetraHelmGal::nodalDisplacement(Scalar disp[4][3], Scalar pel[4], CoordSet&cs, double w, double rho)
{
  // u = Grad(p) / (w^2.rho)
  // gradient is constant
  Node nd[4];
  for(int i=0;i<4;i++) nd[i] = cs.getNode(nn[i]);

  double coord[4][3];
  for(int i=0;i<4;i++) {
    coord[i][0] = nd[i].x;
    coord[i][1] = nd[i].y;
    coord[i][2] = nd[i].z;
  }

  double jac[3][3];
                                                                                                                           
  //Jacobian
  // J_ij = dx_i/dxi_j
  jac[0][0] = coord[1][0] - coord[0][0]; 
  jac[0][1] = coord[2][0] - coord[0][0]; 
  jac[0][2] = coord[3][0] - coord[0][0]; 
  jac[1][0] = coord[1][1] - coord[0][1]; 
  jac[1][1] = coord[2][1] - coord[0][1]; 
  jac[1][2] = coord[3][1] - coord[0][1]; 
  jac[2][0] = coord[1][2] - coord[0][2]; 
  jac[2][1] = coord[2][2] - coord[0][2];
  jac[2][2] = coord[3][2] - coord[0][2]; 

  // compute determinant of jac
  double dOmega = jac[0][0] * (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1]) +
                  jac[1][0] * (jac[0][2] * jac[2][1] - jac[0][1] * jac[2][2]) +
                  jac[2][0] * (jac[0][1] * jac[1][2] - jac[0][2] * jac[1][1]);
                                                                           
  // compute inverse matrix of jac
  // Maple code used
  double t17 = -1.0/dOmega;
                                                                                                                           
  //compute shape function gradients
  // Shape function gradients dN_i/dx_i = dN/dxi * transpose(jInv)
  // Note: 1st index = shape function #
  // 2nd index = direction (0=x, 1=y, 2=z)
  double nGrad[4][3];                                                                                                                         
  nGrad[1][0] =  (-jac[1][1] * jac[2][2] + jac[1][2] * jac[2][1] ) * t17;
  nGrad[1][1] =  ( jac[0][1] * jac[2][2] - jac[0][2] * jac[2][1] ) * t17;
  nGrad[1][2] = -( jac[0][1] * jac[1][2] - jac[0][2] * jac[1][1] ) * t17;
                                                                                                                           
  nGrad[2][0] = -(-jac[1][0] * jac[2][2] + jac[1][2] * jac[2][0] ) * t17;
  nGrad[2][1] = -( jac[0][0] * jac[2][2] - jac[0][2] * jac[2][0] ) * t17;
  nGrad[2][2] =  ( jac[0][0] * jac[1][2] - jac[0][2] * jac[1][0] ) * t17;
                                                                                                                           
  nGrad[3][0] = -( jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0] ) * t17;
  nGrad[3][1] =  ( jac[0][0] * jac[2][1] - jac[0][1] * jac[2][0] ) * t17;
  nGrad[3][2] = -( jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0] ) * t17;
                                                                                                                           
  nGrad[0][0] = -( nGrad[1][0] + nGrad[2][0] + nGrad[3][0] );
  nGrad[0][1] = -( nGrad[1][1] + nGrad[2][1] + nGrad[3][1] );
  nGrad[0][2] = -( nGrad[1][2] + nGrad[2][2] + nGrad[3][2] );

  double c = 1.0/(rho*w*w);
  for(int i=0; i<4; i++){
    disp[i][0] = c*(nGrad[0][0]*pel[0] + nGrad[1][0]*pel[1] + nGrad[2][0]*pel[2] + nGrad[3][0]*pel[3]);
    disp[i][1] = c*(nGrad[0][1]*pel[0] + nGrad[1][1]*pel[1] + nGrad[2][1]*pel[2] + nGrad[3][1]*pel[3]);
    disp[i][2] = c*(nGrad[0][2]*pel[0] + nGrad[1][2]*pel[1] + nGrad[2][2]*pel[2] + nGrad[3][2]*pel[3]);
  }                                                                                                                           

  return dOmega/6.;
}

double Tet::computeGradientP1Function(SVec<double,3> &nodes, double nGrad[4][3])
{
                                                                                                                           
  double jac[3][3];
                                                                                                                           
  //Jacobian
  // J_ij = dx_i/dxi_j
  jac[0][0] = nodes[ nodeNum[1] ][0] - nodes[ nodeNum[0] ][0];
  jac[0][1] = nodes[ nodeNum[2] ][0] - nodes[ nodeNum[0] ][0];
  jac[0][2] = nodes[ nodeNum[3] ][0] - nodes[ nodeNum[0] ][0];
  jac[1][0] = nodes[ nodeNum[1] ][1] - nodes[ nodeNum[0] ][1];
  jac[1][1] = nodes[ nodeNum[2] ][1] - nodes[ nodeNum[0] ][1];
  jac[1][2] = nodes[ nodeNum[3] ][1] - nodes[ nodeNum[0] ][1];
  jac[2][0] = nodes[ nodeNum[1] ][2] - nodes[ nodeNum[0] ][2];
  jac[2][1] = nodes[ nodeNum[2] ][2] - nodes[ nodeNum[0] ][2];
  jac[2][2] = nodes[ nodeNum[3] ][2] - nodes[ nodeNum[0] ][2];
                                                                                                                           
  // compute determinant of jac
  double dOmega = jac[0][0] * (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1]) +
                  jac[1][0] * (jac[0][2] * jac[2][1] - jac[0][1] * jac[2][2]) +
                  jac[2][0] * (jac[0][1] * jac[1][2] - jac[0][2] * jac[1][1]);
                                                                                                                           
  // compute inverse matrix of jac
  // Maple code used
  double t17 = -1.0/dOmega;
                                                                                                                           
  //compute shape function gradients
  nGrad[1][0] =  (-jac[1][1] * jac[2][2] + jac[1][2] * jac[2][1] ) * t17;
  nGrad[1][1] =  ( jac[0][1] * jac[2][2] - jac[0][2] * jac[2][1] ) * t17;
  nGrad[1][2] = -( jac[0][1] * jac[1][2] - jac[0][2] * jac[1][1] ) * t17;
                                                                                                                           
  nGrad[2][0] = -(-jac[1][0] * jac[2][2] + jac[1][2] * jac[2][0] ) * t17;
  nGrad[2][1] = -( jac[0][0] * jac[2][2] - jac[0][2] * jac[2][0] ) * t17;
  nGrad[2][2] =  ( jac[0][0] * jac[1][2] - jac[0][2] * jac[1][0] ) * t17;
                                                                                                                           
  nGrad[3][0] = -( jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0] ) * t17;
  nGrad[3][1] =  ( jac[0][0] * jac[2][1] - jac[0][1] * jac[2][0] ) * t17;
  nGrad[3][2] = -( jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0] ) * t17;
                                                                                                                           
  // Shape function gradients dN_i/dx_i = dN/dxi * transpose(jInv)
  // Note: 1st index = shape function #
  // 2nd index = direction (0=x, 1=y, 2=z)
                                                                                                                           
  nGrad[0][0] = -( nGrad[1][0] + nGrad[2][0] + nGrad[3][0] );
  nGrad[0][1] = -( nGrad[1][1] + nGrad[2][1] + nGrad[3][1] );
  nGrad[0][2] = -( nGrad[1][2] + nGrad[2][2] + nGrad[3][2] );
                                                                                                                           
  return sixth * dOmega;
                                                                                                                           
}

*/


