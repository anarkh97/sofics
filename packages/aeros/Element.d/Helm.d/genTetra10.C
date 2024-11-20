#include <cstdio>
/*
n[1] := (1-r-s-t)*(2*(1-r-s-t)-1)
n[2] := r*(2r-1)
n[3] := s*(2s-1)
n[4] := t*(2t-1)
n[5] := 4r*(1-r-s-t)
n[6] := 4*r*s
n[7] := 4*s*(1-r-s-t)
n[8] := 4*t*(1-r-s-t)
n[9] := 4*r*t
n[10] := 4*s*t


dNdr[i_] := D[n[i],r]
dNds[i_] := D[n[i],s]
dNdt[i_] := D[n[i],t]


For[i=1, i<=10, i++,
  Write["p2m",dNdr[i]]; ];
For[i=1, i<=10, i++,
  Write["p2m",dNds[i]]; ];
For[i=1, i<=10, i++,
  Write["p2m",dNdt[i]]; ];
*/

int main() {

double dNdr[10][10];
double dNds[10][10];
double dNdt[10][10];

int i,j;

double r[10],s[10],t[10];

r[0] = 0.0; s[0] = 0.0; t[0] = 0.0;
r[1] = 1.0; s[1] = 0.0; t[1] = 0.0;
r[2] = 0.0; s[2] = 1.0; t[2] = 0.0;
r[3] = 0.0; s[3] = 0.0; t[3] = 1.0;
r[4] = 0.5; s[4] = 0.0; t[4] = 0.0;
r[5] = 0.5; s[5] = 0.5; t[5] = 0.0;
r[6] = 0.0; s[6] = 0.5; t[6] = 0.0;
r[7] = 0.0; s[7] = 0.0; t[7] = 0.5;
r[8] = 0.5; s[8] = 0.0; t[8] = 0.5;
r[9] = 0.0; s[9] = 0.5; t[9] = 0.5;

for(i=0;i<10;i++) {
  
  dNdr[0][i] = 1.0 - 4.0*(1.0 - r[i] - s[i] - t[i]);
  dNdr[1][i] = -1.0 + 4.0*r[i];
  dNdr[2][i] = 0.0;
  dNdr[3][i] = 0.0;
  dNdr[4][i] = -4.0*r[i] + 4.0*(1.0 - r[i] - s[i] - t[i]);
  dNdr[5][i] = 4.0*s[i];
  dNdr[6][i] = -4.0*s[i];
  dNdr[7][i] = -4.0*t[i];
  dNdr[8][i] = 4.0*t[i];
  dNdr[9][i] = 0.0;


  dNds[0][i] = 1.0 - 4.0*(1.0 - r[i] - s[i] - t[i]);
  dNds[1][i] = 0.0;
  dNds[2][i] = -1.0 + 4.0*s[i];
  dNds[3][i] = 0.0;
  dNds[4][i] = -4.0*r[i];
  dNds[5][i] = 4.0*r[i];
  dNds[6][i] = -4.0*s[i] + 4.0*(1.0 - r[i] - s[i] - t[i]);
  dNds[7][i] = -4.0*t[i];
  dNds[8][i] = 0.0;
  dNds[9][i] = 4.0*t[i];


  dNdt[0][i] = 1.0 - 4.0*(1.0 - r[i] - s[i] - t[i]);
  dNdt[1][i] = 0.0;
  dNdt[2][i] = 0.0;
  dNdt[3][i] = -1.0 + 4.0*t[i];
  dNdt[4][i] = -4.0*r[i];
  dNdt[5][i] = 0.0;
  dNdt[6][i] = -4.0*s[i];
  dNdt[7][i] = 4.0*(1.0 - r[i] - s[i] - t[i]) - 4.0*t[i];
  dNdt[8][i] = 4.0*r[i];
  dNdt[9][i] = 4.0*s[i];
}

printf("{\n");
for(j=0;j<10;j++) {
  printf("{\n");
  for(i=0;i<10;i++) {
    printf("{ %.16f, %.16f, %.16f}",dNdr[i][j],dNds[i][j],dNdt[i][j]);
    if (i==9) printf("\n");
    else printf(",\n");
  }
  printf("}");
  if (j==9) printf("\n");
  else printf(",\n");
}
printf("}\n");
}
