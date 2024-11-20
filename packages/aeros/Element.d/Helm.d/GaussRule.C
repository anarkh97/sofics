#include "GaussRule.h"


GaussRule::GaussRule(int o) {
 order = o;
 pool = 0;
/* if (order>4)
  ngauss = 7;
 else
  ngauss = 4; 
*/
ngauss = 7;
}


void GaussRule::init() {
 wgauss = pool;
 xigauss = pool+ngauss;
 if  (ngauss==4) {
   xigauss[0] = -0.861136311594053;
   xigauss[1] = -0.339981043584856;
   xigauss[2] = 0.339981043584856;
   xigauss[3] = 0.861136311594053;
   wgauss[0] = 0.347854845137454;
   wgauss[1] = 0.652145154862546;
   wgauss[2] = 0.652145154862546;
   wgauss[3] = 0.347854845137454;
 } else {
   xigauss[0] = 0.949107912342758524526190;
   xigauss[1] = -0.949107912342758524526190;
   xigauss[2] = 0.741531185599394439863865;
   xigauss[3] = -0.741531185599394439863865;
   xigauss[4] = 0.405845151377397166906607;
   xigauss[5] = -0.405845151377397166906607;
   xigauss[6] = 0;
   wgauss[0] = 0.12948496616886969327061;
   wgauss[1] = 0.12948496616886969327061;
   wgauss[2] = 0.27970539148927666790147;
   wgauss[3] = 0.27970539148927666790147;
   wgauss[4] = 0.38183005050511894495037;
   wgauss[5] = 0.38183005050511894495037;
   wgauss[6] = 0.41795918367346938775510;
 }
}
