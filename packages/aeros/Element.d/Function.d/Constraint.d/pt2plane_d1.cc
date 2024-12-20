#include <math.h>

void pt2plane_d1 (double x[12], double dd_dx[12])
{
  double d[1];
  double t1;
  double t10;
  double t102;
  double t11;
  double t114;
  double t12;
  double t124;
  double t13;
  double t136;
  double t15;
  double t16;
  double t18;
  double t19;
  double t2;
  double t20;
  double t21;
  double t22;
  double t23;
  double t24;
  double t26;
  double t28;
  double t29;
  double t3;
  double t30;
  double t32;
  double t34;
  double t35;
  double t36;
  double t37;
  double t38;
  double t39;
  double t4;
  double t40;
  double t41;
  double t43;
  double t49;
  double t5;
  double t56;
  double t58;
  double t6;
  double t68;
  double t7;
  double t70;
  double t8;
  double t80;
  double t9;
  double t92;
  t1 = x[7];
  t2 = x[4];
  t3 = t1 - t2;
  t4 = x[11];
  t5 = x[5];
  t6 = t4 - t5;
  t7 = t3 * t6;
  t8 = x[8];
  t9 = t8 - t5;
  t10 = x[10];
  t11 = t10 - t2;
  t12 = t9 * t11;
  t13 = t7 - t12;
  t15 = x[3];
  t16 = x[0] - t15;
  t18 = x[6];
  t19 = t18 - t15;
  t20 = t19 * t6;
  t21 = x[9];
  t22 = t21 - t15;
  t23 = t9 * t22;
  t24 = -t20 + t23;
  t26 = x[1] - t2;
  t28 = t19 * t11;
  t29 = t3 * t22;
  t30 = t28 - t29;
  t32 = x[2] - t5;
  t34 = t13 * t16 + t24 * t26 + t30 * t32;
  t35 = t13 * t13;
  t36 = t24 * t24;
  t37 = t30 * t30;
  t38 = t35 + t36 + t37;
  t39 = sqrt(t38);
  t40 = 1.0 / t39;
  d[0] = t34 * t40;
  dd_dx[0] = t13 * t40;
  dd_dx[1] = t24 * t40;
  dd_dx[2] = t30 * t40;
  t41 = t4 - t8;
  t43 = -t10 + t1;
  t49 = t34 / t39 / t38;
  dd_dx[3] = (t26 * t41 + t32 * t43 + t12 - t7) * t40 - t49 * (2.0 * t24 * t41 + 2.0 * t30 * t43) / 2.0;
  t56 = -t41;
  t58 = -t18 + t21;
  dd_dx[4] = (t16 * t56 + t32 * t58 + t20 - t23) * t40 - t49 * (2.0 * t13 * t56 + 2.0 * t30 * t58) / 2.0;
  t68 = -t43;
  t70 = -t58;
  dd_dx[5] = (t16 * t68 + t26 * t70 - t28 + t29) * t40 - t49 * (2.0 * t13 * t68 + 2.0 * t24 * t70) / 2.0;
  t80 = -t6;
  dd_dx[6] = (t11 * t32 + t26 * t80) * t40 - t49 * (2.0 * t11 * t30 + 2.0 * t24 * t80) / 2.0;
  t92 = -t22;
  dd_dx[7] = (t16 * t6 + t32 * t92) * t40 - t49 * (2.0 * t13 * t6 + 2.0 * t30 * t92) / 2.0;
  t102 = -t11;
  dd_dx[8] = (t102 * t16 + t22 * t26) * t40 - t49 * (2.0 * t102 * t13 + 2.0 * t22 * t24) / 2.0;
  t114 = -t3;
  dd_dx[9] = (t114 * t32 + t26 * t9) * t40 - t49 * (2.0 * t114 * t30 + 2.0 * t24 * t9) / 2.0;
  t124 = -t9;
  dd_dx[10] = (t124 * t16 + t19 * t32) * t40 - t49 * (2.0 * t124 * t13 + 2.0 * t19 * t30) / 2.0;
  t136 = -t19;
  dd_dx[11] = (t136 * t26 + t16 * t3) * t40 - t49 * (2.0 * t13 * t3 + 2.0 * t136 * t24) / 2.0;
  return;
}
