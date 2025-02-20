// Integrals over the reference domain; shape functions as in p2.math
double p2_dNdxiN[3][10][10] = 
{
{ 
{ 
1/120.0, 
-1/360.0, 
-1/360.0, 
-1/360.0, 
1/90.0, 
-1/90.0, 
-1/90.0, 
1/90.0, 
1/90.0, 
-1/90.0, 
},
{ 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
},
{ 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
},
{ 
1/360.0, 
1/360.0, 
1/360.0, 
-1/120.0, 
1/90.0, 
1/90.0, 
-1/90.0, 
-1/90.0, 
1/90.0, 
-1/90.0, 
},
{ 
-1/90.0, 
0.0, 
-1/90.0, 
-1/90.0, 
2/45.0, 
2/45.0, 
1/45.0, 
1/45.0, 
1/45.0, 
2/45.0, 
},
{ 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
},
{ 
1/90.0, 
1/90.0, 
0.0, 
1/90.0, 
-1/45.0, 
-2/45.0, 
-2/45.0, 
-1/45.0, 
-2/45.0, 
-1/45.0, 
},
{ 
-1/90.0, 
0.0, 
0.0, 
1/90.0, 
-1/45.0, 
0.0, 
1/45.0, 
0.0, 
-1/45.0, 
1/45.0, 
},
{ 
-1/90.0, 
-1/90.0, 
0.0, 
-1/90.0, 
1/45.0, 
2/45.0, 
2/45.0, 
1/45.0, 
2/45.0, 
1/45.0, 
},
{ 
1/90.0, 
0.0, 
1/90.0, 
1/90.0, 
-2/45.0, 
-2/45.0, 
-1/45.0, 
-1/45.0, 
-1/45.0, 
-2/45.0, 
},
 },
{ 
{ 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
},
{ 
-1/360.0, 
1/120.0, 
-1/360.0, 
-1/360.0, 
1/90.0, 
1/90.0, 
-1/90.0, 
-1/90.0, 
-1/90.0, 
1/90.0, 
},
{ 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
},
{ 
1/360.0, 
1/360.0, 
1/360.0, 
-1/120.0, 
1/90.0, 
1/90.0, 
-1/90.0, 
-1/90.0, 
1/90.0, 
-1/90.0, 
},
{ 
0.0, 
-1/90.0, 
-1/90.0, 
-1/90.0, 
2/45.0, 
1/45.0, 
1/45.0, 
2/45.0, 
2/45.0, 
1/45.0, 
},
{ 
-1/90.0, 
-1/90.0, 
0.0, 
-1/90.0, 
1/45.0, 
2/45.0, 
2/45.0, 
1/45.0, 
2/45.0, 
1/45.0, 
},
{ 
1/90.0, 
1/90.0, 
0.0, 
1/90.0, 
-1/45.0, 
-2/45.0, 
-2/45.0, 
-1/45.0, 
-2/45.0, 
-1/45.0, 
},
{ 
0.0, 
1/90.0, 
1/90.0, 
1/90.0, 
-2/45.0, 
-1/45.0, 
-1/45.0, 
-2/45.0, 
-2/45.0, 
-1/45.0, 
},
{ 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
},
{ 
0.0, 
-1/90.0, 
0.0, 
1/90.0, 
-1/45.0, 
-1/45.0, 
1/45.0, 
1/45.0, 
0.0, 
0.0, 
},
 },
{ 
{ 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
},
{ 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
},
{ 
-1/360.0, 
-1/360.0, 
1/120.0, 
-1/360.0, 
-1/90.0, 
1/90.0, 
1/90.0, 
-1/90.0, 
1/90.0, 
-1/90.0, 
},
{ 
1/360.0, 
1/360.0, 
1/360.0, 
-1/120.0, 
1/90.0, 
1/90.0, 
-1/90.0, 
-1/90.0, 
1/90.0, 
-1/90.0, 
},
{ 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
0.0, 
},
{ 
-1/90.0, 
0.0, 
-1/90.0, 
-1/90.0, 
2/45.0, 
2/45.0, 
1/45.0, 
1/45.0, 
1/45.0, 
2/45.0, 
},
{ 
0.0, 
0.0, 
-1/90.0, 
1/90.0, 
0.0, 
-1/45.0, 
0.0, 
1/45.0, 
-1/45.0, 
1/45.0, 
},
{ 
0.0, 
1/90.0, 
1/90.0, 
1/90.0, 
-2/45.0, 
-1/45.0, 
-1/45.0, 
-2/45.0, 
-2/45.0, 
-1/45.0, 
},
{ 
0.0, 
-1/90.0, 
-1/90.0, 
-1/90.0, 
2/45.0, 
1/45.0, 
1/45.0, 
2/45.0, 
2/45.0, 
1/45.0, 
},
{ 
1/90.0, 
0.0, 
1/90.0, 
1/90.0, 
-2/45.0, 
-2/45.0, 
-1/45.0, 
-1/45.0, 
-1/45.0, 
-2/45.0, 
},
 }};
