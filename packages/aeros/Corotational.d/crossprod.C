#ifndef NO_INLINE_COROT_UTILS
void crossprod(double a[3], double b[3], double c[3])
/***************************************************************
 *     Compute cross product of a and b and place result in c
 ***************************************************************/
{
	c[0] =  a[1]*b[2] - b[1]*a[2];
	c[1] = -a[0]*b[2] + b[0]*a[2];
	c[2] =  a[0]*b[1] - b[0]*a[1];
}
#endif
void dcrossprod(double  a[3], double  b[3], double  c[3],
                double da[3], double db[3], double dc[3])
/***************************************************************
 *     Compute cross product of a and b and place result in c
                       and its derivatives
 ***************************************************************/
{
	c[0] =  a[1]*b[2] - b[1]*a[2];
	c[1] = -a[0]*b[2] + b[0]*a[2];
	c[2] =  a[0]*b[1] - b[0]*a[1];

        dc[0] =  da[1]*b[2] + a[1]*db[2] - db[1]*a[2] - b[1]*da[2];
        dc[1] = -da[0]*b[2] - a[0]*db[2] + db[0]*a[2] + b[0]*da[2];
        dc[2] =  da[0]*b[1] + a[0]*db[1] - db[0]*a[1] - b[0]*da[1];
}
