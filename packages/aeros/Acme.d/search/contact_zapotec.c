#include <stdio.h>
#include <math.h>
#include <stdlib.h>

struct tet {		/* pointers to coordinates for 4 corner points. */
   double *corner[4];
};

#define FALSE 0
#define TRUE  1

static void intersect_line(double *, double *, int, double, double *);

static double intersect_3d_recurse(double **, double **, int *, int *, int *,
                                   int *, double *);
extern double tet_volume(double **); 
/*extern void slice_tet(double **, int, double, double (*)[], struct tet *,*/
extern void slice_tet(double **, int, double, double p[4][3], struct tet *,
                      int *, struct tet *, int *);

double intersect_3d(double *t0, double *t1, double *t2, double *t3,
                    double *xmesh, double *ymesh, double *zmesh, double *areas)

/* Calculate the volume of intersection between a set of rectilinear
   blocks and a tetrahedron */

/* double t0[3], t1[3], t2[3], t3[3]; vertices of tetrahedron */
/* double *xmesh, *ymesh, *zmesh; x, y and z locations of grid lines */
/* double *areas;	area of intersection for all nx*ny*nz cells */

 /* ordered x-major, then y then z */
{
  int       global_mesh_min[3];       /* starting indices of mesh to explore*/
  int       global_mesh_max[3];       /* ending indices of mesh to explore  */
  double   *mesh[3];                  /* 3 vectors of mesh spacing          */
  double   *corner[4];                /* array pointing to corner points    */
  int       i, j;                     /* loop counters                      */

  corner[0] = t0;
  corner[1] = t1;
  corner[2] = t2;
  corner[3] = t3;

  /* Subtract one to convert from Fortran to C numbering. */
  global_mesh_min[0] = 0;
  global_mesh_min[1] = 0;
  global_mesh_min[2] = 0;

  global_mesh_max[0] = 1;
  global_mesh_max[1] = 1;
  global_mesh_max[2] = 1;

  mesh[0] = xmesh;
  mesh[1] = ymesh;
  mesh[2] = zmesh;

  /* Zero the return values */
  j = 1;
  for (i = 0; i < j; ++i) {
    areas[i] = 0;
  }

  return (intersect_3d_recurse(corner, mesh, 
                               global_mesh_min, global_mesh_max, 
                               global_mesh_min, global_mesh_max, areas));
}


static double intersect_3d_recurse(
double   *corner[4],		/* vertices of tetrahedron */
double   *mesh[3],		/* x, y and z locations of grid lines */
int       local_mesh_min[3],	/* start of grid cells I'm computing */
int       local_mesh_max[3],	/* end of grid cells I'm computing */
int       global_mesh_min[3],	/* start of full set of grid cells */
int       global_mesh_max[3],	/* end of full set of grid cells */
double   *areas 		/* area of intersection for all nx*ny*nz cells*/
)
 /* ordered x-major, then y then z */
{
  double    tmin, tmax;       /* bounding box around tetrahedron */
  double   *cmesh;            /* selects between x, y, and z-mesh arrays */
  double    total_volume;     /* total volume of intersected tet */
  double    points[4][3];     /* space for storing intersection points */
  struct tet high_tet[3];     /* number of subtets above cut */
  struct tet low_tet[3];      /* number of subtets below cut */
  int       local_submesh_min[3];     /* start of grid cells for recursion */
  int       local_submesh_max[3];     /* end of grid cells for recursion */
  int       index[3];         /* offsets into array of volume results */
  int       order[3];         /* order to check dimensions */
  int       iorder;           /* index into iorder */
  int       delta[3];         /* span of box in each dimension */
  int       lo, hi;           /* indices of grid lines bounding tet */
  int       cut_dir;          /* direction I'm cutting mesh */
  int       mid;              /* index into middle of mesh array */
  int       done;             /* flag for completion of cutting loop */
  int       nlow, nhigh;      /* number of tets above/below cut */
  int       i, j;             /* loop counters */

  /* Find a plane which slices the tet. Divide the pieces of the tet into
     smaller tets. Recurse on the pieces as appropriate. */

  done = FALSE;
  total_volume = 0;

  /* To improve performance & numerics, cut tet in longest direction. */
  /* As a quick stand-in, use # grid lines in bounding box for now */
  for (i = 0; i < 3; ++i) {
    local_submesh_min[i] = local_mesh_min[i];
    local_submesh_max[i] = local_mesh_max[i];
    delta[i] = local_submesh_max[i] - local_submesh_min[i];
  }
  if (delta[0] <= delta[1] && delta[0] <= delta[2])
    order[0] = 0;
  else if (delta[1] < delta[0] && delta[1] <= delta[2])
    order[0] = 1;
  else
    order[0] = 2;

  if (delta[0] > delta[1] && delta[0] > delta[2])
    order[2] = 0;
  else if (delta[1] >= delta[0] && delta[1] > delta[2])
    order[2] = 1;
  else
    order[2] = 2;

  order[1] = 3 - order[0] - order[2];

  /* Compute quantities needed for binary search */
  iorder = 0;
  cut_dir = order[iorder];
  cmesh = mesh[cut_dir];
  tmax = corner[0][cut_dir];
  tmin = corner[0][cut_dir];
  for (j = 1; j < 4; ++j) {
    if (corner[j][cut_dir] < tmin)
      tmin = corner[j][cut_dir];
    if (corner[j][cut_dir] > tmax)
      tmax = corner[j][cut_dir];
  }
  lo = local_mesh_min[cut_dir];
  hi = local_mesh_max[cut_dir];
  if (tmin >= cmesh[hi] || tmax <= cmesh[lo]) {       /* Zero intersection */
    done = TRUE;
    iorder = 3;
  }

  /* Use a binary search to find a grid plane which slices the tet. */
  while (!done) {
    mid = (hi + lo) / 2;
    if (cmesh[mid] < tmax && cmesh[mid] > tmin) {
      /* Found a plane which cuts tetrahedron */
      done = TRUE;
    }
    else {
      if (cmesh[mid] <= tmin) {
        lo = mid;
        /*local_mesh_min[cut_dir] = mid;*/
      }
      else {
        hi = mid;
        /*local_mesh_max[cut_dir] = mid;*/
      }

      if (hi <= lo + 1) {
        /* Still need to try the hi or lo plane */
        mid = hi;
        if (cmesh[mid] < tmax && cmesh[mid] > tmin) {
          /* Found a plane which cuts tetrahedron */
          done = TRUE;
        }
        else {
          mid = lo;
          if (cmesh[mid] < tmax && cmesh[mid] > tmin) {
            /* Found a plane which cuts tetrahedron */
            done = TRUE;
          }
        }
        if (!done) {
          index[cut_dir] = lo - global_mesh_min[cut_dir];
          if (++iorder < 3) { /* Try the next direction */
            cut_dir = order[iorder];
            cmesh = mesh[cut_dir];
            tmax = corner[0][cut_dir];
            tmin = corner[0][cut_dir];
            for (j = 1; j < 4; ++j) {
              if (corner[j][cut_dir] < tmin)
                  tmin = corner[j][cut_dir];
              if (corner[j][cut_dir] > tmax)
                  tmax = corner[j][cut_dir];
            }
            lo = local_mesh_min[cut_dir];
            hi = local_mesh_max[cut_dir];
            if (tmin >= cmesh[hi] || tmax <= cmesh[lo]) {
              /* Zero intersection */
              done = TRUE;
              iorder = 3;
            }
          }
          else {      /* No plane cuts tet. */
            done = TRUE;
            total_volume = tet_volume(corner);
            i = (index[2]*(global_mesh_max[1] - global_mesh_min[1])
               + index[1])*(global_mesh_max[0] - global_mesh_min[0])
               + index[0];

            areas[i] += total_volume;

          }
        }
      }
    }
  }

  if (iorder < 3) {           /* Found a plane which cuts tetrahedron */
    slice_tet(corner, cut_dir, cmesh[mid], points, high_tet, &nhigh,
              low_tet, &nlow);
    if (mid < global_mesh_max[cut_dir]) {   /* Recurse above cut */
      for (i = 0; i < nhigh; ++i) {
        for (j = 0; j < 3; ++j) {
          local_submesh_min[j] = local_mesh_min[j];
          local_submesh_max[j] = local_mesh_max[j];
        }
        local_submesh_min[cut_dir] = mid;
        total_volume += intersect_3d_recurse(high_tet[i].corner, mesh,
                              local_submesh_min, local_submesh_max,
                              global_mesh_min, global_mesh_max, areas);
      }
    }
    if (mid > global_mesh_min[cut_dir]) {   /* Recurse below cut */
      for (i = 0; i < nlow; ++i) {
        for (j = 0; j < 3; ++j) {
          local_submesh_min[j] = local_mesh_min[j];
          local_submesh_max[j] = local_mesh_max[j];
        }
        local_submesh_max[cut_dir] = mid;
        total_volume += intersect_3d_recurse(low_tet[i].corner, mesh,
                              local_submesh_min, local_submesh_max,
                              global_mesh_min, global_mesh_max, areas);
      }
    }
  }

  return (total_volume);
}


void check_bound(
double   *corner[4],		/* vertices of tetrahedron */
double   *mesh[3],		/* x, y and z locations of grid lines */
int       local_mesh_min[3],	/* start of grid cells I'm computing */
int       local_mesh_max[3] 	/* end of grid cells I'm computing */
)
{
  int i, j;

  for (i = 0; i < 4; ++i) {
    for (j = 0; j < 3; ++j) {
      if (corner[i][j] < mesh[j][local_mesh_min[j]]) {
        printf("point too low\n");
      }
      if (corner[i][j] > mesh[j][local_mesh_max[j]]) {
        printf("point too low\n");
      }
    }
  }
}


/* Calculate the volume of a tetrahedron */

double tet_volume(double *corner[4])
{
   double v1[3], v2[3], v3[3];   /* vectors moved to origin */
   double volume;                /* computed and returned volume */
   int i;                        /* loop counter */

   for (i = 0; i < 3; ++i) {
      v1[i] = corner[1][i] - corner[0][i];
      v2[i] = corner[2][i] - corner[0][i];
      v3[i] = corner[3][i] - corner[0][i];
   }

   volume = (v1[0] * (v2[1] * v3[2] - v2[2] * v3[1]) -
             v1[1] * (v2[0] * v3[2] - v2[2] * v3[0]) +
             v1[2] * (v2[0] * v3[1] - v2[1] * v3[0])) / 6.0;

   if (volume < 0) volume = - volume;

   return(volume);
}
#define dist_sq(a, b, c, d, e, f) (fabs((a)-(d)) + fabs((b)-(e)) + fabs((c)-(f)))

/* Given a cutting plane, divide tet into set of smaller tets */

void slice_tet(
double *corner[4],		/* tetrahedron being divided */
int       cut,			/* x/y/z  = 0/1/2 direction of cut */
double    cval,			/* value of cut */
double    point[4][3],		/* intersections of plane with tets */
struct tet high_tet[3],	/* array of sub-tets from region above plane */
int      *pnhigh,		/* number of sub-tets above plane */
struct tet low_tet[3],		/* array of sub-tets from region above plane */
int      *pnlow 		/* number of sub-tets below plane */
)
{
    double    dist1, dist2;	/* distances between points */
    int       npoint;		/* counts number of intersection points */
    int       nplus, nminus, nzero;	/* number of tet points above/below/at cut */
    int       side[4];		/* side of the cut for each vertex */
    int       error;		/* checks for unanticipated case */
    int       zeros[2];		/* vertices on the cut */
    int       pluses[3];	/* vertices above the cut */
    int       minuses[3];	/* vertices below the cut */
    int       i;		/* loop counters */

    /* First figure out which side of the cut each vertex is on */
    nplus = nminus = nzero = 0;
    for (i = 0; i < 4; ++i) {
        if (corner[i][cut] > cval) {
	    side[i] = 1;
	    pluses[nplus++] = i;
        }
        else if (corner[i][cut] < cval) {
	    side[i] = -1;
	    minuses[nminus++] = i;
	}
        else {
	    side[i] = 0;
	    zeros[nzero++] = i;
        }
    }
    /* Check for consistency of point totals & +/- values */
    /* nplus & nminus should be at least 1 */
    /* nzero = 2  => 3 points */
    /* nzero = 1  => 3 points */
    /* 3/1 distribution => 3 points */
    /* 2/2 distribution => 4 points */
    error = FALSE;
    if (nplus >= 4 || nminus >= 4 || nzero > 2 || nplus + nminus + nzero != 4) {
	printf("Unexpected case in slice_tet: nplus = %d, nminus = %d\n", nplus, nminus);
	error = TRUE;
    }

    if (error) {
	printf("Cut dimension = %d, value = %.9g\n", cut, cval);
	printf("Corner values = %.9g %.9g %.9g %.9g\n",
	    corner[0][cut], corner[1][cut], corner[2][cut], corner[3][cut]);
	printf("Sides = %d %d %d %d\n", side[0], side[1], side[2], side[3]);
	exit(0);
    }


    /* Cases: (# tets <= 2n-7: stare from vtx to opposing faces )
       1 tet vertex & 3 points (easy) Only 1 tet
       2 tet vertices & 3 points Only 2 tets.
	 pick non-tet point for origin.
	 point + 3 tets, 3 points + 1 tet vtx
       2 tet vertices & 4 points
	  3 faces have 4 points, 2 faces have 3 => identical to below
       3 tet vertices & 3 points
	  Only 3 tets 'cause 2 opposing faces from view of a point
	  are degenerate. */

    if (nplus == 1) {
	/* upper side is single tet. */
        npoint = 0;
	for (i = 0; i < nminus; ++i) {
	    intersect_line(corner[pluses[0]], corner[minuses[i]], cut, cval, point[npoint]);
	    npoint++;
	}

	high_tet[0].corner[0] = corner[pluses[0]];
	high_tet[0].corner[1] = point[0];

	if (nzero > 1) high_tet[0].corner[2] = corner[zeros[1]];
	else high_tet[0].corner[2] = point[1];

	if (nzero > 0) high_tet[0].corner[3] = corner[zeros[0]];
	else high_tet[0].corner[3] = point[2];

	*pnhigh = 1;

	if (nminus == 1) {	/* nzero must be 2. */
	    /* Only one low tet */
	    low_tet[0].corner[0] = corner[minuses[0]];
	    low_tet[0].corner[1] = point[0];

	    if (nzero > 1) low_tet[0].corner[2] = corner[zeros[1]];
	    else low_tet[0].corner[2] = point[1];

	    if (nzero > 0) low_tet[0].corner[3] = corner[zeros[0]];
	    else low_tet[0].corner[3] = point[2];

	    *pnlow = 1;
	}
	else if (nminus == 2) {		/* nzero must be 1. */
	    /* Two low tets */
	    dist1 = dist_sq(corner[minuses[0]][0], corner[minuses[0]][1], 
	        corner[minuses[0]][2], point[1][0], point[1][1], point[1][2]);
	    dist2 = dist_sq(corner[minuses[1]][0], corner[minuses[1]][1], 
	        corner[minuses[1]][2], point[0][0], point[0][1], point[0][2]);

	    low_tet[0].corner[0] = corner[zeros[0]];
	    low_tet[0].corner[1] = point[0];
	    low_tet[0].corner[2] = point[1];

	    low_tet[1].corner[0] = corner[zeros[0]];
	    low_tet[1].corner[1] = corner[minuses[1]];
	    low_tet[1].corner[2] = corner[minuses[0]];

	    if (dist1 < dist2) {
	        low_tet[0].corner[3] = corner[minuses[0]];
	        low_tet[1].corner[3] = point[1];
	    }
	    else {
	        low_tet[0].corner[3] = corner[minuses[1]];
	        low_tet[1].corner[3] = point[0];
	    }

	    *pnlow = 2;
	}
	else {		/* nminus = 3 */
	    /* Three low tets */
	    low_tet[0].corner[0] = corner[minuses[0]];
	    low_tet[0].corner[1] = corner[minuses[1]];
	    low_tet[0].corner[2] = corner[minuses[2]];
	    low_tet[0].corner[3] = point[0];

	    low_tet[1].corner[0] = point[0];
	    low_tet[1].corner[1] = point[1];
	    low_tet[1].corner[2] = point[2];
	    low_tet[1].corner[3] = corner[minuses[1]];

	    low_tet[2].corner[0] = corner[minuses[1]];
	    low_tet[2].corner[1] = corner[minuses[2]];
	    low_tet[2].corner[2] = point[0];
	    low_tet[2].corner[3] = point[2];


	    *pnlow = 3;
	}
    }
    else if (nplus == 2) {
	if (nzero == 1) {	/* Only two tets */
	    /* 2 points */
	    intersect_line(corner[pluses[0]], corner[minuses[0]], cut, cval, point[0]);
	    intersect_line(corner[pluses[1]], corner[minuses[0]], cut, cval, point[1]);

	    dist1 = dist_sq(corner[pluses[0]][0], corner[pluses[0]][1], 
	        corner[pluses[0]][2], point[1][0], point[1][1], point[1][2]);
	    dist2 = dist_sq(corner[pluses[1]][0], corner[pluses[1]][1], 
	        corner[pluses[1]][2], point[0][0], point[0][1], point[0][2]);

	    high_tet[0].corner[0] = corner[zeros[0]];
	    high_tet[0].corner[1] = point[0];
	    high_tet[0].corner[2] = point[1];

	    high_tet[1].corner[0] = corner[zeros[0]];
	    high_tet[1].corner[1] = corner[pluses[1]];
	    high_tet[1].corner[2] = corner[pluses[0]];

	    if (dist1 < dist2) {
	        high_tet[0].corner[3] = corner[pluses[0]];
	        high_tet[1].corner[3] = point[1];
	    }
	    else {
	        high_tet[0].corner[3] = corner[pluses[1]];
	        high_tet[1].corner[3] = point[0];
	    }

	    *pnhigh = 2;

	    /* nminus = 1, only one low tet */

	    low_tet[0].corner[0] = corner[minuses[0]];
	    low_tet[0].corner[1] = point[0];
	    low_tet[0].corner[2] = point[1];
	    low_tet[0].corner[3] = corner[zeros[0]];

	    *pnlow = 1;
	}
	else if (nzero == 0) {	/* Three tets */
	    /* 4 points */
	    intersect_line(corner[pluses[0]], corner[minuses[0]], cut, cval, point[0]);
	    intersect_line(corner[pluses[0]], corner[minuses[1]], cut, cval, point[1]);
	    intersect_line(corner[pluses[1]], corner[minuses[0]], cut, cval, point[2]);
	    intersect_line(corner[pluses[1]], corner[minuses[1]], cut, cval, point[3]);

	    high_tet[0].corner[0] = corner[pluses[0]];
	    high_tet[0].corner[1] = corner[pluses[1]];
	    high_tet[0].corner[2] = point[0];
	    high_tet[0].corner[3] = point[1];

	    high_tet[1].corner[0] = corner[pluses[1]];
	    high_tet[1].corner[1] = point[2];
	    high_tet[1].corner[2] = point[3];
	    high_tet[1].corner[3] = point[0];

	    high_tet[2].corner[0] = point[0];
	    high_tet[2].corner[1] = point[1];
	    high_tet[2].corner[2] = corner[pluses[1]];
	    high_tet[2].corner[3] = point[3];

	    *pnhigh = 3;

	    low_tet[0].corner[0] = corner[minuses[0]];
	    low_tet[0].corner[1] = corner[minuses[1]];
	    low_tet[0].corner[2] = point[0];
	    low_tet[0].corner[3] = point[2];

	    low_tet[1].corner[0] = corner[minuses[1]];
	    low_tet[1].corner[1] = point[1];
	    low_tet[1].corner[2] = point[3];
	    low_tet[1].corner[3] = point[0];

	    low_tet[2].corner[0] = point[0];
	    low_tet[2].corner[1] = point[2];
	    low_tet[2].corner[2] = corner[minuses[1]];
	    low_tet[2].corner[3] = point[3];

	    *pnlow = 3;
	}
    }
    else if (nplus == 3) {	/* Three tets */
	/* should be 3 points */
	intersect_line(corner[pluses[0]], corner[minuses[0]], cut, cval, point[0]);
	intersect_line(corner[pluses[1]], corner[minuses[0]], cut, cval, point[1]);
	intersect_line(corner[pluses[2]], corner[minuses[0]], cut, cval, point[2]);

	high_tet[0].corner[0] = corner[pluses[0]];
	high_tet[0].corner[1] = corner[pluses[1]];
	high_tet[0].corner[2] = corner[pluses[2]];
	high_tet[0].corner[3] = point[0];

	high_tet[1].corner[0] = point[0];
	high_tet[1].corner[1] = point[1];
	high_tet[1].corner[2] = point[2];
	high_tet[1].corner[3] = corner[pluses[1]];

	high_tet[2].corner[0] = corner[pluses[1]];
	high_tet[2].corner[1] = corner[pluses[2]];
	high_tet[2].corner[2] = point[0];
	high_tet[2].corner[3] = point[2];

	*pnhigh = 3;

	/* Only a single low tet */
	low_tet[0].corner[0] = corner[minuses[0]];
	low_tet[0].corner[1] = point[0];
	low_tet[0].corner[2] = point[1];
	low_tet[0].corner[3] = point[2];

	*pnlow = 1;
    }
}


static void intersect_line(
double p0[3],			/* segment end point */
double p1[3],			/* segment end point */
int cut,			/* dimension being sliced */
double cval,			/* value at which slice occurs */
double point[3] 		/* returned intersection point */
)
{
    double ratio;		/* fractional distance along segment */
    int i;			/* loop counter */

    /* Choose to take differences to keep ratio from being small. */
    if (p0[cut] - cval >= cval - p1[cut]) {
	ratio = (p0[cut] - cval) / (p0[cut] - p1[cut]);
        for (i = 0; i < 3; ++i) {
	    if (i != cut) point[i] = p0[i] - ratio * (p0[i] - p1[i]);
	}
    }
    else {
	ratio = (cval - p1[cut]) / (p0[cut] - p1[cut]);
        for (i = 0; i < 3; ++i) {
	    if (i != cut) point[i] = p1[i] + ratio * (p0[i] - p1[i]);
	}
    }
    point[cut] = cval;
}
