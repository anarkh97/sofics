
//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  Michel Lesoinne                     Denis Vanderstraeten            //
//  June 1996                           January 1993                    //
//                                                                      //
//  Department of Aerospace Eng         Unite de Mecanique Appliquee    //
//  College of Engineering              Universite Catholique de Louvain//
//  University of Colorado              Louvain-la-Neuve - Belgium      //
//  Boulder - USA                                                       //
//                                                                      //
//									//
//  Note : The original code is due to Charbel Farhat			//
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdio>
#include <Utils.d/Connectivity.h>
#include <Dec.d/Decomp.d/Decomp.h>

void nomask( Connectivity* mn, int* maskno, int nelem ) {

int i;
for (i = 0; i < mn->num(nelem); i++)  maskno[(*mn)[nelem][i]]--;
}

// ----------------------------------------------------------------------

struct DDLink {
    int previous, next;
};

class IList {
    int nn;
    int *n;
    int *map;
   public:
     IList(int maxSize);
     ~IList() { delete[]n; }
     int operator[](int i) { return n[i]; }
     void add(int);
     void remove(int);
     int size() { return nn; }
};

IList::IList(int maxSize)
{
 nn = 0;
 n = new int[2*maxSize];
 map = n +maxSize;
}

void
IList::add(int x)
{
 //map[xx] = nn;
 n[nn++] = x;
}

void
IList::remove(int x)
{
 map[n[nn]] = map[x];
 n[map[x]] = n[nn];
 map[x] = -1;
 nn--;
}

void Decomposition::greedyBuild(int nsub, Connectivity *etoe,
                                Connectivity *ntoe, Connectivity *eton) {

	int i,minw,minw2;
	int node,node2,ntot,nelem, ptr1, ptr2;
	int *maskno, *elsub, *me, *color;
	int *newelsub;

	int exact_num = etoe->numNonZeroP();
	int numele    = etoe->csize();
	int numnod    = ntoe->csize();
	elsub = new int [nsub+1];
	newelsub = new int[nsub+2];
	me = new int[exact_num+1];
	color = new int[numele];
	maskno = new int[numnod];

//  initialize newelsub as a pointer into me 

	newelsub[0] = 0;
	int div = exact_num / nsub;
	int rem = exact_num % nsub;

// Compute an a priori balanced distribution of element numbers
	for( i=0; i < nsub; i++ )
		elsub[i] = (i < rem) ? div+1 : div;

	for( i=0; i < nsub; i++ )
		newelsub[i+1] = newelsub[i] + elsub[i];

	for (i = 0; i < numele; i++) color[i] = -1;

	for (i = 0; i < numnod; i++) maskno[i] = ntoe->num(i);

	ptr1 = 0;
	ptr2 = 0;
	int numsub = 0;
	ntot = 0;


	while (ntot < elsub[numsub] && numsub < nsub) {
// Locate a starting (interface if possible) node
// The part of code is in O(Numnod). It could be done in
// O(I) (I = interface size) by handling properly the list of
// interface nodes.

		minw  = minw2 = 32000000;
		node  = -1;
		for (i = 0; i < numnod; i++) {
			if(maskno[i] == 0) continue;
			if(maskno[i] == ntoe->num(i)) {
				if(maskno[i] < minw2) {
					minw2 = maskno[i];
					node2 = i;
				}
			} else {
				if(maskno[i] < minw) {
					minw = maskno[i];
					node = i;
				}
			}
		}

		if (node < 0) node = node2;

// Initialize the list of elements atached to the starting node

		ptr1 = ptr2;
		for (i = 0; i < ntoe->num(node) && ntot < elsub[numsub]; i++) {

			/* find an element attached to "node" */
			nelem = (*ntoe)[node][i];
			if(color[nelem] >= 0) continue;
			color[nelem] = numsub;
			me[ptr2] = nelem;
			ptr2++;
			ntot++;

			// reduce mask value for all nodes connected to this element
			nomask( eton, maskno, nelem);

		}

/* RECURSIVELY ADD TO LIST NEW ADJACENT ELEMENTS */


		while (ptr2 > ptr1 && ntot<elsub[numsub]){
			for (i = 0; i<etoe->num(me[ptr1]) && ntot<elsub[numsub];i++) {
				nelem = (*etoe)[me[ptr1]][i];
				if(color[nelem] >= 0) continue;
				color[nelem] = numsub;
				me[ptr2] = nelem;
				ntot++;
				ptr2++;
				nomask( eton, maskno, nelem);
			}
			ptr1++;
		}

		if(ntot == elsub[numsub]) {
			numsub++;
			ntot = 0;
		}
	}

	for (i = 0; i < numele; i++){
		if (etoe->num(i) && color[i] < 0) {me[ptr2++] = i; color[i] = nsub;}

	}

	pele = newelsub;
	eln = me;

// delete [] maskno;
// delete [] color;
	delete [] elsub;
}


void
Decomposition::outputDump(FILE *outFile, int type)
{
 //fprintf(outFile,"Decomposition %s for %s\n", Name(), esname);
 fprintf(outFile,"Decomposition %s for %s\n", "dec", esname);
 fprintf(outFile, " %d\n",nsub);
 int isub, iEle;
 for(isub=0; isub < nsub; ++isub) {
   fprintf(outFile, " %d\n",pele[isub+1]-pele[isub]);
   for(iEle = pele[isub]; iEle < pele[isub+1]; ++iEle)
     fprintf(outFile, "%d\n",eln[iEle]+1);
 }
}


// ----------------------------------------------------------------------

void Decomposition::cpuGreedyBuild(int nsub, Connectivity *etoe,
                      Connectivity *ntoe, Connectivity *eton,
                      long *sizes, int numSubdomains) {

int i,minw,minw2;
int node,node2,ntot,nelem, ptr1, ptr2;
int *maskno, *elsub, *me, *color;
int *newelsub;

// Compute total profile size
int totalSize = 0;
for(i=0; i<numSubdomains; ++i)
  totalSize += sizes[i];

double percent = 0.95;

int exact_num = etoe->numNonZeroP();
int numele    = etoe->csize();
int numnod    = ntoe->csize();
elsub = new int [nsub+1];
newelsub = new int[nsub+2];
me = new int[exact_num+1];
color = new int[numele];
maskno = new int[numnod];

//  initialize newelsub as a pointer into me 

newelsub[0] = 0;

int div = totalSize / nsub;
int rem = totalSize % nsub;

fprintf(stderr,"Per Cluster %d Remainder %d\n",div,rem);

int div1 = exact_num / nsub;
int rem1 = exact_num % nsub;

// Compute an a priori balanced distribution of element numbers
for( i=0; i < nsub; i++ ) {
   elsub[i] = (i < rem1) ? div1+1 : div1;
}

for( i=0; i < nsub; i++ ) {
   newelsub[i+1] = newelsub[i] + elsub[i];
}

for( i=0; i < nsub; i++ ) {
   elsub[i] = (i < rem) ? div+1 : div;
}

for (i = 0; i < numele; i++) color[i] = -1;

for (i = 0; i < numnod; i++) maskno[i] = ntoe->num(i);

ptr1 = 0;
ptr2 = 0;
int numsub = 0;
ntot = 0;

while ( ntot <= percent*elsub[numsub] && numsub < nsub ) {

// Locate a starting (interface if possible) node
// The part of code is in O(Numnod). It could be done in
// O(I) (I = interface size) by handling properly the list of
// interface nodes.

   minw  = minw2 = 32000000;
   node  = -1;
   for (i = 0; i < numnod; i++) {
     if(maskno[i] == 0) continue;
     if(maskno[i] == ntoe->num(i)) {
       if(maskno[i] < minw2) {
        minw2 = maskno[i];
        node2 = i;
       }
     } else {
       if(maskno[i] < minw) {
           minw = maskno[i];
           node = i;
       }
     }
   }
 
   if (node < 0) node = node2;

   // Initialize the list of elements atached to the starting node

   ptr1 = ptr2;
   for (i = 0; i < ntoe->num(node) && ntot < percent*elsub[numsub]; i++) {
 
      /* find an element attached to "node" */
      nelem = (*ntoe)[node][i];
      if(color[nelem] >= 0) continue;
      color[nelem] = numsub;
      me[ptr2] = nelem;
      ptr2++;
      ntot += sizes[nelem]; // MODIFICATION
      // reduce mask value for all nodes connected to this element 
      nomask( eton, maskno, nelem);
   }
 
/* RECURSIVELY ADD TO LIST NEW ADJACENT ELEMENTS */

  while (ptr2 > ptr1 && ntot < percent*elsub[numsub] ){
   for (i = 0; i<etoe->num(me[ptr1]) && ntot<percent*elsub[numsub];i++) {
        nelem = (*etoe)[me[ptr1]][i];
        if(color[nelem] >= 0) continue;
        color[nelem] = numsub;
        me[ptr2] = nelem;
        ntot += sizes[nelem]; // MODIFICATION
        ptr2++;
        nomask( eton, maskno, nelem);
    }
    ptr1++;
  }
 
  // MODIFICATION
  if( ntot >= percent*elsub[numsub] && ntot <= elsub[numsub]/percent) {
    numsub++;
    ntot = 0;
  } 
}

for (i = 0; i < numele; i++){ 
     if (etoe->num(i) && color[i] < 0) 
        {me[ptr2++] = i; color[i] = nsub;}

}

   pele = newelsub;
   eln  = me;

// delete [] maskno;
// delete [] color;

delete [] elsub;

}

