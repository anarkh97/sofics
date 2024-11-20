#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <complex>
#include <cmath>

int main() {

	char *filename = (char *) "out";
	char buf[256];
	double k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16;
	int toto = scanf("%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &k1, &k2, &k3, &k4, &k5, &k6, &k7, &k8, &k9,
	                 &k10, &k11, &k12, &k13, &k14, &k15, &k16);

	int i;
	FILE *fin = fopen("pdep", "r");
	int dep[16];
	int indep[16];
	int ndep;
	int nindep;
	toto = fscanf(fin, "%d", &nindep);
	for (i = 0; i < nindep; i++) toto = fscanf(fin, "%d", indep + i);
	for (i = 0; i < nindep; i++) indep[i]--;
	toto = fscanf(fin, "%d", &ndep);
	for (i = 0; i < ndep; i++) toto = fscanf(fin, "%d", dep + i);
	for (i = 0; i < ndep; i++) dep[i]--;

	double *pdep = new double[(nindep + 1) * ndep];
	for (i = 0; i < (nindep + 1) * ndep; i++) toto = fscanf(fin, "%lf", pdep + i);

	strcpy(buf, filename);
	strcat(buf, ".sgi");

	FILE *nodes = fopen(buf, "r");
	if (nodes == NULL) {
		fprintf(stderr, "Could not open the node file.\n");
		return -1;
	}

	int j;
	int nno = 0;
	char buffer[1024];
	char *toto2 = fgets(buffer, 256, nodes);
	while (!feof(nodes)) {
		int dummy;
		char *p = fgets(buffer, 1024, nodes);
		int st = sscanf(buffer, "%d", &dummy);
		if (st == 0) break;
		if (p == NULL) break;
		nno++;
	}
	fseek(nodes, 0, 0);
	toto2 = fgets(buffer, 256, nodes);
	double (*no)[3] = new double[nno][3];
	for (i = 0; i < nno; i++) {
		int dummy;
		toto2 = fgets(buffer, 1024, nodes);
		int st = sscanf(buffer, "%d%lf%lf%lf", &dummy, &no[i][0], &no[i][1], &no[i][2]);
		if (st == 0) break;
	}

	fprintf(stderr, "Read %d nodes.\n", nno);


	toto2 = fgets(buffer, 256, nodes);
	toto2 = fgets(buffer, 256, nodes);

	int nel = 0;
	int offset = ftell(nodes);
	while (!feof(nodes)) {
		int dummy;
		char *p = fgets(buffer, 1024, nodes);
		int st = sscanf(buffer, "%d", &dummy);
		if (st == 0) break;
		if (p == NULL) break;
		nel++;
	}

	fseek(nodes, offset, 0);
	int (*el)[4] = new int[nel][4];
	for (i = 0; i < nel; i++) {
		int dummy;
		toto2 = fgets(buffer, 1024, nodes);
		int st = sscanf(buffer, "%d%d%d%d%d%d", &dummy, &dummy, &el[i][0], &el[i][1], &el[i][2], &el[i][3]);
		if (st == 0) break;
		el[i][0]--;
		el[i][1]--;
		el[i][2]--;
		el[i][3]--;
	}

	fprintf(stderr, "Read %d elements.\n", nel);


	toto2 = fgets(buffer, 256, nodes);
	toto2 = fgets(buffer, 256, nodes);

	int nsubel = 0;
	offset = ftell(nodes);
	while (!feof(nodes)) {
		int dummy;
		char *p = fgets(buffer, 1024, nodes);
		int st = sscanf(buffer, "%d", &dummy);
		if (st == 0) break;
		if (p == NULL) break;
		nsubel++;
	}

	fseek(nodes, offset, 0);
	int (*subel)[3] = new int[nsubel][3];
	for (i = 0; i < nsubel; i++) {
		int dummy;
		toto2 = fgets(buffer, 1024, nodes);
		int st = sscanf(buffer, "%d%d%d%d%d", &dummy, &dummy, &subel[i][0], &subel[i][1], &subel[i][2]);
		if (st == 0) break;
		subel[i][0]--;
		subel[i][1]--;
		subel[i][2]--;
	}

	fprintf(stderr, "Read %d elements.\n", nsubel);


	toto2 = fgets(buffer, 256, nodes);

	int nextel = 0;
	offset = ftell(nodes);
	while (!feof(nodes)) {
		int dummy;
		char *p = fgets(buffer, 1024, nodes);
		int st = sscanf(buffer, "%d", &dummy);
		if (st == 0) break;
		if (p == NULL) break;
		nextel++;
	}

	fseek(nodes, offset, 0);
	int (*extel)[3] = new int[nextel][3];
	for (i = 0; i < nextel; i++) {
		int dummy;
		toto2 = fgets(buffer, 1024, nodes);
		int st = sscanf(buffer, "%d%d%d%d%d", &dummy, &dummy, &extel[i][0], &extel[i][1], &extel[i][2]);
		if (st == 0) break;
		extel[i][0]--;
		extel[i][1]--;
		extel[i][2]--;
	}

	fprintf(stderr, "Read %d elements.\n", nextel);


	toto2 = fgets(buffer, 256, nodes);

	int ncylel = 0;
	offset = ftell(nodes);
	while (!feof(nodes)) {
		int dummy;
		char *p = fgets(buffer, 1024, nodes);
		int st = sscanf(buffer, "%d", &dummy);
		if (st == 0) break;
		if (p == NULL) break;
		ncylel++;
	}

	fseek(nodes, offset, 0);
	int (*cylel)[3] = new int[ncylel][3];
	for (i = 0; i < ncylel; i++) {
		int dummy;
		toto2 = fgets(buffer, 1024, nodes);
		int st = sscanf(buffer, "%d%d%d%d%d", &dummy, &dummy, &cylel[i][0], &cylel[i][1], &cylel[i][2]);
		if (st == 0) break;
		cylel[i][0]--;
		cylel[i][1]--;
		cylel[i][2]--;
	}

	fprintf(stderr, "Read %d elements.\n", ncylel);


	toto2 = fgets(buffer, 256, nodes);

	int ntailel = 0;
	offset = ftell(nodes);
	while (!feof(nodes)) {
		int dummy;
		char *p = fgets(buffer, 1024, nodes);
		int st = sscanf(buffer, "%d", &dummy);
		if (st == 0) break;
		if (p == NULL) break;
		ntailel++;
	}

	fseek(nodes, offset, 0);
	int (*tailel)[3] = new int[ntailel][3];
	for (i = 0; i < ntailel; i++) {
		int dummy;
		toto2 = fgets(buffer, 1024, nodes);
		int st = sscanf(buffer, "%d%d%d%d%d", &dummy, &dummy, &tailel[i][0], &tailel[i][1], &tailel[i][2]);
		if (st == 0) break;
		tailel[i][0]--;
		tailel[i][1]--;
		tailel[i][2]--;
	}

	fprintf(stderr, "Read %d elements.\n", ntailel);


	toto2 = fgets(buffer, 256, nodes);

	int nheadel = 0;
	offset = ftell(nodes);
	while (!feof(nodes)) {
		int dummy;
		char *p = fgets(buffer, 1024, nodes);
		int st = sscanf(buffer, "%d", &dummy);
		if (st == 0) break;
		if (p == NULL) break;
		nheadel++;
	}

	fseek(nodes, offset, 0);
	int (*headel)[3] = new int[nheadel][3];
	for (i = 0; i < nheadel; i++) {
		int dummy;
		toto2 = fgets(buffer, 1024, nodes);
		int st = sscanf(buffer, "%d%d%d%d%d", &dummy, &dummy, &headel[i][0], &headel[i][1], &headel[i][2]);
		if (st == 0) break;
		headel[i][0]--;
		headel[i][1]--;
		headel[i][2]--;
	}

	fprintf(stderr, "Read %d elements.\n", nheadel);


	toto2 = fgets(buffer, 256, nodes);

	int ntowerel = 0;
	offset = ftell(nodes);
	while (!feof(nodes)) {
		int dummy;
		char *p = fgets(buffer, 1024, nodes);
		int st = sscanf(buffer, "%d", &dummy);
		if (st == 0) break;
		if (p == NULL) break;
		ntowerel++;
	}

	fseek(nodes, offset, 0);
	int (*towerel)[3] = new int[ntowerel][3];
	for (i = 0; i < ntowerel; i++) {
		int dummy;
		toto2 = fgets(buffer, 1024, nodes);
		int st = sscanf(buffer, "%d%d%d%d%d", &dummy, &dummy, &towerel[i][0], &towerel[i][1], &towerel[i][2]);
		if (st == 0) break;
		towerel[i][0]--;
		towerel[i][1]--;
		towerel[i][2]--;
	}

	fprintf(stderr, "Read %d elements.\n", ntowerel);


	toto2 = fgets(buffer, 256, nodes);

	int nwingsel = 0;
	offset = ftell(nodes);
	while (!feof(nodes)) {
		int dummy;
		char *p = fgets(buffer, 1024, nodes);
		int st = sscanf(buffer, "%d", &dummy);
		if (st == 0) break;
		if (p == NULL) break;
		nwingsel++;
	}

	fseek(nodes, offset, 0);
	int (*wingsel)[3] = new int[nwingsel][3];
	for (i = 0; i < nwingsel; i++) {
		int dummy;
		toto2 = fgets(buffer, 1024, nodes);
		int st = sscanf(buffer, "%d%d%d%d%d", &dummy, &dummy, &wingsel[i][0], &wingsel[i][1], &wingsel[i][2]);
		if (st == 0) break;
		wingsel[i][0]--;
		wingsel[i][1]--;
		wingsel[i][2]--;
	}

	fprintf(stderr, "Read %d elements.\n", nwingsel);


	toto2 = fgets(buffer, 256, nodes);

	int nttopel = 0;
	offset = ftell(nodes);
	while (!feof(nodes)) {
		int dummy;
		char *p = fgets(buffer, 1024, nodes);
		int st = sscanf(buffer, "%d", &dummy);
		if (st == 0) break;
		if (p == NULL) break;
		nttopel++;
	}

	fseek(nodes, offset, 0);
	int (*ttopel)[3] = new int[nttopel][3];
	for (i = 0; i < nttopel; i++) {
		int dummy;
		toto2 = fgets(buffer, 1024, nodes);
		int st = sscanf(buffer, "%d%d%d%d%d", &dummy, &dummy, &ttopel[i][0], &ttopel[i][1], &ttopel[i][2]);
		if (st == 0) break;
		ttopel[i][0]--;
		ttopel[i][1]--;
		ttopel[i][2]--;
	}

	fprintf(stderr, "Read %d elements.\n", nttopel);


	int (*ndunion)[7] = new int[nno][7];
	for (i = 0; i < nno; i++) {
		ndunion[i][0] = ndunion[i][1] = ndunion[i][2] = ndunion[i][3] = ndunion[i][4] = ndunion[i][5] = ndunion[i][6] = 0;
	}
	for (i = 0; i < nsubel; i++) {
		ndunion[subel[i][0]][6]++;
		ndunion[subel[i][1]][6]++;
		ndunion[subel[i][2]][6]++;
	}
	for (i = 0; i < ncylel; i++) {
		ndunion[cylel[i][0]][0]++;
		ndunion[cylel[i][1]][0]++;
		ndunion[cylel[i][2]][0]++;
	}
	for (i = 0; i < ntailel; i++) {
		ndunion[tailel[i][0]][1]++;
		ndunion[tailel[i][1]][1]++;
		ndunion[tailel[i][2]][1]++;
	}
	for (i = 0; i < nheadel; i++) {
		ndunion[headel[i][0]][2]++;
		ndunion[headel[i][1]][2]++;
		ndunion[headel[i][2]][2]++;
	}
	for (i = 0; i < ntowerel; i++) {
		ndunion[towerel[i][0]][3]++;
		ndunion[towerel[i][1]][3]++;
		ndunion[towerel[i][2]][3]++;
	}
	for (i = 0; i < nttopel; i++) {
		ndunion[ttopel[i][0]][4]++;
		ndunion[ttopel[i][1]][4]++;
		ndunion[ttopel[i][2]][4]++;
	}
	for (i = 0; i < nwingsel; i++) {
		ndunion[wingsel[i][0]][5]++;
		ndunion[wingsel[i][1]][5]++;
		ndunion[wingsel[i][2]][5]++;
	}


	for (i = 0; i < nextel; i++) {
		int tmp;
		tmp = extel[i][1];
		extel[i][1] = extel[i][2];
		extel[i][2] = tmp;
	}

	for (i = 0; i < nsubel; i++) {
		int tmp;
		tmp = subel[i][1];
		subel[i][1] = subel[i][2];
		subel[i][2] = tmp;
	}


	int nsn = 0;
	j = 0;
	for (i = 0; i < nno; i++) {
		if (ndunion[i][6] > 0) nsn++;
		if (ndunion[i][0] > 0 || ndunion[i][1] > 0 || ndunion[i][2] > 0 ||
		    ndunion[i][3] > 0 || ndunion[i][4] > 0 || ndunion[i][5] > 0)
			j++;
	}
	if (nsn != j) fprintf(stderr, "sub neq cyl .+ tail .+ head .+ tower %d %d\n", nsn, j);

	int *subnd = new int[nsn];
	nsn = 0;
	for (i = 0; i < nno; i++) {
		if (ndunion[i][6] > 0) subnd[nsn++] = i;
	}


	double (*deriv)[6][16][3] = new double[nsn][6][16][3];
	int l, k;
	for (i = 0; i < nsn; i++)
		for (j = 0; j < 6; j++)
			for (k = 0; k < 16; k++)
				for (l = 0; l < 3; l++)
					deriv[i][j][k][l] = 0.0;

	for (i = 0; i < nsn; i++) {
		double x = no[subnd[i]][0];
		double y = no[subnd[i]][1];
		double z = no[subnd[i]][2];
		if (ndunion[subnd[i]][0] > 0) {
			double e = y * y / k2 / k2 + z * z / k3 / k3 - 1;
//     fprintf(stderr,"%d %d %f\n",0,subnd[i]+1,e);
			deriv[i][0][2 - 1][1] = y / k2;
			deriv[i][0][3 - 1][2] = z / k3;
		}
		if (ndunion[subnd[i]][1] > 0) {
			double e = (x + k1) * (x + k1) / k16 / k16 + y * y / k2 / k2 + z * z / k3 / k3 - 1.0;
//     fprintf(stderr,"%d %d %f\n",1,subnd[i]+1,e);
			deriv[i][1][1 - 1][0] = -1.0;
			deriv[i][1][16 - 1][0] = (x + k1) / k16;
			deriv[i][1][2 - 1][1] = y / k2;
			deriv[i][1][3 - 1][2] = z / k3;
		}
		if (ndunion[subnd[i]][2] > 0) {
			double e = (x - k1) * (x - k1) / k15 / k15 + y * y / k2 / k2 + z * z / k3 / k3 - 1.0;
//     fprintf(stderr,"%d %d %f\n",2,subnd[i]+1,e);
			deriv[i][2][1 - 1][0] = 1.0;
			deriv[i][2][15 - 1][0] = (x - k1) / k15;
			deriv[i][2][2 - 1][1] = y / k2;
			deriv[i][2][3 - 1][2] = z / k3;
		}
		if (ndunion[subnd[i]][3] > 0) {
			double e = (x - k8) * (x - k8) / k4 / k4 + y * y / k5 / k5 + z * z / k6 / k6 - 1.0;
//     fprintf(stderr,"%d %d %f\n",3,subnd[i]+1,e);
			deriv[i][3][4 - 1][0] = (x - k8) / k4;
			deriv[i][3][5 - 1][1] = y / k5;
			deriv[i][3][6 - 1][2] = z / k6;
			deriv[i][3][8 - 1][0] = 1.0;
		}
		if (ndunion[subnd[i]][4] > 0) {
			double e = y - k7;
//     fprintf(stderr,"%d %d %f\n",4,subnd[i]+1,e);
			deriv[i][4][7 - 1][1] = 1.0;
			deriv[i][4][8 - 1][0] = 1.0;
		}
		if (ndunion[subnd[i]][5] > 0) {
			double e = (x - k8) * (x - k8) / k9 / k9 + (y - k12) * (y - k12) / k10 / k10 + z * z / k11 / k11 - 1.0;
//     fprintf(stderr,"%d %d %f\n",5,subnd[i]+1,e);
			deriv[i][5][8 - 1][0] = 1.0;
			deriv[i][5][9 - 1][0] = (x - k8) / k9;
			deriv[i][5][10 - 1][1] = (y - k12) / k10;
			deriv[i][5][11 - 1][2] = z / k11;
			deriv[i][5][12 - 1][1] = 1.0;
		}
	}

	FILE *fem = fopen("fem", "w");
	fprintf(fem, "NODES\n");
	for (i = 0; i < nno; i++) {
		fprintf(fem, "%d %.6f %.6f %.6f\n", i + 1, no[i][0], no[i][1], no[i][2]);
	}
	fprintf(fem, "TOPO\n");
	for (i = 0; i < nel; i++) {
		fprintf(fem, "%d 40 %d %d %d %d\n", i + 1, el[i][0] + 1, el[i][1] + 1, el[i][2] + 1, el[i][3] + 1);
	}
	fprintf(fem, "ATTRIB\n1 %d 1\n", nel);
	fprintf(fem, "HARBO\n");
	for (i = 0; i < nextel; i++) {
		fprintf(fem, "%d 3 %d %d %d\n", i + 1, extel[i][0] + 1, extel[i][1] + 1, extel[i][2] + 1);
	}
	fprintf(fem, "HSCBO\n");
	for (i = 0; i < nsubel; i++) {
		fprintf(fem, "%d 3 %d %d %d\n", i + 1, subel[i][0] + 1, subel[i][1] + 1, subel[i][2] + 1);
	}
	fprintf(fem, "HDIR\n");
	for (i = 0; i < nsn; i++) {
		fprintf(fem, "%d 8 1.0 0.0\n", subnd[i] + 1);
	}


	fprintf(fem, "FUNSCATTER %d %d\n", nsn, nindep);
	for (i = 0; i < nsn; i++) {
		for (int jj = 0; jj < nindep; jj++) {
			j = indep[jj];
			int c = 0;
			double d[3] = {0, 0, 0};
			for (l = 0; l < 6; l++) {
				if (ndunion[subnd[i]][l] > 0) {
					d[0] += deriv[i][l][j][0];
					d[1] += deriv[i][l][j][1];
					d[2] += deriv[i][l][j][2];
					for (int kk = 0; kk < ndep; kk++) {
						int k = dep[kk];
						d[0] += deriv[i][l][k][0] * pdep[(nindep + 1) * kk + jj];
						d[1] += deriv[i][l][k][1] * pdep[(nindep + 1) * kk + jj];
						d[2] += deriv[i][l][k][2] * pdep[(nindep + 1) * kk + jj];
					}
					c++;
				}
			}
			d[0] /= c;
			d[1] /= c;
			d[2] /= c;
			fprintf(fem, "%.3f %.3f %.3f ", d[0], d[1], d[2]);
		}
		fprintf(fem, "\n");
	}

	fclose(nodes);
	return 0;
}
