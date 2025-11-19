#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "random.h"
#include <assert.h>
#include <stdbool.h>
#include <windows.h>
#include <omp.h>

// 3 dimensional Ising model, # of nearest neighbors == 6

#define J		(1e-3*1.60219e-19)	
#define mu_B	(9.2741e-24)		//Bohr magneton
#define mu_0	(1.2566e-6)			//Vacuum permiability
#define kB		(1.38062e-23)		//Boltzmann factor

#define N		25 // 2d: 100, 3d: 25
int NN = N + 1;
double N_sq = ((double)N * (double)N);
double N_cb = (((double)N * (double)N) * (double)N);

#define nmax	10000
#define npat	10

#define hrng	250.0
#define hstp	0.1
#define nearone 0.999999999999999999


// ---------------------------------------------------------
char*** alloc_3d() {

	char*** ptr_dum = (char***)malloc(sizeof(char**) * (N + 2));
	if (ptr_dum == NULL) { fprintf(stderr, "malloc error\n"); return NULL; }
	for (int i = 0; i < N + 2; i++) {
		
		(ptr_dum)[i] = (char**)malloc(sizeof(char*) * (N + 2));
		if ((ptr_dum)[i] == NULL) { fprintf(stderr, "malloc error\n"); return NULL; }
		for (int j = 0; j < N + 2; j++) {

			(ptr_dum)[i][j] = (char*)malloc(sizeof(char) * (N + 2));
			if ((ptr_dum)[i][j] == NULL) { fprintf(stderr, "malloc error\n"); return NULL; }
		}
	}
	return ptr_dum;
}
int free_3d(char*** s_3d) {
	for (int i = 0; i < N + 2; i++) {
		for (int j = 0; i < N + 2; i++) {
			free(s_3d[i][j]);
		}
		free(s_3d[i]);
	}
	free(s_3d);
	return 0;
}

int initialize_3d(char*** s_3d) {
	for (int i = 0; i < (N + 2); i++) {
		for (int j = 0; j < (N + 2); j++) {
			for (int k = 0; k < (N + 2); k++) {
				if (i == 0 || i == NN || j == 0 || j == NN || k == 0 || k == NN)
					{ s_3d[i][j][k] = 0; }
				else
					{ s_3d[i][j][k] = 1; }
			}
		}
	}
	return 0;
}
int iter_3d(char*** s_3d, double T, int accept_more, double r, bool periodic, bool rand_parity_order){

	double delta_E[13];
	for(int l=0;l<13;l++){delta_E[l] = (J*((double)l-6));}
	// arr[l] = J*(l-6)
	double prob[2][13] = { 0, };
	double dum;
	
	for (int c = 0; c <= 1; c++) {
	for (int l = 0; l <= 12; l++) {
		if		(accept_more == 2) {
				dum = 2.0 * (double)(2 * (double)c - 1) * delta_E[l] / (kB * T);
				if (dum < 0) { dum = 0.0; }
		}
		else if (accept_more == 1) { dum = ((double)(2 * (double)c - 1) * delta_E[l] + 6.0 * J) / (kB * T); }
		else if (accept_more == 0) { dum = ((double)(2 * (double)c - 1) * delta_E[l] + 12.0 * J) / (kB * T); }
		else { dum = 1.0; }

		prob[c][l] = r * exp(-1.0 * dum);
	}
	}
	// prob[i][j] = Accptance ratio for flipping s[a]==(2i-1) with
	// sum_<a,b> s[b] = 2j-6.
	
	int flipped = 0;
	int i;
	int order[2] = { 0, 0 };
	int idx, parity;

	if (genrand64_real1() < 0.5 && rand_parity_order) { order[0] = 1; order[1] = 0; }
	else											  { order[0] = 0; order[1] = 1; }

	for (idx = 0; idx <= 1; idx++) {
		parity = order[idx];
	#pragma omp parallel shared(s_3d,parity) reduction(+:flipped)
	#pragma omp for	schedule(static)
		for (i = 1; i <= N; i++) {

			for (int j = 1; j <= N; j++) {
			for (int k = 1; k <= N; k++) {
			if ((i + j + k) % 2 == parity) {
				int count = 0;

				if (periodic) {
					count += s_3d[(i == N) ? 1 : i + 1][j][k];
					count += s_3d[(i == 1) ? N : i - 1][j][k];
					count += s_3d[i][(j == N) ? 1 : j + 1][k];
					count += s_3d[i][(j == 1) ? N : j - 1][k];
					count += s_3d[i][j][(k == N) ? 1 : k + 1];
					count += s_3d[i][j][(k == 1) ? N : k - 1];
				}
				else {
					count += s_3d[i + 1][j][k];
					count += s_3d[i - 1][j][k];
					count += s_3d[i][j + 1][k];
					count += s_3d[i][j - 1][k];
					count += s_3d[i][j][k + 1];
					count += s_3d[i][j][k - 1];
				}

				double p = prob[(s_3d[i][j][k] + 1) / 2][(count + 6)];
				
				if (genrand64_real1() < p) {
					s_3d[i][j][k] *= -1;
					flipped += 1;
				}

			}	
			}
			}
		}
	}

	return flipped;
} // iter function returns the number of flipped spins.

double get_avg_3d(char*** s_3d){

	double m = 0;
	int i = 0;
	#pragma omp parallel shared(s_3d) reduction(+:m)
	#pragma omp for	schedule(static)
		for (i = 1; i <= N; i++) {
			for (int j = 1; j <= N; j++) { for (int k = 1; k <= N; k++) { m += s_3d[i][j][k]; } } 
		}
	m /= N_cb;
	return m;
}
// ---------------------------------------------------------
char** alloc_2d() {

	char** ptr_dum = (char**)malloc(sizeof(char*) * (N + 2));
	if (ptr_dum == NULL) { fprintf(stderr, "malloc error\n"); return NULL; }
		for (int i = 0; i < N + 2; i++) {
			(ptr_dum)[i] = (char*)malloc(sizeof(char) * (N + 2));
			if ((ptr_dum)[i] == NULL) { fprintf(stderr, "malloc error\n"); return NULL; }
		}
	return ptr_dum;
}
int free_2d(char** s_2d) {
	for (int i = 0; i < N+2; i++) {
		free(s_2d[i]);
	}
	free(s_2d);
	return 0;
}
void copy_config_2d(char** s_2d_original, char** s_2d_copy) {

	int i = 0;
#pragma omp parallel shared(s_2d_original, s_2d_copy)
#pragma omp for	schedule(static)
	for (i = 0; i < (N+2); i++) {
		for (int j = 0; j < (N+2); j++) { s_2d_copy[i][j] = s_2d_original[i][j]; }
	}

}

int initialize_2d(char** s_2d) {
	for (int i = 0; i < (N + 2); i++) {
		for (int j = 0; j < (N + 2); j++) {
			if (i == 0 || i == NN || j == 0 || j == NN)
				 { s_2d[i][j] = 0; }
			else 
				 { s_2d[i][j] = 1; }
		}
	}
	return 0;
}
int iter_2d(char** s_2d, double T, int accept_more , double r, bool periodic , bool rand_parity_order) {

	double delta_E[9];
	for (int l = 0; l < 9; l++) { delta_E[l] = (J * ((double)l - 4)); }
	// arr[l] = J*(2l-4) i.e {l,arr[l]} = {0,-4J},{1,-3J},...... {8,4J}
	double prob[2][9] = { 0, };
	double dum;

	for (int c = 0; c <= 1; c++) {
	for (int l = 0; l <= 8; l++) {
		if		(accept_more == 2) { dum = 2.0 * (double)(2 * (double)c - 1) * delta_E[l] / (kB * T); 
								 if (dum < 0) { dum = 0.0; }
		}
		else if (accept_more == 1) { dum = ((double)(2 * (double)c - 1) * delta_E[l] + 4.0 * J) / (kB * T); }
		else if (accept_more == 0) { dum = ((double)(2 * (double)c - 1) * delta_E[l] + 8.0 * J) / (kB * T); }
		else { dum = 1.0; }

		prob[c][l] = r * exp(-1.0 * dum);
	}
	}
	// prob[i][j] = Accptance ratio for flipping s[a]==(2i-1) with
	// sum_<a,b> s[b] = 2j-6.

	int flipped = 0;
	int i;
	int order[2] = { 0, 0};
	int idx,parity;

	if (genrand64_real1() < 0.5 && rand_parity_order) { order[0] = 1; order[1] = 0; }
	else											  { order[0] = 0; order[1] = 1; }
	
	for (idx = 0; idx <= 1; idx++) {
		parity = order[idx];
	#pragma omp parallel shared(s_2d,parity) reduction(+:flipped)
	#pragma omp for	schedule(static)
		for (i = 1; i <= N; i++) {

			for (int j = 1; j <= N; j++) {
			if ((i + j) % 2 == parity) {
				int count = 0;

				if (periodic) {
					count += s_2d[(i == N) ? 1 : i + 1][j];
					count += s_2d[(i == 1) ? N : i - 1][j];
					count += s_2d[i][(j == N) ? 1 : j + 1];
					count += s_2d[i][(j == 1) ? N : j - 1];
				}
				else {
					count += s_2d[i + 1][j];
					count += s_2d[i - 1][j];
					count += s_2d[i][j + 1];
					count += s_2d[i][j - 1];
				}

				double p = prob[(s_2d[i][j] + 1) / 2][(count + 4)];
				
				if (genrand64_real1() < p) {
					s_2d[i][j] *= -1;
					flipped += 1;
				}

			}	
			}
		}
	}

	return flipped;
}

double get_avg_2d(char** s_2d) {

	double m = 0;
	int i = 0;
#pragma omp parallel shared(s_2d) reduction(+:m)
#pragma omp for	schedule(static)
	for (i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) { m += s_2d[i][j]; }
	}
	m /= N_sq;
	return m;
}
void get_parity_avg_2d(char** s_2d, double* m_p) {

	double m = 0.0;
	int i = 0;
	int parity;
	for (parity = 0; parity <= 1; parity++) {
		m = 0.0;
	#pragma omp parallel shared(s_2d,parity) reduction(+:m)
	#pragma omp for	schedule(static)
		for (i = 1; i <= N; i++) {
			for (int j = 1; j <= N; j++) {
					if ((i + j) % 2 == parity) { m += s_2d[i][j]; }
			}
		}
		m /= N_sq;
		m *= 2.0;
		m_p[parity] = m;
	}
}
// ---------------------------------------------------------
int newwrite_d(FILE* fp, double d, bool frontspace) {
	if (fp != NULL) { 
		if (frontspace) { fprintf(fp, " %lf", d); return 0; }
		else { fprintf(fp, "%lf", d); return 0; }
	}
	else { return 1; }
}
int newwrite_i(FILE* fp, int i) {
	if (fp != NULL) { fprintf(fp, "%d ", i); return 0; }
	else { return 1; }
}
int newline(FILE* fp) {
	if (fp != NULL) { fprintf(fp, "\n"); return 0; }
	else { return 1; }
}
// ----------------------------------------------------------

//   CurieT: Averaging for small random magnetic field.
int main(){
	
	unsigned long long init[4] = { 0x12345ULL, 0x23456ULL, 0x34567ULL, 0x45678ULL }, length = 4;
	init_by_array64(init, length);
	int num_fp = 8;
	
	double T;
	int i, l, n, t;
	int repeat_per_T = 100;
	int iter_per_datum = 50000;
	double inv, r_i = 1.0, r_f = 0.05;

	int option = 3;
	// option 3: 3d
	// option 2: 2d

	if (option == 3) {	
		num_fp = 1;

		FILE** fp_fp_metro = (FILE **)malloc(sizeof(FILE*) * num_fp);
		char filename[20];
		double* m_p = (double*)malloc(sizeof(double) * num_fp);

		if (fp_fp_metro) {
			for (i = 0; i < num_fp; i++) {
				sprintf_s(filename, 20, "3d_metro_%d.txt", i);
				if (fopen_s(&(fp_fp_metro[i]), filename, "w") != 0) { return 1; }
			}
		}

		char**** s_3d_ptr = (char****)malloc(sizeof(char***) * num_fp);
		if (s_3d_ptr) {
			for (i = 0; i < num_fp; i++) {
				s_3d_ptr[i] = alloc_3d();
			}
		}

		repeat_per_T = 100;
		iter_per_datum = 50000;
		inv = 1.0 / (double)iter_per_datum;
		r_i = 1.0; r_f = 0.01;

		double Tmpt[36] = { 49.91,49.92,49.93,49.94,49.95,49.96,49.97,49.98,49.99,50.01,50.02,50.03,50.04,50.05,50.06,50.07,50.08,50.09,
						   50.11,50.12,50.13,50.14,50.15,50.16,50.17,50.18,50.19,50.21,50.22,50.23,50.24,50.25,50.26,50.27,50.28,50.29 };

		if (m_p && fp_fp_metro && s_3d_ptr) {

			for (t = 15; t <= 35; t++) {
		//		T = 49.50 + 0.1* (double)t;
				T = Tmpt[t];
				printf("		      %.5f", T);
				for (l = 1; l <= repeat_per_T; l++) {
					
					for (i = 0; i < num_fp; i++) {
						newwrite_d(fp_fp_metro[i], T, false);
						initialize_3d(s_3d_ptr[i]);
					}

					for (n = 1; n <= iter_per_datum; n++) {

						iter_3d(s_3d_ptr[0], T, 2, r_i - (r_i - r_f) * pow(((double)n * inv), 0.5), false, true);
				//		iter_3d(s_3d_ptr[1], T, 2, r_i - (r_i - r_f) * pow(((double)n * inv), 0.5), false, true);

						for (i = 0; i < num_fp; i++) {
							m_p[i] += get_avg_3d(s_3d_ptr[i]);
						}

						if (n % 1000 == 0) {
							printf("\r%05d | %5.3f", n,m_p[0]/1000.0);
							for (i = 0; i < num_fp; i++) {
								m_p[i] /= 1000.0;
								newwrite_d(fp_fp_metro[i], m_p[i], true);
								m_p[i] = 0.0;
							}
						}

					}
					printf("\r		%03d", l);
					for (i = 0; i < num_fp; i++) {
						newline(fp_fp_metro[i]);
					}
				}
				printf("\n");
			}
			for (i = 0; i < num_fp; i++) {
				if (fp_fp_metro[i] != NULL) { fclose(fp_fp_metro[i]); }
			}

			for (i = 0; i < num_fp; i++) {
				free_3d(s_3d_ptr[i]);
			}
			free(s_3d_ptr);
			free(fp_fp_metro);
			free(m_p);
		}
	}

	//----------------------------------------------------------------------------------------------------------------------

	if (option == 2) {
		num_fp = 8;

		FILE** fp_fp_metro = (FILE**)malloc(sizeof(FILE*) * num_fp);
		char filename[20];
		double* m_p = (double *)malloc(sizeof(double) * num_fp);

		if (fp_fp_metro) {
			for (i = 0; i < num_fp; i++) {
				sprintf_s(filename, 20, "2d_metro_%d.txt", i);
				if (fopen_s(&(fp_fp_metro[i]), filename, "w") != 0) { return 1; }
			}
		}

		char*** s_2d_ptr = (char ***)malloc(sizeof(char**) * num_fp);
		if (s_2d_ptr) {
			for (i = 0; i < num_fp; i++) {
				s_2d_ptr[i] = alloc_2d();
			}
		}

	//	double T[15] = {26.25,26.30,26.31,26.315,26.32,26.325,26.33,26.331,26.332,26.333,26.3331,26.3332,26.3333,26.3334,26.3335};

		iter_per_datum = 50000;
		inv = 1.0 / (double)iter_per_datum;
		r_i = 1.0; r_f = 0.05;

		if (m_p && fp_fp_metro && s_2d_ptr) {

			for (t = 0; t <= 0; t++) {
				T = 26.3335 + 0.01 * (double)t;
				for (l = 1; l <= 50; l++) {

					printf("				%.5f", T);
					for (i = 0; i < num_fp; i++) {
						newwrite_d(fp_fp_metro[i], T,false);
						initialize_2d(s_2d_ptr[i]);
					}

					for (n = 1; n <= iter_per_datum; n++) {

						iter_2d(s_2d_ptr[0], T, 2, r_i - (r_i - r_f) * pow(((double)n * inv), 0.5), false, true);
						iter_2d(s_2d_ptr[1], T, 2, 0.3, false, true);
						iter_2d(s_2d_ptr[2], T, 2, 0.5, false, true);
						iter_2d(s_2d_ptr[3], T, 2, 0.7, false, true);
						iter_2d(s_2d_ptr[4], T, 2, 1.0, false, true);
						iter_2d(s_2d_ptr[5], T, 2, 1.0, false, false);
						iter_2d(s_2d_ptr[6], T, 2, 1.0, true, false);
						iter_2d(s_2d_ptr[7], T, 2, 1.0, true, true);

						for (i = 0; i < num_fp; i++) {
							m_p[i] += get_avg_2d(s_2d_ptr[i]);
						}

						if (n % 1000 == 0) {
							printf("\r%d", n);
							for (i = 0; i < num_fp; i++) {
								m_p[i] /= 1000.0;
								newwrite_d(fp_fp_metro[i], m_p[i],true);
								m_p[i]  = 0.0;
							}
						}

					}
					printf("\r			%d", l);
					for (i = 0; i < num_fp; i++) {
						newline(fp_fp_metro[i]);
					}
				}
				printf("\n");
			}
			for (i = 0; i < num_fp; i++) {
				if (fp_fp_metro[i] != NULL) { fclose(fp_fp_metro[i]); }
			}

			for (i = 0; i < num_fp; i++) {
				free_2d(s_2d_ptr[i]);
			}
			free(s_2d_ptr);
			free(fp_fp_metro);
			free(m_p);
		}
	}

	return 0;
}
