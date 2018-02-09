#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>

#define BUF_SIZE 4096

#define USE_DOUBLES

#ifdef USE_FLOATS
#define REAL float
#define FORMAT_82 "%4d  %f,%f  %f,%f  %f,%f"
#define FORMAT_SIGMA "%4d  %f,%f  %f,%f  %f,%f  %f,%f  %f,%f  %f,%f  %f,%f  %f,%f  %f,%f"
#endif

#ifdef USE_DOUBLES
#define REAL double
#define FORMAT_82 "%4d  %lf,%lf  %lf,%lf  %lf,%lf"
#define FORMAT_SIGMA "%4d  %lf,%lf  %lf,%lf  %lf,%lf  %lf,%lf  %lf,%lf  %lf,%lf  %lf,%lf  %lf,%lf  %lf,%lf"
#endif

#define EPS 0.001

static char buffer[BUF_SIZE];


int jump_to(int len, const char* text, char buffer[len], FILE* file)
{
	while (fgets(buffer, len, file) != NULL)
	{
		if (strstr(buffer, text))
		{
			return 0;
		}
	}
	return 1;
}

void get_dimension_data(int* n,
					   int* nbe,
					   int* ni,
					   FILE* file
					   )
{
	#define MALHA_STR    " N\332MERO DE ELEMENTOS DA MALHA="
	#define CONTORNO_STR " N\332MERO DE ELEMENTOS DE CONTORNO="
	#define N_INTERNO_STR  " N\332MERO DE PONTOS INTERNOS="

	rewind(file);

	while (fgets(buffer, BUF_SIZE, file) != NULL)
	{
		if (strstr(buffer, MALHA_STR))
		{
			sscanf(buffer, MALHA_STR "%5d\n", n);
		}
		else if (strstr(buffer, CONTORNO_STR))
		{
			sscanf(buffer, CONTORNO_STR "%5d", nbe);
		}
		else if (strstr(buffer, N_INTERNO_STR))
		{
			sscanf(buffer, N_INTERNO_STR "%5d", ni);
		}

	}
	rewind(file);
}

void replace_D_with_E(int bufsize, char buffer[bufsize])
{
	unsigned int i = 0;
	while (i < bufsize && buffer[i] != '\0')
	{	if (buffer[i] == 'd' || buffer[i] == 'D')
			buffer[i] = 'E';
		++i;
	}
}

void get_nos_contorno(int nbe, REAL complex nos_contorno[nbe][3], FILE* file)
{
	REAL reals[3][2];
	int index, i, j;

	rewind(file);

	#define NOS_CONTORNO_STR " N\323S DO CONTORNO"

	jump_to(BUF_SIZE, NOS_CONTORNO_STR, buffer, file);

	/*Descarta as tres próximas linhas*/
	for (i = 0; i < 3; ++i)
	{
		fgets(buffer, BUF_SIZE, file);
	}

	/*Lê de fato os resultados.*/
	for (i = 0; i < nbe; ++i)
	{
		if (fgets(buffer, BUF_SIZE, file) == NULL)
			return;
		replace_D_with_E(BUF_SIZE, buffer);
		sscanf(buffer, FORMAT_82, &index, &reals[0][0], &reals[0][1], &reals[1][0], &reals[1][1], &reals[2][0], &reals[2][1]);
		
		for (j = 0; j < 3; ++j)
			nos_contorno[index-1][j] = reals[j][0] + I*reals[j][1];
	}
}

void get_tractions(int nbe, REAL complex tractions[nbe][3], FILE* file)
{
	REAL reals[3][2];
	int index, i, j;

	rewind(file);

	#define TRACTION_STR "TRACTION X"
	
	jump_to(BUF_SIZE, TRACTION_STR, buffer, file);

	/*Descarta a próxima linha*/
	fgets(buffer, BUF_SIZE, file);

	/*Lê de fato os resultados.*/
	for (i = 0; i < nbe; ++i)
	{
		if (fgets(buffer, BUF_SIZE, file) == NULL)
			return;
		replace_D_with_E(BUF_SIZE, buffer);
		sscanf(buffer, FORMAT_82, &index, &reals[0][0], &reals[0][1], &reals[1][0], &reals[1][1], &reals[2][0], &reals[2][1]);
		
		for (j = 0; j < 3; ++j)
			tractions[index-1][j] = reals[j][0] + I*reals[j][1];
	}
	rewind(file);
}

void get_deslocamentos_internos(int ni, REAL complex deslocamentos[ni][3], FILE* file)
{
	REAL reals[3][2];
	int index, i, j;

	rewind(file);

	#define INTERNO_STR "  PONTOS INTERNOS"
	
	jump_to(BUF_SIZE, INTERNO_STR, buffer, file);

	for (i = 0; i < 2; ++i)
		fgets(buffer, BUF_SIZE, file);


	/*Lê de fato os resultados.*/
	for (i = 0; i < ni; ++i)
	{
		if (fgets(buffer, BUF_SIZE, file) == NULL)
			return;
		replace_D_with_E(BUF_SIZE, buffer);
		sscanf(buffer, FORMAT_82, &index, &reals[0][0], &reals[0][1], &reals[1][0], &reals[1][1], &reals[2][0], &reals[2][1]);
		
		for (j = 0; j < 3; ++j)
			deslocamentos[index-1][j] = reals[j][0] + I*reals[j][1];
	}
	rewind(file);
}

void get_sigmas_internos(int ni, REAL complex sigmas[ni][3][3], FILE* file)
{
	REAL reals[3][3][2];
	int index, i, j, k;

	rewind(file);

	#define SIGMA_STR "SIGMAXX"

	jump_to(BUF_SIZE, SIGMA_STR, buffer, file);

	/*Lê de fato os resultados.*/
	for (i = 0; i < ni; ++i)
	{
		if (fgets(buffer, BUF_SIZE, file) == NULL)
			return;
		replace_D_with_E(BUF_SIZE, buffer);
		sscanf(buffer, FORMAT_SIGMA, &index, 
				&reals[0][0][0], &reals[0][0][1], 	
				&reals[0][1][0], &reals[0][1][1], 
				&reals[0][2][0], &reals[0][2][1], 
				&reals[1][0][0], &reals[1][0][1], 
				&reals[1][1][0], &reals[1][1][1], 
				&reals[1][2][0], &reals[1][2][1], 
				&reals[2][0][0], &reals[2][0][1], 
				&reals[2][1][0], &reals[2][1][1], 
				&reals[2][2][0], &reals[2][2][1] 
		);	
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
			sigmas[index-1][j][k] = reals[j][k][0] + I*reals[j][k][1];
	}
	rewind(file);
}

double norm_1_delta(int m, int n, REAL complex A[m][n], REAL complex B[m][n])
{
	int i, j;
	double acc = 0;
	double acc_max = 0;

	for (i = 0; i < m; ++i)
	{
		acc = 0;
		for (j = 0; j < n; ++j)
			acc += cabs(A[i][j] - B[i][j]);
	
		if (acc > acc_max)
			acc_max = acc;
	}
	return acc_max;
}

#define ASSERT(acc) if (acc > EPS) return false
bool compare_nos_contorno(int nbe, REAL complex nos_contorno[nbe][3], REAL complex nos_contorno_sol[nbe][3])
{
	double acc = norm_1_delta(nbe, 3, nos_contorno, nos_contorno_sol);

	printf("%E\n",acc);
	ASSERT(acc);

	return true;
}

bool compare_tractions(int nbe, REAL complex tractions[nbe][3], REAL complex tractions_sol[nbe][3])
{
	double acc = norm_1_delta(nbe, 3, tractions, tractions_sol);
	
	printf("%E\n",acc);
	ASSERT(acc);
	return true;
}

bool compare_deslocamentos_internos(int ni, REAL complex deslocamentos[ni][3], REAL complex deslocamentos_sol[ni][3])
{
	double acc = norm_1_delta(ni, 3, deslocamentos, deslocamentos_sol);
	
	printf("%E\n",acc);
	ASSERT(acc);
	return true;
}

bool compare_sigmas_internos(int ni, REAL complex sigma[ni][3][3], REAL complex sigma_sol[ni][3][3])
{
	double acc = 0;
	double acc_max = 0;
	int i;


	for (i = 0; i < ni; ++i)
		acc += norm_1_delta(3, 3, sigma[i], sigma_sol[i]);
	if (acc > acc_max)
		acc_max = acc;

	printf("%E\n",acc);
	ASSERT(acc);
	
	return true;
}

bool has_nans(int m, int n, REAL complex A[m][n])
{
	int i, j;

	for (i = 0; i < m; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			if (isnan(creal(A[i][j])) || isnan(cimag(A[i][j])))
			{
				return true;
			}
		}
	}
	return false;
}
int main(int argc, char* argv[])
{
	FILE* file, *file_sol;
	int n, nbe, ni;

	if (argc != 3)
	{
		fputs("Error: Wrong number of arguments\n", stderr);
		return 500;
	}

	file = fopen(argv[1], "r");
	file_sol = fopen(argv[2], "r");
	
	if (!file || !file_sol)	
	{
		fputs("Error: Invalid filepath\n", stderr);
		return 501;
	}
	get_dimension_data(&n, &nbe, &ni, file);

	REAL complex nos_contorno[nbe][3], tractions[nbe][3], 
		         deslocamentos[ni][3], sigmas[ni][3][3];

	REAL complex nos_contorno_sol[nbe][3], tractions_sol[nbe][3], 
		         deslocamentos_sol[ni][3], sigmas_sol[ni][3][3];

	bool validity = true;

	get_nos_contorno(nbe, nos_contorno, file);
	get_tractions(nbe, tractions, file);
	get_deslocamentos_internos(ni, deslocamentos, file);
	get_sigmas_internos(ni, sigmas, file);
	
	validity &= !has_nans(nbe, 3, nos_contorno);
	validity &= !has_nans(nbe, 3, tractions);
	validity &= !has_nans(ni,  3, deslocamentos);
	validity &= !has_nans(ni,  9, sigmas);
	
	if (validity == 0)
		printf("NaN in calculations = %d\n", (unsigned int) validity);

	get_nos_contorno(nbe, nos_contorno_sol, file_sol);
	get_tractions(nbe, tractions_sol, file_sol);
	get_deslocamentos_internos(ni, deslocamentos_sol, file_sol);
	get_sigmas_internos(ni, sigmas_sol, file_sol);

	validity &= compare_nos_contorno(nbe, nos_contorno, nos_contorno_sol);
	validity &= compare_tractions(nbe, tractions, tractions_sol);
	validity &= compare_deslocamentos_internos(ni, deslocamentos, deslocamentos_sol);
	validity &= compare_sigmas_internos(ni, sigmas, sigmas_sol);


	if (validity)
		return 0;
	else return 1;
}
