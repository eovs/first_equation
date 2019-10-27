#include <iostream>

#include <stdio.h>

#include "data_structures.h"
#include "equations.h"
#include "check_girth.h"


#define MATRIX_HB "HB.TXT"
#define MATRIX_HC "HC.TXT"

#define MIN_GIRTH 6
#define MAX_GIRTH 32
#define Q_MOD  16
#define TAILBITE_LENGTH 7//32

//#define CYCL 4

using namespace std;

#if 01
#define NROWS 3
#define NCOLS 4
int HB[NROWS][NCOLS] = 
{
	0,  7, -1,  2,
	10,  16,  5, -1,
	8, -1,  17,  3
};

int HC[NROWS][NCOLS] = 
{
	1,  1, -1,  1,
	1,  1,  1, -1,
	1, -1,  1,  1
};
#else
#define NROWS 4
#define NCOLS 8
int HB[NROWS][NCOLS] = 
{
	0, -1, -1, 31, 10, -1,  0, 16,   
	0,  0, -1, -1, 23, 12,  1, 14,
	-1,  0,  0,  0,  5, 27, -1,  3,
	-1, -1,  0, 31, -1, 30,  6,  5  
};

int HC[NROWS][NCOLS] = 
{
	6, -1, -1, 10, 14, -1,  7,  8,   
	6, 14, -1, -1, 11, 10,  4,  4,
	-1, 14,  9,  8, 12,  6, -1,  6,
	-1, -1,  9, 10, -1,  6, 11,  2  
};
#endif

int main( int argc, char **argv )
{
	matrix<int> current_HB;
	matrix<int> current_HC;
	int q_mod;
	int rows, columns;
	int M;
	int starting_length;
	FILE *fp;

	M               = TAILBITE_LENGTH;
	q_mod           = Q_MOD;

	//q_mod = 2;   // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

#if 01
	if( argc != 3 && argc != 5 )
	{
		printf("correct run for binary codes:\n");
		printf("  first_equation HB M\n");
		printf("correct run for non-binary codes:\n");
		printf("  first_equation HB HC M Q\n");
		printf("where HB, HC - matrices\n");
		printf("      M - lifting\n");
		printf("      Q - 2^m\n");
		printf("\n");
	
		system("pause");
		return 1;
	}

	M     = argc == 3 ? atoi(argv[2]) : atoi(argv[3]); 
	if( argc == 3 )
		q_mod = 2;
	else
		q_mod = atoi(argv[4]);

	if( q_mod > 2 )
	{
		printf("parameters:\n");
		printf("HB - %s\n", argv[1]);
		printf("HC - %s\n", argv[2]);
		printf("M  - %d\n", M);
		printf("Q  - %d\n", q_mod);
		printf("\n");
	}
	else
	{
		printf("parameters:\n");
		printf("HB - %s\n", argv[1]);
		printf("M  - %d\n", M);
		printf("\n");
	}

	fp = fopen(argv[1], "rt");
	if( fp == NULL )
		return 1;

	fscanf_s( fp, "%d", &rows );
	fscanf_s( fp, "%d", &columns );

	current_HB = matrix<int>(rows, columns);
	current_HC = matrix<int>(rows, columns);

	for( int i = 0; i < rows; i++ )
	{
		for( int j = 0; j < columns; j++ )
		{
			int tmp;
			fscanf_s( fp, "%d", &tmp );
			current_HB(i, j) = tmp;
		}
	}
	fclose( fp );

	if( q_mod > 2 )
	{
		fp = fopen(argv[2], "rt");
		if( fp == NULL )
			return 1;

		fscanf_s( fp, "%d", &rows );
		fscanf_s( fp, "%d", &columns );

		for( int i = 0; i < rows; i++ )
		{
			for( int j = 0; j < columns; j++ )
			{
				int tmp;
				fscanf_s( fp, "%d", &tmp );
				current_HC(i, j) = tmp;
			}
		}
		fclose( fp );
	}

	starting_length = current_HB.n_cols();

#else
	starting_length = NCOLS;

	current_HB = matrix<int>(NROWS, NCOLS);
	current_HC = matrix<int>(NROWS, NCOLS);
	rows    = current_HB.n_rows();
	columns = current_HB.n_cols();

	for( int i = 0; i < rows; i++ )
		for( int j = 0; j < columns; j++ )
			current_HB(i, j) = HB[i][j];

	for( int i = 0; i < rows; i++ )
		for( int j = 0; j < columns; j++ )
			current_HC(i, j) = HC[i][j];
#endif

#ifdef CYCL
	for( int i = 0; i < rows; i++ )
		for( int j = 0; j < columns; j++ )
			current_HB(i, j) = current_HB(i, j) >= 0 ? current_HB(i,j) + CYCL : -1;
#endif

	for( int girth = MIN_GIRTH; girth < MAX_GIRTH; girth += 2 )
	{
		int flag; 
		GIRTH_ATTRIBUTE girth_attribute = check_girth( current_HB, current_HC, M, q_mod, girth );

		if( girth_attribute.equations < 0 )
			flag = 1; 
		else
			flag = q_mod > 2 ? flag = girth_attribute.badHC : girth_attribute.badHB;

		if( flag == 0 )
			printf("girth >= %2d, ", girth); 
		else
			printf("girth  = %2d, ", girth-2); 

		if( girth_attribute.equations < 0 )
		{
			printf("-- BALANCED_CYCLE\n");
		}
		else
		{
			printf("number of eqs: %d", girth_attribute.equations);
			if( q_mod == 2 )
				printf("   bad HB eqs %3d\n", girth_attribute.badHB);
			else
				printf("   bad HB eqs  %3d, bad HC eqs  %3d\n", girth_attribute.badHB, girth_attribute.badHC);
		}

		if( flag )
			break;
	}
	system("pause");
	return 0;
}