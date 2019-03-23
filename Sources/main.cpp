#include <iostream>

#include <stdio.h>

#include "data_structures.h"
#include "equations.h"
#include "find_girth.h"

#define START_GIRTH 4
#define Q_MOD 16
#define TAILBITE_LENGTH 32

using namespace std;


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


int main( void )
{
	matrix<int> current_HB;
	matrix<int> current_HC;
	int q_mod;
	int rows, columns;
	int M;
	int starting_length;
	GIRTH_ATTRIBUTE girth_attribute;
	FILE *fp;

	M               = TAILBITE_LENGTH;
	q_mod           = Q_MOD;


#if 1
	fp = fopen("HB.txt", "rt");
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

	fp = fopen("HC.txt", "rt");
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
	for( int i = 0; i < current_HC.n_rows(); i++ )
	{
		for( int j = 0; j < current_HC.n_cols(); j++ )
		{
			if( current_HB(i, j) > -1  )
				current_HC(i, j) = current_HC(i, j) % M;
		}
	}


	if( q_mod > 2 )
	{
		for( int i = 0; i < current_HC.n_rows(); i++ )
		{
			for( int j = 0; j < current_HC.n_cols(); j++ )
			{
				if( current_HB(i, j) > -1 && current_HC(i, j) > -1 )
					current_HC(i, j) = current_HC(i, j) % (q_mod-1);
				else
					current_HC(i, j) = -1;
			}
		}
	}

	girth_attribute = find_girth( current_HB, current_HC, M, q_mod, START_GIRTH );

/*	cout << "girth is " << girth_attribute.girth << endl << endl;
	cout << "girth << : " << "equations " << girth_attribute.equations
		 << ", bad HB eqs " << girth_attribute.badHB
		 << ", bad HC eqs " << girth_attribute.badHC
		 << endl;*/

	system("pause");
	return 0;
}