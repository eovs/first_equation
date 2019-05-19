#include <iostream>
#include "data_structures.h"
#include "equations.h"
#include "check_girth.h"

#ifndef SKIP_MEX
#include <mex.h>
#define SKIP_COUT
#endif

#define SKIP_COUT

using namespace std;


int gf2( int m, int p, short *gf_log, short *gf_alog )
{
	// v is field array of gf(2^m) over prim polynom p
	int i;
	int m2 = 1 << m;
	int element;
/*	
	if( nargin == 1 )
	{
		P=[3,7,11,19,37,67,137, 285, 529, 1033];
		p=P(m);
	}
*/

	element = 1;
	gf_log[0] = -1; 
	gf_alog[0] = 1;

	for( i = 1; i < m2; i++ )
	{
		gf_alog[i-1] = element;
		gf_log[element] = i-1;

		element <<= 1;
		if( element >= m2 )
			element ^= p;

		if( element == 1 && i != (m2-1))
		{
			printf("p is not primitive\n");
			return -1;
		}
	}

	return 0;
}

static int prim_polinomial[10][100] =
{
	{03, 0},
	{07, 0},
	{015, 0},
	{023,037, 0},  
	{045, 075, 067, 0},   
	{0103,0147,0155, 0}, 
	{0203, 0211, 0217, 0235, 0367, 0277, 0325,  0313, 0345, 0},   
	{0435, 0551, 0453, 0545, 0543, 0537, 0703, 0747, 0}, 
	{
		01021, 01131, 01461, 01423, 01055, 01167, 01541, 01333, 
			01605, 01751, 01743, 01617, 01553, 01157, 01715, 01563, 
			01713, 01175, 01725, 01225, 01275, 01773, 01425, 01267, 0
	},
	{
		02011, 02415, 03771, 02157, 03515, 02773, 02033, 02443, 
			02461, 03023, 03543, 02745, 02431, 03177, 03525, 02617, 
			03471, 03323, 03507, 03623, 02707, 02327, 03265, 02055, 
			03575, 03171, 02047, 03025, 03337, 03211, 0
		}
};	

#if 01
GIRTH_ATTRIBUTE check_girth( matrix<int> &HB_org, matrix<int> &HC_org, int M, int q_mod, int girth )
{
	GIRTH_ATTRIBUTE gattr;
	short *gf2log;
	short *gf2alog;
	int p;
	int q_bits;

	int starting_length = HB_org.n_cols();
	int nrows = HB_org.n_rows();
	int ncols = HB_org.n_cols();

	matrix< int >HM(nrows, ncols);

	matrix< int >HB(nrows, ncols);
	for( int i = 0; i < HB.n_rows(); i++ )
		for( int j = 0; j < HB.n_cols(); j++ )
			HB(i, j) = HB_org(i, j);

	matrix< int >HC(nrows, ncols);
	for( int i = 0; i < HB.n_rows(); i++ )
		for( int j = 0; j < HB.n_cols(); j++ )
			HC(i, j) = HC_org(i, j);

	gattr.girth     = girth;
	gattr.equations = 0;
	gattr.badHB     = 0;
	gattr.badHC     = 0;


	for( int i = 0; i < HB.n_rows(); i++ )
	{
		for( int j = 0; j < HB.n_cols(); j++ )
		{
			if( HB(i, j) > -1  )
				HB(i, j) = HB(i, j) % M;
		}
	}

	if( q_mod > 2 )
	{
		for( int i = 0; i < HB.n_rows(); i++ )
		{
			for( int j = 0; j < HB.n_cols(); j++ )
			{
				if( HB(i, j) > -1  )
				{
					if( HC(i, j) == 0 )
						HC(i, j) = q_mod-1;
					else
						HC(i, j) = HC(i, j) % q_mod;
				}
				else
					HC(i, j) = -1;

			}
		}

		gf2log  = (short*)calloc( q_mod, sizeof(short) );
		gf2alog = (short*)calloc( q_mod, sizeof(short) );

		switch( q_mod )
		{
		case    2: q_bits =  1; break;
		case    4: q_bits =  2; break;
		case    8: q_bits =  3; break;
		case   16: q_bits =  4; break;
		case   32: q_bits =  5; break;
		case   64: q_bits =  6; break;
		case  128: q_bits =  7; break;
		case  256: q_bits =  8; break;
		case  512: q_bits =  9; break;
		case 1024: q_bits = 10; break;
		default:   q_bits =  0;
		}

		p = prim_polinomial[q_bits-1][0];
		gf2( q_bits, p, gf2log, gf2alog );

		// transform natural representation to power representation 
		for( int i = 0; i < HB.n_rows(); i++ )
		{
			for( int j = 0; j < HB.n_cols(); j++ )
			{
				if( HC(i, j) > -1  )
					HC(i, j) = gf2log[HC(i, j)];
			}
		}

		free( gf2log );
		free( gf2alog );
	}


	for( int i = 0; i < nrows; i++ )
	{
		for( int j = 0; j < ncols; j++ )
		{
			switch( HB(i, j) )
			{ 
			case -1: HM(i, j) = 0; break;
			case  0: HM(i, j) = 2; break;
			default: HM(i, j) = 1;
			}
		}
	}	

	// 1.1. Base matrix
	graph< vector_bag > first_equations_base;
	int first_equations_base_edges = 0;
	for( int i = 0; i < starting_length; i++ ) 
	{
		for( int j = 0; j < nrows; j++ ) 
		{
			if( HM(j, i) > 0 ) 
			{
				vector_bag fw(first_equations_base_edges++);
				first_equations_base.add_bidi_edge(i, j + starting_length, fw, -fw);
			}
		}
	}

#ifndef SKIP_COUT
	cout << " Tanner graph size: nt=" << first_equations_base.n_edges() / 2
		<< ", rt=" << first_equations_base.n_vertices() << endl;
	cout << "Generating first equations" << endl;
#endif
	// 1.2. Equations
	equation_builder first_equations(first_equations_base, girth);
	if (first_equations.tree_girth() < girth) 
	{
#ifndef SKIP_COUT
		cout << "problem: required girth " << girth
			<< ", found only " << first_equations.tree_girth() << endl;
#endif
		return gattr;
	}
#ifndef SKIP_COUT
	cout << "ok: g=" << first_equations.tree_girth()
		<< ", eqs=" << first_equations.n_equations()
		<< ", nodes=" << first_equations.n_vertices() << endl;
#endif

	voltage_check_result check;
	vector< int > a(first_equations_base_edges);
	vector< int > coef(first_equations_base_edges);    // used for non-binary codes only 

	for( int k = 0, i = 0; i < starting_length; i++ ) 
	{
		for( int j = 0; j < nrows; j++ ) 
		{
			//			if (current_HM(j, i) > 0) 
			if( HB(j, i) > -1 ) 
			{
				a[k]    = HB(j, i);
				coef[k] = HC(j, i);
				k++;
			}
		}
	}

	if( q_mod == 2 )
		check = first_equations.check_voltages(a, M);
	else
		check = first_equations.check_voltages_and_coefs( a, M, coef, q_mod-1 );

	gattr.equations = check.equations;
	gattr.badHB     = check.badHB;
	gattr.badHC     = check.badHC;
	
	return gattr;
}

#endif

#ifndef SKIP_MEX
void unpackMatrix(double m[], int height, int width, double *result[] )
{
	int i,j;
	for( i = 0; i < height; ++i)
	{
		for( j = 0; j < width; ++j)
		{
			result[i][j] = m[j * height + i];
		}
	}
}

void unpackRow(double m[], int height, int width, double result[] )
{
	int j;
	for( j = 0; j < width; ++j)
	{
		result[j] = m[j];
	}
}


void unpackMatrix_double2int(double m[], int height, int width, int *result[] )
{
	int i,j;
	for( i = 0; i < height; ++i)
	{
		for( j = 0; j < width; ++j)
		{
			result[i][j] = (int)m[j * height + i];
		}
	}
}
void unpackMatrix_double2short(double m[], int height, int width, short *result[] )
{
	int i,j;
	for( i = 0; i < height; ++i)
	{
		for( j = 0; j < width; ++j)
		{
			result[i][j] = (short)m[j * height + i];
		}
	}
}

void unpackMatrix(double m[], int height, int width, matrix<int> &result )
{
	int i,j;
	for( i = 0; i < height; ++i)
	{
		for( j = 0; j < width; ++j)
		{
			result(i, j) = (short)m[j * height + i];
		}
	}
}
#endif  //SKIP_MEX

#ifndef SKIP_MEX
// %
// % decoder interface:
// % 1. SETUP ( before simulation )
// % ========
// % result = decoder( decoder_type, c, HB, HC, M );
// %    decodec_type - defined decoder:
// %        0 - bp
// %        1 - sum prod
// %		2 - advanced sum prod
// %		3 - min-sum
// %		4 - int min-sum
// %		5 - int advanced sum prod
// %		6 - fht sum prod
// %        7 - int fht sum prod
// %    c  - size of symbol (in bits)
// %    HB - circulant matrix
// %    HC - symbol matrix
// %    M  - matrix extension

//
// %  result == 1, if OK
// %
// % 2. DECODE
// % =========
// %  [it, soft, decword] = decoder( chan_out, maxiter );
// %	 chan_out -  channel output
// %	 maxiter  -  max number of iterations
// %	 
// %     it       -  number of made iterations 
// %     soft     -  decoded codeword ( soft decision )
// %     decword  -  decoded codeword ( hard decision )
// %
// %
// % 3. CLEAR ( free memory allocations )
// % ========
// %  decoder();     
// % 


void mexFunction(int nOut, mxArray *pOut[], int nInp, const mxArray *pInp[])
{
	static int nh, mh, M, q_mod;
	static int girth; 
	static matrix<int> HB;
	static matrix<int> HC;

	mh = (int)mxGetM(pInp[0]);
	nh = (int)mxGetN(pInp[0]);

	HB = matrix<int>(mh, nh);
	HC = matrix<int>(mh, nh);

	M     = (int)mxGetPr(pInp[2])[0];
	q_mod = (int)mxGetPr(pInp[3])[0];
	girth = (int)mxGetPr(pInp[4])[0];

	unpackMatrix(mxGetPr(pInp[0]), mh, nh, HB );
    if( q_mod > 2 )
        unpackMatrix(mxGetPr(pInp[1]), mh, nh, HC );

    GIRTH_ATTRIBUTE girth_attribute = check_girth( HB, HC, M, q_mod, girth );
	
    pOut[0] =  mxCreateDoubleScalar((double)girth_attribute.equations);
	pOut[1] =  mxCreateDoubleScalar((double)girth_attribute.badHB);
	pOut[2] =  mxCreateDoubleScalar((double)girth_attribute.badHC);
}
#endif
