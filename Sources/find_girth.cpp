#include <iostream>
#include "data_structures.h"
#include "equations.h"
#include "find_girth.h"

using namespace std;

GIRTH_ATTRIBUTE find_girth( matrix<int> &HB, matrix<int> &HC, int M, int q_mod, int start_girth )
{
	GIRTH_ATTRIBUTE gattr;

	int girth = start_girth;
	int starting_length = HB.n_cols();
	int nrows = HB.n_rows();
	int ncols = HB.n_cols();

	matrix< int >HM(nrows, ncols);

	gattr.girth     = girth;
	gattr.equations = 0;
	gattr.badHB     = 0;
	gattr.badHC     = 0;

	while( 1 )
	{
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

		cout << " Tanner graph size: nt=" << first_equations_base.n_edges() / 2
			<< ", rt=" << first_equations_base.n_vertices() << endl;
		cout << "Generating first equations" << endl;
		// 1.2. Equations
		equation_builder first_equations(first_equations_base, girth);
		if (first_equations.tree_girth() < girth) 
		{
			cout << "problem: required girth " << girth
				<< ", found only " << first_equations.tree_girth() << endl;
			return gattr;
		}
		cout << "ok: g=" << first_equations.tree_girth()
			<< ", eqs=" << first_equations.n_equations()
			<< ", nodes=" << first_equations.n_vertices() << endl;

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
			check = first_equations.check_voltages_and_coefs( a, M, coef );

		if( check.first_failed_equation == -1 )
		{
			cout << "----- girth " << girth << ": OK\n\n";
			gattr.girth = girth;
			girth += 2;
		}
		else
		{
			cout << "----- girth " << girth << ": BAD -- eqs: " << check.equations << ", badHB " << check.badHB << ", badHC " << check.badHC << "\n\n";
			gattr.equations = check.equations;
			gattr.badHB = check.badHB;
			gattr.badHC = check.badHC;
			break;
		}

		//break;
	}
	return gattr;
}
