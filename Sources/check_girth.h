#ifndef _CHECK_GIRTH_H_
#define _CHECK_GIRTH_H_

typedef struct  
{
	int girth;
	int equations;
	int badHB;
	int badHC;

}GIRTH_ATTRIBUTE;

GIRTH_ATTRIBUTE check_girth( matrix<int> &HB, matrix<int> &HC, int M, int q_mod, int girth );

#endif //_CHECK_GIRTH_H_