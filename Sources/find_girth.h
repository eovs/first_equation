#ifndef _FIND_GIRTH_H_
#define _FIND_GIRTH_H_

typedef struct  
{
	int girth;
	int equations;
	int badHB;
	int badHC;

}GIRTH_ATTRIBUTE;

GIRTH_ATTRIBUTE find_girth( matrix<int> &HB, matrix<int> &HC, int M, int q_mod, int start_girth );

#endif //_FIND_GIRTH_H_