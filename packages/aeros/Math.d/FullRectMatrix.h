#ifndef	_FULLRECTMATRIX_H_
#define	_FULLRECTMATRIX_H_

#include <cstdio>

class FullRectMatrix {

private:
	
	int	sizem;
	int	sizen;
	int     myval;

	double*	value;

public:
        FullRectMatrix(int, int, double *l=0);
 
	FullRectMatrix(); 
	
	~FullRectMatrix();
	
        FullRectMatrix operator *= (double v);

	double *operator[] (int row);
        int numRow() const { return sizem; }
        int numCol() const { return sizen; }
	void zero();

        double* data() { return value; }
};

inline
double *
FullRectMatrix::operator[](int row)
 { return value+sizen*row; }

#endif
