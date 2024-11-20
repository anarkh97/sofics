#ifndef INT_FULLM_H_
#define INT_FULLM_H_

// Integer Full Matrix class

class IntFullM {
 protected:
   int nrow, ncolumn;
   int *v;
 public:
   IntFullM(); // Creates an empty matrix
   IntFullM(int _nr);
   IntFullM(int _nr, int _nc);

   ~IntFullM();

   void zero();

   int numRow() const { return nrow;    }
   int numCol() const { return ncolumn; }

   int *operator[](int i);

   void print();
};

inline int *
IntFullM::operator[](int i)
{ return v+i*ncolumn; }

#endif
