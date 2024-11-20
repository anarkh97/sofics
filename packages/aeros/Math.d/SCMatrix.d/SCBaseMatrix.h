#ifndef SCBASEMATRIX_H_
#define SCBASEMATRIX_H_

#include <mpi.h>
#include <string>

#ifdef SCARRAYS_DEV
#include "scpblas.h"
#include "scblacs.h"
#else
#include <Math.d/SCMatrix.d/scpblas.h>
#include <Math.d/SCMatrix.d/scblacs.h>
#endif


class SCBaseMatrix {
    public:
        SCBaseMatrix(int context, int m, int n, int mb, int nb, MPI_Comm comm);
        SCBaseMatrix(std::string filename, int context, int mb, int nb, MPI_Comm comm);
        ~SCBaseMatrix();

        int getNumberOfRows()       {return _m;}
        int getNumberOfRowsLocal()  {return _mlocal;}
        int getNumberOfCols()       {return _n;}
        int getNumberOfColsLocal()  {return _nlocal;}
        int getRowBlockingFactor()  {return _mb;}
        int getColBlockingFactor()  {return _nb;}
        int getContext()            {return _context;}
        int getFirstProcessRow()    {return _rsrc;}
        int getFirstProcessCol()    {return _csrc;} // this was set to rsrc????
        int getNumberOfProcsRow()   {return _mprow;}
        int getNumberOfProcsCol()   {return _npcol;}
        int getSizelocal()          {return _sizelocal;}
        int getLeadingLocalDim()    {return _lld;}
        int * getDesc()             {return _desc;}
        void printDescription();
        void printSummary();
        void setScope(char scope);
        void setScope();
        void setSrc(int rsrc=0, int csrc=0) {_rsrc=rsrc; _csrc=csrc;}
        void setTopology(char top=' ') {_top=top;}
        int getMatrixSize(std::string filename);
        int getLworkPxgeqrf(int ia=1, int ja=1, int rsrc=0, int csrc=0);
        int getLworkPxormqr(SCBaseMatrix& C, int ia=1, int ja=1, int rsrca=0,
                            int csrca=0, int ic=1, int jc=1, int rsrcc=0, int csrcc=0);
        int getLworkPdlapiv(char rowcol);
        char getPivroc();
        double getWallTime();
        int getLocalOffset(int ig, int jg);
        int getGlobalRowIndex(int iloc) {return _FORTRAN(indxl2g)(&iloc, &_mb, &_myrow, &ZERO, &_mprow);}
        int getGlobalColIndex(int jloc) {return _FORTRAN(indxl2g)(&jloc, &_nb, &_mycol, &ZERO, &_npcol);}
        int getLocalRowIndex(int iglo)  {return _FORTRAN(indxg2l)(&iglo, &_mb, &ZERO,   &ZERO, &_mprow);}
        int getLocalColIndex(int jglo)  {return _FORTRAN(indxg2l)(&jglo, &_nb, &ZERO,   &ZERO, &_npcol);}
        int getRowProc(int iglo)        {return _FORTRAN(indxg2p)(&iglo, &_mb, &ZERO,   &ZERO, &_mprow);}
        int getColProc(int jglo)        {return _FORTRAN(indxg2p)(&jglo, &_mb, &ZERO,   &ZERO, &_npcol);}
        int getProc(int ig, int jg);
        void setRowColComms();
        MPI_Comm getMpiComm() {return _comm;}
        int distributeVector(int    *vec, int n); // Convenience routine for distributing local arrays.
        int distributeVector(double *vec, int n); // Convenience routine for distributing local arrays.

        static int ZERO;


    protected:
        int _m;           // Global number of rows for this Scalapack matrix
        int _n;           // Global number of cols
        int _mlocal;      // Local number of rows
        int _nlocal;      // Local number of columns
        int _mb;          // Row blocking factor
        int _nb;          // Column blocking factor
        int _rsrc;        // First row processor. Defaults to 0.
        int _csrc;        // First column processor. Defaults to 0.
        int _lld;         // Local leading dimension of A
        int _mprow;       // Number of row processors in the grid
        int _npcol;       // Nuber of column processors in the grid
        int _myrow;       // Global index to first local row
        int _mycol;       // Global index to first local column
        int _context;     // PBlas Context
        int _nprocs;      // Number of processors
        int _mypid;       // My MPI process ID
        int _sizelocal;   // Size of local matrix
        char _scope;      // Scope of matrix. Changes depending on desired commuication pattern
        char _top;        // Communication topology.
        int _desc[DLEN_]; // Scalapack descriptor for matrix _A (_A will be defined in derived class)
        bool _row_col_comm_set;
        MPI_Comm _comm;
        MPI_Comm _row_comm;
        MPI_Comm _col_comm;

        void init();
};

#endif // SCBASEMATRIX_H_
