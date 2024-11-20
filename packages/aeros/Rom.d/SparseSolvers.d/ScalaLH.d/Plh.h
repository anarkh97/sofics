#ifndef PLH_H_
#define PLH_H_

#include <vector>
#include <cstring>
#include <Eigen/Core>
#include <ctime>

#ifdef NNLS_DEV
#include "SCDoubleMatrix.h"
#else
#include "Math.d/SCMatrix.d/SCDoubleMatrix.h"
#endif

#define TIME_MAIN_LOOP 0
#define TIME_NEXT_VECTOR 1
#define TIME_GRADF 2
#define TIME_UPDATEQR 3
#define TIME_PDORMQR 4
#define TIME_HOUSE 5
#define TIME_UPDATEQTB 6
#define TIME_SOLVER 7
#define TIME_PDTRSV 8
#define TIME_UPDATEX 9
#define TIME_MOVEFROMPTOZ 10
#define TIME_PERMUTE_A 11
#define TIME_GETMAX 12
#define TIME_INIT_NNLS 13
#define TIME_INITLS_ZERO 14
#define TIME_INITLS_SET 15
#define TIME_INITLS_SIZE 16
#define TIME_INITLS_COPY 17
#define TIME_COPYREDIST 18
#define TIME_INVERT_MASK 19
#define TIME_COPYXQRTOX 20
#define TIME_WRITE_RESIDUAL 21
#define TIME_RESIDUAL2NORM 22
#define TIME_ITER 23
#define TIME_LDCHECK 24
#define TIME_PZCHECK 25
#define TIME_MCOPYQTOA 26
#define TIME_PDORMQR_GRADF 27
#define TIME_MULT_GRADF 28
#define TIME_LOADMATRIX 29
#define TIME_LOADRHS 30
#define TIME_REJECT_VECTOR 31
#define TIME_GET_SOLUTION 32
#define TIME_DOWNDATE 33
#define TIME_COLUMNSCALING 34
#define TIME_COLUMNNORMS 35
#define N_TIMES 36

#define MBA_DEFAULT 16
#define NBA_DEFAULT 16
#define MBQ_DEFAULT 16
#define NBQ_DEFAULT 16

#define MAXITE_DEFAULT 10.0

#define HEADER_INCR 30


class Plh {

    public:
        Plh();
        Plh(int m, int n);
        Plh(const std::vector< Eigen::Map<Eigen::MatrixXd> >& A);
        ~Plh();

        int setMatrix(SCDoubleMatrix &A);
        int setMatrixColumn(int j, double *col);
        int setMatrixRow(int j, double *row);
        int setRHS(double *rhs);
        int solve();
        int close();

        int writeMatrix(std::string filename, bool compact=false);
        int writeRhs(std::string filename, bool compact=false);
        int writeX(std::string filename, bool compact=false);
        double getResidualNorm() {return _rnorm2;};
        double getComputeTime() {return getTime(TIME_MAIN_LOOP);}
        double getDistributeMatrixTime() {return getTime(TIME_LOADMATRIX);}
        double getDownDateTime() {return getTime(TIME_DOWNDATE);}
        void initMatrix(double *A);
        void initRhs(double *b);
        void setResidualIncrement(int incr) {_residualIncr = incr;};
        void loadMatrix(const std::vector< Eigen::Map<Eigen::MatrixXd> >& em);
        void loadRhs(const Eigen::Ref<const Eigen::VectorXd> &eb);
        void summary();

        void setABlockSize(int mb, int nb) {_mb=mb; _nb=nb;}
        void setAProcGrid(int mprow, int npcol) {_mprow=mprow; _npcol=npcol;}
        void setQBlockSize(int mbq, int nbq) {_mbq=mbq; _nbq=nbq;}
        void setQProcGrid(int mprowQR, int npcolQR) {_mprowQR=mprowQR; _npcolQR=npcolQR;}
        void setRtol(double rtol) {_rtol=rtol;}
        void setMaxNP(int maxNP) {_maxNP=maxNP;}
        void setVerbose(int verbose) {_verbose=verbose;}
        void setContext(int context, int mprow, int npcol, MPI_Comm comm);

        void init();
        void hotStart() { _sizeHS = hotIndices.size(); _hotInd = 0; 
                          if(_sizeHS > 0) _hotStart = true; };
        void getQtoA(SCIntMatrix &QtoAHot) { _QtoA->copy(QtoAHot); };
        void init(const std::vector< Eigen::Map<Eigen::MatrixXd> >& eigenMatrix,
                  const Eigen::Ref<const Eigen::VectorXd>& eigenRhs);
        void setMatrixSize(int m, int n);
        void setMatrixSize(const std::vector< Eigen::Map<Eigen::MatrixXd> >& eigenMatrix, const Eigen::Ref<const Eigen::VectorXd>& eigenRhs);
        void singleContext();
        void writeSolution(bool compact=false);
        void printTimes(bool debug=false);
        void setResidualFileName(std::string fname, int incr=1);
        Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,1> getSolution();
        void getSolution(double *x);
        void write(std::string filename, Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,1> & x);
        void testCommunicators();
        void setDownDateMask();
        void setMaxIterRatio( double maxite ) {_max_iter=maxite*_n;}
        void setMaxIter( double max_iter ) {_max_iter=max_iter;}
        void setColumnScaling();
        void setPolytopeFacesPursuit();
        void setOrthogonalMatchingPursuit();
        void setConstraint(int zeroedCol) {_constraint = zeroedCol; }
        int getStatus() {return _status;}
        int getContext() {return _context;}
        int getIter() {return _iter;}
        SCDoubleMatrix& getRhsVector() {return *_b;}
        SCDoubleMatrix& getSolutionVector() {return *_x;}

    private:
        // A context
        SCDoubleMatrix * _A;
        SCDoubleMatrix * _x;
        SCDoubleMatrix * _b;
        SCDoubleMatrix * _w;
        SCDoubleMatrix * _wmask;
        SCDoubleMatrix * _colnorms;
        SCDoubleMatrix * _Atb;

        // QR context
        SCDoubleMatrix * _Q;
        SCDoubleMatrix * _Qtb;
        SCDoubleMatrix * _rhs;
        SCDoubleMatrix * _workm;
        SCDoubleMatrix * _bQR;
        SCDoubleMatrix * _zQR;
        SCDoubleMatrix * _xQR;
        SCDoubleMatrix * _wQR;
        SCDoubleMatrix * _rQR;
        SCDoubleMatrix * _workmQR;

        // data structures for Polytope Face Pursuit
        SCDoubleMatrix * _vertex;
        SCDoubleMatrix * _oneVecQR; 
        SCDoubleMatrix * _trslvQR; 
        SCDoubleMatrix * _trslv;
        SCDoubleMatrix * _Atv; 

        double * _work_qr;

        SCIntMatrix * _QtoA;
        SCIntMatrix * _QtoAHot; // pointer to hotStart data
        std::vector<int> hotIndices; 

        int _hotInd; 
        int _sizeHS; 

        // Comm groups. Needed when MPI routines are more handy than BLACS routines.
        MPI_Comm _comm;
        MPI_Comm _row_comm;
        MPI_Comm _col_comm;
        MPI_Comm _row_commQR;
        MPI_Comm _col_commQR;

        bool _contextInitialized;
        bool _matrixInitialized;
        bool _initializedWithEigen;
        bool _ddmask;
        bool _col_scaling;
        bool _PFP; 
        bool _OMP; 
        bool _hotStart;
        int _residualIncr;
        std::string _iterstring;
        std::string _residualFileName;
        FILE * _residualFilePtr;
        int _nEigenSubDomains;
        int * _eigenSubDomainSize;
        int * _eigenColsPerProc;
        int * _eigenSubdomainStart;
        int * _eigenSubDomainsPerProc;

        DoubleInt _wmax;
        DoubleInt _zmin;

        double _rnorm2;
        double _alpha;
        double _rtol;
        int _verbose;

        // Timings
        double _wallclock[N_TIMES];
        double _wallclock_total[N_TIMES];

        int _status;
        int _m;
        int _n;
        int _mb;
        int _nb;
        int _mbq;
        int _nbq;
        int _mprow;
        int _npcol;
        int _myrow;
        int _mycol;
        int _context;
        int _mprowQR;
        int _npcolQR;
        int _myrowQR;
        int _mycolQR;
        int _contextQR;
        int _nprocs;
        int _mypid;
        int _iter;
        int _iter_total;
        int _max_iter;
        int _subiter;
        int _lwork_qr;
        int _nP;
        int _nZ;
        int _ialpha;
        int _maxNP;
        int _constraint;

        void initDefaults();
        void initWithSize();
        void initWithEigen(const std::vector< Eigen::Map<Eigen::MatrixXd> >& A);
        int initplh();
        int gradf();
        double residual2Norm();
        void iteration_output();
        void sub_iteration_output(int iqr);
        void header();
        std::string intToString(int i);

        int updateX();
        int computeX();
        bool nextVector();
        bool rejectVector();
        int updateQ(int Acol, int Qcol);
        int getAlpha();
        void writeSet(std::string filename);
        void writeSet();
        void writeHotSet();
        void copyToHotInd();
        int distributeVector(int context, char scope, char top, double *vec, int n, int proc_coord=0);
        int distributeVector(int context, char scope, char top, int    *vec, int n, int proc_coord=0);
        int updateQROld(int iqr);
        bool updateQR(int iqr);
        bool updateQtb(int iqr=0);
        void updateVertex();
        void solveR();
        int moveFromPToZSwap();
        int moveFromPToZShift();
        int copyxQR2x();
        int mcopyQtoA(SCDoubleMatrix * xQR, SCDoubleMatrix * x);
        double getWallTime();
        void writeResidual();
        std::vector<int> factor(int n);
        void defaultProcGrid();
        void startTime(int i) {_wallclock[i] = -getWallTime();}
        void stopTime(int i) {_wallclock[i] += getWallTime(); _wallclock_total[i] += _wallclock[i];}
        double getTime(int i) { return _wallclock_total[i]; }
        int getEigenProc(int j);
        void wMaskByColNorm();
        void loadMaxTimes();
        void computeResidual();
};

#endif // PLH_H_
