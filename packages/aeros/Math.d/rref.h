#ifndef __TOREDUCEDROWECHELONFORM__
#define __TOREDUCEDROWECHELONFORM__

/**
 * Row echelon form with full pivoting. Note that columns and
 * rows might be swapped due to the full pivoting, returned in rowmap and colmap if not NULL
 *  
 * @param M            the matr to reduce
 * @param reduced       if true, the reduced row echelon form is returned, that
 *                      is, zeros above the diagonal, too
 * @param rowmap        the row mapping to reestablish the original row ordering. 
 * @param colmap        the column mapping to reestablish the original column ordering.
 * @return              the rank of the matr
 */

//this may be a good choice, especially hooked with a sparse qr factorization. 
#define QR_ROW_ECHELON
#ifdef QR_ROW_ECHELON
#include <Eigen/Dense>
#endif

#include <limits>
#include <Utils.d/MyComplex.h>

template <typename T, typename Matrix>
int rowEchelon(Matrix& M, bool reduced, int* rowmap, int* colmap, int optc = 0, double tol = 10.0, bool usePrescribedThreshold = false) 
{
  // M is an augmented matrix [A;b]
#ifdef QR_ROW_ECHELON
  int p = M.rows();
  int n = M.cols()-1;
  Eigen::Block< Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > A = M.block(0,0,p,n);
  Eigen::ColPivHouseholderQR<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > decA(A); // A*P=Q*R
  const Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> &P = decA.colsPermutation();
  const typename Eigen::ColPivHouseholderQR< Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> >::HouseholderSequenceType &Q = decA.householderQ();
  //const typename Eigen::FullPivHouseholderQR< Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> >::MatrixQReturnType &Q = decA.matrixQ();
  //NOTE: if you want to use the FullPivHouseholderQR then make sure rowmap is set properly

  A = A*P;
  M = Q.transpose()*M; // now M is the row echelon form of the augmented matrix: [R;Q^T*b]
  for(int i=0; i<n; ++i) colmap[i] = P.indices()[i];
  if(usePrescribedThreshold)
    decA.setThreshold(tol*std::numeric_limits<T>::epsilon());
  int r = decA.rank();

  if(reduced) {
    Eigen::Block< Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > R11 = M.block(0,0,r,r); 
    // R11 is upper triangular and non-singular
    R11.template triangularView<Eigen::Upper>().solveInPlace(M.block(0,r,r,n-r+1));
    R11 = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>::Identity(r,r);
    // now M is the reduced row echelon form of the augmented matrix: [R11^{-1}*R12;R11^{-1}*Q^T*b]
  }
  return r;
#else
  int rows = M.rows();
  int cols = M.cols();
  int pivs = std::min(rows, cols-optc); // PJSA: don't pivot last column if optc == 1
        
  //find pivot row/column
  int prow = -1;
  int pcol = -1;
  T pval = -std::numeric_limits<T>::max();
  for (int row = 0; row < rows; row++) {
    for (int col = 0; col < cols-optc; col++) { // PJSA: don't pivot last column if optc == 1
      T val = std::abs(M(row, col));
      if (val > pval) {
        pval = val;
        prow = row;
        pcol = col;
      }
    }
  }
        
  //precondition (each iteration): prow/pcol/pval are set
  for (int pivot = 0; pivot < pivs; pivot++) {
    if (pval <= tol*std::numeric_limits<T>::epsilon()) return pivot;
          
    //swap rows / columns
    if (prow != pivot) {
      M.row(prow).swap(M.row(pivot));
      if (rowmap != NULL) {
        int tmp = rowmap[prow]; 
        rowmap[prow] = rowmap[pivot]; 
        rowmap[pivot] = tmp;
      }
    }

    if (pcol != pivot) {
      M.col(pcol).swap(M.col(pivot));
      if (colmap != NULL) {
        int tmp = colmap[pcol];
        colmap[pcol] = colmap[pivot];
        colmap[pivot] = tmp;
      }
    }
    
    //divide pivot row
    pval = M(pivot, pivot);
    for (int col = pivot + 1; col < cols; col++) {
      M(pivot,col) /= pval;
    }
    M(pivot,pivot) = 1.0;

    //subtract pivot row from other rows
    //find next pivot at the same time
    pval = -std::numeric_limits<T>::max();
    for (int row = pivot + 1; row < rows; row++) {
      T rpiv = M(row, pivot);                          
      M(row, pivot) = 0.0;
      for (int col = pivot + 1; col < cols; col++) {
        T val = M(row, col);
        T sub = M(pivot, col);
        val -= sub * rpiv;
        M(row, col) = val;
        // is this our new pivot?
        if (row < rows && col < cols-optc) { // PJSA: don't pivot last column if optc == 1
          val = std::abs(val);
          if (val > pval) {
            pval = val;
            prow = row;
            pcol = col;
          }
        }
      }
    }
    if (reduced) {
      //subtract pivot from rows above pivot row, too
      for (int row = 0; row < pivot; row++) {
        T rpiv = M(row, pivot);                          
        M(row, pivot) = 0.0;
        for (int col = pivot + 1; col < cols; col++) {
          T val = M(row, col);
          T sub = M(pivot, col);
          M(row, col) = val - sub * rpiv;
        }
      }
    }
  }
  return pivs;
#endif
}

template <typename T, typename Matrix>
void ToReducedRowEchelonForm(Matrix& m, int *rowmap) 
{
  int lead = 0;
  int rowCount = m.rows();
  int colCount = m.cols();

  int i;
  T lv;
 
  for(int r = 0; r < rowCount; r++) {
    if(lead >= colCount)
      return;
    i = r;
    //while(m(i,lead) == 0) {
    while(std::abs<T>(m(i,lead)) <= std::numeric_limits<T>::epsilon()) {
      i++;
      if(i == rowCount) {
        i = r;
        lead++;
        if(lead == colCount)
          return;
      }
    }

    m.row(i).swap(m.row(r));
    if (rowmap != NULL) {
      int tmp = rowmap[i];
      rowmap[i] = rowmap[r];
      rowmap[r] = tmp;
    }

    lv = m(r,lead);
    m.row(r) /= lv;
    for(i = 0; i < rowCount; i++) {
      if(i != r) {
        lv = m(i, lead);
        m.row(i) -= m.row(r)*lv;
      }
    }
    lead++;
  }
}

#endif
