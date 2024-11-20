// ------------------------------------------------------------
// HB - 04-15-05 
// ------------------------------------------------------------
#include <cstdio>
#include <cstdlib>

// HB (04-15-05): check for zero/small jacobian & constant sign of the jacobian
//                over the el.
// Inputs : J    : value of the jacobian
//          jSign: sign of the previous jacobian (+1 if > 0; otherwise -1)
//                 assume jSign = 0 at first call (i.e. no previous jacobian)
//          atol : absolute tolerance for checking zero/small jacobian.
//                 default value is 0.0
//          stop : boolean indicating if the routine stop the program execution 
//                 if the case a zero/small jacobian or non constant sign of the 
//                 jacobian over the el. is encountered.
//                 default value is true.
//          mssg : message to be added to the warning/error message.
//                 default value is 0.
//          file : where the warning/error message is displayed/written.
//                 default value is stderr.
//          
// Outputs: J    : absolute value of the jacobian
//          jSign: sign of the current jacobian
//          error: result of the check: 
//                    error = 0 -> no problem
//                          = 1 -> zero/small jacobian
//                          = 2 -> non constant jacobian's sign over the el.
//
int  
checkJacobian(double *J, int *jSign, int elId, const char* mssg= 0, double atol = 0.0, bool stop=true, FILE* file=stderr)
{ 
  // order of the check if optimized for positive jacobian
  int error = 0;
  if(*J>atol){ // check for constant jacobian sign over the el.
    if(*jSign<0){
      error = 2;
      if(stop){
        fprintf(file, " *** FATAL ERROR: Sign changing Jacobian determinant. %s, element %6d\n", mssg,elId);
        exit(-1);
      } else
        fprintf(file, " *** WARNING: Sign changing Jacobian determinant. %s, element %6d\n", mssg,elId);
    }
    *jSign= +1;
    return error;
  } else {
    if(*J<atol) {
      fprintf(file, " *** WARNING: Negative Jacobian determinant J = %e. %s, element %6d\n", *J, mssg,elId);
      if(*jSign>0){
        error = 2;
        if(stop){
          fprintf(file, " *** FATAL ERROR: Sign changing Jacobian determinant. %s, element %6d\n", mssg,elId);
          exit(-1);
        } else
          fprintf(file, " *** WARNING: Sign changing Jacobian determinant. %s, element %6d\n", mssg,elId);
      }
      *J    = -*J;
      *jSign= -1;
      return error;
    } else { // null/small jacobian
      error = 1;
      if(stop){
        fprintf(file, " *** FATAL ERROR: Zero Jacobian determinant. %s, element %6d\n", mssg,elId);
        exit(-1);
      } else
        fprintf(file, " *** WARNING: Zero Jacobian determinant. %s, element %6d\n", mssg,elId);
    }
  }
  return error;
}

/* previous order of the check
 if(fabs(*J)<atol){ // check for null/small jacobian
    error = 1;
    if(stop){
      fprintf(file, " *** FATAL ERROR: Zero Jacobian determinant. %s\n", mssg);
      exit(-1);
    } else
      fprintf(file, " *** WARNING: Zero Jacobian determinant. %s\n", mssg);
  }
  if(*J<0.0){ // check for constant jacobian sign over the el.
    if(*jSign>0){
      error = 2;
      if(stop){
        fprintf(file, " *** FATAL ERROR: Sign changing Jacobian determinant. %s\n", mssg);
        exit(-1);
      } else
        fprintf(file, " *** WARNING: Sign changing Jacobian determinant. %s\n", mssg);
    }
    *J    = -*J;
    *jSign= -1;
  } else {
    if(*jSign<0){
      error = 2;
      if(stop){
        fprintf(file, " *** FATAL ERROR: Sign changing Jacobian determinant. %s\n", mssg);
        exit(-1);
      } else
        fprintf(file, " *** WARNING: Sign changing Jacobian determinant. %s\n", mssg);
    }
    *jSign  = +1;
  }
  return error;
*/

void  
printShapeFct3D(double* Shape, double (*DShape)[3], int nnodes, char* mssg= 0, FILE* file=stderr)
{
  for(int i=0; i<nnodes; i++){
    fprintf(file," %s, shape[%2d] = %e, DShape[%2d] = %e  %e  %e\n",mssg,i,Shape[i],i,DShape[i][0],DShape[i][1],DShape[i][2]);
  } 
} 
