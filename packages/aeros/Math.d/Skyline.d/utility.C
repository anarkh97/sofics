// Auxiliary definitions
const int BCFIXED = 2;

void
makeUnconstrainedNum( int numdofs, const int *bc, int *unconstrainedNum, int &numUncon )
{
  numUncon  = 0;
  int i;
  for(i=0; i<numdofs; ++i) {
     if( bc != 0 && bc[i] == BCFIXED )
       unconstrainedNum[i] = -1;         // To mark constrained dofs.
     else
       unconstrainedNum[i] = numUncon++; // To renumber unconstrained dofs.
   }
}
