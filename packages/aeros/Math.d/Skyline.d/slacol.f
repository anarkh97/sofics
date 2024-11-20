      SUBROUTINE SLACOL(LD,LACOL,NEQ,MAXLAC,AVELAC)
C 
      INTEGER LD(*),LACOL(*),NEQ,MAXLAC,AVELAC
C
C-----Local Declarations
C
      INTEGER I,J,ROW,MAX
C
C     LACOL stores the last active column at each step
C
      DO 500 I = 1, NEQ
      LACOL(I) = 0
500   CONTINUE
C
      LACOL(1)   = 1
      DO 600 J   = 2, NEQ
      ROW        = J + 1 - (LD(J) - LD(J - 1))
      LACOL(ROW) = J
600   CONTINUE
C
      MAX    = 0
      AVELAC = 0
      DO 700 I = 1, NEQ
      MAX      = MAX0(LACOL(I),MAX)
      LACOL(I) = MAX
      AVELAC   = AVELAC + LACOL(I) - I
700   CONTINUE
      AVELAC = AVELAC/NEQ
C
      MAXLAC = 0
      DO 800 I = 1, NEQ
      MAXLAC   = MAX0(MAXLAC,(LACOL(I) - I))
800   CONTINUE
C
C
      RETURN
      END
