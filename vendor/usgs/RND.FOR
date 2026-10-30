      INTEGER FUNCTION RND*2(B)
      REAL*4 A,B,C
      INTEGER*2 I,J,K
      IF(ABS(B).GT.32767.0)THEN
        RND=32767
        GOTO 10
      END IF
      A=AINT(B)
      I=NINT(A)
      K=IABS(I)
      C=ABS(B-A)*10.
      IF(C.EQ.5.0)THEN
        J=IMOD(K,2)
        IF(J.NE.0)K=K+1
      END IF
      IF(C.GT.5.0)K=K+1
      IF(I.LT.0)THEN
        RND=-K
      ELSE
        RND=K
      END IF
10    RETURN
      END

