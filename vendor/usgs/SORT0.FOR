      SUBROUTINE SORT0(IX,NPTS,IXB,IXL,N)
      INTEGER*2 DV,IX,NPTS,JUMP,LDONE,N,M,ITEMP,IXL,IXB
      DIMENSION IX(1)
      DV=32767
      JUMP=NPTS
    5 JUMP=JUMP/2
   10 LDONE=1
      DO 20 M=1,(NPTS-JUMP)
        N=M+JUMP
        IF(IX(M).LE.IX(N))GOTO 20
        ITEMP=IX(M)
        IX(M)=IX(N)
        IX(N)=ITEMP
        LDONE=0
   20 CONTINUE
      IF(LDONE.EQ.0)GOTO 10
      IF(JUMP.GT.1)GOTO 5
   30 IXL=IX(1)
      N=NPTS
   40 IF(N.LE.0.OR.IX(N).NE.DV)GOTO 50
      N=N-1
      GOTO 40
   50 IXB=IX(N)
      RETURN
      END
