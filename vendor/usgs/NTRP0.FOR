	SUBROUTINE NTRP0(VN,E,NM)
	INTEGER*2 E,I,NR,K,NM
	REAL*4 VN,A
C
      INTEGER*2 NS
	REAL*4 XX,YY,Y,X,B1,DEL
      COMMON/INDIA/XX(124),YY(124),Y(124,2),X(124,2),B1(124,2),
     *DEL(124,2),NS(2)
C
      I=0
	IF(NS(E).EQ.(NM+4))RETURN
      NR=NS(E)-1
    2 IF(I.GT.NR)GOTO 15
      I=I+1
    5 B=X(I+1,E)-X(I,E)
      IF(B.LT.VN)GOTO 2
      A=Y(I+1,E)-Y(I,E)
      K=NR+2
   10 X(K,E)=X(K-1,E)
      Y(K,E)=Y(K-1,E)
      B1(K,E)=B1(K-1,E)
      K=K-1
      IF(K.GT.(I+1))GOTO 10
      X(I+1,E)=X(I,E)+B*0.5
      Y(I+1,E)=Y(I,E)+A*0.5
      NR=NR+1
	IF(NR.LT.(NM+3))GOTO 5
   15 CONTINUE
      NS(E)=NR+1
      IF(NS(E).LT.5)STOP '  NTRP0  15+2'
      IF(NS(E).LT.5)STOP 'NTRP0 15+2'
      RETURN
      END
