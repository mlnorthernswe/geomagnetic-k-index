	SUBROUTINE SPLN2(E,NM)
	INTEGER*2 E,N,N1,NM
	REAL*4 V,A,Q1,QN,C,DVS
	DIMENSION V(124),A(124)
C
	CHARACTER*1 DEBUG
	INTEGER*2 IKM,IXY,IGS,ISP,ICSP,IGK
	COMMON/FOX/DEBUG,IKM,IXY,IGS,ISP,ICSP,IGK
C
      INTEGER*2 NS
	REAL*4 XX,YY,U,X,DEL,S
	COMMON/INDIA/XX(124),YY(124),U(124,2),X(124,2),S(124,2),
     *DEL(124,2),NS(2)
C
C  UNEQUALLY SPACED DATA
C  THIS PROGRAM COMPUTES THE SECOND DERIVATIVE, S(X), OF THE CUBIC
C  SPLINE, SPLIN(X), WHICH INTERPOLATES THE NPTS OF ARBITRARILY
C  SPACED DATA(X,U).
C
      IF(DEBUG.EQ.'Y'.AND.ISP.EQ.1)CALL WRITEDAT(NS(E),U(1,E),X(1,E))
	DVS=6.0                     !Interpolate between X(I,E) if dead
					    ! values span more than DVS hours
	CALL NTRP0(DVS,E,NM)
      IF(DEBUG.EQ.'Y'.AND.ISP.EQ.2)CALL WRITEDAT(NS(E),U(1,E),X(1,E))
      N=NS(E)
      Q1=S(1,E)
	QN=S(N-1,E)
	IF(DEBUG.EQ.'Y'.AND.ISP.EQ.4)WRITE(8,*)'  Q1,QN:',Q1,QN
	IF(N.GT.(NM+4))STOP '  SPLN2 0+5'
      IF(N.LT.5)STOP  '  SPLN2 0+6'
C
C  COMPUTE DEL AND V
C
      DEL(2,E)=X(2,E)-X(1,E)
      IF(DEL(2,E).EQ.0.0)STOP  '  SPLN2 0+7'
      IF(DEL(2,E).EQ.Q1)Q1=0.0
      V(1)=6.0*(((U(2,E)-U(1,E))/DEL(2,E))-Q1)
      N1=N-1
      DO 30 I=2,N1
         DEL(I+1,E)=X(I+1,E)-X(I,E)
         C1=U(I-1,E)/DEL(I,E)
         C2=U(I,E)*((1.0/DEL(I,E))+(1.0/DEL(I+1,E)))
         C3=U(I+1,E)/DEL(I+1,E)
         V(I)=(C1-C2+C3)*6.0
   30 CONTINUE
      V(N)=(QN+(U(N1,E)-U(N,E))/DEL(N,E))*6.0
C
C  GAUSSIAN ELIMINATION AND AUGMENTATION
C
      A(1)=2.0*DEL(2,E)
      A(2)=1.5*DEL(2,E)+2.0*DEL(3,E)
      V(2)=V(2)-0.5*V(1)
      DO 40 I=3,N1
         C=DEL(I,E)/A(I-1)
         A(I)=2.0*(DEL(I,E)+DEL(I+1,E))-C*DEL(I,E)
         V(I)=V(I)-C*V(I-1)
   40 CONTINUE
      C=DEL(N,E)/A(N1)
      A(N)=2.0*DEL(N,E)-C*DEL(N,E)
      V(N)=V(N)-C*V(N1)
C
C  BACK SUBSTITUTION
      S(N,E)=V(N)/A(N)
      DO 50 J=1,N1
         I=N-J
         S(I,E)=(V(I)-DEL(I+1,E)*S(I+1,E))/A(I)
   50 CONTINUE
      IF(DEBUG.EQ.'Y')THEN
	IF(ISP.EQ.3)CALL WRITEDAT(NS(E),DEL(1,E),S(1,E))
      ENDIF
   60 RETURN
      END
