	SUBROUTINE XYSET(E,NM,NPPM,NPTS)
	INTEGER*2 H,E,DV,I,J,NM,NPPM
C
	INTEGER*2 MHV,IDM
	REAL*4 XD
	COMMON/BAKER/MHV(360,2),IDM(3,2),XD
C
	CHARACTER*1 DEBUG
	INTEGER*2 IKM,IXY,IGS,ISP,ICSP,IGK
	COMMON/FOX/DEBUG,IKM,IXY,IGS,ISP,ICSP,IGK
C
	INTEGER*2 NS
	REAL*4 X,Y,YHAT,Z,B1,DEL
	COMMON/INDIA/X(124),Y(124),YHAT(124,2),Z(124,2),B1(124,2),
     *DEL(124,2),NS(2)
C
      DV=32767
	J=1
	A=FLOAT(NPPM)/60.
	B=0.5*A
	IF(MHV(NM-1,E).NE.DV.AND.MHV(NM,E).NE.DV)GOTO 502
	H=NM-2                      !Find value within 1 thru nm-2
	DO I=1,NM-2
	  IF(MHV(H,E).NE.DV)GOTO 30
	  H=H-1
	ENDDO
	H=NM                        !No value, seach nm thru 2*nm
	DO I=1,NM
	  H=H+1
	  IF(MHV(H,E).NE.DV)GOTO 40
	ENDDO
	X(J)=-11.5                  !No value, use daily mean
      Y(J)=IDM(1,E)
      GOTO 50
   30 X(J)=A*FLOAT(H)-(24.0+B)
      Y(J)=MHV(H,E)
      GOTO 50
   40 X(J)=-11.5                  !Found value between nm and 2*nm
	Y(J)=MHV(H,E)               ! but none between 1 and nm
   50 CONTINUE
      J=J+1
  502 IF(MHV(NM-1,E).EQ.DV)GOTO 505
	Y(J)=MHV(NM-1,E)
	X(J)=-3*B
      J=J+1
  505 IF(MHV(NM,E).EQ.DV)GOTO 520
	Y(J)=MHV(NM,E)
	X(J)=-B
      J=J+1
  520 DO 530 H=NM+1,2*NM
	  IF(MHV(H,E).EQ.DV)GOTO 530
	  Y(J)=MHV(H,E)
	  X(J)=A*FLOAT(H)-(24.0+B)
	  J=J+1
  530 CONTINUE
	K=2*NM+1
	IF(MHV(K,E).EQ.DV)GOTO 535
	Y(J)=MHV(K,E)
	X(J)=24.0+B
      J=J+1
  535 IF(MHV(K+1,E).EQ.DV)GOTO 540
	Y(J)=MHV(K+1,E)
	X(J)=24.0+3*B
      J=J+1
  540 IF(MHV(K,E).NE.DV.AND.MHV(K+1,E).NE.DV)GOTO 549
	H=2*NM+3
	DO 110 I=1,NM-2
      IF(MHV(H,E).NE.DV)GOTO 130
      H=H+1
  110 CONTINUE
	H=2*NM
	DO I=1,NM
	  IF(MHV(H,E).NE.DV)GOTO 140
	  H=H-1
	ENDDO
      X(J)=35.5
      Y(J)=IDM(3,E)
      GOTO 150
  130 X(J)=A*FLOAT(H)-(24.0+B)
      Y(J)=MHV(H,E)
      GOTO 150
  140 X(J)=24.0+5*B
	Y(J)=MHV(H,E)
  150 CONTINUE
      GOTO 550
  549 J=J-1
  550 IF(J.GT.2)GOTO 560
C***  GS REQUIRES MIN OF 3 PTS
      J=3
      X(1)=-11.5
      X(2)=12.0
      X(3)=35.5
	DO I=1,3
	  Y(I)=IDM(3,E)
	ENDDO
  560 NPTS=J
	IF(DEBUG.EQ.'Y'.AND.IXY.EQ.1)CALL WRITEDAT(NPTS,Y(1),X(1))
      RETURN
      END
C-----------------------------------------------------------------
	SUBROUTINE WRITEDAT(NPTS,Y,X)
	INTEGER*2 NPTS
	REAL*4 X,Y
	DIMENSION Y(1),X(1)
	I=6
	J=1
	DO WHILE(J.LE.NPTS)
	  IF(I.GT.5)THEN
	    WRITE(8,1000)Y(J),X(J)
 1000       FORMAT(1X,F7.1,F6.2)
	    I=1
	  ELSE
	    WRITE(8,1010)Y(J),X(J)
 1010       FORMAT('&',F7.1,F6.2)
	    I=I+1
	  ENDIF
	  J=J+1
	END DO
	RETURN
	END
