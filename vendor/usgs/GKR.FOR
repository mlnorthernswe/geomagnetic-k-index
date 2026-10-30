	SUBROUTINE GKR(E)
C  Purpose:  Compute component amplitude, SR contribution, 3-hr
C            range, and estimated trucation of amplitude due to
C            minute average.
c      Use statistics compiled for variance of 1-min data
c	about the SR curve and the linear fit of 1-min data minus
c	the SR curve to determine whether the SR curve should be
c     adjusted.  This is consistent with the basic rules for
c     adopting the SR curve by Mayaud.
C
	INTEGER*2 E,IDV,MDV,KS,K,I,M,J,MAX,MIN,NGM,L,MPR
	REAL*4 DV,DVM,YMAX,YMIN,YM,R,C,DM,CSPLN2
C
	INTEGER*2 MHV,IDM
	REAL*4 XD
	COMMON/BAKER/MHV(360,2),IDM(3,2),XD
C
	INTEGER*2 IX,IDUM,YJD
	COMMON/CHARLY/IX(1440,2),IDUM,YJD
C
	CHARACTER*1 DEBUG
	INTEGER*2 IKM,IXY,IGS,ISP,ICSP,IGK
	COMMON/FOX/DEBUG,IKM,IXY,IGS,ISP,ICSP,IGK
C
	INTEGER*2 MP,ICC,I3R
	REAL*4 DEL,SESRC
	COMMON/GOLF/DEL(8,2),SESRC(8,2),MP(8,2),ICC(8,2),I3R(8,2)
C
	CHARACTER*1 SNH
	INTEGER*2 IDT,NBH,NEH,NTR
	REAL*4 SEP,SEPA
	COMMON/HOTEL/IDT,SEP(8,2),SEPA(8,2),SNH,NBH,NEH,NTR(2)
C
      DV=32767.
      IDV=32767
      DVM=-32767.
      MDV=-32767
      K=1
      I=0
      M=1
      DO 70 J=1,8
					!Do least square fit of curve to
	  KS=2                  ! minute data.  Statistics will be
	NI=0                    ! used for adjusting amplitude.
	SX=0.
	  SY=0.
	  SXY=0.
        SXX=0.
        SYY=0.
        B=0.
        KS=2
        YMAX=DVM
        YMIN=DV
        DO L=1,180
          K=(J-1)*180+L
          IF(IX(K,E).NE.IDV)THEN
            NI=NI+1
            SX=SX+FLOAT(L)
            SXX=SXX+FLOAT(L)*FLOAT(L)
            DM=FLOAT(K)/60.
            YHAT1=CSPLN2(KS,DM,E)
            D=0.1*(FLOAT(IX(K,E))-YHAT1)	!tenth to units
            IF(E.EQ.2)D=D*XD			!min to nT
            SY=SY+D				!Sum of residuals
            SYY=SYY+D*D				! and their squares
            SXY=SXY+D*FLOAT(L)
            IF(D.GT.YMAX)YMAX=D
            IF(D.LT.YMIN)YMIN=D
          END IF
        END DO
        IF(NI.GT.0)THEN
          B=(FLOAT(NI)*SXY-SX*SY)/(FLOAT(NI)*SXX-SX*SX)  !Slope of line
          D=(YMAX-YMIN)
        END IF
        IF(NI.GT.2)THEN
          V1=SYY-(SY*SY/FLOAT(NI))	!Variance of residuals (1-min 
					! minus SR curve) about their mean
          V2=B*(SXY-(SX*SY/FLOAT(NI)))	!Variance due to regression 
          VA=(V1-V2)/FLOAT(NI-2)	!Variance about regression
          IF(VA.LE.0.)VA=999.
          VR=SYY/FLOAT(NI-2)		!Variance of 1-min values
					! about regression
        END IF
      KS=2
      DEL(J,E)=DV
	MP(J,E)=IDV
	ICC(J,E)=IDV
	I3R(J,E)=IDV
      MAX=MDV
      MIN=IDV
      YMAX=DVM
      YMIN=DV
      NGM=0
      DO 20 L=1,180
	I=I+1                          !Counter for minute of the day
      IF(IX(I,E).EQ.IDV)GOTO 20
	NGM=NGM+1                      !Number of good minutes in 3-hrs
	DM=FLOATI(I)                   !TIME IN MINUTES OF DAY
      DM=DM/60.
C*** VARIATION IN TENTH-MIN OR TENTH-NT
	YHAT1=CSPLN2(KS,DM,E)          !Ordinate of component SR curve
	YM=IX(I,E)                     !Ordinate of component
	R=(YM-YHAT1)*0.1               !Deviation of component from SR in nT
	C=YHAT1*0.1
      IF(E.EQ.2)R=R*XD
      IF(E.EQ.2)C=C*XD
      YM=YM*0.1
      IF(E.EQ.2)YM=YM*XD
      M=ININT(YM)
      IF(M.GT.MAX)MAX=M
      IF(M.LT.MIN)MIN=M
C*** FIND MAX DIVIATION AND TIME
      IF(R.LE.YMAX)GOTO 13
      YMAX=R
      MAXT=I
      MAXP=ININT(C)
C*** FIND MIN DEVIATION AND TIME
   13 IF(R.GE.YMIN)GOTO 20
      YMIN=R
      MINT=I
      MINP=IFIX(C)
   20 CONTINUE
	IF((180-NGM).GE.60)GOTO 70    !If number of good minutes in
						!3-hr period is greater than 60
						!then compute the amplitude.
      DEL(J,E)=YMAX-YMIN		!AMPLITUDE
      C=DEL(J,E)
	ICC(J,E)=MAXP-MINP            !3-HR SR CONTRIBUTION
      ZZ=FLOAT(MAX)-FLOAT(MIN)
      IF(ZZ.LT.DV)I3R(J,E)=MAX-MIN	!3-HR RANGE
      IF(DEBUG.EQ.'Y')WRITE(8,1005)YJD,J,E,I3R(J,E),DEL(J,E)
				!Amplitude greater than 3-hour range?
      IF(DEL(J,E).GT.FLOAT(I3R(J,E)))THEN
        DEL(J,E)=FLOAT(I3R(J,E))
      ELSE
			!Compute confidence limits for amplitude
        JJ=MAXT/180				!Compute V of max
        NB=MAXT-(JJ*180)
        VXB=(FLOATJ(NB)-90.)*(FLOATJ(NB)-90.)
        VX=SXX -(SX*SX/FLOAT(NI))
        VB=(VR/FLOAT(NI))+VXB*VR/VX
        NL=MINT-(JJ*180)
        VXL=(FLOATJ(NL)-90.)*(FLOATJ(NL)-90.)
        VL=(VR/FLOAT(NI))+VR*VXL/VX
        V2=VB+VL				!Combined V of max/min
        CL1=0.
        IF(V2.GT.0.0)CL1=SQRT(V2)*SEP(J,E) !The standard error of amplitude
	  SESRC(J,E)=CL1/SEP(J,E)           ! times standard normal deviate.
					! The normal standard deviate is input from
					! file *.run so that values may be varied
					! to accommodate nightime intervals.
        V2=FLOATI(I3R(J,E))-DEL(J,E)	! If the difference between range
        IF(V2.GT.CL1.AND.DEL(J,E).GT.CL1)THEN	! and amplitude or the
					! amplitude is significantly
					! greater than zero use the amp-
					! litude otherwize use the range.
          IF(VA.NE.999.)THEN		!Compute adjustment to the SR
            D=B*FLOAT(NB)-B*FLOAT(NL)	! curve using linear
        				! least squares fit to 1-min
            VB=(VA/FLOAT(NI))+VA*VXB/VX	! data minus SR curve.  Compute
            VL=(VA/FLOAT(NI))+VA*VXL/VX	! ??% confidence limits for the
            V2=VB+VL			! adjustment.
            V1=9999.
            IF(V2.GT.0)CL2=SQRT(V2)*SEPA(J,E)
            V3=ABS(D)
	    IF(DEBUG.EQ.'Y')WRITE(8,1003)CL1,ICC(J,E),V3,CL2
		IF(V3.GE.CL2.AND.(ABS(D).LE.0.95*DEL(J,E)))THEN  !If the
              DEL(J,E)=DEL(J,E)-ABS(D)	! adjustment is significantly
					! greater than zero and the
	      IF(DEBUG.EQ.'Y')WRITE(8,1000)VA,DEL(J,E)
					! adjustment is 95% of the
					! amplitude apply the adjust-
					! ment to the amplitude.
            ELSE
	      IF(DEBUG.EQ.'Y')WRITE(8,1001)VA
		END IF
          ENDIF
        ELSE
          DEL(J,E)=FLOATI(I3R(J,E))
        END IF
      END IF
      IF(DEBUG.EQ.'Y'.AND.IGK.EQ.1)WRITE(8,1004)DEL(J,E)
	MP(J,E)=MPR(E,MAXT,MINT)            !ESTIMATE POSSIBLE TRUNCATION
					!DUE TO 1-MIN AVERAGE
      IF(DEL(J,E).GE.DV.OR.DEL(J,E).LT.0.0)DEL(J,E)=DV
   70 CONTINUE
      RETURN
 1000 FORMAT('&',F7.2,F5.1)
 1001 FORMAT('&',F7.2)
 1003 FORMAT('&',',L1=',F5.2,',SRC=',
     *I5,',A2=',F5.1,',L2=',F6.2)
 1004 FORMAT('&',',DEL='F5.1)
 1005 FORMAT(1X,I6,2I2,'R=',I4,',A=',F5.1)
	END
