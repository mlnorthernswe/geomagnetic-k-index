	PROGRAM USGSDK
C  Purpose:  Compute 3-hour amplitude, SE of SR curve, possible error
c	due to 1-minute mean, SR curve contribution, and 3-hour range
c     then output parameters to disk.  The file is used by DKTLB.EXE
c	to generate a K-index table.
c
c  Input:  The path to a file of component 1-minute values in NGIC
c          format is obtained at runtime.
c       :  File *.KAY, contains runtime parameters, is obtained by
c          using program SELMVS.  The filename is input at runtime.
c       :  Input runtime parameters in file *.RUN
C       :  File OBSY.DAT contains observatory parameters that are
c          used by this program and should be in the current working
C          directory.
c
c  Output:  File *.DIG, used in K computation
C        :  File *.CRV, to be used with plot routine
c        :  File *.DKA, digital K-index table
c
c  Subroutines: KAYPRM,KMEAN,XYSET,GS,SPLN2,CSPLN2,GKR,KTAB,MPR,
C               MOVDAT,NTRP0,CURVK,KTABLE,RKV
C  D:\GENERAL\GENERAL.LIB: EXNAM,RND,OINFO,GYJD,UPPER,GETFLN
c
c  Author:  L. R. Wilson 10/29/92
c  Notes:  Original code was used on a Digital Equipment Corp.
c          MicroVax computer.  The code has been modified for
c          use on an IBM-compatible PC.  The program was compiled
c          and linked using Lahey 77 FORTRAN.
	CHARACTER*1 PANS,WCRV,IEL
	CHARACTER*2 IDENT
	CHARACTER*9 IDATE
	CHARACTER*12 NAM,NAM2,NAM3,NAM4,NAM5
	CHARACTER*31 NAM1
	INTEGER*2 X,E,D,DV,NPPM,NM,NCP,IB,IE,NPTS,IER,IY,IZS,
     *INDX,GYJD
	REAL*4 CF,FMIN9
	DIMENSION X(1440,2),CF(2),IY(1986),IZS(1986),INDX(9),IEL(2)
C
      INTEGER*2 IX
      COMMON/ABLE/IX(5762)
C
      INTEGER*2 MHV,IDM
      REAL*4 XD
      COMMON/BAKER/MHV(360,2),IDM(3,2),XD
C
	INTEGER*2 LX
	COMMON/CHARLY/LX(2882)
C
      CHARACTER*2 ID
	INTEGER*2 MHI,MLO,RNG,MENDS,IDA,MQD,TH,VH,RND
	REAL*4 DX,SV,SND,AV,SE,MODE,MEDIAN
	COMMON/DOG/MHI(2),MLO(2),ID,DX,SV(2),SND(2,2),RNG(2),
     *MENDS(2),IDA(10),MQD(10),AV(2,2),SE(2,2),MODE(2,2),
     *MEDIAN(2,2),TH(10),VH(10)
C
	INTEGER*2 IYJD,MAG,CR1,MR,ICR2,J3R
	COMMON/EASY/IYJD,MAG(8),CR1(8),MR(8),ICR2(8),J3R(8)
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
	INTEGER*2 NS,IDUM2
	REAL*4 XX,YY,Z,B1,DELT
	COMMON/INDIA/XX(124),YY(124),IDUM2(496),Z(124,2),B1(124,2),
     *DELT(124,2),NS(2)
C
	CHARACTER*26 APARM
	INTEGER*2 IPARM(4),NPARM
	REAL*4 FPARM(4)
	COMMON/OBSYDAT/FPARM,IPARM,APARM,NPARM
C
	EQUIVALENCE(X(1,1),IX(1)),(IY(1),MHI(1)),(IZS(1),IDUM2(1))
	DATA INDX/25,50,100,200,350,600,1000,1650,2500/
	DATA IEL/'H','D'/
C
	OPEN(UNIT=*,CARRIAGE CONTROL='FORTRAN')
      DV=32767
	WRITE (*,1000)
	READ (*,1010) NAM1
	CALL UPPER(NAM1,31)
	CALL GETFLN(NAM1,NAM,31)
	OPEN (4,FILE=NAM1,FORM='UNFORMATTED',STATUS='OLD',
     *ACCESS='TRANSPARENT',IOSTAT=IER,ERR=5)
	GOTO 6
    5 PRINT,'  IOSTAT ERROR CODE= ',MOD(IER,256),NAM1
	GOTO 999
    6 CALL EXNAM(NAM,NAM5,'.CLN')
      OPEN(7,FILE=NAM5,STATUS='UNKNOWN',
     *CARRIAGE CONTROL='FORTRAN',FORM='FORMATTED')
      WRITE(7,1011)NAM5
 1011 FORMAT(2X,A10)
      WRITE(*,1018)
 1018 FORMAT(1H&,'Beginning Date of Input File [i.e., MM/DD/YY]> ')
	READ(*,1005)IDATE
 1005 FORMAT(A)
	IB=GYJD(IDATE)
	WRITE(*,1017)
 1017 FORMAT(1H&,'Ending Date of Input File [i.e., MM/DD/YY]> ')
	READ(*,1005)IDATE
	IE=GYJD(IDATE)
	CALL EXNAM(NAM,NAM2,'.DIG')
	OPEN (3,FILE=NAM2,STATUS='UNKNOWN',ERR=999)
	CALL EXNAM(NAM,NAM4,'.RUN')
	OPEN(1,FILE=NAM4,STATUS='OLD',ERR=888)
	READ(1,1019)IDT                 !type of data
	READ(1,1016)NPPM                !number of points per mean
	IF((NPPM.LT.12) .OR. (MOD(1440,NPPM).NE.0))THEN
	  CLOSE(1)
	  PRINT *,'  Bad NPPM value, see *.RUN file'
	ENDIF
	NM=1440/NPPM
	READ(1,1015)NCP
	IF(NCP.LT.2)THEN
	  CLOSE(1)
	  PRINT *,'  NCP must be greater than 1, see *.RUN file'
	  GOTO 888
	ELSE
	  NCP=NCP-1             !Number of coordinate pairs used
					! to fit line segments equals NCP+1
	ENDIF
	DO I=1,2
	  READ(1,1144)(SEP(J,I),J=1,8)  !Standard normal deviate SND
						  ! for SR curve
	  READ(1,1144)(SEPA(J,I),J=1,8) !SND for degree of smoothing
      END DO
	READ(1,1145)SNH                 !Smooth nightime hours
	IF(SNH.EQ.'Y' .OR. SNH.EQ.'y')THEN
	  READ(1,1146)NBH,NEH           !Beginning and ending time (hours)
	  READ(1,1148)NTR(1),NTR(2)     !Nightime hourly range in nT
	  DO I=1,2
	    A=FLOAT(NTR(I))*10.*FLOAT(NPPM)/60.
	    NTR(I)=RND(A)
	  ENDDO
	ELSE
	  READ(1,1150)
	  READ(1,1150)
	END IF
	READ(1,1149)DEBUG
	IF(DEBUG.EQ.'y')DEBUG='Y'
	IF(DEBUG.EQ.'Y')THEN
	  READ(1,1050)IKM,IXY,IGS,ISP,ICSP,IGK
	ELSE
	  READ(1,1150)
	ENDIF
	READ(1,1013)WCRV
	IF(WCRV.EQ.'y')WCRV='Y'
	IF(WCRV.EQ.'Y')THEN
	  CALL EXNAM(NAM ,NAM3,'.CRV')
	  OPEN(2,FILE=NAM3,STATUS='UNKNOWN',FORM='UNFORMATTED',
     *  ACCESS='TRANSPARENT')
	ENDIF
	READ(1,1014)PANS
	CLOSE(1)
	IDENT=NAM(1:2)
	CALL OINFO(8,IDENT,'EAMP:')
	FMIN9=FPARM(1)
	CALL CURVK(FMIN9,INDX)        !Get K-index parameters
C
	CALL KAYPRM             !GET OBSY K PARAMETERS
	WRITE(7,1049)
	DO I=1,2
	  WRITE(7,1051)IEL(I),MHI(I),MLO(I),MENDS(I),RNG(I)
	ENDDO
	WRITE(7,1060)(MQD(I),I=1,10)
	WRITE(7,1061)(IDA(I),I=1,10)
	CF(1)=1.                !GENERAL CONVERSION VALUE
      CF(2)=1.
	IF(IDT.EQ.1)THEN        !DATA IS IN MV
	  CALL OINFO(8,IDENT,'OSV:')   !GET OMIS SCALE VALUE
	  SV(1)=FPARM(1)
	  SV(2)=FPARM(2)
        CF(1)=SV(1)*4.0		!SV IS IN nT/MV
        CF(2)=SV(2)*4.0		!SV IS IN MIN/MV
      END IF
      IF(IDT.EQ.2)THEN		!DATA IS IN NT
        CF(2)=1.0/DX		!CONVERT TO MIN
      END IF
      XD=DX
      D=0
	IEX=0
	IF(DEBUG.EQ.'y' .OR. DEBUG.EQ.'Y')THEN
	  OPEN(UNIT=8,FILE='TK0.SES',
     *  CARRIAGE CONTROL='FORTRAN',STATUS='UNKNOWN')
	ENDIF
C
    9 D=0
   10 READ(4,END=14,ERR=12,IOSTAT=IER)IX
      GOTO 15
   12 WRITE(*,1025)IER
   14 DO  E=1,2
	  DO  I=(2*NM+1),3*NM
          MHV(I,E)=DV
        END DO
      END DO
      IEX=1
      GOTO 30
   15 IF(IX(5762).LT.IB)GOTO 9
	WRITE(7,1030)IX(5762)
	CALL WRITEDAY(IX(5762))
	D=D+1
	IF(IDT.NE.0)THEN              !CONVERT IF NOT GAMMA DATA
        DO I=1,2
           DO J=1,1440
              IF(X(J,I).NE.DV) X(J,I)=X(J,I)*CF(I)
           END DO
        END DO
      END IF
	IF(IX(5762).GT.IE) IEX=1
	CALL KMEAN(D,NPPM,NM)          !Select MHV using confidence
					! limits input in subroutine
					!  KAYPRM, modified 7/27/90
      IF(D.LT.2)THEN
        DO E=1,2
	     DO I=(NM+1),2*NM
             MHV(I,E)=DV
           END DO
        END DO
        IYJD=IX(5762)-1
	  IEX=0
	  IF(WCRV.EQ.'Y')WRITE(2,ERR=888)IY
        GOTO 40
      END IF
   30 IYJD=IYJD+1
      DO E=1,2
	  CALL XYSET(E,NM,NPPM,NPTS)    !Generate coordinate pairs
				!(MHV,HOUR) for each day for interval -11
				!to 35.5 or -1.5 to 25.5 (depending on missing
				!MHVs) to assure smooth SR curve tran-
				!sition from day to day
	  CALL GS(E,NM,NPPM,NPTS,NCP)   !Fit coordinate pairs with line
				!segments to estimate SR curve
	  CALL SPLN2(E,NM)         !Generate cubic spline parameters from
				!points of line intercepts generated in
				!subroutine GS.
        CALL GKR(E)		!Get component amplitude, SR contribu-
				!tion, range, and estimate trucation due to
				!1-minute mean.
      END DO
					!Output spline parameters to '*.crv' file
	IF(WCRV.EQ.'Y')WRITE(2,ERR=888)IZS
	CALL KTAB               !Write *.DIG file of 3-hr amplitude, SE
				!of SR curve, SR curve contribution, and
				!possible error due to 1-minute mean
				!truncation.
   40 CALL MOVDAT(NM)       !Shift day of data in MHV array to
				!allow space for next day of MHVs
	IF(IEX.EQ.0)GOTO 10     !Return for next day
	CLOSE(2)
	CLOSE(3)
	CALL KTABLE(NAM,PANS,FMIN9,INDX)
  888 CLOSE(4)
      CLOSE(3)
	CLOSE(2)
	CLOSE(7)
      CLOSE(8)
  999 STOP
 1000 FORMAT(' Enter Path to File of NGIC 1-minute Values : ')
 1010 FORMAT(A31)
 1013 FORMAT(18X,A1)
 1014 FORMAT(13X,A1)
 1015 FORMAT(28X,I2)
 1016 FORMAT(26X,I3)
 1019 FORMAT(10X,I2)
 1025 FORMAT('  I/O ERROR=',I6)
 1030 FORMAT(I6)
 1049 FORMAT('       HI    LO   END  RNG')
 1050 FORMAT(6I2)
 1051 FORMAT(2X,A1,3I6,I5)
 1060 FORMAT(2X,'Magnitude',10I6)
 1061 FORMAT(2X,'Day      ',10I6)
 1144 FORMAT(11X,8F6.2)
 1145 FORMAT(15X,A1)
 1146 FORMAT(28X,2I3)
 1148 FORMAT(21X,2I4)
 1149 FORMAT(7X,A1)
 1150 FORMAT()
	END
