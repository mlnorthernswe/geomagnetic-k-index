	PROGRAM SELMVS
C  Purpose:  Use a specified number of points (NPPM) to compute
c           mean values (MHV) and range of 1-minute values (RHR)
c           to generate statistics to be used in estimating regular
c           variation (SR) curves.  The statistics include H and D
c           frequency distributions of MHV and RHR; their mean, mode,
C               and median; and their confidence limits.
C
C  Author:  L. R. Wilson 10/16/92
c
C  Input:  A file of coded 1-minute component variation data from
c          the user's working directory
c       :  The file 'OBSY.DAT' contains parameters used in program
c          must be in the user's working directory
c       :  The following user input
C       -  The beginning and ending dates for the input file
c       -  Standarda normal deviate (SND) for computing con-
C          fidence limits for RHR.
C       -  SND for computing confidence limits for MHVs.
C       -  Data type: gamma, omis, or syn.
C       -  Omis scale value to convert mv to nT or minutes.
C       -  Previous year's annual mean (PYAM)
C       -  Lower limit of K=5 for observatory is used to select
C          quiet day RHR.
c
C  Output:  File of RHR and MHV statistics with '.kay' extension,
c           note that the program will over-write an existing file
c           with the same name.
C
C  Subroutines: MSTAT,SORT2
C  GENERAL.LIB: EXNAM,WRITEDAY,GYJD
c
c  Comments:  Program and subroutines were compiled and linked
c             using Lahey F77L
C
	INTEGER*2 X,H,E,D,DV,RHR,IX,IDT,IDUM,IB,IL,MHV,IY,MENDS,
     *MHI,MLO,IRNG,LLK5,MDR,MDRS,IBD,IED,GYJD,NM,NPPM,NHRS,IT
	REAL*4 MODE,MED,SV
	CHARACTER*1 IEL,CHA
	CHARACTER*2 ID
	CHARACTER*8 IDATE
	CHARACTER*12 NAM,NAM3
	CHARACTER*26 APARM
	CHARACTER*30 NAM1
	INTEGER*2 IPARM(4),NPARM
	REAL*4 FPARM(4)
	DIMENSION X(1440,2),MHV(120,31,2),IY(120),RHR(120,31,2),
     *PYAM(4),SV(3),MENDS(2),MHI(2),MLO(2),IRNG(2),MDR(31),SNDHR(2),
     *SNDMHV(2),IEL(2),LLK5(2),CF(2),MDRS(31),IT(3720)
      COMMON/ABLE/IX(5762)
	COMMON/EASY/IDUM,SE,IB,IL,AVRG,MODE,MED
	COMMON/OBSYDAT/FPARM,IPARM,APARM,NPARM
	EQUIVALENCE(X(1,1),IX(1))
      DATA IEL/'H','D'/
C
	OPEN(UNIT=*,CARRIAGE CONTROL='FORTRAN')
      DV=32767
      MDV=-32767
	SNDHR(1)=1.96
      SNDHR(2)=1.96
      SNDMHV(1)=2.576
      SNDMHV(2)=2.576 
	WRITE(6,1100)
 1100 FORMAT ('  Enter NGIC Input Filename : ')
	READ (*,1110) NAM1
 1110 FORMAT (A)
	CALL UPPER(NAM1,30)
	CALL GETFLN(NAM1,NAM,30)
	OPEN (UNIT=4,FILE=NAM1,FORM='UNFORMATTED',
     *ACCESS='TRANSPARENT',STATUS='OLD',ERR=998)
	WRITE(*,1035)
 1035 FORMAT(1H&,' Beginning Date of the Input File',
     *' [i.e., MM/DD/YY]> ')
	READ(*,1110)IDATE
	IBD=GYJD(IDATE)
	WRITE(*,1036)
 1036 FORMAT(1H&,' Ending Date of the Input File',
     *' [i.e., MM/DD/YY]> ')
	READ(*,1110)IDATE
	IED=GYJD(IDATE)
	WRITE(*,1130)
 1130 FORMAT('&',' Enter Number of Minute Values per Mean',/,
     *'   Must be MOD 1440, and Greater Than or Equal to 12',/,
     *'   [i.e., 12,20,24,30,... ]> ')
	READ(*,*)NPPM
	WRITE(*,1019)
 1019 FORMAT(1H&,' SND of R1MVs is 1.96, change(Y/[N]) ? ')
	READ(*,1080)CHA
      IF(CHA.EQ.'Y' .OR. CHA.EQ.'y')THEN
	  WRITE(*,1020)
 1020   FORMAT(1H&,' Enter standard deviate for H and D : ')
	  READ(*,*) SNDHR(1),SNDHR(2)
      END IF
	WRITE(*,1021)
 1021 FORMAT(1H&,' SND of MVs is 2.576, change(Y/[N]) ? ')
	READ(*,1080)CHA
 1080 FORMAT(A1)
	IF(CHA.EQ.'Y' .OR. CHA.EQ.'y')THEN
	WRITE(*,1020)
	  READ(*,*) SNDMHV(1),SNDMHV(2)
      END IF
C
c      WRITE(*,1022)
c 1022 FORMAT(1H&,' Data Type :  0=DCP, 1=OMIS = ')
c      READ (*,*) IDT
	IDT=0
	ID=NAM(1:2)
	CALL OINFO(2,ID,'PYAM:')
	DO I=1,NPARM
	  PYAM(I)=FPARM(I)
	ENDDO
      XD=PYAM(1)/3438.
      CF(1)=1.0
      CF(2)=1.0
      SV(1)=1.0                 !GAMMA SCALE VALUE
      SV(2)=1.  
      IF(IDT.EQ.1)THEN          !DATA IS TENTH-MV/4
	  CALL OINFO(2,ID,'OSV:')  !GET OMIS SCALE VALUE
	  DO I=1,NPARM
	    SV(I)=FPARM(I)
	  ENDDO
	CF(1)=SV(1)*4.0         !SV IS IN nT/MV
	CF(2)=SV(2)*4.0         !SV IS MIN/MV
      END IF
	WRITE(*,1023)
 1023 FORMAT(1H&,' Lower Limit of K=5 for this Observatory = ')
	READ(*,*)CLLK5
	A=CLLK5*10.*FLOAT(NPPM)/60.
	LLK5(1)=NINT(A)                             !in tnT
	A=A*34380./PYAM(1)
	LLK5(2)=NINT(A)                             !in tenth-min
C
	CALL EXNAM(NAM,NAM3,'.KAY')
	OPEN(UNIT=1,FILE=NAM3,STATUS='UNKNOWN',
     *CARRIAGE CONTROL='FORTRAN',FORM='FORMATTED')
	WRITE(1,1040)NAM1
 1040 FORMAT(3X,A)
	WRITE(1,1042)NPPM
 1042 FORMAT(' Number of Minute Values per Mean =',I3)
      WRITE(1,1050)XD,SV(1),IEL(1),SV(2),IEL(2)
 1050 FORMAT(' XD,SV:',F8.4,F6.2,A1,2X,F6.2,A1)
      WRITE(1,1059)
 1059 FORMAT(13X,'MVs',16X,'RANGE')
      WRITE(1,1060)SNDMHV(1),IEL(1),SNDMHV(2),IEL(2),
     1SNDHR(1),IEL(1),SNDHR(2),IEL(2)
 1060 FORMAT(' SND',4(F7.3,A1,2X))
C
      D=0
	NM=1440/NPPM
   10 READ(4,END=34,ERR=14,IOSTAT=IER)IX
      GOTO 12
   14 WRITE(*,1105)IER
 1105 FORMAT('  IOSTAT=',I6)
      GOTO 999
   12 IF(IX(5762).LT.IBD)GOTO 10
	IF(IX(5762).GT.IED)GOTO 35
	D=D+1
	CALL WRITEDAY(IX(5762))
C
      IF(IDT.NE.0)THEN                  !IF NOT GAMMA DATA, CONVERT
	DO I=1,2
	  DO J=1,1440
	    IF(X(J,I).NE.DV)THEN
	      A=X(J,I)
	      A=A*CF(I)
	      X(J,I)=NINT(A)
	    END IF
	  END DO
	END DO
      END IF
	DO 32 E=1,2
	   I=0
	 NSD=0
	 SD=0.
	 IF(E.EQ.1)THEN
	   MDR(D)=DV
	   SMDR=0.
	   NMDR=0
	 END IF
	   DO 30 H=1,NM
	    SH=0.
	    NSH=0
	    MHV(H,D,E)=DV
	    MAX=MDV
	    MIN=DV
	    RHR(H,D,E)=DV
		DO 20 L=1,NPPM
		   I=I+1
	       IF(X(I,E).EQ.DV)GOTO 20
	       IF(X(I,E).GT.MAX)MAX=X(I,E)
	       IF(X(I,E).LT.MIN)MIN=X(I,E)
	       SH=SH+X(I,E)
	       NSH=NSH+1
   20       CONTINUE
		IF(NSH.LT.(NPPM/2))GOTO 30
		A=SH/FLOAT(NSH)
		MHV(H,D,E)=NINT(A)
	    SD=SD+MHV(H,D,E)
	    ZZ=FLOAT(MAX)-FLOAT(MIN)
	    IF(ZZ.LT.32767.)THEN
	      RHR(H,D,E)=MAX-MIN
	    ELSE
	      RHR(H,D,E)=32767
	    END IF
	    IF(E.EQ.1)THEN
	      SMDR=SMDR+RHR(H,D,E)
	      NMDR=NMDR+1
	    END IF
	    IF(RHR(H,D,E).GT.LLK5(E))RHR(H,D,E)=DV
   30    CONTINUE
	   IF(NSD.LT.(NM/2))GOTO 32
   32 CONTINUE
	IF(NMDR.GE.(NM/2))THEN
	CC=SMDR/FLOAT(NMDR)
	MDR(D)=NINT(CC)
      END IF
      GOTO 10
   34 CONTINUE
   35 NDS=D
      DO I=1,31                         !SETUP MDRS SEQUENCE
	MDRS(I)=I
      END DO
      CALL SORT2(MDR,MDRS,NDS,IXB,IXL,N) !SORT QUIET DAYS IN ASCENDING
					 ! ORDER                
	NHRS=NDS*NM
	DO 60 E=1,2
	 WRITE(1,1070)IEL(E)
 1070    FORMAT(1X,A1,2X,'MV FREQUENCY DISTRIBUTION')
	   K=1
	   DO J=1,NDS
		DO I=1,NM
		  IT(K)=MHV(I,J,E)
		  K=K+1
		ENDDO
	   ENDDO
	   CALL MSTAT(IT(1),NHRS,1,10,1,0,0,IY)
	 B=SNDMHV(E)*SE
	 A=AVRG+B
	 MHI(E)=NINT(A)
	 A=A-2.0*B
	 MLO(E)=NINT(A)
	 MENDS(E)=NINT(MODE)
	 WRITE(1,1071)IEL(E)
 1071    FORMAT(1X,A1,2X,'RANGE FREQUENCY DISTRIBUTION')
	   K=1
	   DO J=1,NDS
		DO I=1,NM
		  IT(K)=RHR(I,J,E)
		  K=K+1
		ENDDO
	   ENDDO
	   CALL MSTAT(IT(1),NHRS,1,10,1,0,0,IY)
	 B=SNDHR(E)*SE
	 A=AVRG+B
	 IRNG(E)=NINT(A)
	 WRITE(1,1009)IEL(E)
 1009    FORMAT(1X,A1,4X,'MAX',3X,'MIN   IST  RNG')
	 WRITE(1,1010)MHI(E),MLO(E),MENDS(E),IRNG(E)
 1010    FORMAT(3X,3I6,I5)
   60 CONTINUE
      WRITE(1,1140)(MDR(J),J=1,10)
 1140 FORMAT(5X,10I6)
      WRITE(1,1140)(MDRS(J),J=1,10)
	CLOSE(UNIT=1)
	GOTO 999
  998 WRITE(*,1150)NAM1
 1150 FORMAT(/,'  COULD NOT OPEN FILE ',A)
  999 CLOSE(4)
      STOP
      END
