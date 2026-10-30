      SUBROUTINE MSTAT(IY,NPTS,KN,NTV,LST,INTB,INTL,IZ)
C
C       IY:  ARRAY OF DATA
C       NPTS:  NUMBER OF ELEMENTS IN IY
C       KN:  OFFSET/STEP IN IY ELEMENTS
C       NTV:  NUMBER OF INTERVALS FOR DISTRIBUTION
C       LST:  0-> DO NOT LIST DISTRIBUTION, SE, MEAN, MODE, AND MEDIAN
C             1-> DO LIST
C       INTB:  NUMBER OF INTERVALS > MODE INTERVAL FOR HIGH CLEANING
C              LIMITS.  ZERO SETS LIMIT TO MAXIMUM
C       INTL:  NUMBER OF INTERVALS < MODE INTERVAL FOR LOW CLEANING
C              LIMITS.  ZERO SETS LIMIT AT MINIMUM
C       IZ:  ARRAY TO BE CLEANED USING INTB/INTL LIMITS.  ASSUMES
C            IZ ELEMENTS CORRESPOND TO IY.
C
c  Purpose:  define the mean value distribution by generating a
c            histogram, compute the average mean value and the
c            standard deviation of the mean values about the average,
c            and compute the mode and median of the distribution
c
c  Author:  L. R. Wilson 10/16/92
c
c  Subroutines:  SORT0
C
	INTEGER*2 RND,IY,IX,IH,IZ,IHA,IDUM,IB,IL,IDV,MDV,
     *ICL,J,I,KPTS,INT,NTV,KN,LST,INTB,INTL,IXL,IR,IT,
     *IHAA,MAXF,NPTS,IXB,NIX
      REAL*4 MODE,MEDIAN,SG,SFD,SFD2,HINT,A,B,C
	DIMENSION IY(1),IX(3720),IH(10),IZ(1),IHA(10)
      COMMON/EASY/IDUM,SE,IB,IL,AVRG,MODE,MEDIAN
C
	IDV=32767
      MDV=-32767
      ICL=0
      IF(INTB.NE.0.OR.INTL.NE.0)ICL=1
C *** SAVE IY
    1 J=1
      DO 2 I=1,NPTS,KN
      IX(J)=IY(I)
      J=J+1
    2 CONTINUE
      KPTS=J-1
      DV=32767
      CALL SORT0(IX,KPTS,IXB,IXL,NIX)	!ELEMENT 1 SMALLEST, ELEMENT NIX
					! LARGEST
      INT=IXB-IXL			!RANGE OF VALUES
      I=INT/NTV				!INTERVAL WIDTH
      IR=INT-I*NTV			!REMAINDER
      IF(IR.NE.0)INT=INT+NTV-IR		!EXPAND RANGE TO INCLUDE REMAINDER
      INT=INT/NTV			!RECOMPUTE INTERVAL
      IT=IXL+INT			!TOP OF LEAST INTERVAL
C
      J=1				!COMPILE FREQUENCY DISTRIBUTION
      DO 4 I=1,NTV
      IHAA=0
    4 IH(I)=0
      DO 10 I=1,KPTS
    5    IF(IX(I).EQ.DV)GOTO 10
         IF(IX(I).GE.IT)GOTO 8
         IH(J)=IH(J)+1
         IF(I.LT.KPTS)GOTO 10
    8    IHAA=IHAA+IH(J)
         IHA(J)=IHAA
         J=J+1
         IT=IT+INT
         IF(I.LT.KPTS)GOTO 5
   10 CONTINUE
      IF(LST.EQ.1)WRITE(1,1000)(IH(I),I=1,NTV)
 1000 FORMAT(2X,10I6)
      SG=0.
      SFD=0.
      SFD2=0.
      MAXF=0
      HINT=INT-1
      HINT=HINT*0.5
      A=IXL
      A=A+HINT
      IV=RND(A)
      IF(LST.EQ.1)WRITE(1,1005)
 1005 FORMAT('  ')
      DO 20 I=1,NTV
      IF(LST.EQ.1)WRITE(1,1010)IV
 1010 FORMAT('&',I6)
      IF(IH(I).LE.MAXF)GOTO 12
      MAXF=IH(I)
      MF=I
      L1=IXL+(I-1)*INT
   12 IV=IV+INT
      A=I-5
      B=IH(I)
      FD=A*B
      FD2=A*FD
      SG=SG+B
      SFD=SFD+FD
      SFD2=SFD2+FD2
   20 CONTINUE
      IF(ICL.EQ.0)GOTO 100
      K=MF+INTB
      IF(K.GT.NTV)K=MF
      N=IHA(K)
      IB=IX(N)
      K=MF-INTL-1
      IF(K)21,21,22
   21 N=1
      IL=IX(N)
      GOTO 23
   22 N=IHA(K)+1
      IL=IX(N)
   23 IF(INTB.EQ.0)IB=IDV
      IF(INTL.EQ.0)IL=MDV
      DO 98 I=1,NPTS,KN
      IF(IY(I).EQ.IDV)GOTO 98
      IF(IY(I).LE.IB.AND.IY(I).GE.IL)GOTO 98
      IY(I)=IDV
      IZ(I)=IDV
   98 CONTINUE
      ICL=0
      IF(LST.EQ.1)WRITE(1,1020)
 1020 FORMAT('  ')
      GOTO 1
  100 A=IXL
      A=A+HINT
      BINT=INT
      A=A+4.0*BINT
      AVRG=A+(SFD/SG)*BINT
      SE=(SG*SFD2-SFD*SFD)/(SG*(SG-1.0))
      SE=SQRT(SE)*BINT
      IF(MF.GT.1)GOTO 24
      A=MAXF
      B=MAXF-IH(MF+1)
      GOTO 28
   24 IF(MF.GE.NTV)GOTO 27
      A=MAXF-IH(MF-1)
      B=MAXF-IH(MF+1)
      GOTO 28
   27 A=MAXF-IH(MF-1)
      B=MAXF
   28 C=L1
      MODE=C+(A/(A+B))*BINT
      I=MOD(NIX,2)
      IF(I.EQ.0)GOTO 30
      I=NIX/2
      MEDIAN=IX(I)
      GOTO 35
   30 I=NIX/2
      A=IX(I)+IX(I+1)
      MEDIAN=A/2.0
   35 CONTINUE
      IF(LST.EQ.1)WRITE(1,1030)AVRG,SE,MODE,MEDIAN
 1030 FORMAT('  AV,SE,MODE,MEDIAN:',F8.1,F14.2,2F8.1)
      IF(LST.EQ.1)WRITE(1,1020)
      RETURN
      END
