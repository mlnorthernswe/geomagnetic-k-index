	SUBROUTINE KTAB
C  Purpose:  Select maximum  H or D component amplitude and write
C            the amplitude and corresponding standard error of the
C            SR curve, possible error in amplitude due to minute
C            average (carries sign to indicate component choice H
C            or D), the magnitude of the SR contribution, and the
C            3-hr range to the *.dig file.
c
      INTEGER*2 RND,I,NK
      REAL*4 DV
C
      INTEGER*2 IX
      COMMON/ABLE/IX(5762)
C
      INTEGER*2 LX,YJD
      COMMON/CHARLY/LX(2882),YJD
C
	INTEGER*2 IYJD,MAG,MR,ICR2,J3R
	REAL*4 CR1
	COMMON/EASY/IYJD,MAG(8),CR1(8),MR(8),ICR2(8),J3R(8)
C
      INTEGER*2 MP,IC2,I3R
      REAL*4 DEL,SESRC
      COMMON/GOLF/DEL(8,2),SESRC(8,2),MP(8,2),IC2(8,2),I3R(8,2)
C
      DV=32767.
      NK=0
      DO 60 I=1,8
      IF(DEL(I,1).EQ.DV.OR.DEL(I,2).EQ.DV)GOTO 15
      GOTO 20
   15 MAG(I)=9999
	MR(I)=9999
	ICR2(I)=9999
	J3R(I)=9999
      GOTO 60
   20 IF(DEL(I,1).GT.DEL(I,2))GOTO 25
      MAG(I)=RND(DEL(I,2))			!AMPLITUDE
      CR1(I)=SESRC(I,2)				!SE OF SR CURVE
      ICR2(I)=IC2(I,2)				!SR CONTRIBUTION
      MR(I)=-MP(I,2)				!POSSIBLE ERROR DUE TO 
						!1-MIN AVERAGE
	J3R(I)=I3R(I,2)                     !3-hr range
      GOTO 29
   25 MAG(I)=RND(DEL(I,1))
      MR(I)=MP(I,1)
      CR1(I)=SESRC(I,1)
      ICR2(I)=IC2(I,1)
      J3R(I)=I3R(I,1)
   29 NK=NK+1
   60 CONTINUE
      WRITE(3,1000) IYJD,(MAG(I),CR1(I),MR(I),ICR2(I),J3R(I),
     *I=1,8)
 1000 FORMAT(I6,4(I6,F6.1,3I6),/,6X,4(I6,F6.1,3I6))
      RETURN
      END
