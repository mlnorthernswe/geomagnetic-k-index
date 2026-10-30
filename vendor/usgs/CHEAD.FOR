      SUBROUTINE CHEAD(ID,MO,IDAY,IYR,EL,HEAD)
      CHARACTER*(*) ID,HEAD,EL
      CHARACTER*3 JMO,IMO
      DIMENSION JMO(12)
      DATA JMO/'JAN','FEB','MAR','APR','MAY','JUN','JUL',
     *'AUG','SEP','OCT','NOV','DEC'/
C
      HEAD(1:3)=ID(1:3)
      HEAD(4:5)='  '
      J=IDAY/10
      K=IDAY-(J*10)
      HEAD(6:6)=CHAR(J+48)
      HEAD(7:7)=CHAR(K+48)
      HEAD(8:8)='-'
      IMO=JMO(MO)
      HEAD(9:11)=IMO(1:3)
      HEAD(12:12)='-'
      J=IYR/10
      K=IYR-(J*10)
      HEAD(13:13)=CHAR(J+48)
      HEAD(14:14)=CHAR(K+48)
      HEAD(15:17)='   '
      HEAD(18:19)=EL(1:2)
      RETURN
      END

