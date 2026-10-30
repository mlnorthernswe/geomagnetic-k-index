      REAL FUNCTION PSPLIN*4(K,Y,E)
      INTEGER*2 E,IX,IW,NS,K
      REAL*4 U,X,S,DEL,Y
      DIMENSION X(28,2),U(28,2),S(28,2),DEL(28,2),NS(2)
      COMMON/ABLE/IX(5762),IW(2882)
      EQUIVALENCE(U(1),IW(1)),(X(1),IW(113)),
     1(S(1),IW(225)),(DEL(1),IW(337)),(NS(1),IW(449))
C
C      TYPE *,'  K,Y,E,NS(E):',K,Y,E,NS(E)
C      TYPE *,(X(I,E),U(I,E),S(I,E),DEL(I,E),I=1,NS(E))
    2 IF(Y.LE.X(K,E))GOTO 4
      K=K+1
      GOTO 2
    4 K1=K-1
      FF1=S(K1,E)*(X(K,E)-Y)**3
      FF2=S(K,E)*(Y-X(K1,E))**3
      FF3=1.0/(6.0*DEL(K,E))
      F1=(FF1+FF2)*FF3
      F2=(Y-X(K1,E))*((U(K,E)/DEL(K,E))-(S(K,E)*DEL(K,E))/6.0)
      F3=(X(K,E)-Y)*((U(K1,E)/DEL(K,E))-(S(K1,E)*DEL(K,E))/6.0)
      PSPLIN=F1+F2+F3
      RETURN
      END
