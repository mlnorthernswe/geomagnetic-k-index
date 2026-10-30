      REAL*4 FUNCTION RKV(V,FMIN9,INDEX,J)
      INTEGER*2 RND,INDEX,J,IV
      REAL*4 FMIN9,V,SV,VT,VK
      DIMENSION INDEX(1)
      SV=V
      IV=V
      IF(IV.GE.INDEX(1))GOTO 22
      VT=INDEX(1)
      VK=V/VT
      GOTO 30
   22 V=V/FMIN9
      V=ALOG10(V)-1.39794
      VK=1.003+3.232717*V+0.190954*V*V*V
   30 VK=VK*10.
      J=RND(VK)
      VK=J*0.1
      IF(VK.GE.10)VK=9.9     !REPRESENTS MAXIMUM DIGITAL K
      J=VK+1
      RKV=VK
      V=SV
      RETURN
      END
