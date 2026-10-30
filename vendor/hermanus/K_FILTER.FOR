C**************************************************************************
C*                                                                        *
C*       ESTIMATION OF THE QUIET DAY PATTERN USING THE LRNS METHOD        *
C*                                                                        *
C*    For a mathematical description on the LRNS method (i.e. LINEAR      *
C*    PHASE, ROBUST NON-LINEAR SMOOTHING METHOD) refer to the article:    *
C*    Hattingh, M., L. Loubser & D. Nagtegaal - "Computer K-index         *
C*    estimation by a new linear phase, robust smoothing method",         *
C*    Geophys. J. Int., 99, 533-547, 1989.                                *
C*                                                                        *
C*    Developed by M. Hattingh                                            *
C*    FORTRAN-77 coding by M. M. Van Dyk                                  *
C*                                                                        *
C*    Hermanus Magnetic Observatory                                       *
C*    C S I R                                                             *
C*    P.O. Box 32,  HERMANUS 7200,  South Africa.                         *
C*    Telephone:  + 27-283-21196      Facsimilie:  + 27-283-22039         *
C*                                                                        *
C**************************************************************************
C*                                                                        *
C*                     HOW  TO  USE  THE  SUBROUTINE                      *
C*                                                                        *
C*      This subroutine needs only one day's digital data                 *
C*      (i.e. 1440 data points of 1-minute data) in nT.                   *
C*                                                                        *
C*      Take a few quiet days of your observatory, determine "q start"    *
C*      the approximate starting time of a quiet day pattern and          *
C*      "q end" the approximate end time of a quiet day pattern           *
C*      (values are not critical)                                         *
C*                                                                        *
C*      When the quiet day pattern occurs in the middle of a UT day       *
C*      set the "q index" = "M".                                          *
C*                                                                        *
C*      When the quiet day pattern occurs at the end of a UT day          *
C*      set the "q index" = "E".                                          *
C*                                                                        *
C*      Set the "q factor" = reference amplitude value below which        *
C*      K = 0 for your observatory.                                       *
C*                                                                        *
C*      Vary the smoothing factor "mu" to fit a typical quiet day,        *
C*      normally between 0.08 and 0.2 on the standard world data set.     *
C*      Try this value of mu on a selection of quiet, medium and          *
C*      disturbed days and fine tune to get the best estimate of the      *
C*      hand values taken, eg over a year's data.                         *
C*                                                                        *
C*      Note!  Dramatic estimates of the hand values can be achieved      *
C*      if only the mu factor is slightly varied in a man/computer        *
C*      interactive way, but to use it without human intervention,        *
C*      a fixed value of mu should be used for all the data.              *
C*      This eliminates subjectiveness.                                   *
C*                                                                        *
C**************************************************************************
C*                                                                        *
C*                      Parameters for Hermanus                           *
C*                                                                        *
C*  Components  q-index   q-factor   q-start   q-end   h-absolute   mu    *
C*                                                                        *
C*     D & H       M        3.0        300     1009      11115     0.17   *
C*                                                                        *
C**************************************************************************
C
      subroutine kfilter
     +           ( signal,raw,smooth,mu,qindex,qfactor,qstart,qend )
C
C**************************************************************************
C                                                                         *
C     signal - signal to be smoothed (free of spikes and steps)           *
C        raw - copy of original signal                                    *
C     smooth - smoothed signal (estimate of quiet day or Sq pattern)      *
C         mu - smoothing factor (between 0.08 and 0.20)                   *
C     qindex - M(middle) - quiet day pattern in middle of a UT day        *
C            - E(end)    - quiet day pattern at  ends  of a UT day        *
C    qfactor - reference amplitude value below which K=0 for observatory  *
C     qstart - start of quiet day pattern in UT                           *
C       qend - end of quiet day pattern in UT                             *
C     u2fact - ratio of L2/L1 smoothing factor (0.08)                     *
C          x - reference signal (initial guide of the adaptive filter)    *
C          w - data adaptive filter weight                                *
C         sw - starting weight                                            *
C          y - forward filter output                                      *
C         yy - reverse filter output                                      *
C         n2 - number of data points per day (1440 for 1-minute data)     *
C    cascade - order of filter (for K-index determination = 2)            *
C                                                                         *
C**************************************************************************
C
       implicit none
         real*4 qfactor,sss,sw,mu,u,u2,u2fact,w
         real*4 raw(*),signal(*),smooth(*),error
         real*4 sum, x(1440), y(1440), yy(1440)
         real*4 pmax,pmin,difrange
        integer i,cascade,iter,lz,n1,n2,qstart,qend
      character qindex*1
         real*8 ra,random                                               
C
C*** Initialise.
C
      cascade = 2      
      n1 = 1
      n2 = 1440
      pmax = 10.0E-06
      pmin = 10.0E+06
      sw = 0.0
      u2fact = 0.08
      do 10  i = n1, n2
         y(i) = 0.0
         yy(i) = 0.0
         ra = random()                           ! random number
         x(i) = 1.0 + 0.005*ra
   10 continue
C
C*** Subtract mean.
C
      sum = 0.0
      do 20  i = n1, n2
         sum = sum + signal(i)
   20 continue
      sum = sum / (n2 - n1 + 1)                  ! find average  
      do 30  i = n1, n2
         signal(i) = signal(i) - sum             ! subtract average
   30 continue
      do 40  i = 1, n2
         raw(i) = signal(i)
         smooth(i) = signal(i)
   40 continue
      u = mu
      u2 = u * u2fact
      do 50  i = n1, n2-2
         sss = signal(i+2) - signal(i)
         if( sss.gt.pmax )  pmax = sss
         if( sss.lt.pmin )  pmin = sss
   50 continue
      difrange = pmax - pmin
      if( difrange.lt.qfactor )  then
         u2 = u2 * 3.0                           ! skip L1 - norm  
      else   
         if( qindex.eq.'M' )  then               ! quiet day pattern in middle
            do 60  i = n1, n2
               sw = sw + signal(i)
   60       continue
            sw = sw / n2                         ! average
         else
            do 70  i = n1, 4                     ! quiet day pattern at end
               sw = sw + signal(i)
   70       continue
            sw = sw / 4
         endif
C
C*** Start smoothing operation (L1-NORM : Rough smoothing)
C
         do 110  iter = 1, cascade         
            w = sw
C
C*** Forward filter operation
C           do 80  i = n1, n2
               y(i) = y(i) + w*x(i)              ! filter ref signal   
               error = signal(i) - y(i)
               lz = ifix( sign(1.0, error) )     ! determines sign of error
               if( error.eq.0.0 )  lz = 0
               if( abs(error).lt.1.0E-3 )  lz = 0
               if( qindex.eq.'M' )  then         ! update w for 'M' case
                  if( i.le.qstart .or. i.ge.qend )  then
                     w = w + 2.0*u/3.0*lz*x(i)   ! w outside Sq pattern
                  else
                     w = w + 2.0*u*lz*x(i)       ! w inside Sq pattern
                  endif
               else                              ! update w for 'E' case
                  if( i.gt.qstart .and. i.lt.qend )  then
                     w = w + 2.0*u/3.0*lz*x(i)
                  else
                     w = w + 2.0*u*lz*x(i)
                  endif
               endif
   80       continue 
C
C*** Reverse filter operation
C
            if( qindex.eq.'E' )  w = signal(n2)  
            do 90  i = n2, n1, -1 
               yy(i) = yy(i) + w*x(i)
               error = signal(i) - yy(i)
               lz = ifix( sign(1.0,error) )
               if( error.eq.0.0 )  lz = 0
               if( abs(error).lt.1.0E-3 )  lz = 0
               if( qindex.eq.'M' )  then
                  if( i.le.qstart .or. i.ge.qend )  then
                     w = w + 2.0*u/3.0*lz*x(i)
                  else
                     w = w + 2.0*u*lz*x(i)
                  endif
               else
                  if( i.gt.qstart .and. i.lt.qend )  then
                     w = w + 2.0*u/3.0*lz*x(i)
                  else
                     w = w + 2.0*u*lz*x(i)
                  endif
               endif
   90       continue
C
C*** Average of Forward & Reverse process to remove Phase shift.
C
            sw = yy(n1)
            do 100  i = n1, n2
               signal(i) = signal(i) - (y(i) + yy(i))/2.0
               y(i) = 0.0
               yy(i) = 0.0
  100       continue   
            sw = signal(n1)
  110    continue
         do 120  i = n1, n2
            smooth(i) = raw(i) - signal(i)
            signal(i) = raw(i) - signal(i)
  120    continue
      endif
C
C*** L2-NORM Section (Final smoothing)
C
      do 160  iter = 1, cascade
         w = signal(n1)
         do 130  i = n1, n2
            y(i) = y(i) + w*x(i)
            error = signal(i) - y(i)
            if( qindex.eq.'M' )  then
               if( i.le.qstart .or. i.ge.qend )  then
                  w = w + 2.0*u2/2.0*error*x(i)
               else
                  w = w + 2.0*u2*error*x(i)
               endif
            else
               if( i.gt.qstart .and. i.lt.qend )  then
                  w = w + 2.0*u2/2.0*error*x(i)
               else
                  w = w + 2.0*u2*error*x(i)
               endif
            endif
  130    continue
         if( qindex.eq.'E' )  w = signal(n2)
         do 140  i = n2, n1, -1
            yy(i) = yy(i) + w*x(i)
            error = signal(i) - yy(i)
            if( qindex.eq.'M' )  then
               if( i.le.qstart .or. i.ge.qend )  then
                  w = w + 2.0*u2/2.0*error*x(i)
               else
                  w = w + 2.0*u2*error*x(i)
               endif
            else
               if( i.gt.qstart .and. i.lt.qend )  then
                  w = w + 2.0*u2/2.0*error*x(i)
               else
                  w = w + 2.0*u2*error*x(i)
               endif
            endif
  140    continue
         do 150  i = n1, n2
            signal(i) = signal(i) - (y(i) + yy(i))/2.0  
            y(i) = 0.0
            yy(i) = 0.0
  150    continue
  160 continue
C
      do  i = n1, n2
         smooth(i) = smooth(i) - signal(i)       ! final estimate of Sq pattern
      enddo
C
      return
      end
