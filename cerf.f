C=======================================================================  
      FUNCTION CERF(Z,EMACH)  
      IMPLICIT NONE 
C     ------------------------------------------------------------------  
C     GIVEN COMPLEX ARGUMENT Z, PROVIDES THE FADDEEVA 
C     COMPLEX ERROR FUNCTION (Eq. (7.2.3) of \ct{Ol}): 
!
!     W(Z)=EXP(-Z**2)*(1.0-ERF(-I*Z))  ... 7.1.3 AS; 7.2.3 Ol
!
C     THE  EVALUATION  ALWAYS  TAKES   PLACE  IN  THE  FIRST   QUADRANT.  
C     ONE  OF  THREE METHODS  IS  EXPLOYED  DEPENDING ON THE SIZE OF THE 
C     ARGUMENT (A POWER SERIES,A RECURRENCE BASED ON CONTINUED FRACTIONS  
C     THEORY, OR AN ASYMPTOTIC SERIES). 
C     EMACH IS THE MACHINE ACCURACY  
!
!     THE ARGUMENT IS TRANSLATED TO THE FIRST QUADRANT FROM  
!     THE NN_TH QUADRANT, BEFORE THE METHOD FOR THE FUNCTION  
!     EVALUATION IS CHOSEN
!     SYMMETRY RELATIONS ARE NOW USED TO TRANSFORM THE FUNCTION  
!     BACK TO QUADRANT NN 
!
!     AS =ABRAMOWITZ AND STEGUN HANDBOOK OF MATHEMATICAL FUNCTIONS
!     Ol = Olver NIST HANDBOOK OF MATHEMATICAL FUNCTIONS at https://dlmf.nist.gov
C     ------------------------------------------------------------------  
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      REAL*8, intent(in) :: EMACH
      COMPLEX*16, intent(in) :: Z  
C  
C ..  LOCAL SCALARS  ..  
C  
      INTEGER    NN,N 
      REAL*8     ABSZ,ABTERM,API,EPS,FACT,FACTD,FACTN,PI  
      REAL*8     Q,RTPI,TEST,X,Y,YY  
      COMPLEX*16 ZZ,CONE,CI,CZERO,SUM,ZZS,XZZS,CER,CERF   
      COMPLEX*16 H1,H2,H3,U1,U2,U3,TERM1,TERM2  
C  
C ..  INTRINSIC FUNCTIONS  ..  
C  
*     INTRINSIC DCMPLX,EXP,DCONJG  
C  
C ..  DATA STATEMENTS  ..  
C  
      DATA PI/3.14159265358979D0/  
      DATA CONE/(1.D0,0.D0)/,CI/(0.D0,1.D0)/,CZERO/(0.D0,0.D0)/  
C     ------------------------------------------------------------------  
C  
      EPS=5.0D0*EMACH  
      API=1.0D0/PI 

! a switch for Z==0 
      IF(ABS(Z))2,1,2  
   1  CERF=CONE  
      GOTO 29  
!
!     THE ARGUMENT Z/=0 IS TRANSLATED FROM THE NN_TH QUADRANT INTO ZZ
!     IN THE FIRST QUADRANT, BEFORE THE METHOD FOR THE FUNCTION  
C     EVALUATION IS CHOSEN  
C  
   2  X=DBLE(Z)  
      Y=DIMAG(Z)  
      YY=Y  
      IF(Y)6,3,3  
   3  IF(X)5,4,4  
   4  ZZ=Z  
      NN=1  
      GOTO 9  
   5  ZZ=DCMPLX(-X,Y)  
      NN=2  
      GOTO 9  
   6  YY=-Y  
      IF(X)7,8,8  
   7  ZZ=-Z  
      NN=3  
      GOTO 9  
   8  ZZ=DCMPLX(X,-Y)  
      NN=4  
   9  ZZS=ZZ*ZZ  
      XZZS=EXP(-ZZS)  
      ABSZ=ABS(ZZ)  
      IF(ABSZ-10.0D0)10,10,23  
  10  IF(YY-1.0D0)11,12,12  
  11  IF(ABSZ-4.0D0)13,18,18  
  12  IF(ABSZ-1.0D0)13,18,18  
C  
! POWER SERIES (SEE 7.1.5 AS, p. 297; 7.6.1 Ol) 
! In the series, n-th term for -iz is a z**2*(2n-1)/(n*(2n+1))
! multiple of the (n-1)-th term
! The overall prefactor of the series is -2*i/\sqrt{\pi} for (-iz)
! argument
!
! One wonders why not the power series 7.1.8 AS, p. 297; 7.6.3 Ol
! is being used?? The latter would lead directly 
! to the Faddeeva error function w(z)??? TODO
C  
  13  Q=1.0D0          !provides n!
      FACTN=-1.0D0     !provides (2n-1) 
      FACTD=1.0D0      !provides (2n+1) 
      TERM1=ZZ         !==z**2  ... the first term of the sum for n=0
      SUM=ZZ           !==z**2  ... the first term of the sum for n=0 

  14  DO N=1,5  
      FACTN=FACTN+2.0D0  
      FACTD=FACTD+2.0D0  
      FACT=FACTN/(Q*FACTD)  !=(2n-1)/(n*(2n+1))

      TERM1=FACT*ZZS*TERM1  !generates nth term from the (n-1)th term
      SUM=SUM+TERM1         !summation

      Q=Q+1.0D0
      enddo
  
      ABTERM=ABS(TERM1)  
      IF(ABTERM-EPS)17,16,16  
  16  IF(Q-100.0D0)14,17,17  

  17  FACT=2.0D0*SQRT(API)  
      SUM=FACT*CI*SUM       !adds the missing prefactor 2*i/\sqrt{\pi} for 
                            !the (-iz) argument to provide -erf(-iz)
      CER=XZZS+XZZS*SUM     !generates Faddeeva error function w(z)
                            !according to 7.1.3 AS; 7.2.3 Ol
      GOTO 24  
C  
C     CONTINUED FRACTION THEORY (W(Z) IS RELATED TO THE LIMITING  
C     VALUE OF U(N,Z)/H(N,Z), WHERE U AND H OBEY THE SAME  
C     RECURRENCE RELATION IN N. SEE FADDEEVA AND TERENTIEV  
C     (TABLES OF VALUES OF W(Z) FOR COMPLEX ARGUMENTS, PERGAMON  
C     N.Y. 1961)  
!
!  7.9 Ol
C  
  18  TERM2=DCMPLX(1.D6,0.0D0)  
      Q=1.0D0  
      H1=CONE  
      H2=2.0D0*ZZ  
      U1=CZERO  
      RTPI=2.0D0*SQRT(PI)  
      U2=DCMPLX(RTPI,0.0D0)  
  19  TERM1=TERM2  
      DO 20 N=1,5  
      H3=H2*ZZ-Q*H1  
      U3=U2*ZZ-Q*U1  
      H1=H2  
      H2=2.0D0*H3  
      U1=U2  
      U2=2.0D0*U3  
  20  Q=Q+1.0D0  
      TERM2=U3/H3  
      TEST=ABS((TERM2-TERM1)/TERM1)  
      IF(TEST-EPS)22,21,21  
  21  IF(Q-60.0D0)19,19,13  
  22  CER=API*CI*TERM2  
      GOTO 24  
C 
!
!      7.12 Ol 
C     ASYMPTOTIC SERIES: SEE ABRAMOWITZ AND STEGUN, P328  
C  
  23  CER=0.5124242D0/(ZZS-0.2752551D0)+0.05176536D0/(ZZS-2.724745D0)  
      CER=CI*ZZ*CER  
C  
C     SYMMETRY RELATIONS ARE NOW USED TO TRANSFORM THE FUNCTION  
C     BACK TO QUADRANT NN  
C  
  24  GOTO(28,26,27,25),NN  
  25  CER=2.0D0*XZZS-CER  
  26  CERF=DCONJG(CER)  
      GOTO 29  
  27  CERF=2.0D0*XZZS-CER  
      GOTO 29  
  28  CERF=CER  
  29  RETURN  
      END  