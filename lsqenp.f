        SUBROUTINE LSQENP( NOBS, NPAR, NVAR, Y, X, B, NFIX, IFIX,
     *  IDVT, ICON, MAXITER, IPRNT,
     *  FCODE, PCODE )
*
*  Performs weighted least-squares fit of models to data
*
*  Input:
*       NOBS    = Number of observations
*       NPAR    = Number of parameters
*       NVAR    = Number of independent variables
*       Y       = Dependent variable
*       X       = Independent variables
*       B       = Parameter values
*       NFIX    = Number of fixed parameters
*       IFIX    = Index of fixed parameters
*       IDVT    = 0 for analytic derivatives via PCODE
*                 1 for numerical derivatives
*       ICON    = 0
*       MAXITER = Maximum allowed number of iterations
*       IPRNT   = Print option
*       FCODE   = External routine to compute predicted data
*       PCODE   = External routine to compute parameter derivatives
*  Output:
*       B       = Revised parameter values
*
*  See also common blocks for confidence limits
*
*  Origin: Central computing centre at Caltech?
*
C REVISED 6-12-70 TO CALCULATE SUM(( F -Y(I))**2) IN DBLE PRECISION
C REVISED 11-20-72 OPTION TO OMIT PRINTING OBSERVED AND CALCULATED
C      VALUES AFTER CONVERSION.  SET IPRNT LESS THAN 0.
C REVISED 9-18-73
C LSQENP AND CONLIM REVISED FEB., 1975
C LSQENP REVISED JUNE 24, 1976
* NOV 1984 Common blocks for weights, confidence limits  KDH at IOA.
* JUN 1985 FCODE and PCODE made external arguments  KDH.
* Mar 1987 KDH @ STScI Added common/lsqiout/ for output unit number
*
        PARAMETER (MAXPARM = 50)  ! Maximum allowed number of parameters
        PARAMETER (MAXDATA = 500) ! Maximum allowed number of data points
        REAL*4 X(MAXDATA, 1), Y(NOBS), B(NPAR)
        CHARACTER*1 IBCH, IOCH, IPCH, IXCH, IYCH, IP1, IP2
        INTEGER*4 IFIX(NFIX)
        REAL*4 BS(MAXPARM), DB(MAXPARM), BA(MAXPARM), G(MAXPARM)
        REAL*4 SA(MAXPARM), P(MAXPARM), A(MAXPARM, MAXPARM + 1)
        REAL*8 PHD, XLL, DTG, GTG
        EXTERNAL FCODE, PCODE

        COMMON/LSQIOUT/ IOUT    ! Output unit number
        COMMON /LSQPLT/ IFP, YMIN, YMAX
        COMMON /PARAM/FF, T, E, TAU, AL, GAMCR, DEL, ZETA
        COMMON /WEITS/WTS(MAXDATA)
        COMMON /AFCLSQ/ SE, PHI, PHIZ, WS, XL, IFSS2, IFSS3, IWS6,
     *  IOBS, JPAR, JJ, IBK2, IBKA, IBKM, IPLOT, BS, SA, A

        DATA    IFP,    YMIN,   YMAX
     *  /       0,      0.,     0.      /
        DATA
     *  FF,     T,      E,      TAU,    AL,     GAMCR,  DEL,    ZETA
     *  /
     *  4.,     2.,     1.E-6,  5.E-6,  1.E-3,  45.,    1.E-5,  1.E-30
     *  /
        DATA WTS/MAXDATA*1.0/

        IWS4 = MAXITER
        IWS6 = ICON
        IPLOT = IFP
        IF(IFP .EQ. 1) SPRD = YMAX - YMIN
C     ------------------------------------------------------------------
C                       IWHER = 1 MEANS GET P S AND F
C                       IWHER GREATER THAN 1 MEANS GET F ONLY
   23   IWHER = 0
        GOTO 43

   31   CALL FCODE(Y, X, B, F, IOBS)
        IF(IWHER.NE.1) GOTO 39
   34   IF(IFSS2.NE.0) GOTO 39

   35   CALL PCODE(P, X, B, F, IOBS)
   39   IF(IWHER.EQ.0) GOTO 43
C                1  2   3   4
   42   GOTO(106, 363, 114, 125), IWHER
   43   ITCT = 0
        IF(IPLOT .LE. 0) GOTO 55
   49   IBCH = ' '
        IOCH = 'O'
        IPCH = 'P'
        IXCH = 'X'
        IYCH = 'Y'
   55 IF(IPRNT.GT.0) THEN
        DO I = 1, NFIX
          IF(IFIX(I).LE.0) THEN
            WRITE(IOUT, 402)
  402       FORMAT( ' BAD DATA, SUBSCRIPTS FOR UNUSED BS = 0', 
     &      /// )
            STOP
          END IF
        END DO
      END IF
        XKDB = 1.
C     ..................................................................
C                       START THE CALCULATION OF THE PTP MATRIX
        IBKA = 1
        IF( IPRNT.GT.0 .AND. IOUT.GT.0 ) WRITE(IOUT, 383)
     *  NOBS, NPAR, NFIX, NVAR, IFP,
     *  GAMCR, DEL, FF, T, E, TAU, AL, ZETA
  383   FORMAT( /' NOBS =', I3, ' NPAR =', I3, ' NFIX =', I3,
     *  ' NVAR =', I3, ' IFP =', I3, 'GAMMA CRIT =', E10.3,
     *  'DEL =', E10.3,
     *  /, ' FF = ', E10.3, ' T =', E10.3, ' E =', E10.3,
     *  ' TAU =', E10.3, ' XL =', E10.3, ' ZETA =', E10.3, /)

   82 DO I = 1, NPAR
        G(I) = 0.
      DO J = 1, NPAR
        A(I, J) = 0.
      END DO
      END DO
      IF( IBKA.EQ.1 ) THEN
        IFSS3 = IPRNT
        IFSS2 = IDVT
      ELSE
        IFSS3 = 1
      END IF
   92   IF(IFSS3 .GT. 0.AND. IOUT.GT.0)
     #          WRITE(IOUT, 384) (B(J), J = 1, NPAR)
  384   FORMAT( /, ' PARAMETERS ', 5E18.8, /(12X, 5E18.8) )
        IF(IFSS3 .LE. 0 .OR. IPRNT.LT. 0) GOTO 99

      IF(IPLOT .LE. 0) THEN
        IF( IOUT.GT.0 ) WRITE(IOUT, 386)
  386   FORMAT( //T12, 'X(I,1)', T31, 'Y OBS.',
     *  T49, 'Y PRED.', T68, 'DIFF', / )
      ELSE
        WS = YMIN + SPRD
        IF( IOUT.GT.0 ) WRITE(IOUT, 382) YMIN, WS
  382   FORMAT( 1X, 1E9.2, 86X, 1E9.2, /1X, ' + ',99X, ' + ' )
      END IF
   99   IOBS = 1
        PHD = 0.
        IF(IFSS2.EQ.0) GOTO 104
        GOTO 112

  103   IF(IFSS2.EQ.1) GOTO 112
  104   IWHER = 1
C                       GET P S AND F
        GOTO 31

  106 IF(NFIX.LE.0) GOTO 132
  107 DO I = 1, NFIX
        IWS = IFIX(I)
        P(IWS) = 0.
      END DO
        GOTO 132
C     ..................................................................
C                       THIS IS THE ESTIMATED P S ROUTINE
  112   IWHER = 3
        GOTO 31
  114   FWS = F
        JPAR = 1
  116   IF(NFIX.LE.0) GOTO 120
  117 DO I = 1, NFIX
        IF(JPAR.EQ.IFIX(I)) GOTO 128
      END DO
  120   DBW = B(J)*DEL
        IF(B(JPAR) .EQ. 0.) DBW = DEL
        TWS = B(JPAR)
        B(JPAR) = B(JPAR) + DBW
        PRINT *, 'NUMERICAL DERIVATIVE FOR PARAMETER', JPAR
        IWHER = 4
        GOTO 31
  125   B(JPAR) = TWS
        P(JPAR) = (F - FWS)/DBW
        GOTO 129
  128   P(JPAR) = 0.
  129   JPAR = JPAR + 1
        IF(JPAR.EQ.NPAR) GOTO 116
  131   F = FWS
C               END OF ESTIMATED P S ROUTINE
C     ..................................................................
C               NOW, USE THE P S TO MAKE PARTIALS MATRIX
  132 DO J = 1, NPAR
c       G(J) = G(J) + (Y(IOBS) - F)*P(J) ! jhw bug fixed 29 Apr 92
        G(J) = G(J) + wts(iobs)*(Y(IOBS) - F)*P(J)
      DO I = J, NPAR
        A(I, J) = A(I, J) + P(I)*P(J)
        A(J, I) = A(I, J)
      END DO
      END DO
        IF(IPLOT .LE. 0) GOTO 184
  138   IF(IFSS3.LE.0) GOTO 188
C                       PLOTTING Y(IOBS), F
  139   IO = (Y(IOBS) - YMIN)*100./SPRD
        IPP = (F - YMIN)*100./SPRD
        IF(IO.EQ.IPP) GOTO 148
        IF(IO.GT.IPP) GOTO 153
C                       Y(IOBS) OUT FIRST
  143   IP1 = IOCH
        IP2 = IPCH
        I1 = IO
        I2 = IPP
        GOTO 157
C                       ONLY ONE CHARACTER
  148   IP1 = IYCH
        IP2 = IBCH
        I1 = IO
        I2 = IPP
        GOTO 157
C                       F OUT FIRST
  153   IP1 = IPCH
        IP2 = IOCH
        I1 = IPP
        I2 = IO
C                       ZERO PLOTS IN THE LEFT HAND COLUMN, SO I1 IS ITS
C                       OWN BLANK COUNTER
C                       OVERFLOWS PLOT X IN COLUMN 102
C                       UNDERFLOWS ALSO PLOT X IN COLUMN ZERO
  157   IF(I2.LE.101) GOTO 165
  158   I2 = 101
        IP2 = IXCH
        IF(I1.LT.101) GOTO 165
  161   I1 = 101
        IP1 = IXCH
        IP2 = IBCH
        GOTO 171
  165   IF(I1.GE.0) GOTO 171
  166   I1 = 0
        IP1 = IXCH
        IF(I2.GT.0) GOTO 171
  169   I2 = 1
        IP2 = IBCH
  171   I1M1 = I1
        I1M2 = I2 - I1 - 1
        IF(I1M1.GT.0) GOTO 179
  174   IF(I1M2.GT.0) GOTO 177
  175   IF( IOUT.GT.0 ) WRITE(IOUT, 404) IP1, IP2
        GOTO 188
  177   IF( IOUT.GT.0 ) WRITE(IOUT, 404) IP1, (IBCH, I = 1, I1M2), IP2
        GOTO 188
  179   IF(I1M2.GT.0) GOTO 182
  180   IF( IOUT.GT.0 ) WRITE(IOUT, 404) (IBCH, I = 1, I1M1), IP1, IP2
        GOTO 188
  182   IF( IOUT.GT.0 ) WRITE(IOUT, 404) (IBCH, I = 1, I1M1), IP1,
     *                  (IBCH, I = 1, I1M2), IP2
  404   FORMAT( ' ', 110A1  )
        GOTO 188

  184   WS = Y(IOBS) - F
        IF(IFSS3 .LE. 0 .OR. IPRNT.LT. 0) GOTO 188
  187   IF( IOUT.GT.0 ) WRITE(IOUT, 401) X(IOBS, 1), Y(IOBS), F, WS
  401   FORMAT( 5X, 6E18.8, /59X, 2E18.8)
  188   WS = WTS(IOBS)*(Y(IOBS) - F)
        PHD = PHD + WS*WS*1.0D0
        IOBS = IOBS + 1
        IF(IOBS.LE.NOBS) GOTO 103
        PHI = PHD

      IF(NFIX.GT.0) THEN
  193 DO J = 1, NFIX
        IWS = IFIX(J)
      DO I = 1, NPAR
        A(IWS, I) = 0.
        A(I, IWS) = 0.
      END DO
        A(IWS, IWS) = 1.
      END DO
      END IF
*  insert 1
  199 IF(IBKA .NE. 1) THEN
        IBKS = 1
        CALL CONLIM( NOBS, NPAR, NVAR, Y, X, B, NFIX, IFIX, IBKS, IBD)
        GOTO(82, 314, 38, 359, 38), IBD
C                       SAVE SQUARE ROOTS OF DIAGONAL ELEMENTS
      END IF
      DO I = 1, NPAR
        SA(I) = SQRT (A(I, I))
      END DO
      DO I = 1, NPAR
      DO J = 1, NPAR
        WS = SA(I)*SA(J)
      IF(WS.GT.0.) THEN
        A(I, J) = A(I, J)/WS
      ELSE
        A(I, J) = 0.
      END IF
      END DO
      IF(SA(I).GT.0.) THEN
        G(I) = G(I)/SA(I)
      ELSE
        G(I) = 0.
      END IF
      END DO
      DO I = 1, NPAR
        A(I, I) = 1.
      END DO
        PHIZ = PHI
C                       WE NOW HAVE PHI ZERO
C     ..................................................................
        IF(ITCT.GT.0) GOTO 230
C                       FIRST ITERATION
  225   XL = AL
        ITCT = 1
      DO J = 1, NPAR
        BS(J) = B(J)
      END DO
C                       BS(J) CORRESPONDS TO PHIZ
  230   IBK1 = 1
        WS = NOBS - NPAR + NFIX
        SE = SQRT(PHIZ/WS)
        IF(IPRNT .LE. 0) GOTO 310
        IF(IFSS3.GT.0) GOTO 239
        IF(IFSS2.EQ.0) GOTO 237
  235   IF( IOUT.GT.0 ) WRITE(IOUT, 387) PHIZ, SE, XLL, GAMMA, XL
  387   FORMAT( /13X,' PHI',14X,' S E',11X,' LENGTH', 6X,' GAMMA ',6X,
     1  ' LAMBDA', 6X, 'ESTIMATED PARTIALS USED  ',/5X,2E18.8,3E13.3)
        GOTO 310
  237   IF( IOUT.GT.0 ) WRITE(IOUT, 388) PHIZ, SE, XLL, GAMMA, XL
  388   FORMAT( /13X,' PHI',14X,' S E',11X,' LENGTH',6X,' GAMMA ',6X,
     1  ' LAMBDA', 6X,'ANALYTIC PARTIALS USED  ',/5X,2E18.8,3E13.3)
        GOTO 310
  239   IF(IFSS2.EQ.0) GOTO 242
  240   IF( IOUT.GT.0 ) WRITE(IOUT, 379) PHIZ, SE, XL
  379   FORMAT( /13X, ' PHI' ,14X, ' S E', 9X, ' LAMBDA', 6X,
     1  ' ESTIMATED PARTIALS USED ', /5X, 2E18.8, E13.3 )
        GOTO 310
  242   IF( IOUT.GT.0 ) WRITE(IOUT, 385) PHIZ, SE, XL
  385   FORMAT( /13X, ' PHI', 14X, ' S E', 9X, ' LAMBDA', 6X,
     1  ' ANALYTIC PARTIALS USED  ', /5X, 2E18.8, E13.3 )
        GOTO 310
  244   PHIL = PHI
C                       WE NOW HAVE PHI LAMBDA
      DO J = 1, NPAR
        IF(ABS(DB(J)/(ABS(B(J)) + TAU)).GE.E) GOTO 251
      END DO
        IF( IOUT.GT.0 ) WRITE(IOUT, '(A)' ) ' EPSILON TEST '
        IBS = 4
        GOTO 371
  251   IF(IWS4.EQ.0) GOTO 257
        IF(IWS4.EQ.1) GOTO 255
        IWS4 = IWS4 - 1
        GOTO 257
  255   IF( IOUT.GT.0 ) WRITE(IOUT, '(A)' ) ' FORCE OFF'
        GOTO 371
  257   XKDB = 1.
        IF(PHIL.GT.PHIZ) GOTO 281
  259   XLS = XL
      DO 262 J = 1, NPAR
        BA(J) = B(J)
  262   B(J) = BS(J)
        IF(XL.GT..00000001) GOTO 268
  264 DO 266 J = 1, NPAR
        B(J) = BA(J)
  266   BS(J) = B(J)
        GOTO 82
  268   XL = XL/10.
        IBK1 = 2
        GOTO 310
  271   PHL4 = PHI
C                       WE NOW HAVE PHI(LAMBDA/10)
        IF(PHL4.GT.PHIZ) GOTO 276
  273 DO 274 J = 1, NPAR
  274   BS(J) = B(J)
        GOTO 82
  276   XL = XLS
      DO 279 J = 1, NPAR
        BS(J) = BA(J)
  279   B(J) = BA(J)
        GOTO 82
  281   IBK1 = 4
        XLS = XL
        XL = XL/10.
      DO 285 J = 1, NPAR
  285   B(J) = BS(J)
        GOTO 310
  287   IF(PHI.LE.PHIZ) GOTO 296
  288   XL = XLS
        IBK1 = 3
  290   XL = XL*10.
  291 DO 292 J = 1, NPAR
  292   B(J) = BS(J)
        GOTO 310
  294   PHIT4 = PHI
C                       WE NOW HAVE PHI(10*LAMBDA)
        IF(PHIT4.GT.PHIZ) GOTO 299
  296 DO 297 J = 1, NPAR
  297   BS(J) = B(J)
        GOTO 82
  299   IF(GAMMA.GE.GAMCR) GOTO 290
  300   XKDB = XKDB/2.
      DO J = 1, NPAR
        IF(ABS(DB(J)/(ABS(B(J)) + TAU)).GE.E) GOTO 291
      END DO
  305   B(J) = BS(J)
        IF( IOUT.GT.0 ) WRITE(IOUT, '(A)' ) ' GAMMA EPSILON TEST'
        IBS = 4
        GOTO 371
C
C     ..................................................................
C                       SET UP FOR MATRIX INVERSION
  310 DO I = 1, NPAR
        A(I, I) = A(I, I) + XL
      END DO
C                       GET INVERSE OF A AND SOLVE FOR DB (J)S
        IBKM = 1
C     ..................................................................
C                       THIS IS THE MATRIX INVERSION ROUTINE
C                       NPAR IS THE SIZE OF THE MATRIX
  314   CALL GJR(A, NPAR, ZETA, MSING)
        GOTO(316, 38), MSING
C     INSERT 2
  316   IF(IBKM .EQ. 1) GOTO 321
        IBKS = 2
        CALL CONLIM( NOBS, NPAR, NVAR, Y, X, B, NFIX, IFIX, IBKS, IBD)
        GOTO( 82, 314, 38, 359, 38), IBD
C               END OF MATRIX INVERSION, SOLVE FOR DB(J)

*  Apply correction to parameters

  321 DO 325 I = 1, NPAR
        DB(I) = 0.
      if(nfix.gt.0) then
      do J=1,nfix
        if(ifix(J).eq.I) goto 325
      end do
      end if
      DO J = 1, NPAR
        DB(I) = A(I, J) * G(J) + DB(I)
      END DO
        DB(I) = XKDB * DB(I)
  325 END DO
        XLL = 0.0
        DTG = 0.
        GTG = 0.
      DO J = 1, NPAR
        XLL = XLL + DB(J)*DB(J)
        DTG = DTG + DB(J)*G(J)
        GTG = GTG + G(J)*G(J)
        IF( SA(J).GT.0. ) DB(J) = DB(J) / SA(J)
        B(J) = B(J) + DB(J)
      END DO

        KIP = NPAR - NFIX
        IF(KIP.EQ.1) GOTO 350
      if( xll*gtg .gt.0. ) then
        CGAM = DTG/DSQRT(XLL*GTG)
      else
        print *, char(7), '** LSQENP FLOATING DIVIDE BY 0 AVERTED.'
        print *, ' GTG = ', gtg, ' XLL = ', xll
        CGAM = DTG
      end if
        JGAM = 1
        IF(CGAM.GT..0) GOTO 342
  340   CGAM = ABS(CGAM)
        JGAM = 2
  342 IF(CGAM .GT. 1.0) CGAM = 1.0
        GAMMA = 57.2957795*(1.5707288 + CGAM*(
     *  - 0.2121144 + CGAM*(0.074261 - CGAM*.0187293)))
     *  *SQRT(1. - CGAM)
        GOTO(351, 344), JGAM
  344   GAMMA = 180. - GAMMA
        IF(XL.LT.1.0) GOTO 351
  346   IF( IOUT.GT.0 ) WRITE(IOUT, '(A,5X,2E13.3)' )
     #          ' GAMMA LAMBDA TEST', XL, GAMMA
        IBS = 4
        GOTO 371
  350   GAMMA = 0.
  351   XLL = DSQRT(XLL)
        IBK2 = 1
        GOTO 359

  354   IF(IFSS3.LE.0) GOTO 358
        IF( IOUT.GT.0 ) WRITE(IOUT, 380) (DB(J), J = 1, NPAR)
  380   FORMAT( /,' INCREMENTS ', 5E18.8, /(12X, 5E18.8) )
        IF( IOUT.GT.0 ) WRITE(IOUT, 381) PHI, XL, GAMMA, XLL
  381   FORMAT( 13X, ' PHI', 10X, ' LAMBDA', 6X, ' GAMMA ', 
     &  6X, ' LENGTH', /, 5X, E18.8, 3E13.3)
  358   GOTO(244, 271, 294, 287), IBK1
C
C     ..................................................................
C                       CALCULATE PHI
  359   IOBS = 1
        PHD = 0.
        IWHER = 2
        GOTO 31

  363   PHD = PHD + ((Y(IOBS) - F)*1.0D0*WTS(IOBS))**2
        IOBS = IOBS + 1
        IF(IOBS.LE.NOBS) GOTO 31
        PHI = PHD
C     INSERT 3
        IF(IBK2 .EQ. 1) GOTO 354
        IBKS = 3
        CALL CONLIM( NOBS, NPAR, NVAR, Y, X, B, NFIX, IFIX, IBKS, IBD)
        GOTO( 82, 314, 38, 359, 38), IBD

  371   IBKS = 4
        CALL CONLIM( NOBS, NPAR, NVAR, Y, X, B, NFIX, IFIX, IBKS, IBD)
        GOTO( 82, 314, 38, 359, 38), IBD

   38   RETURN
        END

        SUBROUTINE CONLIM(NOBS,NPAR,NVAR,Y,X,B,NFIX,IFIX,IBKS,IBD)
*
*  Computes and reports confidence limits
*
        PARAMETER (MAXPARM = 50)  ! Maximum allowed number of parameters
        PARAMETER (MAXDATA = 500) ! Maximum allowed number of data points
        PARAMETER (IOUT=6)        ! Output unit number
        REAL*4 BSAVE(MAXPARM), SA(MAXPARM), A(MAXPARM, MAXPARM + 1)
        REAL*4 X(MAXDATA, 1), Y(NOBS), B(NPAR)
        INTEGER*4 IFIX(NFIX)
        COMMON /PARAM/FF, T, E, TAU, AL, GAMCR, DEL, ZETA
        COMMON /AFCLSQ/ SE, PHI, PHIZ, WS, XL, IFSS2, IFSS3, IWS6,
     *  I, J, JJ, IBK2, IBKA, IBKM, IPLOT, BSAVE, SA, A
        COMMON/LSQSTE/ STE(MAXPARM)     ! Standard 1-sigma uncertainty
        COMMON/LSQOPL/ OPL(MAXPARM)     ! One-parameter lower limits
        COMMON/LSQOPU/ OPU(MAXPARM)     ! One-parameter upper limits
        COMMON/LSQSPL/ SPL(MAXPARM)     ! Support plane lower limits
        COMMON/LSQSPU/ SPU(MAXPARM)     ! Support plane upper limits
        COMMON/LSQBU/ BU(MAXPARM)       ! Nonlinear upper limits
        COMMON/LSQBL/ BL(MAXPARM)       ! Nonlinear lower limits

        INDEX=0
        IBKP=0
        DD=0.0
        IBKN=0
        D=0.0
        PHI1=0.0

C     TO INITIATE CONLIM
        IF(IBKS .EQ. 4) GOTO 21
        GOTO(15, 17, 18), IBKS
   15   IBKA1 = IBKA - 1
        GOTO(27, 32), IBKA1
   17   GOTO 43
   18   IBK21 = IBK2 - 1
        J = INDEX
        GOTO(158, 27, 125, 134, 144), IBK21
   21 DO J = 1, NPAR
        B(J) = BSAVE(J)
      END DO

       WRITE(IOUT, 201) NOBS, NPAR, NFIX, NVAR, FF, T, E, TAU
  201   FORMAT( /, ' NOBS =', I3, ' NPAR =', I3,
     *  ' NFIX =', I3, ' NVAR =', I3,
     *  /' FF =', E10.3, ' T =', E10.3,
     *  ' E =', E10.3, ' TAU =', E10.3, / )

        IBKA = 2
C                       THIS WILL PRINT THE Y, YHAT, DELTA Y
        IBD = 1
        GOTO 204
   27   IF(IPLOT .LE. 0) GOTO 32
   28   IBKA = 3
        IPLOT = 0
        IBD = 1
        GOTO 204
   32   WS = NOBS - NPAR + NFIX
        SE = SQRT(PHI/WS)
        PHIZ = PHI
        IF(IFSS2.EQ.0) GOTO 38
36      CONTINUE
       WRITE(IOUT, 189) PHIZ, SE, XL
  189   FORMAT( /13X, ' PHI', 14X, ' S E', 9X, ' LAMBDA', 6X,
     *  ' ESTIMATED PARTIALS USED ', /5X, 2E18.8, E13.3 )
        GOTO 39
38      CONTINUE
       WRITE(IOUT, 190) PHIZ, SE, XL
  190   FORMAT( /13X, ' PHI', 14X, ' S E', 9X, ' LAMBDA', 6X,
     *  ' ANALYTIC PARTIALS USED  ', /5X, 2E18.8, E13.3)
C                       NOW WE HAVE MATRIX A
   39   IBKM = 2
        IBD = 2
        GOTO 204
C                       NOW WE HAVE C = A INVERSE
   43 DO J = 1, NPAR
        IF(A(J, J).LT..0) GOTO 47
        SA(J) = SQRT(A(J, J))
      END DO
        IBOUT = 0
        GOTO 48
   47   IBOUT = 1
   48   KST = -4
       WRITE(IOUT, '(//A)' ) ' PTP INVERSE '
   50   KST = KST + 5
        KEND = KST + 4
        IF(KEND.LT.NPAR) GOTO 54
        KEND = NPAR
54      CONTINUE
      DO I = 1, NPAR
       WRITE(IOUT, '(I5,5E18.8)' ) I, (A(I, J), J = KST, KEND)
      END DO

        IF(KEND.LT.NPAR) GOTO 50
      IF( IBOUT.NE.0 ) THEN
        WRITE(IOUT, '(A)' ) ' NEGATIVE DIAGONAL ELEMENT '
        IBD = 3
        GOTO 204
      END IF

      DO J = 1, NPAR
      DO I = 1, NPAR
        WS = SA(I)*SA(J)
      IF(WS.GT. 0.) THEN
        A(I, J) = A(I, J)/WS
      ELSE
        A(I, J) = 0.
      END IF
      END DO
      END DO

      DO J = 1, NPAR
        A(J, J) = 1.
      END DO
       WRITE(IOUT, '(//A)' ) ' PARAMETER CORRELATION MATRIX '
        KST = -9
   73   KST = KST + 10
        KEND = KST + 9
        KEND = MAX(KEND,NPAR)
      DO I = 1, NPAR
       WRITE(IOUT, '(I8,2X,10F10.4)' ) I, (A(I, J), J = KST, KEND)
      END DO
        IF(KEND.LT.NPAR) GOTO 73
C                       GET T*SE*SQRT(C(I, I))
      DO J = 1, NPAR
        SA(J) =  SE*SA(J)
      END DO
*
*  Report Confidence limits based on quadratic form
*
       WRITE(IOUT, 197)
  197   FORMAT( ////13X, ' STD', 17X, ' ONE - PARAMETER', 21X,
     * ' SUPPORT PLANE', /3X, ' B', 7X, ' ERROR', 12X, ' LOWER', 12X,
     * ' UPPER', 12X, ' LOWER', 12X, ' UPPER' )
        WS = NPAR - NFIX
      DO 98 J = 1, NPAR
      IF( NFIX.GT.0 ) THEN
      DO I = 1, NFIX                    ! Is this parameter fixed ?
      IF(J.EQ.IFIX(I)) THEN
       WRITE(IOUT, '(I5,A)' ) J, ' PARAMETER NOT USED'
        GOTO 98
      END IF
      END DO
      END IF
        HJTD = SQRT(WS*FF)*SA(J)
        STE(J) = SA(J)
        OPL(J) = BSAVE(J) - SA(J)*T
        OPU(J) = BSAVE(J) + SA(J)*T
        SPL(J) = BSAVE(J) - HJTD
        SPU(J) = BSAVE(J) + HJTD
       WRITE(IOUT, '(I5,5E18.8)' )
     *         J, STE(J), OPL(J), OPU(J), SPL(J), SPU(J)
   98 END DO
*
*  Nonlinear confidence limits
*
        IF(IWS6.EQ.1) IBD = 3
        IF(IWS6.EQ.1) GOTO 204
        WS = NPAR - NFIX
        WS1 = NOBS - NPAR + NFIX
        PKN = WS/WS1
        PC = PHIZ*(1. + FF*PKN)

       WRITE(IOUT, 198) PC
  198   FORMAT( ////'  NONLINEAR CONFIDENCE LIMITS ',
     *  //' PHI CRITICAL =', E15.8 )
       WRITE(IOUT, 199)
  199   FORMAT( //'  PARA', 6X, ' LOWER B', 8X, ' LOWER PHI', 10X,
     *  ' UPPER B', 8X, ' UPPER PHI' )

        IFSS3 = 1
        J = 1
  109   IBKP = 1
      DO JJ = 1, NPAR
        B(JJ) = BSAVE(JJ)
      END DO
      IF(NFIX.GT.0) THEN
      DO JJ = 1, NFIX
        IF(J.EQ.IFIX(JJ)) GOTO 173
      END DO
      END IF
        DD = -1.
        IBKN = 1
  119   D = DD
        B(J) = BSAVE(J) + D*SA(J)
        IBK2 = 4
        IBD = 4
        INDEX = J
        GOTO 204

  125   PHI1 = PHI
        IF(PHI1.GE.PC) GOTO 137
  127   D = D + DD
        IF(D/DD.GE.5.) GOTO 177
  129   B(J) = BSAVE(J) + D*SA(J)
        IBK2 = 5
        IBD = 4
        INDEX = J
        GOTO 204

  134   PHID = PHI
        IF(PHID.LT.PC) GOTO 127
        GOTO 146
  137   D = D/2.
        IF(D/DD.LE..001) GOTO 177
  139   B(J) = BSAVE(J) + D*SA(J)
        IBK2 = 6
        IBD = 4
        INDEX = J
        GOTO 204

  144   PHID = PHI
        IF(PHID.GT.PC) GOTO 137
  146   XK1 = PHIZ/D + PHI1/(1. - D) + PHID/(D*(D - 1.))
        XK2 = -(PHIZ*(1. + D)/D + D/(1. - D)*PHI1 + PHID/(D*(D - 1.)))
        XK3 = PHIZ - PC
        BC = (SQRT(XK2*XK2 - 4.*XK1*XK3) - XK2)/(2.*XK1)
        GOTO(151, 153), IBKN
  151   B(J) = BSAVE(J) - SA(J)*BC
        GOTO 154
  153   B(J) = BSAVE(J) + SA(J)*BC
  154   IBK2 = 2
        IBD = 4
        INDEX = J
        GOTO 204
  158   GOTO(159, 164), IBKN
  159   IBKN = 2
        DD = 1.
        BL(J) = B(J)
        PL = PHI
        GOTO 119
  164   BU(J) = B(J)
        PU = PHI
        GOTO(167, 169, 171, 175), IBKP
 167    CONTINUE
       WRITE(IOUT, '(I5,5E18.8)') J, BL(J), PL, BU(J), PU
        GOTO 185
169     CONTINUE
       WRITE(IOUT, '(I5,36X,2E18.8)' ) J, BU(J), PU
        GOTO 185
171     CONTINUE
       WRITE(IOUT, '(I5,5E18.8)' ) J, BL(J), PL
        GOTO 185
173     CONTINUE
       WRITE(IOUT, '(I5,A)' ) J, ' PARAMETER NOT USED '
        GOTO 185
175     CONTINUE
       WRITE(IOUT, '(I5,A)' ) J, ' NONE FOUND '
        GOTO 185

  177   GOTO(178, 180), IBKN
C                       DELETE LOWER PRINT
  178   IBKP = 2
        GOTO 158
  180   GOTO(181, 183), IBKP
C                       DELETE UPPER PRINT
  181   IBKP = 3
        GOTO 158
C                       LOWER IS ALREADY DELETED, SO DELETE BOTH
  183   IBKP = 4
        GOTO 158
  185   J = J + 1
        J1 = J - 1
        IF(J1 .NE. NPAR) GOTO 109
      DO JJ = 1, NPAR
        B(JJ) = BSAVE(JJ)
      END DO
        IBD = 5
  204   RETURN
        END
