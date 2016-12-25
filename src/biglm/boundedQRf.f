

      SUBROUTINE INCLUD(NP, NRBAR, WEIGHT, XROW, YELEM, D,
     +      RBAR, THETAB, SSERR, IER)
C
C     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2
C     Modified from algorithm AS 75.1
C
C     Calling this routine updates d, rbar, thetab and sserr by the
C     inclusion of xrow, yelem with the specified weight.   The number
C     of columns (variables) may exceed the number of rows (cases).
C
C**** WARNING: The elements of XROW are overwritten  ****
C
      INTEGER NP, NRBAR, IER
      DOUBLE PRECISION WEIGHT, XROW(NP), YELEM, D(NP), RBAR(*),
     +    THETAB(NP), SSERR
C
C     Local variables
C
      INTEGER I, K, NEXTR
      DOUBLE PRECISION ZERO, W, Y, XI, DI, WXI, DPI, CBAR, SBAR, XK
C
      DATA ZERO/0.D0/
C
C     Some checks.
C
      IER = 0
      IF (NP .LT. 1) IER = 1
      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 2
      IF (IER .NE. 0) RETURN
C
      W = WEIGHT
      Y = YELEM
      NEXTR = 1
      DO 30 I = 1, NP
C
C     Skip unnecessary transformations.   Test on exact zeroes must be
C     used or stability can be destroyed.
C
      IF (W .EQ. ZERO) RETURN
       XI = XROW(I)
      IF (XI .EQ. ZERO) THEN
       NEXTR = NEXTR + NP - I
      GO TO 30
      END IF
      DI = D(I)
      WXI = W * XI
      DPI = DI + WXI*XI
      CBAR = DI / DPI
      SBAR = WXI / DPI
      W = CBAR * W
      D(I) = DPI
      IF (I .EQ. NP) GO TO 20
      DO 10 K = I+1, NP
      XK = XROW(K)
      XROW(K) = XK - XI * RBAR(NEXTR)
      RBAR(NEXTR) = CBAR * RBAR(NEXTR) + SBAR * XK
      NEXTR = NEXTR + 1
   10   CONTINUE
   20   XK = Y
       Y = XK - XI * THETAB(I)
       THETAB(I) = CBAR * THETAB(I) + SBAR * XK
   30  CONTINUE
C
C     Y * SQRT(W) is now equal to Brown & Durbin's recursive residual.
C
      SSERR = SSERR + W * Y * Y
C
      RETURN
      END
C
C
      SUBROUTINE TOLSET(NP, NRBAR, D, RBAR, TOL, WORK, IER)
C
C     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2
C
C     Sets up array TOL for testing for zeroes in an orthogonal
C     reduction formed using AS75.1.
C
      INTEGER NP, NRBAR, IER
      DOUBLE PRECISION D(NP), RBAR(*), TOL(NP), WORK(NP)
C
C     Local variables.
C
      INTEGER COL, ROW, POS
      DOUBLE PRECISION EPS, SUM, ZERO
C
C     EPS is a machine-dependent constant.   For compilers which use
C     the IEEE format for floating-point numbers, recommended values
C     are 1.E-06 for single precision and 1.D-12 for double precision.
C
c     changed EPS from 10^-12 to 5x10^-10 to try to fix a bug
      DATA EPS/1.D-12/, ZERO/0.D0/
C
C     Some checks.
C
      IER = 0
      IF (NP .LT. 1) IER = 1
      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 2
      IF (IER .NE. 0) RETURN
C
C     Set TOL(I) = sum of absolute values in column I of RBAR after
C     scaling each element by the square root of its row multiplier.
C
      DO 10 ROW = 1, NP
      WORK(ROW) = SQRT(D(ROW))
  10  END DO
      DO 30 COL = 1, NP
      POS = COL - 1
      IF (COL .LE. NP) THEN
      SUM = WORK(COL)
      ELSE
      SUM = ZERO
      END IF
      DO 20 ROW = 1, MIN(COL-1, NP)
      SUM = SUM + ABS(RBAR(POS)) * WORK(ROW)
      POS = POS + NP - ROW - 1
  20  CONTINUE
      TOL(COL) = EPS * SUM
  30  CONTINUE
C
      RETURN
      END

      SUBROUTINE SINGCHK(NP, NRBAR, D, RBAR, THETAB, SSERR, TOL,
     +   LINDEP, WORK, IER)
C
C     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2
C
C     Checks for singularities, reports, and adjusts orthogonal
C     reductions produced by AS75.1.
C
      INTEGER NP, NRBAR, IER
      DOUBLE PRECISION D(NP), RBAR(NRBAR), THETAB(NP), SSERR,
     +      TOL(NP), WORK(NP)
      LOGICAL LINDEP(NP)
C
C     Local variables
C
      DOUBLE PRECISION ZERO, TEMP
      INTEGER COL, POS, ROW, NC2, POS2
C
      DATA ZERO/0.D0/
C
C     Check input parameters
C
      IER = 0
      IF (NP .LT. 1) IER = 1
      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 2
      IF (IER .NE. 0) RETURN
C
      DO 10 COL = 1, NP
      WORK(COL) = SQRT(D(COL))
   10 END DO
C
      DO 40 COL = 1, NP
C
C     Set elements within RBAR to zero if they are less than TOL(COL) in
C     absolute value after being scaled by the square root of their row
C     multiplier.
C
      TEMP = TOL(COL)
      POS = COL - 1
      DO 30 ROW = 1, COL-1
      IF (ABS(RBAR(POS)) * WORK(ROW) .LT. TEMP) RBAR(POS) = ZERO
      POS = POS + NP - ROW - 1
   30 CONTINUE
C
C     If diagonal element is near zero, set it to zero, set appropriate
C     element of LINDEP, and use INCLUD to augment the projections in
C     the lower rows of the orthogonalization.
C
      LINDEP(COL) = .FALSE.
      IF (WORK(COL) .LE. TEMP) THEN
      LINDEP(COL) = .TRUE.
      IER = IER - 1
      IF (COL .LT. NP) THEN
      NC2 = NP - COL
      POS2 = POS + NP - COL + 1
        IF (NC2 .GT. 1) THEN
       CALL INCLUD(NC2, NC2*(NC2-1)/2, D(COL), RBAR(POS+1),
     +            THETAB(COL), D(COL+1), RBAR(POS2), THETAB(COL+1),
     +            SSERR, IER)
        ELSE
          CALL INCLUD(1, 0, D(COL), RBAR(POS+1),
     +            THETAB(COL), D(COL+1), RBAR(1), THETAB(COL+1),
     +            SSERR, IER)
        END IF
      ELSE
      SSERR = SSERR + D(COL) * THETAB(COL)**2
      END IF
      D(COL) = ZERO
      WORK(COL) = ZERO
      THETAB(COL) = ZERO
      END IF
   40 CONTINUE
      RETURN
      END


      SUBROUTINE REGCF(NP, NRBAR, D, RBAR, THETAB, TOL, BETA,
     +     NREQ, IER)
C
C     ALGORITHM AS274  APPL. STATIST. (1992) VOL 41, NO. x
C
C     Modified version of AS75.4 to calculate regression coefficients
C     for the first NREQ variables, given an orthogonal reduction from
C     AS75.1.
C
      INTEGER NP, NRBAR, NREQ, IER
      DOUBLE PRECISION D(NP), RBAR(*), THETAB(NP), TOL(NP),
     +     BETA(NP)
C
C     Local variables
C
      INTEGER I, J, NEXTR
      DOUBLE PRECISION ZERO
C
      DATA ZERO/0.D0/
C
C     Some checks.
C
      IER = 0
      IF (NP .LT. 1) IER = 1
      IF (NRBAR .LT. NP*(NP-1)/2) IER = IER + 2
      IF (NREQ .LT. 1 .OR. NREQ .GT. NP) IER = IER + 4
      IF (IER .NE. 0) RETURN
C
      DO 20 I = NREQ, 1, -1
      IF (SQRT(D(I)) .LT. TOL(I)) THEN
      BETA(I) = ZERO
      D(I) = ZERO
      GO TO 20
      END IF
      BETA(I) = THETAB(I)
      NEXTR = (I-1) * (NP+NP-I)/2 + 1
      DO 10 J = I+1, NREQ
      BETA(I) = BETA(I) - RBAR(NEXTR) * BETA(J)
      NEXTR = NEXTR + 1
   10 CONTINUE
   20 CONTINUE
C
      RETURN
      END
