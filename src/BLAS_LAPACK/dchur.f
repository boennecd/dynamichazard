*  Changes by Benjamin Christoffersen
*  ==================================
*
*  The CHARACTER UPLO and TRANS arguments is replaced by
*  LOGICAL UPPER and DOTRAN so LSAME is not needed
*
      SUBROUTINE DCHUR(UPLO,TRANS,N,M,R,LDR,X,Z,LDZ,Y,RHO,C,S,INFO)
*     .. Scalar Arguments ..
      CHARACTER UPLO,TRANS
*      LOGICAL UPPER,DOTRAN
      INTEGER N,M,LDR,LDZ,INFO
*     .. Array Arguments ..
      DOUBLE PRECISION R(LDR,*)
      DOUBLE PRECISION X(*)
      DOUBLE PRECISION Z(LDZ,*)
      DOUBLE PRECISION Y(*)
      DOUBLE PRECISION RHO(*)
      DOUBLE PRECISION C(*)
      DOUBLE PRECISION S(*)
*
*  Purpose
*  =======
*
*  DCHUR performs a rank 1 update of a Cholesky factor R and
*  optionnally, of associated vectors Z
*
*  If UPLO = 'U', this solves the problems
*      R2'*R2 = R1'*R1 - x*x'
*      R2'*op(Z2) = R1'*op(Z1) - x*y'
*  Or if UPLO='L', the transposed versions
*      R2*R2'=R1*R1' - x*x'
*      R2*op(Z2)=R1*op(Z1) - x*y'
*  where R1 is R on entry, R2 is R on exit, ' means transpose,
*  Z1 is Z on entry, Z2 is Z on exit, op is transpose or identity.
*
*  The problem can also be written if UPLO='U' as finding U such that
*
*                              (R  op(Z))     (R2  op(Z2))
*                         U  * (        )  =  (          ) ,
*                              (X      Y)     ( 0     Y2 )
*
*  Or its transposed form if UPLO='L'
*
*                         (R      X)          (R2      0 )
*                         (        ) * U'  =  (          ) ,
*                         (op(Z)' Y)          (op(Z2)' Y2)
*
*  Where U is a product of Givens plane rotations U=U(N)*...*U(1)
*  And U(I) is a rotation in the plane (I,N+1)
*
*                        (     c(i)      s(i) )
*                  U(I)= (                    ) .
*                        (    -s(i)      c(i) )
*
*  When op(Z) are the factors of a least square problems,
*  and RHO is the norm of residuals of least square solution,
*  this is equivalent to adding a new observation (X,Y),
*  and update the factored solution R,Z,RHO of the problem.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  R is a upper triangular factor of R'*R
*          = 'L':  R is a lower triangular factor of R*R'.
*
*  TRANS    (input) CHARACTER*1
*          = 'N':  Z is a N-BY-M Matrix and op=identity.
*          = 'T':  Z is a M-BY-N Matrix and op=transpose
*
*  N       (input) INTEGER
*          The order of the matrix R.  N >= 0.
*
*  M       (input) INTEGER
*          The number of columns the matrix Z.  M >= 0.
*          In the case M = 0, the array Z,Y,RHO are not referenced
*
*  R       (input/output) DOUBLE PRECISION array, dimension (LDR,N)
*          On entry, the triangular matrix R.  If UPLO = 'U', the leading
*          N-by-N upper triangular part of R contains the upper
*          triangular part of the matrix R, and the strictly lower
*          triangular part of R is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of R contains the lower
*          triangular part of the matrix R, and the strictly upper
*          triangular part of R is not referenced.
*
*          On exit, if INFO = 0, R contains the update Cholesky factor
*
*  LDR     (input) INTEGER
*          The leading dimension of the array R, LDR >= max(1,N).
*
*  X       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, X contains the vector X used to update R.
*          On exit, X is destructed.
*
*  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,k)
*          On entry, Z contains a list of M vectors of dimension N
*          to be updated with R.
*          If TRANS = 'N', Z is considered to be N-by-M (column vectors)
*          If TRANS = 'T', Z is considered to be M-by-N (row vectors)
*          On exit, if INFO = 0, Z contains the updated matrix Z
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.
*          If TRANS = 'N' LDZ >= max(1,N).
*          If TRANS = 'T' LDZ >= max(1,M).
*
*  Y       (input/output) DOUBLE PRECISION array, dimension (M)
*          On entry, Y contains the scalar quantities used to update Z.
*          on exit, Y contains the updated values Y2.
*
*  RHO     (input/output) DOUBLE PRECISION array, dimension (M)
*          on entry, RHO contains the norm of residual vectors.
*          on exit, RHO contains the updated norm of residual vectors.
*
*  C       (output) DOUBLE PRECISION array, dimension (N)
*          C contains the cosines of performed Givens rotations.
*
*  S       (output) DOUBLE PRECISION array, dimension (N)
*          S contains the sines of performed Givens rotations.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the ith residual RHO was negative.
*
*
*  Note:
*  This is a variant of LINPACKs' dchud [1] & [2]
*  but with the following differences taking [3] into account
*  - uses BLAS drot and drotg for Given's rotations
*  - forces positive elements on the diagonal of R
*  - accepts upper or lower Cholesky factor R
*
*  Author: Nicolas Cellier
*
*  References:
*  [1] Linpack dchud http://www.netlib.org/linpack/dchud.f
*  [2] LINPACK User Guide http://www.uploading.com/files/F89HJTGH/linguide.zip.html
*  [3] Low Rank Updates for the Cholesky Decomposition - Matthias Seeger
*      http://infoscience.epfl.ch/record/161468/files/cholupdate.pdf
*
*     .. Local Scalars ..
      LOGICAL UPPER,DOTRAN
      INTEGER I
      DOUBLE PRECISION SC,W
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     .. External Subroutines ..
      EXTERNAL DROT,DROTG
*     .. Intrinsic Functions ..
      INTRINSIC DSQRT
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO=0
      UPPER=LSAME(UPLO,'U')
      DOTRAN=LSAME(TRANS,'T')
      IF( .NOT. UPPER .AND. .NOT. LSAME(UPLO,'L') ) THEN
        INFO=-1
      ELSE IF( .NOT.DOTRAN .AND. .NOT. LSAME(TRANS,'N') ) THEN
        INFO=-2
      ELSE IF( N.LT.0 ) THEN
      ELSE IF( N.LT.0 ) THEN
        INFO=-3
      ELSE IF( M.LT.0 ) THEN
        INFO=-4
      ELSE IF( LDR.LT.N ) THEN
        INFO=-6
      ELSE IF (M.GT.0) THEN
        IF(DOTRAN) THEN
          IF( LDZ.LT.M ) INFO=-9
        ELSE
          IF( LDZ.LT.N ) INFO=-9
        ENDIF
      ENDIF
      IF(INFO.NE.0) RETURN
*
*     Quick return if possible
*
      IF(N.EQ.0) RETURN
*
* ------ FIRST PART: UPDATE R ------
*
      IF(UPPER) THEN
*       Compute givens rotations to annihilate x
        DO 10 I=1,N
          CALL DROTG(R(I,I), X(I), C(I), S(I))
*         Force positive diagonal values
          IF(R(I,I).LT.0.0D0) THEN
            R(I,I)=-R(I,I)
            C(I)=-C(I)
            S(I)=-S(I)
          ENDIF
*         Apply rotation on each row
          IF(I.LT.N) THEN
            CALL DROT(N-I,R(I,I+1),LDR,X(I+1),1,C(I),S(I))
          ENDIF
 10     CONTINUE
      ELSE
*       Compute givens rotations to annihilate x
        DO 20 I=1,N
          CALL DROTG(R(I,I), X(I), C(I), S(I))
*         Force positive diagonal values
          IF(R(I,I).LT.0.0D0) THEN
            R(I,I)=-R(I,I)
            C(I)=-C(I)
            S(I)=-S(I)
          ENDIF
*         Apply rotation on each column
          IF(I.LT.N) THEN
            CALL DROT(N-I,R(I+1,I),1,X(I+1),1,C(I),S(I))
          ENDIF
 20     CONTINUE
      ENDIF
*
* ------  SECOND PART: UPDATE Z AND RHO IF REQUESTED ------
*
      IF(M.GT.0) THEN
        IF(DOTRAN) THEN
*         Update colums of z
          DO 30 I=1,N
            CALL DROT(M,Z(1,I),1,Y(1),1,C(I),S(I));
 30       CONTINUE
        ELSE
*         Update rows of z
          DO 40 I=1,N
            CALL DROT(M,Z(I,1),LDZ,Y(1),1,C(I),S(I));
 40       CONTINUE
        ENDIF
*       Update rho
        DO 50 I=1,M
          W=DABS(Y(I))
          IF(RHO(I).GT.0.0D0) THEN
            SC=RHO(I)+W
            RHO(I)=SC*DSQRT((RHO(I)/SC)**2+(W/SC)**2)
          ELSE
            INFO=I
            GOTO 100
          ENDIF
 50     CONTINUE
      ENDIF
 100  CONTINUE
      RETURN
      END
