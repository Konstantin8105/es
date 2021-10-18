package main

import "github.com/Konstantin8105/f4go/intrinsic"
import "math"

func SSPACE(A []float64, B int, MAXA []int, R int, EIGV int, TT int, W int, AR []float64, BR int, VEC int, D int, RTOLV int, BUP int, BLO int, BUPC int, NN int, NNM int, NWK int, NWM int, NROOT int, RTOL int, NC int, NNC int, NITEM int, IFSS int, IFPR int, NSTIF int, IOUT int) {
	var ART float64
	//C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//C .                                                                   .
	//C .   P R O G R A M                                                   .
	//C .        TO SOLVE FOR THE SMALLEST EIGENVALUES-- ASSUMED .GT. 0 --  .
	//C .        AND CORRESPONDING EIGENVECTORS IN THE GENERALIZED          .
	//C .        EIGENPROBLEM USING THE SUBSPACE ITERATION METHOD           .
	//C .                                                                   .
	//C .  - - INPUT VARIABLES - -                                          .
	//C .        A(NWK)    = STIFFNESS MATRIX IN COMPACTED FORM (ASSUMED    .
	//C .                    POSITIVE DEFINITE)                             .
	//C .        B(NWM)    = MASS MATRIX IN COMPACTED FORM                  .
	//C .        MAXA(NNM) = VECTOR CONTAINING ADDRESSES OF DIAGONAL        .
	//C .                    ELEMENTS OF STIFFNESS MATRIX A                 .
	//C .        R(NN,NC)  = STORAGE FOR EIGENVECTORS                       .
	//C .        EIGV(NC)  = STORAGE FOR EIGENVALUES                        .
	//C .        TT(NN)    = WORKING VECTOR                                 .
	//C .        W(NN)     = WORKING VECTOR                                 .
	//C .        AR(NNC)   = WORKING MATRIX STORING PROJECTION OF K         .
	//C .        BR(NNC)   = WORKING MATRIX STORING PROJECTION OF M         .
	//C .        VEC(NC,NC)= WORKING MATRIX                                 .
	//C .        D(NC)     = WORKING VECTOR                                 .
	//C .        RTOLV(NC) = WORKING VECTOR                                 .
	//C .        BUP(NC)   = WORKING VECTOR                                 .
	//C .        BLO(NC)   = WORKING VECTOR                                 .
	//C .        BUPC(NC)  = WORKING VECTOR                                 .
	//C .        NN        = ORDER OF STIFFNESS AND MASS MATRICES           .
	//C .        NNM       = NN + 1                                         .
	//C .        NWK       = NUMBER OF ELEMENTS BELOW SKYLINE OF            .
	//C .                    STIFFNESS MATRIX                               .
	//C .        NWM       = NUMBER OF ELEMENTS BELOW SKYLINE OF            .
	//C .                    MASS MATRIX                                    .
	//C .                      I. E. NWM=NWK FOR CONSISTENT MASS MATRIX     .
	//C .                            NWM=NN  FOR LUMPED MASS MATRIX         .
	//C .        NROOT     = NUMBER OF REQUIRED EIGENVALUES AND EIGENVECTORS.
	//C .        RTOL      = CONVERGENCE TOLERANCE ON EIGENVALUES           .
	//C .                    ( 1.E-06 OR SMALLER )                          .
	//C .        NC        = NUMBER OF ITERATION VECTORS USED               .
	//C .                    (USUALLY SET TO MIN(2*NROOT, NROOT+8), BUT NC  .
	//C .                    CANNOT BE LARGER THAN THE NUMBER OF MASS       .
	//C .                    DEGREES OF FREEDOM)                            .
	//C .        NNC       = NC*(NC+1)/2 DIMENSION OF STORAGE VECTORS AR,BR .
	//C .        NITEM     = MAXIMUM NUMBER OF SUBSPACE ITERATIONS PERMITTED.
	//C .                    (USUALLY SET TO 16)                            .
	//C .                    THE PARAMETERS NC AND/OR NITEM MUST BE         .
	//C .                    INCREASED IF A SOLUTION HAS NOT CONVERGED      .
	//C .        IFSS      = FLAG FOR STURM SEQUENCE CHECK                  .
	//C .                      EQ.0  NO CHECK                               .
	//C .                      EQ.1  CHECK                                  .
	//C .        IFPR      = FLAG FOR PRINTING DURING ITERATION             .
	//C .                      EQ.0  NO PRINTING                            .
	//C .                      EQ.1  PRINT                                  .
	//C .        NSTIF     = SCRATCH FILE                                   .
	//C .        IOUT      = UNIT USED FOR OUTPUT                           .
	//C .                                                                   .
	//C .  - - OUTPUT - -                                                   .
	//C .        EIGV(NROOT) = EIGENVALUES                                  .
	//C .        R(NN,NROOT) = EIGENVECTORS                                 .
	//C .                                                                   .
	//C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//C .   THIS PROGRAM IS USED IN SINGLE PRECISION ARITHMETIC ON CRAY     .
	//C .   EQUIPMENT AND DOUBLE PRECISION ARITHMETIC ON IBM MACHINES,      .
	//C .   ENGINEERING WORKSTATIONS AND PCS. DEACTIVATE ABOVE LINE FOR     .
	//C .   SINGLE PRECISION ARITHMETIC.                                    .
	//C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//C
	//C     SET TOLERANCE FOR JACOBI ITERATION
	TOLJ = 1.0e-12
	//C
	//C     INITIALIZATION
	//C
	ICONV = 0
	NSCH = 0
	NSMAX = 12
	N1 = NC + 1
	NC1 = NC - 1
	intrinsic.REWIND(NSTIF)
	// Unused by f4go :  WRITE ( NSTIF ) A
	fmt.Println("WRITE SOMETHING")
	for I = 1; I <= NC; I++ {
		D(I) = 0.
	//Label2:
	}
	//C
	//C     ESTABLISH STARTING ITERATION VECTORS
	//C
	ND = NN / NC
	if NWM > NN {
		goto Label4
		//  J = 0
	}
	J = 0
	for I = 1; I <= NN; I++ {
		II = MAXA[I-1]
		R(I, funcint{y:=1;return &y}) = (*B(I))
		if (*B(I)) > 0 {
			J = J + 1
		}
		W(I) = B(I) / A[II-1]
	//Label6:
	}
	if NC <= J {
		goto Label16
		//  WRITE ( IOUT , 1007 )
	}
	// Unused by f4go :  WRITE ( IOUT , 1007 )
	fmt.Println("WRITE SOMETHING")
	goto Label800
	//  4 DO 10 I = 1 , NN
Label4:
	;
	for I = 1; I <= NN; I++ {
		II = MAXA[I-1]
		R(I, funcint{y:=1;return &y}) = (*B(II))
		W(I) = B(II) / A[II-1]
	//Label10:
	}
Label16:
	;
	for J = 2; J <= NC; J++ {
		for I = 1; I <= NN; I++ {
			R(I, J) = 0.
		//Label20:
		}
	}
	//C
	L = NN - ND
	for J = 2; J <= NC; J++ {
		RT = 0.
		for I = 1; I <= L; I++ {
			if (*W(I)) < RT {
				goto Label40
				//  RT = W ( I )
			}
			RT = (*W(I))
			IJ = I
		Label40:
		}
		for I = L; I <= NN; I++ {
			if (*W(I)) <= RT {
				goto Label50
				//  RT = W ( I )
			}
			RT = (*W(I))
			IJ = I
		Label50:
		}
		TT(J) = (*FLOAT(IJ))
		W(IJ) = 0.
		L = L - ND
		R(IJ, J) = 1.
	//Label30:
	}
	//C
	// Unused by f4go :  WRITE ( IOUT , 1008 )
	fmt.Println("WRITE SOMETHING")
	// Unused by f4go :  WRITE ( IOUT , 1002 ) ( TT ( J ) , J = 2 , NC )
	fmt.Println("WRITE SOMETHING")
	//C
	//C     A RANDOM VECTOR IS ADDED TO THE LAST VECTOR
	//C
	PI = 3.141592654e0
	XX = 0.5e0
	for K = 1; K <= NN; K++ {
		XX = math.Pow((PI + XX), 5)
		IX = (*INT(XX))
		XX = XX - FLOAT(IX)
		R(K, NC) = R(K, NC) + XX
	//Label60:
	}
	//C
	//C     FACTORIZE MATRIX A INTO (L)*(D)*(L(T))
	//C
	ISH = 0
	DECOMP((A), (MAXA), NN, ISH, IOUT)
	//C
	//C - - - S T A R T   O F   I T E R A T I O N   L O O P
	//C
	NITE = 0
	TOLJ2 = 1.0e-24
Label100:
	;
	NITE = NITE + 1
	if IFPR == 0 {
		goto Label90
		//  WRITE ( IOUT , 1010 ) NITE
	}
	// Unused by f4go :  WRITE ( IOUT , 1010 ) NITE
	fmt.Println("WRITE SOMETHING")
	//C
	//C     CALCULATE THE PROJECTIONS OF A AND B
	//C
Label90:
	;
	IJ = 0
	for J = 1; J <= NC; J++ {
		for K = 1; K <= NN; K++ {
			TT(K) = (*R(K, J))
		//Label120:
		}
		REDBAK((A), TT, (MAXA), NN)
		for I = J; I <= NC; I++ {
			ART = 0.
			for K = 1; K <= NN; K++ {
				ART = ART + R(K, I)*TT(K)
			//Label140:
			}
			IJ = IJ + 1
			AR[IJ-1] = ART
		//Label130:
		}
		for K = 1; K <= NN; K++ {
			R(K, J) = (*TT(K))
		//Label150:
		}
	//Label110:
	}
	IJ = 0
	for J = 1; J <= NC; J++ {
		MULT(TT, B, R(funcint{y:=1;return &y}, J), (MAXA), NN, NWM)
		for I = J; I <= NC; I++ {
			BRT = 0.
			for K = 1; K <= NN; K++ {
				BRT = BRT + R(K, I)*TT(K)
			//Label190:
			}
			IJ = IJ + 1
			BR(IJ) = BRT
		//Label180:
		}
		if ICONV > 0 {
			goto Label160
			//  DO 200 K = 1 , NN
		}
		for K = 1; K <= NN; K++ {
			R(K, J) = (*TT(K))
		//Label200:
		}
	Label160:
	}
	//C
	//C     SOLVE FOR EIGENSYSTEM OF SUBSPACE OPERATORS
	//C
	if IFPR == 0 {
		goto Label320
		//  IND = 1
	}
	IND = 1
Label210:
	;
	// Unused by f4go :  210 WRITE ( IOUT , 1020 )
	fmt.Println("WRITE SOMETHING")
	II = 1
	for I = 1; I <= NC; I++ {
		ITEMP = II + NC - I
		// Unused by f4go :  WRITE ( IOUT , 1005 ) ( AR ( J ) , J = II , ITEMP )
		fmt.Println("WRITE SOMETHING")
		II = II + N1 - I
	//Label300:
	}
	// Unused by f4go :  WRITE ( IOUT , 1030 )
	fmt.Println("WRITE SOMETHING")
	II = 1
	for I = 1; I <= NC; I++ {
		ITEMP = II + NC - I
		// Unused by f4go :  WRITE ( IOUT , 1005 ) ( BR ( J ) , J = II , ITEMP )
		fmt.Println("WRITE SOMETHING")
		II = II + N1 - I
	//Label310:
	}
	if IND == 2 {
		goto Label350
		//  C
	}
	//C
Label320:
	;
	JACOBI((AR), BR, VEC, EIGV, W, NC, NNC, TOLJ, NSMAX, IFPR, IOUT)
	//C
	if IFPR == 0 {
		goto Label350
		//  WRITE ( IOUT , 1040 )
	}
	// Unused by f4go :  WRITE ( IOUT , 1040 )
	fmt.Println("WRITE SOMETHING")
	IND = 2
	goto Label210
	//  C
	//C
	//C     ARRANGE EIGENVALUES IN ASCENDING ORDER
	//C
Label350:
	;
	IS = 0
	II = 1
	for I = 1; I <= NC1; I++ {
		ITEMP = II + N1 - I
		if (*EIGV(I + 1)) >= (*EIGV(I)) {
			goto Label360
			//  IS = IS + 1
		}
		IS = IS + 1
		EIGVT = (*EIGV(I + 1))
		EIGV(I + 1) = (*EIGV(I))
		EIGV(I) = EIGVT
		BT = (*BR(ITEMP))
		BR(ITEMP) = (*BR(II))
		BR(II) = BT
		for K = 1; K <= NC; K++ {
			RT = (*VEC(K, I+1))
			VEC(K, I+1) = (*VEC(K, I))
			VEC(K, I) = RT
		//Label370:
		}
		II = ITEMP
	Label360:
	}
	if IS > 0 {
		goto Label350
		//  IF ( IFPR .EQ. 0 ) goto 375
	}
	if IFPR == 0 {
		goto Label375
		//  WRITE ( IOUT , 1035 )
	}
	// Unused by f4go :  WRITE ( IOUT , 1035 )
	fmt.Println("WRITE SOMETHING")
	// Unused by f4go :  WRITE ( IOUT , 1006 ) ( EIGV ( I ) , I = 1 , NC )
	fmt.Println("WRITE SOMETHING")
	//C
	//C     CALCULATE B TIMES APPROXIMATE EIGENVECTORS (ICONV.EQ.0)
	//C        OR     FINAL EIGENVECTOR APPROXIMATIONS (ICONV.GT.0)
	//C
Label375:
	;
	for I = 1; I <= NN; I++ {
		for J = 1; J <= NC; J++ {
			TT(J) = (*R(I, J))
		//Label422:
		}
		for K = 1; K <= NC; K++ {
			RT = 0.
			for L = 1; L <= NC; L++ {
				RT = RT + TT(L)*VEC(L, K)
			//Label430:
			}
			R(I, K) = RT
		//Label424:
		}
	//Label420:
	}
	//C
	//C     CALCULATE ERROR BOUNDS AND CHECK FOR CONVERGENCE OF EIGENVALUES
	//C
	for I = 1; I <= NC; I++ {
		VDOT = 0.
		for J = 1; J <= NC; J++ {
			VDOT = VDOT + VEC(I, J)*VEC(I, J)
		//Label382:
		}
		EIGV2 = EIGV(I) * EIGV(I)
		DIF = VDOT - EIGV2
		RDIF = intrinsic.MAX(DIF, TOLJ2*EIGV2) / EIGV2
		RDIF = (intrinsic.SQRT(RDIF))
		RTOLV(I) = RDIF
	//Label380:
	}
	if IFPR == 0 && ICONV == 0 {
		goto Label385
		//  WRITE ( IOUT , 1050 )
	}
	// Unused by f4go :  WRITE ( IOUT , 1050 )
	fmt.Println("WRITE SOMETHING")
	// Unused by f4go :  WRITE ( IOUT , 1005 ) ( RTOLV ( I ) , I = 1 , NC )
	fmt.Println("WRITE SOMETHING")
Label385:
	;
	if ICONV > 0 {
		goto Label500
		//  C
	}
	//C
	for I = 1; I <= NROOT; I++ {
		if (*RTOLV(I)) > RTOL {
			goto Label400
			//  390 CONTINUE
		}
	//Label390:
	}
	// Unused by f4go :  WRITE ( IOUT , 1060 ) RTOL
	fmt.Println("WRITE SOMETHING")
	ICONV = 1
	goto Label100
	//  400 IF ( NITE .LT. NITEM ) goto 100
Label400:
	;
	if NITE < NITEM {
		goto Label100
		//  WRITE ( IOUT , 1070 )
	}
	// Unused by f4go :  WRITE ( IOUT , 1070 )
	fmt.Println("WRITE SOMETHING")
	ICONV = 2
	IFSS = 0
	goto Label100
	//  C
	//C
	//C - - - E N D   O F   I T E R A T I O N   L O O P
	//C
Label500:
	;
	// Unused by f4go :  500 WRITE ( IOUT , 1100 )
	fmt.Println("WRITE SOMETHING")
	// Unused by f4go :  WRITE ( IOUT , 1006 ) ( EIGV ( I ) , I = 1 , NROOT )
	fmt.Println("WRITE SOMETHING")
	// Unused by f4go :  WRITE ( IOUT , 1110 )
	fmt.Println("WRITE SOMETHING")
	for J = 1; J <= NROOT; J++ {
		// Unused by f4go :  WRITE ( IOUT , 1005 ) ( R ( K , J ) , K = 1 , NN )
		fmt.Println("WRITE SOMETHING")
	//Label530:
	}
	//C
	//C     CALCULATE AND PRINT ERROR MEASURES
	//C
	intrinsic.REWIND(NSTIF)
	//  READ ( NSTIF ) A
	//C
	for L = 1; L <= NROOT; L++ {
		RT = (*EIGV(L))
		MULT(TT, (A), R(funcint{y:=1;return &y}, L), (MAXA), NN, NWK)
		VNORM = 0.
		for I = 1; I <= NN; I++ {
			VNORM = VNORM + TT(I)*TT(I)
		//Label590:
		}
		MULT(W, B, R(funcint{y:=1;return &y}, L), (MAXA), NN, NWM)
		WNORM = 0.
		for I = 1; I <= NN; I++ {
			TT(I) = TT(I) - RT*W(I)
			WNORM = WNORM + TT(I)*TT(I)
		//Label600:
		}
		VNORM = (intrinsic.SQRT(VNORM))
		WNORM = (intrinsic.SQRT(WNORM))
		D(L) = WNORM / VNORM
	//Label580:
	}
	// Unused by f4go :  WRITE ( IOUT , 1115 )
	fmt.Println("WRITE SOMETHING")
	// Unused by f4go :  WRITE ( IOUT , 1005 ) ( D ( I ) , I = 1 , NROOT )
	fmt.Println("WRITE SOMETHING")
	//C
	//C     APPLY STURM SEQUENCE CHECK
	//C
	if IFSS == 0 {
		goto Label900
		//  CALL SCHECK ( EIGV , RTOLV , BUP , BLO , BUPC , D , NC , NEI , RTOL , SHIFT , IOUT )
	}
	SCHECK(EIGV, RTOLV, BUP, BLO, BUPC, D, NC, NEI, RTOL, SHIFT, IOUT)
	//C
	// Unused by f4go :  WRITE ( IOUT , 1120 ) SHIFT
	fmt.Println("WRITE SOMETHING")
	//C
	//C     SHIFT MATRIX A
	//C
	intrinsic.REWIND(NSTIF)
	//  READ ( NSTIF ) A
	if NWM > NN {
		goto Label645
		//  DO 640 I = 1 , NN
	}
	for I = 1; I <= NN; I++ {
		II = MAXA[I-1]
		A[II-1] = A[II-1] - B(I)*SHIFT
	//Label640:
	}
	goto Label660
	//  645 DO 650 I = 1 , NWK
Label645:
	;
	for I = 1; I <= NWK; I++ {
		A[I-1] = A[I-1] - B(I)*SHIFT
	//Label650:
	}
	//C
	//C     FACTORIZE SHIFTED MATRIX
	//C
Label660:
	;
	ISH = 1
	DECOMP((A), (MAXA), NN, ISH, IOUT)
	//C
	//C     COUNT NUMBER OF NEGATIVE DIAGONAL ELEMENTS
	//C
	NSCH = 0
	for I = 1; I <= NN; I++ {
		II = MAXA[I-1]
		if A[II-1] < 0. {
			NSCH = NSCH + 1
		}
	//Label664:
	}
	if NSCH == NEI {
		goto Label670
		//  NMIS = NSCH - NEI
	}
	NMIS = NSCH - NEI
	// Unused by f4go :  WRITE ( IOUT , 1130 ) NMIS
	fmt.Println("WRITE SOMETHING")
	goto Label900
	//  670 WRITE ( IOUT , 1140 ) NSCH
Label670:
	;
	// Unused by f4go :  670 WRITE ( IOUT , 1140 ) NSCH
	fmt.Println("WRITE SOMETHING")
	goto Label900
	//  C
	//C
Label800:
	;
	panic("")
Label900:
	;
	return
	//C
//Label1002:
	;
	// Unused by f4go :  1002 FORMAT ( " " , 10 F10 . 0 )
//Label1005:
	;
	// Unused by f4go :  1005 FORMAT ( " " , 12e11 .4 )
//Label1006:
	;
	// Unused by f4go :  1006 FORMAT ( " " , 6e22 .14 )
//Label1007:
	;
	// Unused by f4go :  1007 FORMAT ( // / , " STOP, NC IS LARGER THAN THE NUMBER OF MASS " , "DEGREES OF FREEDOM" )
//Label1008:
	;
	// Unused by f4go :  1008 FORMAT ( // / , " DEGREES OF FREEDOM EXCITED BY UNIT STARTING " , "ITERATION VECTORS" )
//Label1010:
	;
	// Unused by f4go :  1010 FORMAT ( // , " I T E R A T I O N   N U M B E R " , I8 )
//Label1020:
	;
	// Unused by f4go :  1020 FORMAT ( / , " PROJECTION OF A (MATRIX AR)" )
//Label1030:
	;
	// Unused by f4go :  1030 FORMAT ( / , " PROJECTION OF B (MATRIX BR)" )
//Label1035:
	;
	// Unused by f4go :  1035 FORMAT ( / , " EIGENVALUES OF AR-LAMBDA*BR" )
//Label1040:
	;
	// Unused by f4go :  1040 FORMAT ( // , " AR AND BR AFTER JACOBI DIAGONALIZATION" )
//Label1050:
	;
	// Unused by f4go :  1050 FORMAT ( / , " ERROR BOUNDS REACHED ON EIGENVALUES" )
//Label1060:
	;
	// Unused by f4go :  1060 FORMAT ( // / , " CONVERGENCE REACHED FOR RTOL " , E10 . 4 )
//Label1070:
	;
	// Unused by f4go :  1070 FORMAT ( " *** NO CONVERGENCE IN MAXIMUM NUMBER OF ITERATIONS" , " PERMITTED" , / , " WE ACCEPT CURRENT ITERATION VALUES" , / , " THE STURM SEQUENCE CHECK IS NOT PERFORMED" )
//Label1100:
	;
	// Unused by f4go :  1100 FORMAT ( // / , " THE CALCULATED EIGENVALUES ARE" )
//Label1115:
	;
	// Unused by f4go :  1115 FORMAT ( // , " ERROR MEASURES ON THE EIGENVALUES" )
//Label1110:
	;
	// Unused by f4go :  1110 FORMAT ( // , " THE CALCULATED EIGENVECTORS ARE" , / )
//Label1120:
	;
	// Unused by f4go :  1120 FORMAT ( // / , " CHECK APPLIED AT SHIFT " , E22 . 14 )
//Label1130:
	;
	// Unused by f4go :  1130 FORMAT ( // , " THERE ARE " , I8 , " EIGENVALUES MISSING" )
//Label1140:
	;
	// Unused by f4go :  1140 FORMAT ( // , " WE FOUND THE LOWEST " , I8 , " EIGENVALUES" )
	//C
}

func DECOMP(A *[ ]float64, MAXA int, NN int, ISH int, IOUT int) {
	//C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//C .                                                                   .
	//C .   P R O G R A M                                                   .
	//C .        TO CALCULATE (L)*(D)*(L)(T) FACTORIZATION OF               .
	//C .        STIFFNESS MATRIX                                           .
	//C .                                                                   .
	//C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//C
	if NN == 1 {
		goto Label900
		//  C
	}
	//C
	for N = 1; N <= NN; N++ {
		KN = (*MAXA(N))
		KL = KN + 1
		KU = MAXA(N+1) - 1
		KH = KU - KL
		if *KH {
		//Label304:
		}
		, 240 , 210
	//Label210:
		;
		K = N - KH
		IC = 0
		KLT = KU
		for J = 1; J <= KH; J++ {
			IC = IC + 1
			KLT = KLT - 1
			KI = (*MAXA(K))
			ND = MAXA(K+1) - KI - 1
			if *ND {
				, 260 , 270
			}
		//Label260:
		}
	//Label270:
		;
		KK = (*MIN0(IC, ND))
		C = 0.
		for L = 1; L <= KK; L++ {
			C = C + A[KI+L-1]*A[KLT+L-1]
		//Label280:
		}
		A[KLT-1] = A[KLT-1] - C
	//Label260:
		;
		K = K + 1
	//Label240:
		;
		K = N
		B = 0.
		for KK = KL; KK <= KU; KK++ {
			K = K - 1
			KI = (*MAXA(K))
			C = A[KK-1] / A[KI-1]
			if (intrinsic.ABS(C)) < 1.e07 {
				goto Label290
				//  WRITE ( IOUT , 2010 ) N , C
			}
			// Unused by f4go :  WRITE ( IOUT , 2010 ) N , C
			fmt.Println("WRITE SOMETHING")
			goto Label800
			//  290 B = B + C * A ( KK )
		Label290:
			;
			B = B + C*A[KK-1]
			A[KK-1] = C
		//Label300:
		}
		A[KN-1] = A[KN-1] - B
	//Label304:
		;
		if A[KN-1] {
		//Label310:
		}
		, 310 , 200
	//Label310:
		;
		if ISH == 0 {
			goto Label320
			//  IF ( A ( KN ) .EQ. 0. ) A ( KN ) = - 1.e-16
		}
		if A[KN-1] == 0. {
			A[KN-1] = -1.e-16
		}
		goto Label200
		//  320 WRITE ( IOUT , 2000 ) N , A ( KN )
	Label320:
		;
		// Unused by f4go :  320 WRITE ( IOUT , 2000 ) N , A ( KN )
		fmt.Println("WRITE SOMETHING")
		goto Label800
		//  200 CONTINUE
	Label200:
	}
	goto Label900
	//  C
	//C
Label800:
	;
	panic("")
Label900:
	;
	return
	//C
//Label2000:
	;
	// Unused by f4go :  2000 FORMAT ( // " STOP - STIFFNESS MATRIX NOT POSITIVE DEFINITE" , // , " NONPOSITIVE PIVOT FOR EQUATION " , I8 , // , " PIVOT = " , E20 . 12 )
//Label2010:
	;
	// Unused by f4go :  2010 FORMAT ( // " STOP - STURM SEQUENCE CHECK FAILED BECAUSE OF" , " MULTIPLIER GROWTH FOR COLUMN NUMBER " , I8 , // , " MULTIPLIER = " , E20 . 8 )
}

func REDBAK(A *[ ]float64, V int, MAXA int, NN int) {
	//C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//C .                                                                   .
	//C .   P R O G R A M                                                   .
	//C .        TO REDUCE AND BACK-SUBSTITUTE ITERATION VECTORS            .
	//C .                                                                   .
	//C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//C
	//C
	for N = 1; N <= NN; N++ {
		KL = MAXA(N) + 1
		KU = MAXA(N+1) - 1
		if KU - KL {
			, 410 , 410
		}
	//Label400:
	}
//Label410:
	;
	K = N
	C = 0.
	for KK = KL; KK <= KU; KK++ {
		K = K - 1
		C = C + A[KK-1]*V(K)
	//Label420:
	}
	V(N) = V(N) - C
//Label400:
	;
	//C
	for N = 1; N <= NN; N++ {
		K = (*MAXA(N))
		V(N) = V(N) / A[K-1]
	//Label480:
	}
	if NN == 1 {
		goto Label900
		//  N = NN
	}
	N = NN
	for L = 2; L <= NN; L++ {
		KL = MAXA(N) + 1
		KU = MAXA(N+1) - 1
		if KU - KL {
			, 510 , 510
		}
	//Label500:
	}
//Label510:
	;
	K = N
	for KK = KL; KK <= KU; KK++ {
		K = K - 1
		V(K) = V(K) - A[KK-1]*V(N)
	//Label520:
	}
//Label500:
	;
	N = N - 1
	//C
Label900:
	;
	return
}

func MULT(TT *[ ]int, B *[ ]int, RR *[ ]int, MAXA *[ ]int, NN int, NWM int) {
	var AA float64
	//C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//C .                                                                   .
	//C .   P R O G R A M                                                   .
	//C .        TO EVALUATE PRODUCT OF B TIMES RR AND STORE RESULT IN TT   .
	//C .                                                                   .
	//C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//C
	//C
	if NWM > NN {
		goto Label20
		//  DO 10 I = 1 , NN
	}
	for I = 1; I <= NN; I++ {
		TT[I-1] = B[I-1] * RR[I-1]
	//Label10:
	}
	goto Label900
	//  C
	//C
Label20:
	;
	for I = 1; I <= NN; I++ {
		TT[I-1] = 0.
	//Label40:
	}
	for I = 1; I <= NN; I++ {
		KL = MAXA[I-1]
		KU = MAXA[I+1-1] - 1
		II = I + 1
		CC = RR[I-1]
		for KK = KL; KK <= KU; KK++ {
			II = II - 1
			TT[II-1] = TT[II-1] + B[KK-1]*CC
		//Label100:
		}
	}
	if NN == 1 {
		goto Label900
		//  DO 200 I = 2 , NN
	}
	for I = 2; I <= NN; I++ {
		KL = MAXA[I-1] + 1
		KU = MAXA[I+1-1] - 1
		if KU - KL {
			, 210 , 210
		}
	//Label200:
	}
//Label210:
	;
	II = I
	AA = 0.
	for KK = KL; KK <= KU; KK++ {
		II = II - 1
		AA = AA + B[KK-1]*RR[II-1]
	//Label220:
	}
	TT[I-1] = TT[I-1] + AA
//Label200:
	;
	//C
Label900:
	;
	return
}

func SCHECK(EIGV []int, RTOLV []int, BUP []int, BLO []int, BUPC []int, NEIV []int, NC int, NEI int, RTOL int, SHIFT int, IOUT int) {
	//C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//C .                                                                   .
	//C .   P R O G R A M                                                   .
	//C .        TO EVALUATE SHIFT FOR STURM SEQUENCE CHECK                 .
	//C .                                                                   .
	//C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	//C
	//C
	FTOL = 0.01
	//C
	for I = 1; I <= NC; I++ {
		BUP[I-1] = EIGV[I-1] * (1. + FTOL)
		BLO[I-1] = EIGV[I-1] * (1. - FTOL)
	//Label100:
	}
	NROOT = 0
	for I = 1; I <= NC; I++ {
		if RTOLV[I-1] < RTOL {
			NROOT = NROOT + 1
		}
	//Label120:
	}
	if NROOT >= 1 {
		goto Label200
		//  WRITE ( IOUT , 1010 )
	}
	// Unused by f4go :  WRITE ( IOUT , 1010 )
	fmt.Println("WRITE SOMETHING")
	goto Label800
	//  C
	//C
	//C      FIND UPPER BOUNDS ON EIGENVALUE CLUSTERS
	//C
Label200:
	;
	for I = 1; I <= NROOT; I++ {
		NEIV[I-1] = 1
	//Label240:
	}
	if NROOT != 1 {
		goto Label260
		//  BUPC ( 1 ) = BUP ( 1 )
	}
	BUPC[1-1] = BUP[1-1]
	LM = 1
	L = 1
	I = 2
	goto Label295
	//  260 L = 1
Label260:
	;
	L = 1
	I = 2
Label270:
	;
	if BUP[I-1-1] <= BLO[I-1] {
		goto Label280
		//  NEIV ( L ) = NEIV ( L ) + 1
	}
	NEIV[L-1] = NEIV[L-1] + 1
	I = I + 1
	if I <= NROOT {
		goto Label270
		//  280 BUPC ( L ) = BUP ( I - 1 )
	}
Label280:
	;
	BUPC[L-1] = BUP[I-1-1]
	if I > NROOT {
		goto Label290
		//  L = L + 1
	}
	L = L + 1
	I = I + 1
	if I <= NROOT {
		goto Label270
		//  BUPC ( L ) = BUP ( I - 1 )
	}
	BUPC[L-1] = BUP[I-1-1]
Label290:
	;
	LM = L
	if NROOT == NC {
		goto Label300
		//  295 IF ( BUP ( I - 1 ) .LE. BLO ( I ) ) goto 300
	}
Label295:
	;
	if BUP[I-1-1] <= BLO[I-1] {
		goto Label300
		//  IF ( RTOLV ( I ) .GT. RTOL ) goto 300
	}
	if RTOLV[I-1] > RTOL {
		goto Label300
		//  BUPC ( L ) = BUP ( I )
	}
	BUPC[L-1] = BUP[I-1]
	NEIV[L-1] = NEIV[L-1] + 1
	NROOT = NROOT + 1
	if NROOT == NC {
		goto Label300
		//  I = I + 1
	}
	I = I + 1
	goto Label295
	//  C
	//C
	//C      FIND SHIFT
	//C
Label300:
	;
	// Unused by f4go :  300 WRITE ( IOUT , 1020 )
	fmt.Println("WRITE SOMETHING")
	// Unused by f4go :  WRITE ( IOUT , 1005 ) ( BUPC ( I ) , I = 1 , LM )
	fmt.Println("WRITE SOMETHING")
	// Unused by f4go :  WRITE ( IOUT , 1030 )
	fmt.Println("WRITE SOMETHING")
	// Unused by f4go :  WRITE ( IOUT , 1006 ) ( NEIV ( I ) , I = 1 , LM )
	fmt.Println("WRITE SOMETHING")
	LL = LM - 1
	if LM == 1 {
		goto Label310
		//  330 DO 320 I = 1 , LL
	}
Label330:
	;
	for I = 1; I <= LL; I++ {
		NEIV[L-1] = NEIV[L-1] + NEIV[I-1]
	//Label320:
	}
	L = L - 1
	LL = LL - 1
	if L != 1 {
		goto Label330
		//  310 WRITE ( IOUT , 1040 )
	}
Label310:
	;
	// Unused by f4go :  310 WRITE ( IOUT , 1040 )
	fmt.Println("WRITE SOMETHING")
	// Unused by f4go :  WRITE ( IOUT , 1006 ) ( NEIV ( I ) , I = 1 , LM )
	fmt.Println("WRITE SOMETHING")
	L = 0
	for I = 1; I <= LM; I++ {
		L = L + 1
		if NEIV[I-1] >= NROOT {
			goto Label350
			//  340 CONTINUE
		}
	//Label340:
	}
Label350:
	;
	SHIFT = BUPC[L-1]
	NEI = NEIV[L-1]
	goto Label900
	//  C
	//C
Label800:
	;
	panic("")
Label900:
	;
	return
	//C
//Label1005:
	;
	// Unused by f4go :  1005 FORMAT ( " " , 6e22 .14 )
//Label1006:
	;
	// Unused by f4go :  1006 FORMAT ( " " , 6 I22 )
//Label1010:
	;
	// Unused by f4go :  1010 FORMAT ( " *** ERROR ***  SOLUTION STOP IN *SCHECK*" , / , " NO EIGENVALUES FOUND" , / )
//Label1020:
	;
	// Unused by f4go :  1020 FORMAT ( // / , " UPPER BOUNDS ON EIGENVALUE CLUSTERS" )
//Label1030:
	;
	// Unused by f4go :  1030 FORMAT ( // , " NO. OF EIGENVALUES IN EACH CLUSTER" )
//Label1040:
	;
	// Unused by f4go :  1040 FORMAT ( " NO. OF EIGENVALUES LESS THAN UPPER BOUNDS" )
}

func JACOBI(A []float64, B int, X int, EIGV int, D int, N int, NWA int, RTOL int, NSMAX int, IFPR int, IOUT int) {
	var AKK float64
	var AJJ float64
	var AB float64
	var ABCH float64
	var AKKCH float64
	var AJJCH float64
	var AJ float64
	var AK float64
	//C .....................................................................
	//C .                                                                   .
	//C .   P R O G R A M                                                   .
	//C .        TO SOLVE THE GENERALIZED EIGENPROBLEM USING THE            .
	//C .        GENERALIZED JACOBI ITERATION                               .
	//C .....................................................................
	//C
	//C     INITIALIZE EIGENVALUE AND EIGENVECTOR MATRICES
	//C
	N1 = N + 1
	II = 1
	for I = 1; I <= N; I++ {
		if A[II-1] > 0. && B(II) > 0. {
			goto Label4
			//  WRITE ( IOUT , 2020 ) II , A ( II ) , B ( II )
		}
		// Unused by f4go :  WRITE ( IOUT , 2020 ) II , A ( II ) , B ( II )
		fmt.Println("WRITE SOMETHING")
		goto Label800
		//  4 D ( I ) = A ( II ) / B ( II )
	Label4:
		;
		D(I) = A[II-1] / B(II)
		EIGV(I) = (*D(I))
		II = II + N1 - I
	//Label10:
	}
	for I = 1; I <= N; I++ {
		for J = 1; J <= N; J++ {
			X(I, J) = 0.
		//Label20:
		}
		X(I, I) = 1.
	//Label30:
	}
	if N == 1 {
		goto Label900
		//  C
	}
	//C
	//C     INITIALIZE SWEEP COUNTER AND BEGIN ITERATION
	//C
	NSWEEP = 0
	NR = N - 1
Label40:
	;
	NSWEEP = NSWEEP + 1
	if IFPR == 1 {
		// Unused by f4go :  IF ( IFPR .EQ. 1 ) WRITE ( IOUT , 2000 ) NSWEEP
		fmt.Println("WRITE SOMETHING")
	}
	//C
	//C     CHECK IF PRESENT OFF-DIAGONAL ELEMENT IS LARGE ENOUGH TO REQUIRE
	//C     ZEROING
	//C
	EPS = math.Pow((.01), (NSWEEP * 2))
	for J = 1; J <= NR; J++ {
		JP1 = J + 1
		JM1 = J - 1
		LJK = JM1*N - JM1*J/2
		JJ = LJK + J
		for K = JP1; K <= N; K++ {
			KP1 = K + 1
			KM1 = K - 1
			JK = LJK + K
			KK = KM1*N - KM1*K/2 + K
			EPTOLA = (A[JK-1] / A[JJ-1]) * (A[JK-1] / A[KK-1])
			EPTOLB = (B(JK) / B(JJ)) * (B(JK) / B(KK))
			if EPTOLA < EPS && EPTOLB < EPS {
				goto Label210
				//  C
			}
			//C
			//C     IF ZEROING IS REQUIRED, CALCULATE THE ROTATION MATRIX ELEMENTS CA
			//C     AND CG
			//C
			AKK = A[KK-1]*B(JK) - B(KK)*A[JK-1]
			AJJ = A[JJ-1]*B(JK) - B(JJ)*A[JK-1]
			AB = A[JJ-1]*B(KK) - A[KK-1]*B(JJ)
			SCALE = A[KK-1] * B(KK)
			ABCH = AB / SCALE
			AKKCH = AKK / SCALE
			AJJCH = AJJ / SCALE
			CHECK = (ABCH*ABCH + 4.0*AKKCH*AJJCH) / 4.0
			if *CHECK {
			//Label50:
			}
			, 60 , 60
		//Label50:
			;
			// Unused by f4go :  50 WRITE ( IOUT , 2020 ) JJ , A ( JJ ) , B ( JJ )
			fmt.Println("WRITE SOMETHING")
			goto Label800
			//  60 SQCH = SCALE * SQRT ( CHECK )
		//Label60:
			;
			SQCH = SCALE * intrinsic.SQRT(CHECK)
			D1 = AB/2. + SQCH
			D2 = AB/2. - SQCH
			DEN = D1
			if (intrinsic.ABS(D2)) > (intrinsic.ABS(D1)) {
				DEN = D2
			}
			if *DEN {
			//Label80:
			}
			, 70 , 80
		//Label70:
			;
			CA = 0.
			CG = -A[JK-1] / A[KK-1]
			goto Label90
			//  80 CA = AKK / DEN
		//Label80:
			;
			CA = AKK / DEN
			CG = -AJJ / DEN
			//C
			//C     PERFORM THE GENERALIZED ROTATION TO ZERO THE PRESENT OFF-DIAGONAL
			//C     ELEMENT
			//C
		Label90:
			;
			if N - 2 {
			//Label100:
			}
			, 190 , 100
		//Label100:
			;
			if JM1 - 1 {
			//Label130:
			}
			, 110 , 110
		//Label110:
			;
			for I = 1; I <= JM1; I++ {
				IM1 = I - 1
				IJ = IM1*N - IM1*I/2 + J
				IK = IM1*N - IM1*I/2 + K
				AJ = A[IJ-1]
				BJ = (*B(IJ))
				AK = A[IK-1]
				BK = (*B(IK))
				A[IJ-1] = AJ + CG*AK
				B(IJ) = BJ + CG*BK
				A[IK-1] = AK + CA*AJ
				B(IK) = BK + CA*BJ
			//Label120:
			}
		//Label130:
			;
			if KP1 - N {
			//Label140:
			}
			, 140 , 160
		//Label140:
			;
			LJI = JM1*N - JM1*J/2
			LKI = KM1*N - KM1*K/2
			for I = KP1; I <= N; I++ {
				JI = LJI + I
				KI = LKI + I
				AJ = A[JI-1]
				BJ = (*B(JI))
				AK = A[KI-1]
				BK = (*B(KI))
				A[JI-1] = AJ + CG*AK
				B(JI) = BJ + CG*BK
				A[KI-1] = AK + CA*AJ
				B(KI) = BK + CA*BJ
			//Label150:
			}
		//Label160:
			;
			if JP1 - KM1 {
			//Label170:
			}
			, 170 , 190
		//Label170:
			;
			LJI = JM1*N - JM1*J/2
			for I = JP1; I <= KM1; I++ {
				JI = LJI + I
				IM1 = I - 1
				IK = IM1*N - IM1*I/2 + K
				AJ = A[JI-1]
				BJ = (*B(JI))
				AK = A[IK-1]
				BK = (*B(IK))
				A[JI-1] = AJ + CG*AK
				B(JI) = BJ + CG*BK
				A[IK-1] = AK + CA*AJ
				B(IK) = BK + CA*BJ
			//Label180:
			}
		//Label190:
			;
			AK = A[KK-1]
			BK = (*B(KK))
			A[KK-1] = AK + 2.*CA*A[JK-1] + CA*CA*A[JJ-1]
			B(KK) = BK + 2.*CA*B(JK) + CA*CA*B(JJ)
			A[JJ-1] = A[JJ-1] + 2.*CG*A[JK-1] + CG*CG*AK
			B(JJ) = B(JJ) + 2.*CG*B(JK) + CG*CG*BK
			A[JK-1] = 0.
			B(JK) = 0.
			//C
			//C     UPDATE THE EIGENVECTOR MATRIX AFTER EACH ROTATION
			//C
			for I = 1; I <= N; I++ {
				XJ = (*X(I, J))
				XK = (*X(I, K))
				X(I, J) = XJ + CG*XK
				X(I, K) = XK + CA*XJ
			//Label200:
			}
		Label210:
		}
	}
	//C
	//C     UPDATE THE EIGENVALUES AFTER EACH SWEEP
	//C
	II = 1
	for I = 1; I <= N; I++ {
		if A[II-1] > 0. && B(II) > 0. {
			goto Label215
			//  WRITE ( IOUT , 2020 ) II , A ( II ) , B ( II )
		}
		// Unused by f4go :  WRITE ( IOUT , 2020 ) II , A ( II ) , B ( II )
		fmt.Println("WRITE SOMETHING")
		goto Label800
		//  215 EIGV ( I ) = A ( II ) / B ( II )
	Label215:
		;
		EIGV(I) = A[II-1] / B(II)
		II = II + N1 - I
	//Label220:
	}
	if IFPR == 0 {
		goto Label230
		//  WRITE ( IOUT , 2030 )
	}
	// Unused by f4go :  WRITE ( IOUT , 2030 )
	fmt.Println("WRITE SOMETHING")
	// Unused by f4go :  WRITE ( IOUT , 2010 ) ( EIGV ( I ) , I = 1 , N )
	fmt.Println("WRITE SOMETHING")
	//C
	//C     CHECK FOR CONVERGENCE
	//C
Label230:
	;
	for I = 1; I <= N; I++ {
		TOL = RTOL * D(I)
		DIF = (intrinsic.ABS(EIGV(I) - D(I)))
		if DIF > TOL {
			goto Label280
			//  240 CONTINUE
		}
	//Label240:
	}
	//C
	//C     CHECK ALL OFF-DIAGONAL ELEMENTS TO SEE IF ANOTHER SWEEP IS
	//C     REQUIRED
	//C
	EPS = math.Pow(RTOL, 2)
	for J = 1; J <= NR; J++ {
		JM1 = J - 1
		JP1 = J + 1
		LJK = JM1*N - JM1*J/2
		JJ = LJK + J
		for K = JP1; K <= N; K++ {
			KM1 = K - 1
			JK = LJK + K
			KK = KM1*N - KM1*K/2 + K
			EPSA = (A[JK-1] / A[JJ-1]) * (A[JK-1] / A[KK-1])
			EPSB = (B(JK) / B(JJ)) * (B(JK) / B(KK))
			if EPSA < EPS && EPSB < EPS {
				goto Label250
				//  goto 280
			}
			goto Label280
			//  250 CONTINUE
		Label250:
		}
	}
	//C
	//C     SCALE EIGENVECTORS
	//C
Label255:
	;
	II = 1
	for I = 1; I <= N; I++ {
		BB = (intrinsic.SQRT(B(II)))
		for K = 1; K <= N; K++ {
			X(K, I) = X(K, I) / BB
		//Label270:
		}
		II = II + N1 - I
	//Label275:
	}
	goto Label900
	//  C
	//C
	//C     UPDATE  D  MATRIX AND START NEW SWEEP, IF ALLOWED
	//C
Label280:
	;
	for I = 1; I <= N; I++ {
		D(I) = (*EIGV(I))
	//Label290:
	}
	if NSWEEP < NSMAX {
		goto Label40
		//  goto 255
	}
	goto Label255
	//  C
	//C
Label800:
	;
	panic("")
Label900:
	;
	return
	//C
//Label2000:
	;
	// Unused by f4go :  2000 FORMAT ( // , " SWEEP NUMBER IN *JACOBI* = " , I8 )
//Label2010:
	;
	// Unused by f4go :  2010 FORMAT ( " " , 6e20 .12 )
//Label2020:
	;
	// Unused by f4go :  2020 FORMAT ( " *** ERROR *** SOLUTION STOP" , / , " MATRICES NOT POSITIVE DEFINITE" , / , " II = " , I8 , " A(II) = " , E20 . 12 , " B(II) = " , E20 . 12 )
//Label2030:
	;
	// Unused by f4go :  2030 FORMAT ( / , " CURRENT EIGENVALUES IN *JACOBI* ARE" , / )
}
