      PROGRAM MAIN
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)                           

          INTEGER NN, NNM
          INTEGER NWK, NWM
          INTEGER NROOT
          INTEGER NC, NNC
          INTEGER NITEM
          INTEGER NSTIF
          INTEGER IFPR
          INTEGER IOUT
          INTEGER I


          PARAMETER (NN  = 3)
          PARAMETER (NNM = 4)
          PARAMETER (NWK = 5)
          PARAMETER (NWM = 5)
          PARAMETER (NROOT = 2)
          PARAMETER (RTOL = 1e-4)
          PARAMETER (NC  = 4)
          PARAMETER (NNC = 10)
          PARAMETER (NSTIF = 6)
          PARAMETER (NITEM = 16)
          PARAMETER (IFPR = 1)
          PARAMETER (IFSS = 1)
          PARAMETER (IOUT = 6)    

          INTEGER MAXA(NNM)
          DIMENSION A(NWK),B(NWM),R(NN,NC),TT(NN),W(NN),EIGV(NC)
          DIMENSION D(NC),VEC(NC,NC),AR(NNC),BR(NNC),RTOLV(NC),BUP(NC)
          DIMENSION BLO(NC),BUPC(NC)


          DO I=1,NN
            DO J=1,NC
                R(I,J)=0.0
            END DO
          END DO
          DO I=1,NC
            DO J=1,NC
                VEC(I,J)=0.0
            END DO
          END DO
          DATA TT   /NN*0.0/
          DATA W    /NN*0.0/
          DATA EIGV /NC*0.0/
          DATA D    /NC*0.0/
          DATA AR   /NNC*0.0/
          DATA BR   /NNC*0.0/
          DATA RTOLV/NC*0.0/
          DATA BUP  /NC*0.0/
          DATA BLO  /NC*0.0/
          DATA BUPC /NC*0.0/
 
          DATA (A(I), I=1,NWK) / 2,-1, 4,-1, 2/
          DATA B               / 0.5,0, 1,0, 0.5 /
          DATA MAXA            / 1,3,5,5/
C lamda1 = 2  (0.707 0.707 0.707)
C lamda2 = 4  (-1    0     1    )

C         CALL SSPACE (A,B,MAXA,R,EIGV,TT,W,AR,BR,VEC,D,RTOLV,BUP,BLO, 
C    3          BUPC,NN,NNM,NWK,NWM,NROOT,RTOL,NC,NNC,NITEM,IFSS,
C    4          IFPR,NSTIF,IOUT) 

          WRITE(*,*) "DECOMP"
          CALL DECOMP (A,MAXA,NN,ISH,IOUT)
          WRITE(*,*) "A"
          WRITE(*,*)  A
      END


