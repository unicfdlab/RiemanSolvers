C
C CREATED FROM E. Toro report scan.
C Probably contains error due to wrong letters recognition
C
C-------------------------------------------------------------------------------
      PROGRAM DRIVER
      PARAMETER (MD=1000)
      DIMENSION D(0:MD+1),U(0:MD+1),P(0:MD+1),RN(10000),TV(10)
      COMMON/CPGAMMA/GP1,GM1,HGM1,DGAM,G1,G2,G3,G4,G5,G6
      COMMON/STATES/DL,UL,PL,DR,UR,PR
      COMMON/STARSO/US,PS,CL,CR
      COMMON/COVOLU/COVL,COVR,B
      COMMON/GAMTOL/GAMMA,TOL
      DATA NC,TIME,POINTER,TOLTIME/0,0.0,0.5,1.D-06/
      DATA (TV(KT),KT=1,2)/0.0002,0.0004/
      READ(99,*)TUBLEN,M,NOTIST,NOPROF,TOL,CFLCOE,GAMMA,B
      CALL ICDATA(M,TUBLEN,DX,GAMMA,D,U,P)
      CALL VDCK12(RN,NOTIST)
      KT=1
C     COMMENCE TIME STEPPING 
      DO 0001 N=1,NOTIST 
C        REFLECTING BOUNDARY CCMDITIONS APPLIED 
         D(0)  =D(1)
         U(0)  =U(1)
         P(0)  =P(1)
         D(M+1)=D(M)
         U(M+1)=-U(M)
         P(M+1)=P(M)
         CALL CFLCON(B,GAMMA,M,D,U,P,DX,DTMIN)
         DT=CFLCOE*DTMIN 
         TITEST=(TIME+DT)
         IF(TITEST.GT.TV(KT))THEN
            DT=TV(KT)-TIME 
         ENDIF
         TIME=TIME+DT
         RAND=RN(N)
         DTDX=DT/DX
         DXDTL=RAND/DTDX
         DXDTR=(RAND-1.)/DTDX 
C        UPDATE SOLUTION TO NEXT TIME LEVEL 
         DO 0003 I=1,M 
            IF (I.EQ.1) THEN 
C               SOLVE RIEMANN PROBLEM AT THE LEFT BOUNDARY 
               DL=D(I-1)
               UL=U(I-l)
               PL=P(I-1)
               DR=D(I)
               UR=U(I)
               PR=P(I)
               CALL RPCOV 
            ENDIF 
            IF (RAND .LE. POINTER) THEN 
               CALL SAMCOV (D1, U1, P1, DXDTL)
            ENDIF 
C            SOLVE RIEMANN PROBLEM RP(I,I+1) 
             DL=D(I) 
             UL=U(I)
             PL=P(I)
             DR=D(I+1)
             UR=U(I+1)
             PR=P(I+1)
             CALL RPCOV 
             IF (RAND.GT.POINTER) THEN 
                CALL SAMCOV(Dl,Ul,Pl,DXDTR)
             ENDIF 
             D(I)=D1
             U(I)=U1
             P(I)=P1 
 0003    CONTINUE      
C        UPDATING COMPLETED 
         TDIF=ABS(TIME-TV(KT)) 
         IF(TDIF.LE.TOLTIME)THEN
            NC=NC+1 
            CALL OUTPUT(TIME,M,NC,NOPROF,GM1,D,U,P,B) 
            IF(NC.EQ.0)THEN 
               WRITE(6,*)'JOB FINISHED OK'
              STOP 
            ENDIF 
            KT=KT+1 
         ENDIF 
 0001 CONTINUE
C     TIME STEPPING COMPLETED 
      END
C
C-------------------------------------------------------------------------------
C
      SUBROUTINE ICDATA(M,TUBLEN,DX,GAMMA,D,U,P)
      PARAMETER (MD=1000)
      DIMENSION D(0:MD+1),U(0:MD+1),P(0:MD+1)
      COMMON/CPGAMMA/GP1,GM1,HGM1,DGAM,G1,G2,G3,G4,G5,G6
      
      DATA DL0,UL0,PL0/100.0,0.0,100.0E+06/
      DATA DR0,UR0,PR0/1.0,0.0,0.1E+06/
C     PLACE OF DISCONTINUITY
      DATA X0/0.4/
      
      GP1=GAMMA+1.0
      GM1=GAMMA-1.0
      HGM1=0.5*GM1
      HGP1=0.5*GP1
      DGAM=1.0/GAMMA
      G1=HGM1/GAMMA
      G2=HGP1/GAMMA
      G3=1.0/G1
      G4=1.0/HGM1
      G5=2.0/GP1
      G6=GM1/GP1
      DX=TUBLEN/REAL(M)
      DO 1000 I=1,M 
         XP=(REAL(I)-5.0D-1)*DX
         IF(XP.LE.X0)THEN
            D(I)=DL0
            U(I)=UL0
            P(I)=PL0
         ELSE 
            D(I)=DRO
            U(I)=URO
            P(I)=PRO
        ENDIF
1000  CONTINUE
      RETURN
      END 
C
C-------------------------------------------------------------------------------
C
      SUBROUTINE VDCK12(RN,NOTIST)
      PARAMETER  (N1=1000,N2=10000)
      DIMENSION NA(N1),JA(N1),RN(N2)
      DATA K1,K2,NRN0/2,1,100/
      DO 0001 NRN=NRN0,NOTIST+NRN0
         IS=0
         MM=NRN
         DO 0002 I=1,100
         IF(MM.EQ.0)GOTO 8888 
            IS=IS+1
            NA(I)=MOD(MM,K1)
            MM=MM/K1
            KL=K2*NA(I)
            JA(I)=MOD(KL,K1) 
0002     CONTINUE
8888     RANNUM=0.0
         DO 0004 K=1,IS
            RANNUM=RANNUM+REAL(JA(K))/(K1**K)
0004     CONTINUE
         NT=NRN-NRN0+1
         RN(NT)=RANNUM
0001  CONTINUE
      RETURN
      END 
C
C-------------------------------------------------------------------------------
C
      SUBROUTINE CFLCON(B,GAMMA,M,D,U,P,DX,DTMIN)
      PARAMETER (MD=1000)
      DIMENSION D(0:MD+1),U(0:MD+1),P(0:MD+1)
      SMAX=0.
      DO 0001 I=1,M
         DENS=D(I)
         COV=1.0-B*DENS
         A=SQRT(GAMMA*P(I)/(COV*DENS))
         SMUA=ABS(U(I))+A
         IF(SMUA.GT.SMAX)SMAX=SMUA
0001  CONTINUE
      DTMIN=DX/SMAX
      RETURN
      END 
C
C-------------------------------------------------------------------------------
C
      SUBROUTINE RPCOV
      COMMON/STATES/DL,UL,PL,DR,UR,PR
      COMMON/STARSO/US,PS,CL,CR
      COMMON/GAMTOL/GAMMA,TOL
      COMMON/COVOLU/COVL,COVR,B
      COMMON/CPGAMMA/GP1,GM1,HGM1,DGAM,G1,G2,G3,G4,G5,G6
C     SOLVES RIEMANN PROBLEM WITH CONSTANT COVOLUME B
      COVL=1.0-B*DL
      COVR=1.0-B*DR
      CL=SQRT(GAMMA*PL/(COVL*DL))
      CR=SQRT(GAMMA*PR/(COVR*DR))
      DELU=UL-UR
C     GUESSED VALUE FOR PS IS PROVIDED
      CLPLG=CL/PL**G1
      CRPRG=CR/PR**G1
      ABOVE=CL*COVL+CR*COVR+HGM1*DELU
      BELOW=CLPLG*COVL+CRPRG*COVR
      PS=(ABOVE/BELOW)**G3
      PS0 =PS
C     START ITERATION
      DO 0001 IT=1,50
C        LEFT WAVE
         IF(PL.LT.PS)THEN
            S1=SQRT(G5*COVL/DL)
            S2=G6*PL
            S2PS=S2+PS
            DELPLPS=PL-PS
            SQS2PS=1.0/SQRT(S2PS)
            FLEFVAL=S1*DELPLPS*SQS2PS
            FLEFDER=-S1*SQS2PS*(1.0+0.5*DELPLPS/S2PS)
         ELSE
            FLEFVAL=G4*COVL*(CL-CLPLG*PS**G1)
            FLEFDER=-DGAM*COVL*CLPLG*PS**(-G2)
         ENDIF
C        RIGHT WAVE
         IF(PR.LT.PS)THEN
            S1=SQRT(G5*COVR/DR)
            S2=G6*PR
            S2PS=S2+PS
            DELPRPS=PR-PS
            SQS2PS=1.0/SQRT(S2PS)
            FRIGVAL=S1*DELPRPS*SQS2PS
            FRIGDER=-S1*SQS2PS*(1.0+0.5*DELPRPS/S2PS)
         ELSE
            FRIGVAL=G4*COVR*(CR-CRPRG*PS**G1)
            FRIGDER=-DGAM*COVR*CRPRG*PS**(-G2)
         ENDIF
         FUNVAL=FLEFVAL+FRIGVAL+DELU
         FUNDER=FLEFDER+FRIGDER
         PS=PS-FUNVAL/FUNDER
         IF(IT.GT.5)THEN
C           SECANT METHOD
            ABOVE=PS0*FUNVAL-PS*FUNVAL0
            BELOW=FUNVAL-FUNVAL0
            PS=ABOVE/BELOW
         ELSE
C           NEWTON RAPHSON METHOD
         ENDIF
         US=0.5*(FLEFVAL-FRIGVAL+UL+UR)
         TESTPS =ABS((PS-PS0)/PS)
         IF(TESTPS.LE.TOL)GOTO 0002
         IF(PS.LT.TOL)PS=TOL
         PS0=PS
         FUNVAL0=FUNVAL
0001  CONTINUE
      WRITE(6,0003)IT
      STOP
0003  FORMAT('DIVERGENCE IN PSTAR STEP, ITERATION NO. =',I4)
0002  CONTINUE
      RETURN
      END
C
C-------------------------------------------------------------------------------
C
      SUBROUTINE SAMCOV(D,U,P,DXDT)
      COMMON/STATES/DL,UL,PL,DR,UR,PR
      COMMON/STARSO/US,PS,CL,CR
      COMMON/COVOLU/COVL,COVR,B
      COMMON/GAMTOL/GAMMA,TOL
      COMMON/CPGAMMA/GP1,GM1,HGM1,DGAM,G1,G2,G3,G4,G5,G6
      IF(DXDT.GE.US)THEN
C        SAMPLING POINT LIES TO THE RIGHT OF SLIP LINE
         IF(PS.LE.PR)THEN
C           RIGHT WAVE IS A RAREFACTION WAVE
            IF(DXDT.LT.(UR+CR))THEN
               AISEN=(DR/COVR)*(PS/PR)**DGAM
               D3=AISEN/(1.0+B*AISEN)
               COV3=1.0-B*D3
               C3=SQRT(GAMMA*PS/(D3*COV3))
               IF(DXDT.LT.(US+C3))THEN
C                 LEFT OF RIGHT RAREFACTION
                  D=D3
                  U=US
                  P=PS
               ELSE
C                 INSIDE RIGHT RAREFACTION
C                 GUESS VALUE FOR D, MEAN VALUE
                  D=0.5*(DR+D3)
                  RARCON=DXDT-UR
                  CALL RARFAN(DXDT,RARCON,D,C4,P,DR,PR,CR,COVR)
                  U=DXDT-C4
              ENDIF
            ELSE
C              RIGHT OF RIGHT RAREFACTION
               D=DR
               U=UR
               P=PR
            ENDIF
         ELSE
C           RIGHT WAVE IS A SHOCK WAVE
            CONS=0.5*GP1*DR*PR/COVR
            PRERAT=PS/PR
            RMR=SQRT(CONS*(PRERAT+GM1/GP1))
            URS=UR+RMR/DR
            IF(DXDT.GE.URS)THEN
C              RIGHT OF RIGHT SHOCK
               D=DR
               U=UR
               P=PR
            ELSE
C              BEHIND RIGHT SHOCK
               ABOVE=GP1*PRERAT+GM1
               TWIBDR=2.0*B*DR
               BELOW=(GM1+TWIBDR)*PRERAT+GP1-TWIBDR
               D=DR*ABOVE/BELOW
               U=US
               P=PS
            ENDIF
         ENDIF
      ELSE
C        SAMPLING POINT LIES TO THE LEFT OF SLIP LINE
         IF(PS.LE.PL)THEN
C           LEFT WAVE IS A RAREFACTION
            AISEN=(DL/COVL)*(PS/PL)**DGAM
            D3=AISEN/(1.0+B*AISEN)
            COV3 =1.0-B*D3
            C3=SQRT(GAMMA*PS/(D3*COV3))
            IF(DXDT.LT.(US-C3))THEN
               IF(DXDT.LT.(UL-CL))THEN
C                 LEFT OF LEFT RAREFACTION
                  D=DL
                  U=UL
                  P=PL
               ELSE
C                 INSIDE LEFT RAREFACTION
C                 GUESS VALUE FOR D, MEAN VALUE
                  D=0.5*(DL+D3)
                  RARCON=-(DXDT-UL)
                  CALL RARFAN(DXDT,RARCON,D,C4,P,DL,PL,CL,COVL)
                  U=DXDT+C4
               ENDIF
            ELSE
C              RIGHT OF LEFT RAREFACTION
               D=D3
               U=US
               P=PS
            ENDIF
         ELSE
C           LEFT WAVE IS A SHOCK WAVE
            CONS=0.5*GP1*DL*PL/COVL
            PRERAT=PS/PL
            RML=SQRT(CONS*(PRERAT+GM1/GP1))
            ULS=UL-RML/DL
            IF(DXDT.GE.ULS)THEN
C              BEHIND LEFT SHOCK
               ABOVE=GP1*PRERAT+GM1
               TWIBDL=2.0*B*DL
               BELOW= (GM1+TWIBDL)*PRERAT+GP1-TWIBDL
               D=DL*ABOVE/BELOW
               U=US
               P=PS
            ELSE
C              LEFT OF LEFT SHOCK
               D=DL
               U=UL
               P=PL
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END
C
C-------------------------------------------------------------------------------
C
      SUBROUTINE RARFAN(DXDT,RARCON,DF,C4,P,DK,PK,CK,COVK)
      COMMON/COVOLU/COVL,COVR,B
      COMMON/GAMTOL/GAMMA,TOL
      COMMON/CPGAMMA/GP1,GM1,HGM1,DGAM,G1,G2,G3,G4,G5,G6
      Z1=RARCON+2.0*CK*COVK/GM1
      Z2=PK*(COVK/DK)**GAMMA
      ZZ=(Z1*GM1)**2/(GAMMA*Z2)
      DF0=DF
      DO 0001 I=1,100
         COVF=1.0-B*DF
         F1 =GP1-2.0*B*DF
         F2 =COVF**GAMMA
         F3 =F1-2.0
         F4 =DF**GM1
         FVAL=F1*F1*F4-ZZ*F2*COVF
C        NEWTON-RAPHSON ITERATION
         FDER=GP1*(B*ZZ*F2+F1*F3*F4/DF)
         DF =DF-FVAL/FDER
         IF(I.GT.5)THEN
C           SECANT METHOD
            ABOVE=DF0*FVAL-DF*FVAL0
            BELOW=FVAL-FVAL0
            DF=ABOVE/BELOW
         ENDIF
         DETED=ABS((DF-DF0)/DF)
         IF(DETED.LE.TOL)GOTO 0002
         IF(DF.LT.TOL)DF=TOL
         DF0 =DF
         FVAL0=FVAL
0001  CONTINUE
      WRITE(6,0004)I
0004  FORMAT(5X,'DIRVERGENCE INSIDE FAN, NO. OF ITER.=',I5)
      STOP
C     COMPUTE OTHER UNKNOWNS
0002  COV4=1.0-B*DF
      P=Z2*(DF/C0V4)**GAMMA
      C4=SQRT(GAMMA*P/(DF*COV4))
0003  CONTINUE
      RETURN
      END
C
C-------------------------------------------------------------------------------
C
      SUBROUTINE OUTPUT(TIME,M,NC,NOPROF,GM1,D,U,P,B)
      PARAMETER (MD=1000)
      DIMENSION D(0:MD+1),U(0:MD+1),P(0:MD+1)
      DIMENSION TM(20),R1(4,20,MD)
      DATA RMPA/1.0E+06/
      TM(NC)=TIME
      GMCONST=GM1*RMPA
      DO 0001 I=1,M
         R1(1,NC,I)=D(I)
         R1(2,NC,I)=U(I)
         R1(3,NC,I)=P(I)/RMPA
         COV=1.0-B*D(I)
         R1(4,NC,I)=(COV*P(I))/(D(I)*GMCONST)
0001  CONTINUE
      IF(NC.EQ.NOPROF)THEN
         WRITE(1,0004)(TM(J),J=1,NOPROF)
         WRITE(2,0004)(TM(J),J=1,NOPROF)
         WRITE(3,0004)(TM(J),J=1,NOPROF)
         WRITE(4,0004)(TM(J),J=l,NOPROF)
         DO 0002 I=1,M
            WRITE(1,0003)I,(R1(1,J,I),J=1,NOPROF)
            WRITE(2,0003)I,(R1(2,J,I),J=1,NOPROF)
            WRITE(3,0003)I,(R1(3,J,I),J=1,NOPROF)
            WRITE(4,0003)1,(R1(4,J,I),J=1,NOPROF)
0002     CONTINUE
         NC=0
      ENDIF
0003  FORMAT(I4,1X,10(F10.4,1X))
0004  FORMAT(5X,10(F7.4,4X))
      RETURN
      END
C
C     END-OF-FILE
C
