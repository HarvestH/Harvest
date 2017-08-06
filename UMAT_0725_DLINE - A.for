      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
     1 DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,
     2 PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,
     3 COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,
     4 KSPT, KSTEP, KINC)
C 
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
C
      DIMENSION STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS, NTENS),
     1 DDSDDT(NTENS), DRPLDE(NTENS), STRAN(NTENS), DSTRAN(NTENS),
     2 PREDEF(1), DPRED(1), PROPS(NPROPS), COORDS(3), DROT(3, 3),
     3 DFGRD0(3, 3), DFGRD1(3, 3),SIGMA(NTENS),TIME(2)
C
      DIMENSION EELAS(6), EPLAS(6), FLOW(6)
      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,SIX=6.D0,
     1 ENUMAX=.4999D0, NEWTON=10, TOLER=1.0D-6)
C
      REAL*8:: EMOD,ENU,EBULK3,EG2,EG,EG3,ELAM,SJ4,P,THEDA1
      REAL*8:: B,C,ALPHA,THEDA_B,SIGMA_C,SI1,SJ2,SJ3,THEDA
      REAL*8:: F1,F0,R,R1,C1,C2,C3,F2,ZH,DELTA_R1,FM,SLMD,SLMDS,F3,XZS1
      REAL*8:: DSTRESS(NTENS)
      REAL*8:: DF(NTENS)
      REAL*8:: FB,FA(6),FC(6),FD(6,6)
      REAL*8:: D_P(NTENS,NTENS),XZS(NTENS)
      REAL*8:: A1(NTENS),A2(NTENS),A3(NTENS)
      REAL*8:: M
      REAL*8:: ROOT3,HARD,ET
      REAL*8,PARAMETER:: PAI=3.1415926,TOL=1.0D-5
C----------------------------------------------------------------
C   UMAT FOR ISOTROPIC ELASTICITY 
C   CANNOT BE USED FOR PLANE STRESS
C----------------------------------------------------------------
C   PROPS(1) - E
C   PROPS(2) - NU
C   PROPS(3..) - YIELD AND HARDENING DATA
C   CALLS UHARD FOR CURVE OF YIELD STRESS VS. PLASTIC STRAIN
C----------------------------------------------------------------
C 
      IF (NDI/=3) THEN
         WRITE (7, *) 'THIS UMAT MAY ONLY BE USED FOR 3D ELEMENTS'
         CALL XIT
      ENDIF
C-------------------------------------------------------
C   ELASTIC PROPERTIES
      EMOD=PROPS(1)
      ENU=PROPS(2)
      EBULK3=EMOD/(ONE-TWO*ENU)
      EG2=EMOD/(ONE+ENU)
      EG=EG2/TWO
      EG3=THREE*EG
      ELAM=(EBULK3-EG2)/THREE
      ROOT3=SQRT(3.0D0)
	  HARD=0.0D0
	  ET=700.0D0

C 	  WRITE(7,*)"STRESS",STRESS
C 	  WRITE(7,*)"DSTRAN",DSTRAN
C 
C   ELASTIC STIFFNESS
C
      DO K1=1,NTENS
        DO K2=1,NTENS
          DDSDDE(K2, K1)=ZERO
        ENDDO
      ENDDO
      DO K1=1,3
        DO K2=1,3
          DDSDDE(K2, K1)=ELAM
        ENDDO
        DDSDDE(K1, K1)=EG2+ELAM
      ENDDO
      DO K1=NDI+1,NTENS
        DDSDDE(K1, K1)=EG
      ENDDO
      CALL ROTSIG(STATEV( 1       ), DROT, EELAS, 2, NDI, NSHR)
      CALL ROTSIG(STATEV(  NTENS+1), DROT, EPLAS, 2, NDI, NSHR)
C-
C 
C CALCULATE STRESS
C 

C STRNGTH PARAMETER OF MATERIALS
      SIGMA_C=PROPS(3)
      THEDA_B=ATAN(ROOT3/THREE)  
C
C START UPDATING STRESS
C
      DO K1=1,NTENS
        SIGMA(K1)=STRESS(K1)
      ENDDO
C
      DO K1=1,NTENS
        DO K2=1,NTENS
          SIGMA(K2)=SIGMA(K2)+DDSDDE(K2,K1)*DSTRAN(K1)
        ENDDO
      ENDDO
C	  
	  DO K1=1, NTENS 
        EELAS(K1)=STATEV(K1)+DSTRAN(K1)
        EPLAS(K1)=STATEV(K1+NTENS)
      ENDDO
      EQPLAS=STATEV(13)
C	  
      DO K1=1,NTENS
C	DSTRESS	
        DSTRESS(K1)=SIGMA(K1)-STRESS(K1)
      ENDDO	
	  THEDA=ZERO
C-------------------------------------------------------
      SI1=SIGMA(1)+SIGMA(2)+SIGMA(3)
      P=SI1/3.0D0
      SJ2=((SIGMA(1)-SIGMA(2))*(SIGMA(1)-SIGMA(2))+(SIGMA(2)-SIGMA
     1   (3))*(SIGMA(2)-SIGMA(3))+(SIGMA(3)-SIGMA(1))*(SIGMA(3)
     2   -SIGMA(1)))/6.0D0+SIGMA(4)*SIGMA(4)+SIGMA(5)*SIGMA(5)+
     3   SIGMA(6)*SIGMA(6)
      SJ3=(SIGMA(1)-P)*(SIGMA(2)-P)*(SIGMA(3)-P)+2.0D0*SIGMA(4)
     1   *SIGMA(5)*SIGMA(6)-(SIGMA(2)-P)*SIGMA(5)*SIGMA(5)-
     2   (SIGMA(1)-P)*SIGMA(6)*SIGMA(6)-(SIGMA(3)-P)*SIGMA(4)
     3   *SIGMA(4)
      SJ4=SQRT(SJ2*SJ2*SJ2)
      THEDA1=3.0D0*SQRT(3.0D0)*SJ3/((SJ4+1.0D-50)*2.0D0)
      IF(THEDA1.GE.1.0D0) THEN
	    THEDA=0.0D0
      ELSEIF(THEDA1.LE.0.0D0) THEN
		THEDA=1.0D0
	  ELSE
        THEDA=ACOS(THEDA1)/3.0D0
	  ENDIF
C-------------------------------------------------------	
C  CALCULATE EQUIVALENT STRESS FPO and F1
C
      IF((ZERO.LE.THEDA).AND.(THEDA.LE.THEDA_B)) THEN
        FP0=ROOT3*SQRT(SJ2)*COS(THEDA)
      ELSEIF(THEDA_B.LT.THEDA.AND.THEDA.LE.(PAI/THREE)) THEN
        FP0=THREE*SQRT(SJ2)*SIN(THEDA)/TWO
     1     +ROOT3*SQRT(SJ2)*COS(THEDA)/TWO
      ENDIF
	  F1=FP0-SIGMA_C-ET*EQPLAS
C-------------------------------------------------------
      SI1=STRESS(1)+STRESS(2)+STRESS(3)
      P=SI1/3.0D0
      SJ2=((STRESS(1)-STRESS(2))*(STRESS(1)-STRESS(2))+(STRESS(2)-STRESS
     1   (3))*(STRESS(2)-STRESS(3))+(STRESS(3)-STRESS(1))*(STRESS(3)
     2   -STRESS(1)))/6.0D0+STRESS(4)*STRESS(4)+STRESS(5)*STRESS(5)+
     3   STRESS(6)*STRESS(6)
      SJ3=(STRESS(1)-P)*(STRESS(2)-P)*(STRESS(3)-P)+2.0D0*STRESS(4)
     1   *STRESS(5)*STRESS(6)-(STRESS(2)-P)*STRESS(5)*STRESS(5)-
     2   (STRESS(1)-P)*STRESS(6)*STRESS(6)-(STRESS(3)-P)*STRESS(4)
     3   *STRESS(4)
      SJ4=SQRT(SJ2*SJ2*SJ2)
      THEDA1=3.0D0*SQRT(3.0D0)*SJ3/((SJ4+1.0D-50)*2.0D0)
      IF(THEDA1.GE.1.0D0) THEN
	    THEDA=0.0D0
      ELSEIF(THEDA1.LE.0.0D0) THEN
		THEDA=1.0D0
	  ELSE
        THEDA=ACOS(THEDA1)/3.0D0
	  ENDIF
C-------------------------------------------------------	
C
      IF((ZERO.LE.THEDA).AND.(THEDA.LE.THEDA_B)) THEN
        F0=ROOT3*SQRT(SJ2)*COS(THEDA)-SIGMA_C-ET*EQPLAS
      ELSEIF(THEDA_B.LT.THEDA.AND.THEDA.LE.(PAI/THREE)) THEN
        F0=THREE*SQRT(SJ2)*SIN(THEDA)/TWO
     1     +ROOT3*SQRT(SJ2)*COS(THEDA)/TWO-SIGMA_C-ET*EQPLAS
      ENDIF
C----------------------------------------------------------------
C  CALCULATE FACTOR R BASED ON THE STRESS STATE 
C----------------------------------------------------------------
      IF(F1.LE.ZERO .OR. F1.LE.F0)  GOTO 900
	  IF(ZERO.LE.F0 .AND.F0.LE.F1) THEN
	    R=ZERO
      ELSE
	    R1=-F0/(F1-F0)
        DO K1=1,NTENS
	      SIGMA(K1)=STRESS(K1)
        ENDDO
        DO K1=1,NTENS
          SIGMA(K1)=SIGMA(K1)+R1*DSTRESS(K1)
        ENDDO
C-------------------------------------------------------
        SI1=SIGMA(1)+SIGMA(2)+SIGMA(3)
        P=SI1/3.0D0
        SJ2=((SIGMA(1)-SIGMA(2))*(SIGMA(1)-SIGMA(2))+(SIGMA(2)-SIGMA
     1     (3))*(SIGMA(2)-SIGMA(3))+(SIGMA(3)-SIGMA(1))*(SIGMA(3)
     2     -SIGMA(1)))/6.0D0+SIGMA(4)*SIGMA(4)+SIGMA(5)*SIGMA(5)+
     3     SIGMA(6)*SIGMA(6)
        SJ3=(SIGMA(1)-P)*(SIGMA(2)-P)*(SIGMA(3)-P)+2.0D0*SIGMA(4)
     1     *SIGMA(5)*SIGMA(6)-(SIGMA(2)-P)*SIGMA(5)*SIGMA(5)-
     2     (SIGMA(1)-P)*SIGMA(6)*SIGMA(6)-(SIGMA(3)-P)*SIGMA(4)
     3     *SIGMA(4)
        SJ4=SQRT(SJ2*SJ2*SJ2)
        THEDA1=3.0D0*SQRT(3.0D0)*SJ3/((SJ4+1.0D-50)*2.0D0)
        IF(THEDA1.GE.1.0D0) THEN
	      THEDA=0.0D0
		ELSEIF(THEDA1.LE.0.0D0) THEN
		  THEDA=1.0D0
	    ELSE
          THEDA=ACOS(THEDA1)/3.0D0
	    ENDIF
C-------------------------------------------------------
        IF((1*PAI/180 .LE.THEDA).AND.(THEDA.LT.THEDA_B)) THEN
          F2=SQRT(SJ2*3.0D0)*COS(THEDA)-SIGMA_C-ET*EQPLAS
          C1=0.0D0
          C2=ROOT3*COS(THEDA)-ROOT3*SIN(THEDA)/
     1      TAN(3.0D0*THEDA)
          C3=3.0D0*SIN(THEDA)/(2.0D0*SJ2*SIN(3.0D0*THEDA))
C
        ELSEIF(THEDA_B.LT.THEDA.AND.THEDA.LE.59*PAI/180) THEN
          F2=3.0D0*SQRT(SJ2)*SIN(THEDA)/2.0D0
     1     +ROOT3*SQRT(SJ2)*COS(THEDA)/2.0D0-SIGMA_C-ET*EQPLAS
          C1=0.0D0
          C2=ROOT3*COS(THEDA)/2.0D0+3.0D0*SIN(THEDA)/2.0D0+
     1      (3.0D0*COS(THEDA)-ROOT3*SIN(THEDA))/(2.0D0*TAN(3.0D0*THEDA))
          C3=3.0D0*(SIN(THEDA)-ROOT3*COS(THEDA))/(4.0D0*SJ2*
     1     SIN(3.0D0*THEDA))
C
        ELSEIF(0.0D0 .LE. THEDA .AND. THEDA .LT.1.0*PAI/180) THEN
          F2=SQRT(SJ2*3.0D0)*COS(THEDA)-SIGMA_C-ET*EQPLAS
C---THEDA=0,B=1
          C1=0.0D0
          C2=2.0D0/(ROOT3)
          C3=1.0D0/(2.0D0*SJ2)
        ELSEIF(59*PAI/180 .LE.THEDA.AND.THEDA.LT.PAI/3) THEN
          F2=3.0D0*SQRT(SJ2)*SIN(THEDA)/2.0D0
     1     +ROOT3*SQRT(SJ2)*COS(THEDA)/2.0D0-SIGMA_C-ET*EQPLAS
C---THEDA=PAI/3,B=1
          C1=0.0D0
          C2=2.0D0/(ROOT3)
          C3=-1.0D0/(2.0D0*SJ2)
C 
        ELSEIF(THEDA.EQ.THEDA_B) THEN
          F2=SQRT(SJ2*3.0D0)*COS(THEDA)-SIGMA_C-ET*EQPLAS
C---THEDA=THEDA_B
          C1=0.0D0
          C2=ROOT3*COS(THEDA)-(ROOT3*COS(3.0D0*THEDA)/
     1      (3.0D0*COS(THEDA)*COS(THEDA)-SIN(THEDA)*SIN(THEDA)))
          C3=3.0D0/(2.0D0*(SJ2+1.0D-50)*(3.0D0*COS(THEDA)
     1      *COS(THEDA)-SIN(THEDA)*SIN(THEDA)))
        ENDIF
C>>>>>>>>>>>>>>>>CALL FLOWVECTOR(A1,A2,A3,SIGMA,SJ2)>>>>>>
        A1(1)=1.0D0
        A1(2)=1.0D0
	    A1(3)=1.0D0
	    A1(4)=0.0D0
	    A1(5)=0.0D0
	    A1(6)=0.0D0
C	  
        A2(1)=(SIGMA(1)-P)/(2.0D0*SQRT(SJ2))
        A2(2)=(SIGMA(2)-P)/(2.0D0*SQRT(SJ2))		
        A2(3)=(SIGMA(3)-P)/(2.0D0*SQRT(SJ2))
        A2(4)=SIGMA(4)/SQRT(SJ2)
        A2(5)=SIGMA(5)/SQRT(SJ2)
        A2(6)=SIGMA(6)/SQRT(SJ2)
        A3(1)=(SIGMA(2)-P)*(SIGMA(3)-P)-SIGMA(6)*SIGMA(6)+SJ2/3.0D0
        A3(2)=(SIGMA(1)-P)*(SIGMA(3)-P)-SIGMA(5)*SIGMA(5)+SJ2/3.0D0
        A3(3)=(SIGMA(1)-P)*(SIGMA(2)-P)-SIGMA(4)*SIGMA(4)+SJ2/3.0D0
        A3(4)=2.0D0*(SIGMA(5)*SIGMA(6)-(SIGMA(3)-P)*SIGMA(4)) 
        A3(5)=2.0D0*(SIGMA(4)*SIGMA(6)-(SIGMA(2)-P)*SIGMA(5)) 
        A3(6)=2.0D0*(SIGMA(4)*SIGMA(5)-(SIGMA(1)-P)*SIGMA(6))
C-------
        DO K1=1,NTENS
          DF(K1)=C1*A1(K1)+C2*A2(K1)+C3*A3(K1)
        ENDDO
C    CALCULATE DF(I)*DSTRESS(I)
        ZH=0.0D0
        DO K1=1,NTENS
          ZH=ZH+DF(K1)*DSTRESS(K1)
        ENDDO
        DELTA_R1=F2/(ZH+1.0D-50)
        R=R1-DELTA_R1
      ENDIF
C
      DO K1=1,NTENS
        SIGMA(K1)=STRESS(K1)
      ENDDO
      DO K1=1,NTENS
        SIGMA(K1)=SIGMA(K1)+R*DSTRESS(K1)
      ENDDO
C
      M=14
      K=0
      DO WHILE(K.LT.(M-1))
        K=K+1
C>>>>>>>>>>>>>>>>>>>>>CALL INVARIANT>>>>>>>>>>>>
        SI1=SIGMA(1)+SIGMA(2)+SIGMA(3)
        P=SI1/3.0D0
        SJ2=((SIGMA(1)-SIGMA(2))*(SIGMA(1)-SIGMA(2))+(SIGMA(2)-SIGMA
     1     (3))*(SIGMA(2)-SIGMA(3))+(SIGMA(3)-SIGMA(1))*(SIGMA(3)
     2     -SIGMA(1)))/6.0D0+SIGMA(4)*SIGMA(4)+SIGMA(5)*SIGMA(5)+
     3     SIGMA(6)*SIGMA(6)
        SJ3=(SIGMA(1)-P)*(SIGMA(2)-P)*(SIGMA(3)-P)+2.0D0*SIGMA(4)
     1     *SIGMA(5)*SIGMA(6)-(SIGMA(2)-P)*SIGMA(5)*SIGMA(5)-
     2     (SIGMA(1)-P)*SIGMA(6)*SIGMA(6)-(SIGMA(3)-P)*SIGMA(4)
     3     *SIGMA(4)
        SJ4=SQRT(SJ2*SJ2*SJ2)
        THEDA1=3.0D0*SQRT(3.0D0)*SJ3/((SJ4+1.0D-50)*2.0D0)
        IF(THEDA1.GE.1.0D0) THEN
	      THEDA=0.0D0
		ELSEIF(THEDA1.LE.0.0D0) THEN
		  THEDA=1.0D0
	    ELSE
          THEDA=ACOS(THEDA1)/3.0D0
	    ENDIF
C<<------------------------------------------------
        IF((1*PAI/180 .LE.THEDA).AND.(THEDA.LT.THEDA_B)) THEN
          C1=0.0D0
          C2=ROOT3*COS(THEDA)-ROOT3*SIN(THEDA)/
     1      TAN(3.0D0*THEDA)
          C3=3.0D0*SIN(THEDA)/(2.0D0*SJ2*SIN(3.0D0*THEDA))
C
        ELSEIF(THEDA_B.LT.THEDA.AND.THEDA.LE.59*PAI/180) THEN
          C1=0.0D0
          C2=ROOT3*COS(THEDA)/2.0D0+3.0D0*SIN(THEDA)/2.0D0+
     1      (3.0D0*COS(THEDA)-ROOT3*SIN(THEDA))/(2.0D0*TAN(3.0D0*THEDA))
          C3=3.0D0*(SIN(THEDA)-ROOT3*COS(THEDA))/(4.0D0*SJ2*
     1     SIN(3.0D0*THEDA))
C
        ELSEIF(0.0D0 .LE. THEDA .AND. THEDA .LT.1.0*PAI/180) THEN
C---THEDA=0,B=1
          C1=0.0D0
          C2=2.0D0/(ROOT3)
          C3=1.0D0/(2.0D0*SJ2)
        ELSEIF(59*PAI/180 .LE.THEDA.AND.THEDA.LT.PAI/3) THEN
C---THEDA=PAI/3,B=1
          C1=0.0D0
          C2=2.0D0/(ROOT3)
          C3=-1.0D0/(2.0D0*SJ2)
C 
        ELSEIF(THEDA.EQ.THEDA_B) THEN
C---THEDA=THEDA_B
          C1=0.0D0
          C2=ROOT3*COS(THEDA)-(ROOT3*COS(3.0D0*THEDA)/
     1      (3.0D0*COS(THEDA)*COS(THEDA)-SIN(THEDA)*SIN(THEDA)))
          C3=3.0D0/(2.0D0*(SJ2+1.0D-50)*(3.0D0*COS(THEDA)
     1      *COS(THEDA)-SIN(THEDA)*SIN(THEDA)))
        ENDIF
C>>>>>>>>>>>>>>>>CALL FLOWVECTOR(A1,A2,A3,SIGMA,SJ2)>>>>>>
        A1(1)=1.0D0
        A1(2)=1.0D0
	    A1(3)=1.0D0
	    A1(4)=0.0D0
	    A1(5)=0.0D0
	    A1(6)=0.0D0
C
        A2(1)=(SIGMA(1)-P)/(2.0D0*SQRT(SJ2))
        A2(2)=(SIGMA(2)-P)/(2.0D0*SQRT(SJ2))		
        A2(3)=(SIGMA(3)-P)/(2.0D0*SQRT(SJ2))
        A2(4)=SIGMA(4)/SQRT(SJ2)
        A2(5)=SIGMA(5)/SQRT(SJ2)
        A2(6)=SIGMA(6)/SQRT(SJ2)
        A3(1)=(SIGMA(2)-P)*(SIGMA(3)-P)-SIGMA(6)*SIGMA(6)+SJ2/3.0D0
        A3(2)=(SIGMA(1)-P)*(SIGMA(3)-P)-SIGMA(5)*SIGMA(5)+SJ2/3.0D0
        A3(3)=(SIGMA(1)-P)*(SIGMA(2)-P)-SIGMA(4)*SIGMA(4)+SJ2/3.0D0
        A3(4)=2.0D0*(SIGMA(5)*SIGMA(6)-(SIGMA(3)-P)*SIGMA(4)) 
        A3(5)=2.0D0*(SIGMA(4)*SIGMA(6)-(SIGMA(2)-P)*SIGMA(5)) 
        A3(6)=2.0D0*(SIGMA(4)*SIGMA(5)-(SIGMA(1)-P)*SIGMA(6))
C<<<<
        DO K1=1,NTENS
          DF(K1)=C1*A1(K1)+C2*A2(K1)+C3*A3(K1)
        ENDDO
C>>>>>>>>>>CALL CLC_DP>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        FB=0.0D0
	    DO K1=1,6
	      FA(K1)=0.0D0
		  FC(K1)=0.0D0
		ENDDO
		DO K1=1,6
		  DO K2=1,6
	        FA(K1)=FA(K1)+DDSDDE(K1,K2)*DF(K2)
	        FD(K1,K2)=0.0D0
			D_P(K1,K2)=0.0D0
	      ENDDO
	    ENDDO
	    DO K1=1,6
	      FB=FB+DF(K1)*FA(K1)
	    ENDDO
C-
        HARD=EMOD*ET/(EMOD-ET)
		FB=FB+HARD
C-
	    DO K1=1,6
	      DO K2=1,6
		    FC(K1)=FC(K1)+DF(K2)*DDSDDE(K2,K1)
		  ENDDO
	    ENDDO
	    DO K1=1,6
		  DO K2=1,6
	        FD(K1,K2)=FA(K1)*FC(K2)
	      ENDDO
		ENDDO
C		
	    DO K1=1,6
		  DO K2=1,6
	        D_P(K1,K2)=FD(K1,K2)/FB
	      ENDDO
		ENDDO
C<<<<<-------------------------------------------------
        DO K1=1,NTENS
          SIGMA(K1)=SIGMA(K1)+(1-R)*DSTRESS(K1)/M
        ENDDO
        DO K1=1,NTENS
          DO K2=1,NTENS
            SIGMA(K1)=SIGMA(K1)-(1-R)*D_P(K1,K2)*DSTRAN(K2)/M
          ENDDO
        ENDDO
C>>>>>>>>>>>>>CALL STRAIN_P>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        SLMD=0.0D0
        SLMDS=0.0D0
        DO K1=1,6
          DO K2=1,6
            SLMDS=SLMDS+DF(K1)*DDSDDE(K1,K2)*(1-R)*DSTRAN(K2)/M
          ENDDO
        ENDDO
        SLMD=SLMD+SLMDS/(FB)
C		
	    DO K1=1,6
          EPLAS(K1)=EPLAS(K1)+SLMD*DF(K1)
          EELAS(K1)=EELAS(K1)-SLMD*DF(K1)
        ENDDO
C
        SI1=EPLAS(1)+EPLAS(2)+EPLAS(3)
        P=SI1/3.0D0
        SJ2=((EPLAS(1)-EPLAS(2))*(EPLAS(1)-EPLAS(2))+(EPLAS(2)-EPLAS
     1     (3))*(EPLAS(2)-EPLAS(3))+(EPLAS(3)-EPLAS(1))*(EPLAS(3)
     2     -EPLAS(1)))/6.0D0+EPLAS(4)*EPLAS(4)+EPLAS(5)*EPLAS(5)+
     3     EPLAS(6)*EPLAS(6)
        SJ3=(EPLAS(1)-P)*(EPLAS(2)-P)*(EPLAS(3)-P)+2.0D0*EPLAS(4)
     1     *EPLAS(5)*EPLAS(6)-(EPLAS(2)-P)*EPLAS(5)*EPLAS(5)-
     2     (EPLAS(1)-P)*EPLAS(6)*EPLAS(6)-(EPLAS(3)-P)*EPLAS(4)
     3     *EPLAS(4)
        SJ4=SQRT(SJ2*SJ2*SJ2)
C		WRITE(7,*)"(551) EPLAS",EPLAS
C
        IF(0.0D0 .LE.THEDA.AND.THEDA.LE.THEDA_B) THEN
          EQPLAS=SQRT(SJ2*THREE)*COS(THEDA)
        ELSEIF (THEDA_B.LT.THEDA.AND.THEDA.LE.1.047199) THEN
          EQPLAS=THREE*SQRT(SJ2)*SIN(THEDA)/TWO
     1     +ROOT3*SQRT(SJ2)*COS(THEDA)/TWO
        ENDIF
C       WRITE(7,*)"(585) EQPLAS",EQPLAS
C UPDATE STATEV VARIABLES
C
        DO K1=1,6
          STATEV(K1)=EELAS(K1)
          STATEV(K1+6)=EPLAS(K1)
        ENDDO
        STATEV(13)=EQPLAS
C<<<<<
C--------------CALL INVARIANT----------------------------
        SI1=SIGMA(1)+SIGMA(2)+SIGMA(3)
        P=SI1/3.0D0
        SJ2=((SIGMA(1)-SIGMA(2))*(SIGMA(1)-SIGMA(2))+(SIGMA(2)-SIGMA
     1    (3))*(SIGMA(2)-SIGMA(3))+(SIGMA(3)-SIGMA(1))*(SIGMA(3)
     2     -SIGMA(1)))/6.0D0+SIGMA(4)*SIGMA(4)+SIGMA(5)*SIGMA(5)+
     3     SIGMA(6)*SIGMA(6)
        SJ3=(SIGMA(1)-P)*(SIGMA(2)-P)*(SIGMA(3)-P)+2.0D0*SIGMA(4)
     1     *SIGMA(5)*SIGMA(6)-(SIGMA(2)-P)*SIGMA(5)*SIGMA(5)-
     2     (SIGMA(1)-P)*SIGMA(6)*SIGMA(6)-(SIGMA(3)-P)*SIGMA(4)
     3     *SIGMA(4)
        SJ4=SQRT(SJ2*SJ2*SJ2)
        THEDA1=3.0D0*SQRT(3.0D0)*SJ3/((SJ4+1.0D-50)*2.0D0)
        IF(THEDA1.GE.1.0D0) THEN
	      THEDA=0.0D0
		ELSEIF(THEDA1.LE.0.0D0) THEN
		  THEDA=1.0D0
	    ELSE
          THEDA=ACOS(THEDA1)/3.0D0
	    ENDIF
C---------------------------------------------------------------------
        IF((1*PAI/180 .LE.THEDA).AND.(THEDA.LT.THEDA_B)) THEN
          F3=SQRT(SJ2*3.0D0)*COS(THEDA)-SIGMA_C-ET*EQPLAS
          C1=0.0D0
          C2=ROOT3*COS(THEDA)-ROOT3*SIN(THEDA)/
     1      TAN(3.0D0*THEDA)
          C3=3.0D0*SIN(THEDA)/(2.0D0*SJ2*SIN(3.0D0*THEDA))
C
        ELSEIF(THEDA_B.LT.THEDA.AND.THEDA.LE.59*PAI/180) THEN
          F3=3.0D0*SQRT(SJ2)*SIN(THEDA)/2.0D0
     1     +ROOT3*SQRT(SJ2)*COS(THEDA)/2.0D0-SIGMA_C-ET*EQPLAS
          C1=0.0D0
          C2=ROOT3*COS(THEDA)/2.0D0+3.0D0*SIN(THEDA)/2.0D0+
     1      (3.0D0*COS(THEDA)-ROOT3*SIN(THEDA))/(2.0D0*TAN(3.0D0*THEDA))
          C3=3.0D0*(SIN(THEDA)-ROOT3*COS(THEDA))/(4.0D0*SJ2*
     1     SIN(3.0D0*THEDA))
C
        ELSEIF(0.0D0 .LE. THEDA .AND. THEDA .LT.1.0*PAI/180) THEN
          F3=SQRT(SJ2*3.0D0)*COS(THEDA)-SIGMA_C-ET*EQPLAS
C---THEDA=0,B=1
          C1=0.0D0
          C2=2.0D0/(ROOT3)
          C3=1.0D0/(2.0D0*SJ2)
        ELSEIF(59*PAI/180 .LE.THEDA.AND.THEDA.LT.PAI/3) THEN
          F3=3.0D0*SQRT(SJ2)*SIN(THEDA)/2.0D0
     1     +ROOT3*SQRT(SJ2)*COS(THEDA)/2.0D0-SIGMA_C-ET*EQPLAS
C---THEDA=PAI/3,B=1
          C1=0.0D0
          C2=2.0D0/(ROOT3)
          C3=-1.0D0/(2.0D0*SJ2)
C 
        ELSEIF(THEDA.EQ.THEDA_B) THEN
          F3=SQRT(SJ2*3.0D0)*COS(THEDA)-SIGMA_C-ET*EQPLAS
C---THEDA=THEDA_B
          C1=0.0D0
          C2=ROOT3*COS(THEDA)-(ROOT3*COS(3.0D0*THEDA)/
     1      (3.0D0*COS(THEDA)*COS(THEDA)-SIN(THEDA)*SIN(THEDA)))
          C3=3.0D0/(2.0D0*(SJ2+1.0D-50)*(3.0D0*COS(THEDA)
     1      *COS(THEDA)-SIN(THEDA)*SIN(THEDA)))
        ENDIF
C>>>>>>>>>>>>>>>>CALL FLOWVECTOR(A1,A2,A3,SIGMA,SJ2)>>>>>>
        A1(1)=1.0D0
        A1(2)=1.0D0
	    A1(3)=1.0D0
	    A1(4)=0.0D0
	    A1(5)=0.0D0
	    A1(6)=0.0D0
C	  
        A2(1)=(SIGMA(1)-P)/(2.0D0*SQRT(SJ2))
        A2(2)=(SIGMA(2)-P)/(2.0D0*SQRT(SJ2))		
        A2(3)=(SIGMA(3)-P)/(2.0D0*SQRT(SJ2))
        A2(4)=SIGMA(4)/SQRT(SJ2)
        A2(5)=SIGMA(5)/SQRT(SJ2)
        A2(6)=SIGMA(6)/SQRT(SJ2)
        A3(1)=(SIGMA(2)-P)*(SIGMA(3)-P)-SIGMA(6)*SIGMA(6)+SJ2/3.0D0
        A3(2)=(SIGMA(1)-P)*(SIGMA(3)-P)-SIGMA(5)*SIGMA(5)+SJ2/3.0D0
        A3(3)=(SIGMA(1)-P)*(SIGMA(2)-P)-SIGMA(4)*SIGMA(4)+SJ2/3.0D0
        A3(4)=2.0D0*(SIGMA(5)*SIGMA(6)-(SIGMA(3)-P)*SIGMA(4)) 
        A3(5)=2.0D0*(SIGMA(4)*SIGMA(6)-(SIGMA(2)-P)*SIGMA(5)) 
        A3(6)=2.0D0*(SIGMA(4)*SIGMA(5)-(SIGMA(1)-P)*SIGMA(6))
C<<<<<
        DO K1=1,NTENS
          DF(K1)=C1*A1(K1)+C2*A2(K1)+C3*A3(K1)
        ENDDO
C CALCULATE THE MODFYING STRESS INCREMENT XZS
        XZS1=0.0D0
        DO K1=1,NTENS
          XZS1=XZS1-DF(K1)*DF(K1)
        ENDDO
        DO K1=1,NTENS
          XZS(K1)=-F3*DF(K1)/XZS1
        ENDDO
C MODFY THE STRESS
        DO K1=1,6
          SIGMA(K1)=SIGMA(K1)-XZS(K1)
        ENDDO
      ENDDO
C-------------------------------------------------------
C-------------------------------------------------------
C-------------------------------------------------------
C--------------CALL INVARIANT----------------------------
      SI1=SIGMA(1)+SIGMA(2)+SIGMA(3)
      P=SI1/3.0D0
      SJ2=((SIGMA(1)-SIGMA(2))*(SIGMA(1)-SIGMA(2))+(SIGMA(2)-SIGMA
     1   (3))*(SIGMA(2)-SIGMA(3))+(SIGMA(3)-SIGMA(1))*(SIGMA(3)
     2   -SIGMA(1)))/6.0D0+SIGMA(4)*SIGMA(4)+SIGMA(5)*SIGMA(5)+
     3   SIGMA(6)*SIGMA(6)
      SJ3=(SIGMA(1)-P)*(SIGMA(2)-P)*(SIGMA(3)-P)+2.0D0*SIGMA(4)
     1   *SIGMA(5)*SIGMA(6)-(SIGMA(2)-P)*SIGMA(5)*SIGMA(5)-
     2   (SIGMA(1)-P)*SIGMA(6)*SIGMA(6)-(SIGMA(3)-P)*SIGMA(4)
     3   *SIGMA(4)
      SJ4=SQRT(SJ2*SJ2*SJ2)
      THEDA1=3.0D0*SQRT(3.0D0)*SJ3/((SJ4+1.0D-50)*2.0D0)
      IF(THEDA1.GE.1.0D0) THEN
	    THEDA=0.0D0
      ELSEIF(THEDA1.LE.0.0D0) THEN
		THEDA=1.0D0
	  ELSE
        THEDA=ACOS(THEDA1)/3.0D0
	  ENDIF
C----------------------------------------------------------
      IF((1*PAI/180 .LE.THEDA).AND.(THEDA.LT.THEDA_B)) THEN
        F3=SQRT(SJ2*3.0D0)*COS(THEDA)-SIGMA_C-ET*EQPLAS
        C1=0.0D0
        C2=ROOT3*COS(THEDA)-ROOT3*SIN(THEDA)/
     1    TAN(3.0D0*THEDA)
        C3=3.0D0*SIN(THEDA)/(2.0D0*SJ2*SIN(3.0D0*THEDA))
C
      ELSEIF(THEDA_B.LT.THEDA.AND.THEDA.LE.59*PAI/180) THEN
        F3=3.0D0*SQRT(SJ2)*SIN(THEDA)/2.0D0
     1   +ROOT3*SQRT(SJ2)*COS(THEDA)/2.0D0-SIGMA_C-ET*EQPLAS
        C1=0.0D0
        C2=ROOT3*COS(THEDA)/2.0D0+3.0D0*SIN(THEDA)/2.0D0+
     1    (3.0D0*COS(THEDA)-ROOT3*SIN(THEDA))/(2.0D0*TAN(3.0D0*THEDA))
        C3=3.0D0*(SIN(THEDA)-ROOT3*COS(THEDA))/(4.0D0*SJ2*
     1   SIN(3.0D0*THEDA))
C
      ELSEIF(0.0D0 .LE. THEDA .AND. THEDA .LT.1.0*PAI/180) THEN
        F3=SQRT(SJ2*3.0D0)*COS(THEDA)-SIGMA_C-ET*EQPLAS
C---THEDA=0,B=1
        C1=0.0D0
        C2=2.0D0/(ROOT3)
        C3=1.0D0/(2.0D0*SJ2)
      ELSEIF(59*PAI/180 .LE.THEDA.AND.THEDA.LT.PAI/3) THEN
        F3=3.0D0*SQRT(SJ2)*SIN(THEDA)/2.0D0
     1   +ROOT3*SQRT(SJ2)*COS(THEDA)/2.0D0-SIGMA_C-ET*EQPLAS
C---THEDA=PAI/3,B=1
        C1=0.0D0
        C2=2.0D0/(ROOT3)
        C3=-1.0D0/(2.0D0*SJ2)
C 
      ELSEIF(THEDA.EQ.THEDA_B) THEN
        F3=SQRT(SJ2*3.0D0)*COS(THEDA)-SIGMA_C-ET*EQPLAS
C---THEDA=THEDA_B
        C1=0.0D0
        C2=ROOT3*COS(THEDA)-(ROOT3*COS(3.0D0*THEDA)/
     1    (3.0D0*COS(THEDA)*COS(THEDA)-SIN(THEDA)*SIN(THEDA)))
        C3=3.0D0/(2.0D0*(SJ2+1.0D-50)*(3.0D0*COS(THEDA)
     1    *COS(THEDA)-SIN(THEDA)*SIN(THEDA)))
      ENDIF
C>>>>>>>>>>>>>>>>CALL FLOWVECTOR(A1,A2,A3,SIGMA,SJ2)>>>>>>
      A1(1)=1.0D0
      A1(2)=1.0D0
	  A1(3)=1.0D0
	  A1(4)=0.0D0
	  A1(5)=0.0D0
	  A1(6)=0.0D0
C	  
      A2(1)=(SIGMA(1)-P)/(2.0D0*SQRT(SJ2))
      A2(2)=(SIGMA(2)-P)/(2.0D0*SQRT(SJ2))		
      A2(3)=(SIGMA(3)-P)/(2.0D0*SQRT(SJ2))
      A2(4)=SIGMA(4)/SQRT(SJ2)
      A2(5)=SIGMA(5)/SQRT(SJ2)
      A2(6)=SIGMA(6)/SQRT(SJ2)
      A3(1)=(SIGMA(2)-P)*(SIGMA(3)-P)-SIGMA(6)*SIGMA(6)+SJ2/3.0D0
      A3(2)=(SIGMA(1)-P)*(SIGMA(3)-P)-SIGMA(5)*SIGMA(5)+SJ2/3.0D0
      A3(3)=(SIGMA(1)-P)*(SIGMA(2)-P)-SIGMA(4)*SIGMA(4)+SJ2/3.0D0
      A3(4)=2.0D0*(SIGMA(5)*SIGMA(6)-(SIGMA(3)-P)*SIGMA(4)) 
      A3(5)=2.0D0*(SIGMA(4)*SIGMA(6)-(SIGMA(2)-P)*SIGMA(5)) 
      A3(6)=2.0D0*(SIGMA(4)*SIGMA(5)-(SIGMA(1)-P)*SIGMA(6))
C<<<<
      DO K1=1,NTENS
        DF(K1)=C1*A1(K1)+C2*A2(K1)+C3*A3(K1)
      ENDDO
C	  
      XZS1=0.0D0
      DO K1=1,NTENS
        XZS1=XZS1-DF(K1)*DF(K1)
      ENDDO
      DO K1=1,NTENS
        XZS(K1)=-F3*DF(K1)/XZS1
      ENDDO
C  UPDATE THE FINAL STRESS THE LAST TIME
      DO K1=1,NTENS
        SIGMA(K1)=SIGMA(K1)-XZS(K1)
		STRESS(K1)=SIGMA(K1)
      ENDDO
C>>>>>>>>>>>>>CALL STRAIN_P>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        SLMD=0.0D0
        SLMDS=0.0D0
        DO K1=1,6
          DO K2=1,6
            SLMDS=SLMDS+DF(K1)*DDSDDE(K1,K2)*(1-R)*DSTRAN(K2)/M
          ENDDO
        ENDDO
        SLMD=SLMD+SLMDS/(FB)
C
	    DO K1=1,6
          EPLAS(K1)=EPLAS(K1)+SLMD*DF(K1)
          EELAS(K1)=EELAS(K1)-SLMD*DF(K1)
        ENDDO
C
        SI1=EPLAS(1)+EPLAS(2)+EPLAS(3)
        P=SI1/3.0D0
        SJ2=((EPLAS(1)-EPLAS(2))*(EPLAS(1)-EPLAS(2))+(EPLAS(2)-EPLAS
     1     (3))*(EPLAS(2)-EPLAS(3))+(EPLAS(3)-EPLAS(1))*(EPLAS(3)
     2     -EPLAS(1)))/6.0D0+EPLAS(4)*EPLAS(4)+EPLAS(5)*EPLAS(5)+
     3     EPLAS(6)*EPLAS(6)
        SJ3=(EPLAS(1)-P)*(EPLAS(2)-P)*(EPLAS(3)-P)+2.0D0*EPLAS(4)
     1     *EPLAS(5)*EPLAS(6)-(EPLAS(2)-P)*EPLAS(5)*EPLAS(5)-
     2     (EPLAS(1)-P)*EPLAS(6)*EPLAS(6)-(EPLAS(3)-P)*EPLAS(4)
     3     *EPLAS(4)
        SJ4=SQRT(SJ2*SJ2*SJ2)
C		WRITE(7,*)"(551) EPLAS",EPLAS
C
        IF(0.0D0 .LE.THEDA.AND.THEDA.LE.THEDA_B) THEN
          EQPLAS=SQRT(SJ2*THREE)*COS(THEDA)
        ELSEIF (THEDA_B.LT.THEDA.AND.THEDA.LE.1.047199) THEN
          EQPLAS=THREE*SQRT(SJ2)*SIN(THEDA)/TWO
     1     +ROOT3*SQRT(SJ2)*COS(THEDA)/TWO
        ENDIF
C
C UPDATE STATEV VARIABLES
C
        DO K1=1,6
          STATEV(K1)=EELAS(K1)
          STATEV(K1+6)=EPLAS(K1)
        ENDDO
        STATEV(13)=EQPLAS
C<<<<<
C--------------CALL INVARIANT----------------------------
      SI1=SIGMA(1)+SIGMA(2)+SIGMA(3)
      P=SI1/3.0D0
      SJ2=((SIGMA(1)-SIGMA(2))*(SIGMA(1)-SIGMA(2))+(SIGMA(2)-SIGMA
     1   (3))*(SIGMA(2)-SIGMA(3))+(SIGMA(3)-SIGMA(1))*(SIGMA(3)
     2   -SIGMA(1)))/6.0D0+SIGMA(4)*SIGMA(4)+SIGMA(5)*SIGMA(5)+
     3   SIGMA(6)*SIGMA(6)
      SJ3=(SIGMA(1)-P)*(SIGMA(2)-P)*(SIGMA(3)-P)+2.0D0*SIGMA(4)
     1   *SIGMA(5)*SIGMA(6)-(SIGMA(2)-P)*SIGMA(5)*SIGMA(5)-
     2   (SIGMA(1)-P)*SIGMA(6)*SIGMA(6)-(SIGMA(3)-P)*SIGMA(4)
     3   *SIGMA(4)
      SJ4=SQRT(SJ2*SJ2*SJ2)
      THEDA1=3.0D0*SQRT(3.0D0)*SJ3/((SJ4+1.0D-50)*2.0D0)
      IF(THEDA1.GE.1.0D0) THEN
	    THEDA=0.0D0
      ELSEIF(THEDA1.LE.0.0D0) THEN
		THEDA=1.0D0
	  ELSE
        THEDA=ACOS(THEDA1)/3.0D0
	  ENDIF
C-----------------------------------------------------------
      IF((1*PAI/180 .LE.THEDA).AND.(THEDA.LT.THEDA_B)) THEN
        C1=0.0D0
        C2=ROOT3*COS(THEDA)-ROOT3*SIN(THEDA)/
     1    TAN(3.0D0*THEDA)
        C3=3.0D0*SIN(THEDA)/(2.0D0*SJ2*SIN(3.0D0*THEDA))
C
      ELSEIF(THEDA_B.LT.THEDA.AND.THEDA.LE.59*PAI/180) THEN
        C1=0.0D0
        C2=ROOT3*COS(THEDA)/2.0D0+3.0D0*SIN(THEDA)/2.0D0+
     1    (3.0D0*COS(THEDA)-ROOT3*SIN(THEDA))/(2.0D0*TAN(3.0D0*THEDA))
        C3=3.0D0*(SIN(THEDA)-ROOT3*COS(THEDA))/(4.0D0*SJ2*
     1   SIN(3.0D0*THEDA))
C
      ELSEIF(0.0D0 .LE. THEDA .AND. THEDA .LT.1.0*PAI/180) THEN
C---THEDA=0,B=1
        C1=0.0D0
        C2=2.0D0/(ROOT3)
        C3=1.0D0/(2.0D0*SJ2)
      ELSEIF(59*PAI/180 .LE.THEDA.AND.THEDA.LT.PAI/3) THEN
C---THEDA=PAI/3,B=1
        C1=0.0D0
        C2=2.0D0/(ROOT3)
        C3=-1.0D0/(2.0D0*SJ2)
C 
      ELSEIF(THEDA.EQ.THEDA_B) THEN
C---THEDA=THEDA_B
        C1=0.0D0
        C2=ROOT3*COS(THEDA)-(ROOT3*COS(3.0D0*THEDA)/
     1    (3.0D0*COS(THEDA)*COS(THEDA)-SIN(THEDA)*SIN(THEDA)))
        C3=3.0D0/(2.0D0*(SJ2+1.0D-50)*(3.0D0*COS(THEDA)
     1    *COS(THEDA)-SIN(THEDA)*SIN(THEDA)))
      ENDIF
C>>>>>>>>>>>>>>>>CALL FLOWVECTOR(A1,A2,A3,SIGMA,SJ2)>>>>>>
      A1(1)=1.0D0
      A1(2)=1.0D0
	  A1(3)=1.0D0
	  A1(4)=0.0D0
	  A1(5)=0.0D0
	  A1(6)=0.0D0
C	  
      A2(1)=(SIGMA(1)-P)/(2.0D0*SQRT(SJ2))
      A2(2)=(SIGMA(2)-P)/(2.0D0*SQRT(SJ2))		
      A2(3)=(SIGMA(3)-P)/(2.0D0*SQRT(SJ2))
      A2(4)=SIGMA(4)/SQRT(SJ2)
      A2(5)=SIGMA(5)/SQRT(SJ2)
      A2(6)=SIGMA(6)/SQRT(SJ2)
      A3(1)=(SIGMA(2)-P)*(SIGMA(3)-P)-SIGMA(6)*SIGMA(6)+SJ2/3.0D0
      A3(2)=(SIGMA(1)-P)*(SIGMA(3)-P)-SIGMA(5)*SIGMA(5)+SJ2/3.0D0
      A3(3)=(SIGMA(1)-P)*(SIGMA(2)-P)-SIGMA(4)*SIGMA(4)+SJ2/3.0D0
      A3(4)=2.0D0*(SIGMA(5)*SIGMA(6)-(SIGMA(3)-P)*SIGMA(4)) 
      A3(5)=2.0D0*(SIGMA(4)*SIGMA(6)-(SIGMA(2)-P)*SIGMA(5)) 
      A3(6)=2.0D0*(SIGMA(4)*SIGMA(5)-(SIGMA(1)-P)*SIGMA(6))
C<<<<<<-------------------------------------------------------
      DO K1=1,NTENS
        DF(K1)=C1*A1(K1)+C2*A2(K1)+C3*A3(K1)
      ENDDO
C>>>>>>>>>>CALL CLC_DP>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      FB=0.0D0
	  DO K1=1,6
	    FA(K1)=0.0D0
		FC(K1)=0.0D0
	  ENDDO
	  DO K1=1,6
	    DO K2=1,6
		  FA(K1)=FA(K1)+DDSDDE(K1,K2)*DF(K2)
		  FD(K1,K2)=0.0D0
		  D_P(K1,K2)=0.0D0
		ENDDO
	  ENDDO
	  DO K1=1,6
	    FB=FB+DF(K1)*FA(K1)
	  ENDDO
C-
        HARD=EMOD*ET/(EMOD-ET)
		FB=FB+HARD
C-
	  DO K1=1,6
	    DO K2=1,6
		  FC(K1)=FC(K1)+DF(K2)*DDSDDE(K2,K1)
		ENDDO
	  ENDDO
	  DO K1=1,6
	    DO K2=1,6
	      FD(K1,K2)=FA(K1)*FC(K2)
	    ENDDO
	  ENDDO
	  DO K1=1,6
	    DO K2=1,6
	    D_P(K1,K2)=FD(K1,K2)/FB
	    ENDDO
	  ENDDO
C<<<<<<
      DO K1=1,NTENS
        DO K2=1,NTENS
          DDSDDE(K1,K2)=DDSDDE(K1,K2)-D_P(K1,K2)
        ENDDO
      ENDDO
C 
      GOTO 2000
 900  CONTINUE
      DO K1=1,NTENS
        STRESS(K1)=SIGMA(K1)
      ENDDO
2000  CONTINUE
C      WRITE(7,*)"END STRESS=",STRESS
      RETURN
      END