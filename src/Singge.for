************************************************************************
*                                                                      *
*    SUBROTINA QUE CALCULA OS COEFICIENTES DA DIAGONAL DE [G] POR      *
*																	   *	
*  INTEGRAÇÃO NUMÉRICA. O PONTO DE COLOCAÇÃO ESTÁ NO ELEMENTO A SER	   *
*                                                                      *
*                       INTEGRADO                                      *
*                                                                      *
************************************************************************
*
      SUBROUTINE SINGGE(GELEM,CXM,CYM,CZM,ETA,CX,CY,CZ,
     $  N1,N2,N3,N4,NCOX,N,NP,NPG,C1,C2,DELTA,GI,OME)
*
        IMPLICIT REAL (A-H,O-Y)
        IMPLICIT COMPLEX (Z)
        COMMON INP,INQ,IPR,IPS,IPT
        DIMENSION DELTA(3,3)
        DIMENSION HELEM(3,3),GELEM(3,3),GP(3,3)
        DIMENSION CX(NCOX),CY(NCOX),CZ(NCOX)
        DIMENSION CO(4,3),ETA(3)
        DIMENSION GI(NPG), OME(NPG)
*
* INTEGRAÇÃO NUMÉRICA SOBRE TRIÂNGULOS INTERNOS COLAPSANDO OS CANTOS 3-4
* DE UM ELEMENTO QUADRILATERAL
*
        CO(3,1)=CXM
        CO(3,2)=CYM
        CO(3,3)=CZM
        CO(4,1)=CXM
        CO(4,2)=CYM
        CO(4,3)=CZM
*


        HELEM = 0.D0
        GELEM = 0.D0
        
        DO 500 IT=1,4

*           SWITCH IT
            GO TO (20,30,40,50),IT
  20            CO(1,1)=CX(N1)
                CO(1,2)=CY(N1)
                CO(1,3)=CZ(N1)
                CO(2,1)=CX(N2)
                CO(2,2)=CY(N2)
                CO(2,3)=CZ(N2)
            GO TO 60
  30            CO(1,1)=CX(N2)
                CO(1,2)=CY(N2)
                CO(1,3)=CZ(N2)
                CO(2,1)=CX(N3)
                CO(2,2)=CY(N3)
                CO(2,3)=CZ(N3)
            GO TO 60
  40            CO(1,1)=CX(N3)
                CO(1,2)=CY(N3)
                CO(1,3)=CZ(N3)
                CO(2,1)=CX(N4)
                CO(2,2)=CY(N4)
                CO(2,3)=CZ(N4)
            GO TO 60
  50            CO(1,1)=CX(N4)
                CO(1,2)=CY(N4)
                CO(1,3)=CZ(N4)
                CO(2,1)=CX(N1)
                CO(2,2)=CY(N1)
                CO(2,3)=CZ(N1)
   60       CONTINUE
*
* INTEGRAÇÃO SOBRE O ELEMENTO QUADRILATERAL COLAPSADO
*
            CALL NONSINGE_G(GP,CO,CXM,CYM,CZM,ETA,N,NP,NPG,
     $          C1,C2,DELTA,GI,OME)
*
            GELEM = GELEM + GP
*
  500   CONTINUE
*
        RETURN
      END

      SUBROUTINE NONSINGE_G(GELEM,CO,CXP,CYP,CZP,ETA,
     $   N,NP,NPG,C1,C2,DELTA,GI,OME)
        
        IMPLICIT REAL (A-H,O-Y)
        IMPLICIT COMPLEX (Z)
        COMMON INP,INQ,IPR,IPS,IPT
        DIMENSION DELTA(3,3)
        DIMENSION GELEM(3,3),U(3,3)
        DIMENSION GI(NPG),OME(NPG)
        DIMENSION CO(4,3),ETA(3),P(2,4),XJ(2,3),F(4)

        GELEM = 0.0

        DO 400 JG=1,NPG
*
            G2=GI(JG)
            P2=OME(JG)
            SP=1.0+G2
            SM=1.0-G2
            P(1,1)=-0.250*SM
            P(1,2)= 0.250*SM
            P(1,3)= 0.250*SP
            P(1,4)=-0.250*SP
*
            DO 300 IG=1,NPG
*
                G1=GI(IG)
                P1=OME(IG)
                RP=1.0+G1
                RM=1.0-G1
                F(1)=0.250*RM*SM
                F(2)=0.250*RP*SM
                F(3)=0.250*RP*SP
                F(4)=0.250*RM*SP
                P(2,1)=-0.250*RM
                P(2,2)=-0.250*RP
                P(2,3)= 0.250*RP
                P(2,4)= 0.250*RM
*
* CALCULA A RELAÇÃO ENTRE AS COORDENADAS CARTESIANAS E HOMOGÊNEAS
*
                DO 90 I=1,2
                    DO 80 J=1,3
                        TEMP=0.0
                        DO 70 K=1,4
                            TEMP=TEMP+P(I,K)*CO(K,J)
  70                    CONTINUE
                        XJ(I,J)=TEMP
  80                CONTINUE
  90            CONTINUE
*
* CALCULA O JACOBIANO
*
                DET=SQRT((XJ(1,2)*XJ(2,3)-XJ(2,2)*XJ(1,3))**2 +
     $              (XJ(2,1)*XJ(1,3)-XJ(1,1)*XJ(2,3))**2 +
     $              (XJ(1,1)*XJ(2,2)-XJ(2,1)*XJ(1,2))**2  )
                IF (DET < 1.0D-5) THEN

                    WRITE(IPR,1000) DET
                    WRITE(IPR,1100) ((CO(I,J),J=1,3),I=1,4)
 1000   FORMAT(///' NONSING : ERRO, JACOBIANO NULO OU NEGATIVO =',D14.5)
 1100            FORMAT(/1X,'COORDINADAS DOS PONTOS EXTREMOS '/(3D14.5))
                    STOP ' NONSING : ERRO, JACOBIANO NULO OU NEGATIVO'

                ENDIF
* CALCULA AS COORDENADAS DO PONTO DE INTEGRAÇÃO
*
                CXG=0.0
                CYG=0.0
                CZG=0.0
                DO 130 I=1,4
                    CXG=CXG+CO(I,1)*F(I)
                    CYG=CYG+CO(I,2)*F(I)
                    CZG=CZG+CO(I,3)*F(I)
 130            CONTINUE
*
* ACIONA ROTINA QUE CALCULA A SOLUÇÃO FUNDAMENTAL ESTÁTICA 3D
*
                CALL SOLFUNE_G (U,CXP,CYP,CZP,CXG,CYG,CZG,ETA,
     $              N,NP,NPG,C1,C2,DELTA)
*
                P12=P1*P2*DET

                GELEM = GELEM + U*P12
                
*
 300        CONTINUE
 400    CONTINUE
*
        RETURN
      END SUBROUTINE NONSINGE_G

      SUBROUTINE SOLFUNE_G(U,CXP,CYP,CZP,CXG,CYG,CZG,ETA,
     $  N,NP,NPG,C1,C2,DELTA)

        IMPLICIT REAL (A-H,O-Y)
        IMPLICIT COMPLEX (Z)
        COMMON INP,INQ,IPR,IPS,IPT
        DIMENSION U(3,3)
        DIMENSION RD(3)
        DIMENSION ETA(3)
        DIMENSION DELTA(3,3)
*
        R1=CXG-CXP
        R2=CYG-CYP
        R3=CZG-CZP
        R=SQRT(R1**2+R2**2+R3**2)
        RDN=(R1*ETA(1)+R2*ETA(2)+R3*ETA(3))/R
        RD(1)=R1/R
        RD(2)=R2/R
        RD(3)=R3/R
*
        DO 10 J=1,3
             DO 10 I=1,3
                U(I,J)=(C1/R)*(C2*DELTA(I,J)+RD(I)*RD(J))
  10    CONTINUE
*
        RETURN

      END SUBROUTINE SOLFUNE_G
