************************************************************************
*                                                                      *
*  SUBROTINA QUE CALCULA OS COEFICIENTES DE [H] E [G]. O PONTO DE      *
*																	   *
*    COLOCAÇÃO ESTÁ FORA DO CONTORNO, ASSIM NÃO HÁ SINGULARIDADE.      *
*                                                                      *
*         A INTEGRAÇÃO É FEITA USANDO QUADRATURA GAUSSIANA             *
*                                                                      *
************************************************************************
*

      SUBROUTINE NONSINGE(HELEM,CO,CXP,CYP,CZP,ETA,
     $  N,NP,NPG,RNU,RMU,C3,C4,DELTA, GI, OME)
*
        IMPLICIT REAL (A-H,O-Y)
        IMPLICIT COMPLEX (Z)
        COMMON INP,INQ,IPR,IPS,IPT
        DIMENSION DELTA(3,3)
        DIMENSION HELEM(3,3),T(3,3)
        DIMENSION GI(NPG),OME(NPG)
        DIMENSION CO(4,3),ETA(3),P(2,4),XJ(2,3),F(4)

        INTERFACE
            SUBROUTINE SOLFUNE(T,CXP,CYP,CZP,CXG,CYG,CZG,ETA,
     $          N,NP,NPG,RNU,RMU,C3,C4,DELTA)
                IMPLICIT REAL (A-H,O-Y)
                IMPLICIT COMPLEX (Z)
                COMMON INP,INQ,IPR,IPS,IPT
                DIMENSION T(3,3)
                DIMENSION RD(3)
                DIMENSION ETA(3)
                DIMENSION DELTA(3,3)
            END SUBROUTINE SOLFUNE
        END INTERFACE




*
* ZERA AS MATRIZES ELEMENTARES HELEM E GELEM
*
        HELEM = 0.0
*
*
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
                CALL SOLFUNE (T,CXP,CYP,CZP,CXG,CYG,CZG,ETA,
     $              N,NP,NPG,RNU,RMU,C3,C4,DELTA)
*
                P12=P1*P2*DET

                HELEM = HELEM + T*P12
*
 300        CONTINUE
 400    CONTINUE
*
        RETURN
      END
