************************************************************************
*                                                                      *
*  SUBROTINA QUE CALCULA AS SOLUÇÕES FUNDAMENTAIS DINÂMICAS 3D         *
*                                                                      *
************************************************************************
*
      SUBROUTINE SOLFUND(ZU,ZT,CXP,CYP,CZP,CXG,CYG,CZG,RN,ZGE,ZCS,ZCP,
     $  DELTA,PI,FR)
*
        IMPLICIT REAL (A-H,O-Y)
        IMPLICIT COMPLEX (Z)
        COMMON INP,INQ,IPR,IPS,IPT
        DIMENSION ZU(3,3),ZT(3,3)
        DIMENSION RD(3)
        DIMENSION RN(3)
        DIMENSION DELTA(3,3)
*
        R1=CXG-CXP
        R2=CYG-CYP
        R3=CZG-CZP
        R=SQRT(R1**2+R2**2+R3**2)
        DRN=(R1*RN(1)+R2*RN(2)+R3*RN(3))/R
        RD(1)=R1/R
        RD(2)=R2/R
        RD(3)=R3/R
*
* CONSTANTES USADAS NAS SOLUÇÕES FUNDAMENTAIS DINÂMICAS
*
        ZWI=CMPLX(0.0,FR)
        ZC0=(1.0,0.0)/(CMPLX(4.0*PI,0.0)*ZGE)
        ZC1=(ZCP/ZCS)**2
        ZC2=(ZCS/ZCP)**2
        ZKP=-ZWI/ZCP
        ZKS=-ZWI/ZCS
        ZZP=ZKP*R
        ZZS=ZKS*R
        ZEZP=CEXP(ZZP)
        ZEZS=CEXP(ZZS)
        ZP2=ZZP*ZZP
        ZS2=ZZS*ZZS
*
        ZFHI=(1.0+1.0/ZS2-1.0/ZZS)*ZEZS/R-
     @      ZC2*(1.0/ZP2-1.0/ZZP)*ZEZP/R
        ZCAPPA=(1.0+3.0/ZS2-3.0/ZZS)*ZEZS/R-
     @      ZC2*(1.0+3.0/ZP2-3.0/ZZP)*ZEZP/R
        ZFHIDR=(-2.0+ZZS+3.0/ZZS-3.0/ZS2)*ZEZS/R**2-
     @      ZC2*(-1.0+3.0/ZZP-3.0/ZP2)*ZEZP/R**2
        ZCAPPADR=(ZZS-4.0+9.0/ZZS-9.0/ZS2)*ZEZS/R**2-
     @      ZC2*(ZZP-4.0+9.0/ZZP-9.0/ZP2)*ZEZP/R**2
*
        ZAA=ZFHIDR-ZCAPPA/R
        ZBB=4.0*ZCAPPA/R-2.0*ZCAPPADR
        ZCC=(ZC1-2.0)*(ZAA+5.E-1*ZBB-3.0*ZCAPPA/R)-2.0*ZCAPPA/R      
*
        DO 40 J=1,3
            DO 40 I=1,3
                ZU(I,J)=ZC0*(ZFHI*DELTA(I,J)-ZCAPPA*RD(J)*RD(I))
           ZT(I,J)=(1.0/(4.0*PI))*((ZAA*(DRN*DELTA(I,J)+RD(J)*RN(I)))+
     $              RD(I)*RD(J)*DRN*ZBB+RD(I)*RN(J)*ZCC)
 40     CONTINUE
*
        RETURN
      END
