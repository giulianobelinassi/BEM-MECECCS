************************************************************************
*                                                                      *
*  SUBROTINA QUE CALCULA AS SOLUÇÕES FUNDAMENTAIS ESTÁTICAS 3D		 *
*																	 *
************************************************************************
*
      SUBROUTINE SOLFUNE(T,CXP,CYP,CZP,CXG,CYG,CZG,ETA,
     $  N,NP,NPG,RNU,RMU,C3,C4,DELTA)
*
        IMPLICIT REAL (A-H,O-Y)
        IMPLICIT COMPLEX (Z)
        COMMON INP,INQ,IPR,IPS,IPT
        DIMENSION T(3,3)
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
                T(I,J)=(C3/(R*R))*(RDN*(C4*DELTA(I,J)+3.0*RD(I)*RD(J))+
     $              C4*(RD(J)*ETA(I)-RD(I)*ETA(J)))
  10    CONTINUE
*
        RETURN
      END
