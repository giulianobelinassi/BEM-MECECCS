************************************************************************
*                                                                      *
* SUBROTINA QUE CONTÉM OS PONTOS DA QUADRATURA GAUSSIANA E PONDERAÇÕES *
*                                                                      *
*                         (N PONTOS)                                   *
*                                                                      *
************************************************************************
*
      SUBROUTINE GAULEG (X1,X2,X,W,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL X1,X2,X(N),W(N)
      PARAMETER (EPS=3.D-14)
      M=(N+1)/2
      XM=0.50*(X2+X1)
      XL=0.50*(X2-X1)
      DO 12 I=1,M
        Z=DCOS(3.141592654D0*(I-.25)/(N+.5))
1       CONTINUE
          P1=1.0
          P2=0.0
          DO 11 J=1,N
            P3=P2
            P2=P1
            P1=((2.0*J-1.0)*Z*P2-(J-1.0)*P3)/J
11        CONTINUE
          PP=N*(Z*P1-P2)/(Z*Z-1.0)
          Z1=Z
          Z=Z1-P1/PP
        IF(DABS(Z-Z1).GT.EPS)GO TO 1
        X(I)=REAL(XM-XL*Z)
        X(N+1-I)=REAL(XM+XL*Z)
        W(I)=REAL(2.0*XL/((1.-Z*Z)*PP*PP))
        W(N+1-I)=W(I)
12    CONTINUE
      RETURN
      END
