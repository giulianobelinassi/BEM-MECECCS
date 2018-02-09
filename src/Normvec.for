! Subrotina que calcula o produto vetorial usado originalmente em Ghmatecd.for,
! Ghmatece.for e Interec.for, evitando cálculos redundantes.

      FUNCTION NORMVEC(cone, cx, cy, cz, n, np) RESULT(ETAS)
        IMPLICIT NONE
        
        INTEGER, DIMENSION(n,4), INTENT(IN) :: CONE
        REAL, DIMENSION(np), INTENT(IN) :: CX, CY, CZ
        INTEGER, INTENT(IN) :: n, np
        REAL, DIMENSION(:,:), ALLOCATABLE:: ETAS

        INTEGER J, N1,N2,N3,N4, stats
        REAL A, B, C, R

        ALLOCATE(ETAS(3, n), STAT = stats)
        IF (stats /= 0) THEN
            PRINT*, "MEMÓRIA INSUFICIENTE!"
        ENDIF

!$OMP  PARALLEL DO DEFAULT(SHARED)
!$OMP& PRIVATE(N1,N2,N3,N4,J,A,B,C,R)
        DO J=1, N
            N1=CONE(J,1)
            N2=CONE(J,2)
            N3=CONE(J,3)
            N4=CONE(J,4)
            A=(CY(N2)-CY(N1))*(CZ(N3)-CZ(N1)) - 
     $           (CZ(N2)-CZ(N1))*(CY(N3)-CY(N1))
            B=(CZ(N2)-CZ(N1))*(CX(N3)-CX(N1)) - 
     $           (CX(N2)-CX(N1))*(CZ(N3)-CZ(N1))
            C=(CX(N2)-CX(N1))*(CY(N3)-CY(N1)) - 
     $           (CY(N2)-CY(N1))*(CX(N3)-CX(N1))
            R=SQRT(A*A+B*B+C*C)
            ETAS(1,j)=A/R
            ETAS(2,j)=B/R
            ETAS(3,j)=C/R
        ENDDO
!$OMP END PARALLEL DO
      END FUNCTION NORMVEC
