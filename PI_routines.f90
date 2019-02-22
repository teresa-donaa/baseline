MODULE PI_routines
!
USE globals
USE QL_routines
!
! Various routines used to compute PI matrices at runtime
!
IMPLICIT NONE
!
CONTAINS
! 
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE computePIMatricesGiven ( DemandParameters, NashPrices, CoopPrices, &
        PI, NashProfits, CoopProfits, &
        indexNashPrices, indexCoopPrices, NashMarketShares, CoopMarketShares, &
        PricesGrids )
    !
    ! Computes the total common payoff matrix PI
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: DemandParameters(numDemandParameters)
    REAL(8), DIMENSION(numAgents), INTENT(IN) :: NashPrices, CoopPrices
    REAL(8), INTENT(OUT) :: PI(numActions,numAgents)
    REAL(8), DIMENSION(numAgents), INTENT(OUT) :: NashProfits, CoopProfits, &
        NashMarketShares, CoopMarketShares
    INTEGER, DIMENSION(numAgents), INTENT(OUT) :: indexNashPrices, indexCoopPrices
    REAL(8), INTENT(OUT) :: PricesGrids(numPrices,numAgents)
    !
    ! Declaring local variables
    !
    INTEGER :: iPrice, jPrice, h1, h2, i
    !
    ! Beginning execution
    !
    ! Computing PI matrices for 2 agents;
    ! with more than 2 agents there is no clear definition of "equivalent" state
    ! hence the numAgents > 2 case is ruled out here
    !
    IF (numAgents .GT. 2) THEN
        !
        PRINT*, 'Problem in computePIMatricesGiven:'
        PRINT*, 'With more than 2 agents there is no clear definition of "equivalent" state'
        PRINT*, 'Hence the numAgents > 2 case is ruled out here'
        STOP
        !
    END IF
    !
    PI(:,1) = DemandParameters
    DO iPrice = 1, numPrices
        !
        DO jPrice = 1, numPrices
            !
            h1 = (iPrice-1)*numPrices+jPrice
            h2 = (jPrice-1)*numPrices+iPrice
            PI(h1,2) = PI(h2,1)
            !
        END DO
        !
    END DO
    indexNashPrices = 1
    indexCoopPrices = numActions
    NashProfits = PI(1,:)
    CoopProfits = PI(numActions,:)
    NashMarketShares = 0.d0
    CoopMarketShares = 0.d0
    PricesGrids = SPREAD((/ ( DBLE(i), i = 1, numPrices ) /),DIM = 2,NCOPIES = numAgents)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE computePIMatricesGiven
! 
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE computePIMatricesSinghVives ( DemandParameters, NashPrices, CoopPrices, &
        PI, NashProfits, CoopProfits, &
        indexNashPrices, indexCoopPrices, NashMarketShares, CoopMarketShares, &
        PricesGrids )
    !
    ! Computes the Singh&Vives common payoff matrix PI
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: DemandParameters(numDemandParameters)
    REAL(8), DIMENSION(numAgents), INTENT(IN) :: NashPrices, CoopPrices
    REAL(8), INTENT(OUT) :: PI(numActions,numAgents)
    REAL(8), DIMENSION(numAgents), INTENT(OUT) :: NashProfits, CoopProfits, &
        NashMarketShares, CoopMarketShares
    INTEGER, DIMENSION(numAgents), INTENT(OUT) :: indexNashPrices, indexCoopPrices
    REAL(8), INTENT(OUT) :: PricesGrids(numPrices,numAgents)
    !
    ! Declaring local variables
    !
    REAL(8) :: gamma
    !
    ! Beginning execution
    !
    ! Extract demand parameters
    !
    gamma = DemandParameters(1)
    !
    ! Computing PI matrices for numAgents == 2 and numAgents == 3
    !
    IF (numPrices .EQ. 2) THEN
        !
        PI(1,1) = 0.d0
        PI(2,1) = -(1.d0-gamma)/(1.d0-2.d0*gamma)
        PI(3,1) = 2.d0
        PI(4,1) = 1.d0
        !
        PI(1,2) = 0.d0
        PI(2,2) = 2.d0
        PI(3,2) = -(1.d0-gamma)/(1.d0-2.d0*gamma)
        PI(4,2) = 1.d0
        !
        PricesGrids(1,:) = (1.d0-2.d0*gamma)/(2.d0-3.d0*gamma)
        PricesGrids(2,:) = 0.5d0
        !
    ELSE IF (numPrices .EQ. 3) THEN
        !
        PI(1,1) = 0.d0
        PI(2,1) = -gamma**2/(4.d0-12.d0*gamma+8.d0*gamma**2)
        PI(3,1) = -(1.d0-gamma)/(1.d0-2.d0*gamma)
        PI(4,1) = gamma/(1.d0-gamma)
        PI(5,1) = gamma*(4.d0-5.d0*gamma)/(4.d0*(1.d0-gamma)**2)
        PI(6,1) = (-2.d0+6.d0*gamma-5.d0*gamma**2)/(2.d0-6.d0*gamma+4.d0*gamma**2)
        PI(7,1) = 2.d0
        PI(8,1) = (8.d0-24.d0*gamma+17.d0*gamma**2)/(4.d0-12.d0*gamma+8.d0*gamma**2)
        PI(9,1) = 1.d0
        !
        PI(1,2) = 0.d0
        PI(2,2) = gamma/(1.d0-gamma)
        PI(3,2) = 2.d0
        PI(4,2) = -gamma**2/(4.d0-12.d0*gamma+8.d0*gamma**2)
        PI(5,2) = gamma*(4.d0-5.d0*gamma)/(4.d0*(1.d0-gamma)**2)
        PI(6,2) = (8.d0-24.d0*gamma+17.d0*gamma**2)/(4.d0-12.d0*gamma+8.d0*gamma**2)
        PI(7,2) = -(1.d0-gamma)/(1.d0-2.d0*gamma)
        PI(8,2) = (-2.d0+6.d0*gamma-5.d0*gamma**2)/(2.d0-6.d0*gamma+4.d0*gamma**2)
        PI(9,2) = 1.d0
        !
        PricesGrids(1,:) = (1.d0-2.d0*gamma)/(2.d0-3.d0*gamma)
        PricesGrids(2,:) = 0.25d0*(3.d0-1.d0/(1.d0-gamma))
        PricesGrids(3,:) = 0.5d0
        !
    END IF
    indexNashPrices = 1
    indexCoopPrices = numActions
    NashProfits = PI(1,:)
    CoopProfits = PI(numActions,:)
    NashMarketShares = 0.d0
    CoopMarketShares = 0.d0
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE computePIMatricesSinghVives
! 
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE computePIMatricesLogit ( DemandParameters, NashPrices, CoopPrices, &
        PI, NashProfits, CoopProfits, &
        indexNashPrices, indexCoopPrices, NashMarketShares, CoopMarketShares, &
        PricesGrids )
    !
    ! Computes the Logit common payoff matrix PI
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: DemandParameters(numDemandParameters)
    REAL(8), DIMENSION(numAgents), INTENT(IN) :: NashPrices, CoopPrices
    REAL(8), INTENT(OUT) :: PI(numActions,numAgents)
    REAL(8), DIMENSION(numAgents), INTENT(OUT) :: NashProfits, CoopProfits, &
        NashMarketShares, CoopMarketShares
    INTEGER, DIMENSION(numAgents), INTENT(OUT) :: indexNashPrices, indexCoopPrices
    REAL(8), DIMENSION(numPrices,numAgents), INTENT(OUT) :: PricesGrids
    !
    ! Declaring local variables
    !
    REAL(8) :: a0, sigma, extend(2), objFct
    REAL(8), DIMENSION(numAgents) :: a, c, d, stepPrices, prices
    INTEGER :: i, j, iter, iAgent
    !
    ! Beginning execution
    !
    ! Computing PI matrices
    !
    ! Extract demand parameters
    !
    a0 = DemandParameters(1)
    a = DemandParameters(2:1+numAgents)
    c = DemandParameters(2+numAgents:1+2*numAgents)
    sigma = DemandParameters(2+2*numAgents)
    extend = DemandParameters(3+2*numAgents:4+2*numAgents)
    !
    ! 1. Compute repeated Nash profits
    !
    NashMarketShares = logitDemands(a0,a,c,sigma,NashPrices)
    NashProfits = (NashPrices-c)*NashMarketShares
    !
    ! 2. Compute cooperation profits
    !
    CoopMarketShares = logitDemands(a0,a,c,sigma,CoopPrices)
    CoopProfits = (CoopPrices-c)*CoopMarketShares
    !
    ! 3. Compute price grid
    !
    PricesGrids(1,:) = NashPrices-extend(1)*(CoopPrices-NashPrices)
    PricesGrids(numPrices,:) = CoopPrices+extend(2)*(CoopPrices-NashPrices)
    stepPrices = (PricesGrids(numPrices,:)-PricesGrids(1,:))/(numPrices-1)
    DO i = 2, numPrices-1
        !
        PricesGrids(i,:) = PricesGrids(i-1,:)+stepPrices
        !
    END DO
    !
    ! 4. Compute Pi matrices
    !
    DO i = 1, numActions
        !
        DO j = 1, numAgents
            !
            prices(j) = PricesGrids(indexActions(i,j),j)
            !
        END DO
        !
        d = logitDemands(a0,a,c,sigma,prices)
        PI(i,:) = (prices-c)*d
        !
    END DO
    !
    ! 5. With logit demand, the repeated Nash prices do no necessarily belong to the
    !    prices grid. Hence, the indexNashPrices vector is empty. Alternatively, we could
    !    look for the row in PricesGrids that is closest to NashPrices (not implemented yet)
    !
    indexNashPrices = 0
    indexCoopPrices = 0
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE computePIMatricesLogit
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    FUNCTION logitDemands ( a0, a, c, sigma, p ) 
    !
    ! Computes logit demands
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: a0, sigma
    REAL(8), DIMENSION(numAgents), INTENT(IN) :: a, c, p
    !
    ! Declaring function's type
    !
    REAL(8), DIMENSION(numAgents) :: logitDemands
    !
    ! Beginning execution
    !
    logitDemands = EXP((a-p)/sigma)
    logitDemands = logitDemands/(SUM(logitDemands)+EXP(a0/sigma))
    !
    ! Ending execution and returning control
    !
    END FUNCTION logitDemands
! 
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
! End of execution
!
END MODULE PI_routines