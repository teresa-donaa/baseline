MODULE PreShockCycles
!
USE globals
USE QL_routines
!
! Computes Impulse Response analysis
!
IMPLICIT NONE
!
CONTAINS
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE computePSCycles ( iModel )
    !
    ! Computes statistics for one model
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: iModel
    !
    ! Declaring local variable
    !
    INTEGER, PARAMETER :: numPeriods = 500
    INTEGER :: p(depthState,numAgents), pPrime(numAgents), iPeriod, jAgent, iDepth, &
        iGame, optimalStrategy(numStates,numAgents), lastObservedState(depthState,numAgents), &
        i, visitedStates(numStates), pHist(numPeriods,DepthState,numAgents)
    REAL(8), DIMENSION(numPeriods,numAgents) :: Prices, Profits, avgPrices, avgProfits
    INTEGER :: OptimalStrategyVec(lengthStrategies), LastStateVec(lengthStates)
    !
    ! Beginning execution
    !
    PRINT*, 'Computing Pre Shock Price and Profit Cycles'
    !
    ! Initializing output file
    !
    avgPrices = 0.d0
    avgProfits = 0.d0
    IF (iModel .EQ. 1) THEN
        !
        WRITE(100021,1) (iDepth, iDepth = 1, DepthState), (iPeriod, iPeriod = 1, numPeriods), &
            (iPeriod, iPeriod = 1, numPeriods), &
            ((iPeriod, iDepth, iDepth = 1, DepthState), iPeriod = 1, numPeriods)
1       FORMAT(' Model   Game  Agent ', <DepthState>('LastStateDepth', I1, ' '), &
            <numPeriods>('  PricePeriod', I3.3, ' '), &
            <numPeriods>('ProfitsPeriod', I3.3, ' '), &
            <numPeriods>('StatePeriod', I3.3, 'Depth', I1, ' '))
        !
    END IF
    !
    ! Reading strategies and states at convergence from file
    !
    OPEN(UNIT = 998,FILE = FileNameIndexStrategies,STATUS = "OLD")
    DO i = 1, lengthStrategies
        !
        IF (MOD(i,10000) .EQ. 0) PRINT*, 'Read ', i, ' lines of indexStrategies'
        READ(998,21) (indexStrategies(i,iGame), iGame = 1, numGames)
    21  FORMAT(<numGames>(I<lengthFormatActionPrint>,1X))
        !
    END DO
    CLOSE(UNIT = 998)                   ! Close indexStrategies file
    OPEN(UNIT = 999,FILE = FileNameIndexLastState,STATUS = "OLD")     ! Open indexLastState file
    DO iGame = 1, numGames
        !
        READ(999,22) indexLastState(:,iGame)
    22  FORMAT(<lengthStates>(I<lengthFormatActionPrint>,1X))
        !
    END DO
    PRINT*, 'Read indexLastState'
    CLOSE(UNIT = 999)                   ! Close indexLastState file
    !
    ! Beginning loop over games
    !
    DO iGame = 1, numGames        ! Start of loop aver games
        !
        PRINT*, 'iGame = ', iGame
        !
        OptimalStrategyVec = indexStrategies(:,iGame)
        LastStateVec = indexLastState(:,iGame)
        !
        optimalStrategy = RESHAPE(OptimalStrategyVec, (/ numStates,numAgents /) )
        lastObservedState = RESHAPE(LastStateVec, (/ depthState,numAgents /) )
        !
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Pre-shock period analysis
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !
        p = lastObservedState
        pPrime = p(1,:)
        Prices = 0.d0
        Profits = 0.d0
        pHist = 0
        !
        DO iPeriod = 1, numPeriods
            !
            pHist(iPeriod,:,:) = p
            visitedStates(iPeriod) = computeStateNumber(p)
            DO jAgent = 1, numAgents
                !
                Prices(iPeriod,jAgent) = PricesGrids(pPrime(jAgent),jAgent)
                Profits(iPeriod,jAgent) = PI(computeActionNumber(pPrime),jAgent)
                avgPrices(iPeriod,jAgent) = avgPrices(iPeriod,jAgent)+Prices(iPeriod,jAgent)/DBLE(numGames)
                avgProfits(iPeriod,jAgent) = avgProfits(iPeriod,jAgent)+Profits(iPeriod,jAgent)/DBLE(numGames)
                !
            END DO
            pPrime = optimalStrategy(visitedStates(iPeriod),:)
            IF (depthState .GT. 1) p(2:depthState,:) = p(1:depthState-1,:)
            p(1,:) = pPrime
            !
        END DO
        !
        DO jAgent = 1, numAgents
            !
            WRITE(100021,2) iModel, iGame, jAgent, &
                (lastObservedState(iDepth,jAgent), iDepth = 1, DepthState), &
                (Prices(iPeriod,jAgent), iPeriod = 1, numPeriods), &
                (Profits(iPeriod,jAgent), iPeriod = 1, numPeriods), &
                ((pHist(iPeriod,iDepth,jAgent), iDepth = 1, DepthState), iPeriod = 1, numPeriods)
2           FORMAT(<3>(I6, 1X), <DepthState>(I15, 1X), &
                <2*numPeriods>(F16.9, 1X), <numPeriods>(<DepthState>(I20, 1X)))
            !
        END DO
        !
    END DO        ! End of loop over games
    !
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Computing and printing averages 
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    ! Averages of prices and profits
    !
    DO jAgent = 1, numAgents
        !
        WRITE(100021,3) iModel, "  Avg", jAgent, &
            (avgPrices(iPeriod,jAgent), iPeriod = 1, numPeriods), &
            (avgProfits(iPeriod,jAgent), iPeriod = 1, numPeriods)
3       FORMAT(I6, 1X, A6, 1X, I6, 1X, &
            <DepthState>(16X), &
            <2*numAgents*numPeriods>(F16.9, 1X), <numPeriods*DepthState>(21X))
        !
    END DO
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE computePSCycles
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE PreShockCycles
