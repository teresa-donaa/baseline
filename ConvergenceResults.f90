MODULE ConvergenceResults
!
USE globals
USE QL_routines
!
! Computes profit gains and frequency of states of strategies at convergence
!
IMPLICIT NONE
!
CONTAINS
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE ComputeConvResults ( iModel )
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
    INTEGER :: i, j, iGame, iPeriod, iAgent, numPeriods, CycleLength
    INTEGER :: p(DepthState,numAgents), pPrime(numAgents)
    INTEGER :: OptimalStrategyVec(lengthStrategies), LastStateVec(LengthStates)
    INTEGER :: visitedStates(numStates+1), optimalStrategy(numStates,numAgents), &
        LastObservedPrices(DepthState,numAgents)
    REAL(8) :: Profits(numGames,numAgents), visitedProfits(numStates+1,numAgents), AvgProfits(numGames)
    REAL(8), DIMENSION(numAgents) :: meanProfits, seProfit, meanProfitGain, seProfitGain
    REAL(8) :: meanAvgProfit, seAvgProfit, meanAvgProfitGain, seAvgProfitGain
    REAL(8) :: FreqStates(numGames,numStates), meanFreqStates(numStates)
    !
    ! Beginning execution
    !
    PRINT*, 'Computing convergence results (average profits and frequency of prices)'
    !
    ! Initializing variables
    !
    numPeriods = numStates+1        ! If different from numStates, check the dimensions of
                                    ! many of the variables above!!!
    Profits = 0.d0
    FreqStates = 0.d0
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
    22  FORMAT(<LengthStates>(I<lengthFormatActionPrint>,1X))
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
        optimalStrategy = RESHAPE(OptimalStrategyVec, (/ numStates,numAgents /))
        IF (DepthState0 .EQ. 0) THEN
            !
            LastObservedPrices = optimalStrategy
            !
        ELSE IF (DepthState0 .GE. 1) THEN
            !
            LastObservedPrices = RESHAPE(LastStateVec, (/ DepthState,numAgents /))
            !
        END IF
        !
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Convergence analysis
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !
        visitedStates = 0
        visitedProfits = 0.d0
        p = LastObservedPrices
        pPrime = optimalStrategy(computeStateNumber(p),:)
        DO iPeriod = 1, numPeriods
            !
            IF (DepthState .GT. 1) p(2:DepthState,:) = p(1:DepthState-1,:)
            p(1,:) = pPrime
            visitedStates(iPeriod) = computeStateNumber(p)
            DO iAgent = 1, numAgents
                !
                visitedProfits(iPeriod,iAgent) = PI(computeActionNumber(pPrime),iAgent)
                !
            END DO
            !
            ! Check if the state has already been visited
            !
            IF ((iPeriod .GE. 2) .AND. (ANY(visitedStates(:iPeriod-1) .EQ. visitedStates(iPeriod)))) THEN
                !
                CycleLength = iPeriod-MINVAL(MINLOC((visitedStates(:iPeriod-1)-visitedStates(iPeriod))**2))
                Profits(iGame,:) = SUM(visitedProfits(iPeriod-CycleLength+1:iPeriod,:),DIM = 1)/ &
                        DBLE(CycleLength)
                FreqStates(iGame,visitedStates(iPeriod-CycleLength+1:iPeriod)) = 1.d0/DBLE(CycleLength)
                EXIT
                !
            END IF
            pPrime = optimalStrategy(visitedStates(iPeriod),:)
            !
        END DO
        !
    END DO        ! End of loop over games
    !
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Computing averages and descriptive statistics
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    ! Profits
    !
    DO iAgent = 1, numAgents
        !
        meanProfit(iAgent) = SUM(Profits(:,iAgent))/DBLE(numGames)
        seProfit(iAgent) = SQRT((SUM(Profits(:,iAgent)**2)/DBLE(numGames)-meanProfit(iAgent)**2)/DBLE(numGames))
        !
    END DO
    AvgProfits = SUM(Profits,DIM = 2)/DBLE(numAgents)
    meanAvgProfit = SUM(AvgProfits)/DBLE(numGames)
    seAvgProfit = SQRT((SUM(AvgProfits**2)/DBLE(numGames)-meanAvgProfit**2)/DBLE(numGames))
    meanProfitGain = (meanProfit-NashProfits)/(CoopProfits-NashProfits)
    seProfitGain = seProfit/(CoopProfits-NashProfits)
    meanNashProfit = SUM(NashProfits)/numAgents
    meanCoopProfit = SUM(CoopProfits)/numAgents
    meanAvgProfitGain = (meanAvgProfit-meanNashProfit)/(meanCoopProfit-meanNashProfit)
    seAvgProfitGain = seAvgProfit/(meanCoopProfit-meanNashProfit)
    !
    ! States
    !
    DO i = 1, numStates
        !
        meanFreqStates(i) = SUM(freqStates(:,i))/DBLE(numGames)
        !
    END DO
    !
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Printing averages and descriptive statistics
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    IF (iModel .EQ. 1) THEN
        !
        WRITE(100022,1) (i, i = 1, numAgents), (i, i = 1, numExplorationParameters), (i, i = 1, numAgents), &
            (i, i = 1, numDemandParameters), &
            (i, i = 1, numAgents), (i, i = 1, numAgents), &
            (i, i = 1, numAgents), (i, i = 1, numAgents),  &
            (i, i = 1, numAgents), (i, i = 1, numAgents),  &
            ((i, j, j = 1, numPrices), i = 1, numAgents), &
            (i, i, i = 1, numAgents), (i, i, i = 1, numAgents), &
            (labelStates(j), j = 1, numStates)
1       FORMAT('Model ', &
            <numAgents>('    alpha', I1, ' '), &
            <numExplorationParameters>(' MExplPar', I1, ' '), &
            <numAgents>('    delta', I1, ' '), <numDemandParameters>('  DemPar', I2.2, ' '), &
            <numAgents>('NashPrice', I1, ' '), <numAgents>('CoopPrice', I1, ' '), &
            <numAgents>('NashProft', I1, ' '), <numAgents>('CoopProft', I1, ' '), &
            <numAgents>('NashMktSh', I1, ' '), <numAgents>('CoopMktSh', I1, ' '), &
            <numAgents>(<numPrices>('Ag', I1, 'Price', I2.2, ' ')), &
            <numAgents>('  avgProf', I1, 1X, '   seProf', I1, 1X), '   avgProf     seProf ', &
            <numAgents>('avgPrGain', I1, 1X, ' sePrGain', I1, 1X), ' avgPrGain   sePrGain ', &
            <numStates>(A<MAX(10,3+lengthStatesPrint)>, ' ') &
            )
        !
    END IF
    !
    WRITE(100022,2) iModel, &
        alpha, MExpl, delta, DemandParameters, &
        NashPrices, CoopPrices, NashProfits, CoopProfits, NashMarketShares, CoopMarketShares, &
        (PricesGrids(:,i), i = 1, numAgents), &
        (meanProfit(i), seProfit(i), i = 1, numAgents), meanAvgProfit, seAvgProfit, &
        (meanProfitGain(i), seProfitGain(i), i = 1, numAgents), meanAvgProfitGain, seAvgProfitGain, &
        (meanFreqStates(i), i = 1, numStates)
2   FORMAT(I5, 1X, &
        <3*numAgents+numDemandParameters>(F10.3, 1X), &
        <6*numAgents>(F10.7, 1X), &
        <numPrices*numAgents>(F10.5, 1X), &
        <2*(numAgents+1)>(F10.5, 1X), &
        <2*(numAgents+1)>(F10.5, 1X), &
        <numStates>(F<MAX(10,3+lengthStatesPrint)>.6, 1X) &
        )
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE ComputeConvResults
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE ConvergenceResults
