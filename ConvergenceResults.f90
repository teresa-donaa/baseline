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
    INTEGER :: i, j, iGame, iPeriod, iAgent, CycleLength
    INTEGER :: p(DepthState,numAgents), pPrime(numAgents)
    INTEGER :: OptimalStrategyVec(lengthStrategies), LastStateVec(LengthStates)
    INTEGER :: VisitedStates(numPeriods), OptimalStrategy(numStates,numAgents), &
        LastObservedPrices(DepthState,numAgents)
    INTEGER :: pHist(numPeriods,numAgents)
    REAL(8) :: Profits(numGames,numAgents), VisitedProfits(numPeriods,numAgents), AvgProfits(numGames)
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
    Profits = 0.d0
    FreqStates = 0.d0
    !
    ! Reading strategies and states at convergence from file
    !
    OPEN(UNIT = 998,FILE = FileNameIndexStrategies,STATUS = "OLD")
    READ(998,*)     ! Skip 'converged' line
    READ(998,*)     ! Skip 'timeToConvergence' line
20  FORMAT(//)    
    DO i = 1, lengthStrategies
        !
        IF (MOD(i,10000) .EQ. 0) PRINT*, 'Read ', i, ' lines of indexStrategies'
        READ(998,21) (indexStrategies(i,iGame), iGame = 1, numGames)
21      FORMAT(<numGames>(I<lengthFormatActionPrint>,1X))
        !
    END DO
    CLOSE(UNIT = 998)                   ! Close indexStrategies file
    !
    OPEN(UNIT = 999,FILE = FileNameIndexLastState,STATUS = "OLD")     ! Open indexLastState file
    DO iGame = 1, numGames
        !
        READ(999,22) indexLastState(:,iGame)
22      FORMAT(<LengthStates>(I<lengthFormatActionPrint>,1X))
        !
    END DO
    PRINT*, 'Read indexLastState'
    CLOSE(UNIT = 999)                   ! Close indexLastState file
    !
    OPEN(UNIT = 999,FILE = FileNamePriceCycles,STATUS = "REPLACE")        ! Open priceCycles file
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
        OptimalStrategy = RESHAPE(OptimalStrategyVec, (/ numStates,numAgents /))
        IF (DepthState0 .EQ. 0) THEN
            !
            LastObservedPrices = OptimalStrategy
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
        VisitedStates = 0
        VisitedProfits = 0.d0
        pHist = 0
        p = LastObservedPrices
        pPrime = OptimalStrategy(computeStateNumber(p),:)
        DO iPeriod = 1, numPeriods
            !
            IF (DepthState .GT. 1) p(2:DepthState,:) = p(1:DepthState-1,:)
            p(1,:) = pPrime
            pHist(iPeriod,:) = pPrime
            VisitedStates(iPeriod) = computeStateNumber(p)
            DO iAgent = 1, numAgents
                !
                VisitedProfits(iPeriod,iAgent) = PI(computeActionNumber(pPrime),iAgent)
                !
            END DO
            !
            ! Check if the state has already been visited
            !
            IF ((iPeriod .GE. 2) .AND. (ANY(VisitedStates(:iPeriod-1) .EQ. VisitedStates(iPeriod)))) EXIT
            !
            ! Update pPrime and iterate
            !
            pPrime = OptimalStrategy(VisitedStates(iPeriod),:)
            !
        END DO
        !
        CycleLength = iPeriod-MINVAL(MINLOC((VisitedStates(:iPeriod-1)-VisitedStates(iPeriod))**2))
        Profits(iGame,:) = SUM(VisitedProfits(iPeriod-CycleLength+1:iPeriod,:),DIM = 1)/ &
                DBLE(CycleLength)
        FreqStates(iGame,VisitedStates(iPeriod-CycleLength+1:iPeriod)) = 1.d0/DBLE(CycleLength)
        !
        ! Computing and writing price cycles
        !
        pHist(:CycleLength,:) = pHist(iPeriod-CycleLength+1:iPeriod,:)
        pHist(CycleLength+1:,:) = 0.d0
        VisitedStates(:CycleLength) = VisitedStates(iPeriod-CycleLength+1:iPeriod)
        VisitedStates(CycleLength+1:) = 0
        VisitedProfits(:CycleLength,:) = VisitedProfits(iPeriod-CycleLength+1:iPeriod,:)
        VisitedProfits(CycleLength+1:,:) = 0.d0
        WRITE(999,211) CycleLength, &
            VisitedStates(:CycleLength), &
            (pHist(:CycleLength,iAgent), iAgent = 1, numAgents), &
            (VisitedProfits(:CycleLength,iAgent), iAgent = 1, numAgents)
211     FORMAT(I8, 1X, &
            <CycleLength>(I<lengthStatesPrint>,1X), &
            <numAgents>(<CycleLength>(I<lengthFormatActionPrint>,1X)), &
            <numAgents>(<CycleLength>(F8.5, 1X)))
        !
    END DO        ! End of loop over games
    !
    CLOSE(UNIT = 999)                   ! Close priceCycles file
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
        seProfit(iAgent) = SQRT(ABS((SUM(Profits(:,iAgent)**2)/DBLE(numGames)-meanProfit(iAgent)**2)/DBLE(numGames)))
        !
    END DO
    AvgProfits = SUM(Profits,DIM = 2)/DBLE(numAgents)
    meanAvgProfit = SUM(AvgProfits)/DBLE(numGames)
    seAvgProfit = SQRT(ABS((SUM(AvgProfits**2)/DBLE(numGames)-meanAvgProfit**2)/DBLE(numGames)))
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
        WRITE(100022,1) &
            (i, i = 1, numAgents), &
            (i, i = 1, numExplorationParameters), (i, i = 1, numAgents), &
            (i, (j, i, j = 1, 2), i = 1, numAgents), &
            (i, i = 1, numDemandParameters), &
            (i, i = 1, numAgents), (i, i = 1, numAgents), &
            (i, i = 1, numAgents), (i, i = 1, numAgents),  &
            (i, i = 1, numAgents), (i, i = 1, numAgents),  &
            ((i, j, j = 1, numPrices), i = 1, numAgents), &
            (i, i, i = 1, numAgents), (i, i, i = 1, numAgents), &
            (labelStates(j), j = 1, numStates)
1       FORMAT('Model ', &
            <numAgents>('    alpha', I1, ' '), &
            <numExplorationParameters>('     beta', I1, ' '), <numAgents>('    delta', I1, ' '), &
            <numAgents>('typeQini', I1, ' ', <2>('par', I1, 'Qini', I1, ' ')), &
            <numDemandParameters>('  DemPar', I0.2, ' '), &
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
    WRITE(100022,2) codModel, &
        alpha, MExpl, delta, &
        (typeQInitialization(i), parQInitialization(i, :), i = 1, numAgents), &
        DemandParameters, &
        NashPrices, CoopPrices, NashProfits, CoopProfits, NashMarketShares, CoopMarketShares, &
        (PricesGrids(:,i), i = 1, numAgents), &
        (meanProfit(i), seProfit(i), i = 1, numAgents), meanAvgProfit, seAvgProfit, &
        (meanProfitGain(i), seProfitGain(i), i = 1, numAgents), meanAvgProfitGain, seAvgProfitGain, &
        (meanFreqStates(i), i = 1, numStates)
2   FORMAT(I5, 1X, &
        <numAgents>(F10.5, 1X), <numExplorationParameters>(F10.5, 1X), <numAgents>(F10.5, 1X), &
        <numAgents>(A9, 1X, <2>(F9.2, 1X)), &
        <numDemandParameters>(F10.5, 1X), &
        <6*numAgents>(F10.5, 1X), &
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
