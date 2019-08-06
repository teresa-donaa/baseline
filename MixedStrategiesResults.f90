MODULE MixedStrategiesResults
!
USE globals
USE QL_routines
USE generic_routines
!
! Computes mixed strategies results:
! profit gains and trajectories to equilibrium cycle
!
IMPLICIT NONE
!
CONTAINS
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE ComputeMixStratResults ( )
    !
    ! Computes statistics for one model
    !
    IMPLICIT NONE
    !
    ! Declaring local variable
    !
    INTEGER, PARAMETER :: numPeriodsPrint = 25
    INTEGER :: iModel, i, j, iGame, iPeriod, iAgent, jAgent, iCycle, numPriceCycles, iPriceCycle
    INTEGER :: p(DepthState,numAgents), pPrime(numAgents)
    INTEGER :: OptimalStrategyVec(lengthStrategies), CicleLengthVec(numAgents), CycleLength
    INTEGER :: visitedStates(numStates+1), optimalStrategy(numStates,numAgents), &
        LastObservedPrices(DepthState,numAgents)
    INTEGER :: CycleLengthGames(numGames), priceCyclesMat(numStates+1,numAgents), &
        PriceCyclesComb(1000,numAgents)
    INTEGER :: idumRS, ivRS(32), iyRS, idum2RS, u
    REAL(8) :: Profits(numGames,numAgents), visitedProfits(numStates+1,numAgents), AvgProfits(numGames)
    REAL(8), DIMENSION(numAgents) :: meanProfits, seProfit, meanProfitGain, seProfitGain
    REAL(8) :: meanAvgProfit, seAvgProfit, meanAvgProfitGain, seAvgProfitGain
    REAL(8) :: FreqStates(numGames,numStates), meanFreqStates(numStates)
    REAL(8) :: uRandomSampling(numGames,numAgents), price_tmp, profit_tmp
    REAL(8), DIMENSION(numGames,0:numPeriodsPrint,numAgents) :: priceTrajectories, profitTrajectories
    REAL(8), DIMENSION(0:numPeriodsPrint,numAgents) :: &
        AvgPriceTrajectories, seAvgPriceTrajectories, AvgProfitTrajectories, seAvgProfitTrajectories
    !
    ! Beginning execution
    !
    PRINT*, 'Computing mixed strategy results'
    !
    ! Drawing random numbers
    !
    idumRS = -1
    idum2RS = 123456789
    ivRS = 0
    iyRS = 0
    DO iGame = 1, numGAmes
        !
        DO iAgent = 1, numAgents
            !
            uRandomSampling(iGame,iAgent) = ran2(idumRS,ivRS,iyRS,idum2RS)
            !
        END DO
        !
    END DO
    !
    ! Loading strategies and price cycles
    !
    sampledPriceCycles = 0
    DO iAgent = 1, numAgents
        !
        ! Initializing file names
        !
        i = 1+INT(LOG10(DBLE(numModels)))
        WRITE(ModelNumber, "(I0.<i>, A4)") computeMixedStrategies(iAgent), ".txt"
        FileNameIndexStrategies = "indexStrategiesTransposed_" // ModelNumber
        FileNamePriceCycles = "priceCycles_" // ModelNumber
        !
        ! Reading strategies and price cycles at convergence from file
        !
        OPEN(UNIT = 998,FILE = FileNameIndexStrategies,STATUS = "OLD")
        DO i = 1, lengthStrategies
            !
            IF (MOD(i,10000) .EQ. 0) PRINT*, 'Agent: ', iAgent, ': Read ', i, ' lines of indexStrategies'
            READ(998,21) (indexStrategies(i,iGame), iGame = 1, numGames)
        21  FORMAT(<numGames>(I<lengthFormatActionPrint>,1X))
            !
        END DO
        CLOSE(UNIT = 998)                   ! Close indexStrategies file
        !
        priceCycles = 0
        OPEN(UNIT = 999,FILE = FileNamePriceCycles,STATUS = "OLD")     ! Open priceCycles file
        DO iGame = 1, numGames
            !
            READ(999,211) CycleLengthGames(iGame), &
                ((priceCycles((jAgent-1)*CycleLengthGames(iGame)+iCycle,iGame), iCycle = 1, CycleLengthGames(iGame)), jAgent = 1, numAgents)
211         FORMAT(I8, 1X, <numAgents>(<CycleLengthGames(iGame)>(I<lengthFormatActionPrint>,1X)))        
            !
        END DO
        PRINT*, 'Read priceCycles'
        CLOSE(UNIT = 999)                   ! Close priceCycles file
        !
        ! Random sampling of strategies and price cycles for agent 'iAgent'
        !
        DO iGame = 1, numGames
            !
            u = 1+INT(DBLE(numGames)*uRandomSampling(iGame,iAgent))
            sampledIndexStrategies((iAgent-1)*numStates+1:iAgent*numStates,iGame) = &
                indexStrategies((iAgent-1)*numStates+1:iAgent*numStates,u)
            sampledPriceCycles(1,(iGame-1)*numAgents+iAgent) = CycleLengthGames(u)
            sampledPriceCycles(2:CycleLengthGames(u)+1,(iGame-1)*numAgents+iAgent) = &
                priceCycles((iAgent-1)*CycleLengthGames(u)+1:iAgent*CycleLengthGames(u),u)
            !
        END DO
        !
    END DO
    !
    ! Initializing other variables
    !
    Profits = 0.d0
    priceTrajectories = 0.d0
    profitTrajectories = 0.d0
    !
    ! Beginning loop over games
    !
    DO iGame = 1, numGames        ! Start of loop aver games
        !
        PRINT*, 'iGame = ', iGame
        !
        OptimalStrategyVec = sampledIndexStrategies(:,iGame)
        CicleLengthVec = sampledPriceCycles(1,(iGame-1)*numAgents+1:iGame*numAgents)
        PriceCyclesMat = sampledPriceCycles(2:,(iGame-1)*numAgents+1:iGame*numAgents)
        !
        ! Generate the combinations of price cycles
        !
        numPriceCycles = PRODUCT(CicleLengthVec)
        PriceCyclesComb = 0
        CALL generateCombinations(numStates+1,numAgents,CicleLengthVec,PriceCyclesMat,numPriceCycles, &
            PriceCyclesComb(:numPriceCycles,:))
        !
        optimalStrategy = RESHAPE(OptimalStrategyVec, (/ numStates,numAgents /))
        !
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Convergence profits
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !
        DO iPriceCycle = 1, numPriceCycles
            !
            visitedStates = 0
            visitedProfits = 0.d0
            IF (DepthState0 .EQ. 0) THEN
                !
                p = optimalStrategy
                !
            ELSE IF (DepthState0 .GE. 1) THEN
                !
                p = RESHAPE(PriceCyclesComb(iPriceCycle,:), (/ DepthState,numAgents /))
                !
            END IF
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
                IF ((iPeriod .GE. 2) .AND. (ANY(visitedStates(:iPeriod-1) .EQ. visitedStates(iPeriod)))) EXIT
                !
                ! Update pPrime and iterate
                !
                pPrime = optimalStrategy(visitedStates(iPeriod),:)
                !
            END DO
            !
            CycleLength = iPeriod-MINVAL(MINLOC((visitedStates(:iPeriod-1)-visitedStates(iPeriod))**2))
            CALL updateVectorAverage(numAgents,iPriceCycle, &
                SUM(visitedProfits(iPeriod-CycleLength+1:iPeriod,:),DIM = 1)/DBLE(CycleLength),Profits(iGame,:))
            !
        END DO
        !
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Price and profit adjustment trajectories
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !
        DO iPriceCycle = 1, numPriceCycles
            !
            IF (DepthState0 .EQ. 0) THEN
                !
                p = optimalStrategy
                !
            ELSE IF (DepthState0 .GE. 1) THEN
                !
                p = RESHAPE(PriceCyclesComb(iPriceCycle,:), (/ DepthState,numAgents /))
                !
            END IF
            !
            DO iAgent = 1, numAgents
                !               
                price_tmp = PricesGrids(p(1,iAgent),iAgent)
                CALL updateScalarAverage(iPriceCycle,price_tmp,priceTrajectories(iGame,0,iAgent))
                profit_tmp = PI(computeActionNumber(p(1,:)),iAgent)
                CALL updateScalarAverage(iPriceCycle,profit_tmp,profitTrajectories(iGame,0,iAgent))
                !
            END DO
            !
            pPrime = optimalStrategy(computeStateNumber(p),:)
            DO iPeriod = 1, numPeriodsPrint
                !
                IF (DepthState .GT. 1) p(2:DepthState,:) = p(1:DepthState-1,:)
                p(1,:) = pPrime
                DO iAgent = 1, numAgents
                    !
                    price_tmp = PricesGrids(p(1,iAgent),iAgent)
                    CALL updateScalarAverage(iPriceCycle,price_tmp,priceTrajectories(iGame,iPeriod,iAgent))
                    profit_tmp = PI(computeActionNumber(p(1,:)),iAgent)
                    CALL updateScalarAverage(iPriceCycle,profit_tmp,profitTrajectories(iGame,iPeriod,iAgent))
                    !
                END DO
                !
                ! Update pPrime and iterate
                !
                pPrime = optimalStrategy(computeStateNumber(p),:)
                !
            END DO
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
    ! Trajectories
    !
    AvgPriceTrajectories = SUM(priceTrajectories,DIM = 1)/DBLE(numGames)
    seAvgPriceTrajectories = SQRT(ABS((SUM(priceTrajectories**2,DIM = 1)/DBLE(numGames)-AvgPriceTrajectories**2)/DBLE(numGames)))
    AvgProfitTrajectories = SUM(profitTrajectories,DIM = 1)/DBLE(numGames)
    AvgProfitTrajectories = (AvgProfitTrajectories-SPREAD(NashProfits,DIM = 1,Ncopies = numPeriodsPrint+1))/ &
        SPREAD(CoopProfits-NashProfits,DIM = 1,Ncopies = numPeriodsPrint+1)
    seAvgProfitTrajectories = SQRT(ABS((SUM(profitTrajectories**2,DIM = 1)/DBLE(numGames)-AvgProfitTrajectories**2)/DBLE(numGames)))
    seAvgProfitTrajectories = seAvgProfitTrajectories/ &
        SPREAD(CoopProfits-NashProfits,DIM = 1,Ncopies = numPeriodsPrint+1)
    !
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Printing averages and descriptive statistics
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    WRITE(111,1) (i, i = 1, numAgents), &
        (i, i = 1, numAgents), (i, i = 1, numExplorationParameters), (i, i = 1, numAgents), &
        (i, i = 1, numDemandParameters), &
        (i, i = 1, numAgents), (i, i = 1, numAgents), &
        (i, i = 1, numAgents), (i, i = 1, numAgents),  &
        (i, i = 1, numAgents), (i, i = 1, numAgents),  &
        ((i, j, j = 1, numPrices), i = 1, numAgents), &
        (i, i, i = 1, numAgents), (i, i, i = 1, numAgents), &
        ((i, j, j = 0, numPeriodsPrint), i = 1, numAgents), &
        ((i, j, j = 0, numPeriodsPrint), i = 1, numAgents), &
        ((i, j, j = 0, numPeriodsPrint), i = 1, numAgents), &
        ((i, j, j = 0, numPeriodsPrint), i = 1, numAgents)
1   FORMAT(<numAgents>('  Model', I1, 1X), &
        <numAgents>('    alpha', I1, ' '), &
        <numExplorationParameters>('     beta', I1, ' '), &
        <numAgents>('    delta', I1, ' '), <numDemandParameters>('  DemPar', I2.2, ' '), &
        <numAgents>('NashPrice', I1, ' '), <numAgents>('CoopPrice', I1, ' '), &
        <numAgents>('NashProft', I1, ' '), <numAgents>('CoopProft', I1, ' '), &
        <numAgents>('NashMktSh', I1, ' '), <numAgents>('CoopMktSh', I1, ' '), &
        <numAgents>(<numPrices>('Ag', I1, 'Price', I2.2, ' ')), &
        <numAgents>('  avgProf', I1, 1X, '   seProf', I1, 1X), '   avgProf     seProf ', &
        <numAgents>('avgPrGain', I1, 1X, ' sePrGain', I1, 1X), ' avgPrGain   sePrGain ', &
        <numAgents>(<numPeriodsPrint+1>(' avgPriceAg', I1, 'Per', I2.2, ' ')), &
        <numAgents>(<numPeriodsPrint+1>('  sePriceAg', I1, 'Per', I2.2, ' ')), &
        <numAgents>(<numPeriodsPrint+1>('avgProfitAg', I1, 'Per', I2.2, ' ')), &
        <numAgents>(<numPeriodsPrint+1>(' seProfitAg', I1, 'Per', I2.2, ' ')) &
        )
    !
    WRITE(111,2) computeMixedStrategies, &
        alpha, MExpl, delta, DemandParameters, &
        NashPrices, CoopPrices, NashProfits, CoopProfits, NashMarketShares, CoopMarketShares, &
        (PricesGrids(:,i), i = 1, numAgents), &
        (meanProfit(i), seProfit(i), i = 1, numAgents), meanAvgProfit, seAvgProfit, &
        (meanProfitGain(i), seProfitGain(i), i = 1, numAgents), meanAvgProfitGain, seAvgProfitGain, &
        ((AvgPriceTrajectories(j,i), j = 0, numPeriodsPrint), i = 1, numAgents), &
        ((seAvgPriceTrajectories(j,i), j = 0, numPeriodsPrint), i = 1, numAgents), &
        ((AvgProfitTrajectories(j,i), j = 0, numPeriodsPrint), i = 1, numAgents), &
        ((seAvgProfitTrajectories(j,i), j = 0, numPeriodsPrint), i = 1, numAgents)
2   FORMAT(<numAgents>(I8,1X), &
        <3*numAgents+numDemandParameters>(F10.5, 1X), &
        <6*numAgents>(F10.5, 1X), &
        <numPrices*numAgents>(F10.5, 1X), &
        <2*(numAgents+1)>(F10.5, 1X), &
        <2*(numAgents+1)>(F10.5, 1X), &
        <numAgents>(<numPeriodsPrint+1>(F17.5, 1X)), &
        <numAgents>(<numPeriodsPrint+1>(F17.5, 1X)), &
        <numAgents>(<numPeriodsPrint+1>(F17.5, 1X)), &
        <numAgents>(<numPeriodsPrint+1>(F17.5, 1X)) &
        )
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE ComputeMixStratResults
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE MixedStrategiesResults
