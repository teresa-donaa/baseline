MODULE DetailedAnalysis
!
USE globals
USE QL_routines
USE ImpulseResponse
USE QGapToMaximum
!
! Computes disaggregated analysis
!
IMPLICIT NONE
!
CONTAINS
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE ComputeDetailedAnalysis ( iModel )
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
    INTEGER, DIMENSION(numGames) :: CycleLength, converged
    INTEGER, DIMENSION(numShockPeriodsPrint) :: ShockStates
    INTEGER, DIMENSION(numPeriods,numGames) :: CycleStates
    INTEGER, DIMENSION(numAgents,numPeriods,numGames) :: CyclePrices
    INTEGER :: PeriodsLengthPre, ShockLength, PostLength, &
        VisitedStatesPre(numPeriods), VisitedStates(MAX(numShockPeriodsPrint,numPeriods)), &
        VisitedStatesTMP(numPeriods), SameCyclePrePost, &
        p(DepthState,numAgents), pPrime(numAgents), &
        iState, iStatePre, iGame, iAgent, jAgent, iPrice, iPeriod, jPeriod, iPeriodState, &
        indexShockState(LengthStates), i, j, PreCycleLength, QCellCycleLength, IndexDynamicBR(numAgents)
	INTEGER :: OptimalStrategyVec(lengthStrategies), OptimalStrategy(numStates,numAgents)
    INTEGER, DIMENSION(numPeriods,numAgents) :: IndPrePrices
    INTEGER, DIMENSION(numShockPeriodsPrint,numAgents) :: ShockPrices, &
        StaticBRPrices, DynamicBRPrices
    INTEGER, DIMENSION(numAgents) :: flagBRAll, flagBROnPath, flagBROffPath
    INTEGER :: flagEQAll, flagEQOnPath, flagEQOffPath
    !
    REAL(8) :: PIStaticBR, timeToConvergence(numGames)
    REAL(8), DIMENSION(numAgents,numPeriods,numGames) :: CycleProfits
    REAL(8), DIMENSION(numPeriods,numAgents) :: visitedPrices, VisitedProfits, PrePrices, PreProfits
    REAL(8), DIMENSION(numShockPeriodsPrint,numAgents) :: ShockRealPrices, ShockProfits, OptStratQ, DynamicBRQ
    REAL(8), DIMENSION(numAgents) :: DeviationQ, ProfitGains
    REAL(8), DIMENSION(numAgents) :: AvgPrePrices, AvgPreProfits, AvgPostPrices, AvgPostProfits
    REAL(8), DIMENSION(numAgents) :: freqBRAll, freqBROnPath, freqBROffPath
    REAL(8) :: freqEQAll, freqEQOnPath, freqEQOffPath
    REAL(8), DIMENSION(0:numAgents) :: QGapTotGame,QGapOnPathGame,QGapNotOnPathGame,QGapNotBRAllStatesGame, &
            QGapNotBRonPathGame,QGapNotEqAllStatesGame,QGapNotEqonPathGame
    !
    LOGICAL :: FlagReturnedToState
    !
    CHARACTER(len = 30) :: filename
    CHARACTER(len = 3) :: iPriceChar
    !
    ! Beginning execution
    !
    PRINT*, 'Computing Detailed Analysis'
    !
    ! Reading strategies and states at convergence from file
    !
    CALL ReadInfoModel(converged,timeToConvergence, & 
        CycleLength,CycleStates,CyclePrices,CycleProfits,indexStrategies)
    !
    ! Writing header line in global output file
    !
    WRITE(100033,11) &
        (i, i = 1, numAgents), (i, i = 1, numAgents), &
        (i, i = 1, numAgents), (i, i = 1, numAgents), (i, i = 1, numAgents), &
        (i, i = 1, numAgents), (i, i = 1, numAgents), (i, i = 1, numAgents), &
        (i, i = 1, numAgents), (i, i = 1, numAgents), (i, i = 1, numAgents), &
        (i, i = 1, numAgents), (i, i = 1, numAgents), (i, i = 1, numAgents), &
        (i, i = 1, numAgents), (i, i = 1, numAgents), &
        (i, i = 1, numAgents), (i, i = 1, numAgents), &
        (i, i = 1, numAgents), (i, i = 1, numAgents), &
        (i, i = 1, numShockPeriodsPrint), &
        (i, i = 1, numShockPeriodsPrint), &
        (i, i = 1, numShockPeriodsPrint), &
        (i, i = 1, numShockPeriodsPrint), &
        (i, i = 1, numShockPeriodsPrint), &
        (i, i = 1, numShockPeriodsPrint), &
        (i, i = 1, numAgents), (i, i = 1, numAgents)
11  FORMAT('    Game DevTo_Price ', &
        <numAgents>(' NashProfit', I1, ' '), <numAgents>(' CoopProfit', I1, ' '), &
        'PreShock_CycleLength ', &
        <numAgents>('Avg_Pre_Price', I1, ' '), <numAgents>('Avg_Pre_Profit', I1, ' '), <numAgents>('ProfitGain', I1, ' '), &
        'Converged TimeToConvergence PreShock_NumInCycle ', &
        'flagEQAll flagEQOnPath flagEQOffPath ', &
        'freqEQAll freqEQOnPath freqEQOffPath ', &
        <numAgents>('flagBRAll', I1, ' '), <numAgents>('flagBROnPath', I1, ' '), <numAgents>('flagBROffPath', I1, ' '), &
        <numAgents>('freqBRAll', I1, ' '), <numAgents>('freqBROnPath', I1, ' '), <numAgents>('freqBROffPath', I1, ' '), &
        '   QGapTot QGapOnPath QGapNotOnPath QGapNotBRAllStates ', &
            'QGapNotBRonPath QGapNotEqAllStates QGapNotEqonPath ', &
        <numAgents>('   QGapTot', I1, ' '), <numAgents>('QGapOnPath', I1, ' '), <numAgents>('QGapNotOnPath', I1, ' '), &
            <numAgents>('QGapNotBRAllStates', I1, ' '), <numAgents>('QGapNotBRonPath', I1, ' '), &
            <numAgents>('QGapNotEqAllStates', I1, ' '), <numAgents>('QGapNotEqonPath', I1, ' '), &
        <numAgents>('PreShock_Price', I1, ' '), <numAgents>('PreShock_Profit', I1, ' '), &
        'Shock_Agent Obs_Agent  Deviation_Q ShockLength SameCyclePrePost ', &
        <numShockPeriodsPrint>('Shock_Price', I0.3, ' '), &
        <numShockPeriodsPrint>('Shock_Profit', I0.3, ' '), &
        <numShockPeriodsPrint>('StaticBR_Price', I0.3, ' '), &
        <numShockPeriodsPrint>('DynamicBR_Price', I0.3, ' '), &
        <numShockPeriodsPrint>('OptStrat_Q', I0.3, ' '), &
        <numShockPeriodsPrint>('DynamicBR_Q', I0.3, ' '), &
        'PostShock_CycleLength ', <numAgents>('Avg_Post_Price', I1, ' '), <numAgents>('Avg_Post_Profit', I1, ' '))
    !
    ! Beginning loop over games
    !
    DO iGame = 1, numGames        ! Start of loop over games
        !
        PRINT*, 'iGame = ', iGame
        !
        OptimalStrategyVec = indexStrategies(:,iGame)
        OptimalStrategy = RESHAPE(OptimalStrategyVec, (/ numStates,numAgents /) )
        !
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Pre-shock period analysis
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !
        PeriodsLengthPre = CycleLength(iGame)
        VisitedStatesPre = 0
        PrePrices = 0.d0
        PreProfits = 0.d0
        IndPrePrices = 0
        AvgPrePrices = 0.d0
        AvgPreProfits = 0.d0
        !
        DO iPeriod = 1, PeriodsLengthPre
            !
            VisitedStatesPre(iPeriod) = CycleStates(iPeriod,iGame)
            DO iAgent = 1, numAgents
                !
                IndPrePrices(iPeriod,iAgent) = CyclePrices(iAgent,iPeriod,iGame)
                PrePrices(iPeriod,iAgent) = PricesGrids(IndPrePrices(iPeriod,iAgent),iAgent)
                PreProfits(iPeriod,iAgent) = CycleProfits(iAgent,iPeriod,iGame)
                !
            END DO
            !
        END DO
        !
        AvgPrePrices = SUM(PrePrices(:PeriodsLengthPre,:),DIM = 1)/DBLE(PeriodsLengthPre)
        AvgPreProfits = SUM(PreProfits(:PeriodsLengthPre,:),DIM = 1)/DBLE(PeriodsLengthPre)
        !
        ! Compute indicators that depend on the strategy only:
        ! ProfitGain, Statistics on BR and EQ
        !
        ProfitGains = (AvgPreProfits-NashProfits)/(CoopProfits-NashProfits)
        CALL computeEqCheckGame(OptimalStrategy,PeriodsLengthPre,VisitedStatesPre(:PeriodsLengthPre),SlackOnPath,SlackOffPath, &
            freqBRAll,freqBROnPath,freqBROffPath,freqEQAll,freqEQOnPath,freqEQOffPath, &
            flagBRAll,flagBROnPath,flagBROffPath,flagEQAll,flagEQOnPath,flagEQOffPath)
        !
        ! Compute Q gap for the optimal strategy for all agents, in all states and actions
        !
        CALL computeQGapToMaxGame(OptimalStrategy,PeriodsLengthPre,CycleStates(:PeriodsLengthPre,iGame), &
            QGapTotGame,QGapOnPathGame,QGapNotOnPathGame,QGapNotBRAllStatesGame, &
           QGapNotBRonPathGame,QGapNotEqAllStatesGame,QGapNotEqonPathGame)
        !
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! IR analysis with deviation to iPrice
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !
        ! Beginning loop over deviation prices
        !
        DO iPrice = 1, numPrices    ! Start of loop over possible deviations
            !
            DO iStatePre = 1, PeriodsLengthPre      ! Start of loop over pre-shock cycle states
                !
                DO iAgent = 1, numAgents            ! Start of loop over deviating agent 
                    !
                    ShockPrices = 0
                    ShockRealPrices = 0.d0
                    ShockProfits = 0.d0
                    StaticBRPrices = 0
                    DynamicBRPrices = 0
                    StaticBRPrices = 0.d0
                    DynamicBRPrices = 0.d0
                    OptStratQ = 0.d0
                    DynamicBRQ = 0.d0
                    DeviationQ = 0.d0
                    AvgPostPrices = 0.d0
                    AvgPostProfits = 0.d0
                    VisitedStates = 0
                    !
                    ! Find prices and Qs in deviation period n. 1 according to actual price selection:
                    ! iAgent selects each price in turn, other agents stick to the strategy at convergence
                    !
                    pPrime = OptimalStrategy(VisitedStatesPre(iStatePre),:)
                    pPrime(iAgent) = iPrice
                    DO jAgent = 1, numAgents
                        !
                        CALL computeQCell(OptimalStrategy,VisitedStatesPre(iStatePre),pPrime(jAgent),jAgent,delta, &
                            DeviationQ(jAgent),VisitedStatesTMP,PreCycleLength,QCellCycleLength)                                !
                        !
                    END DO
                    !
                    ! Computing individual IRs
                    !
                    CALL computeIndividualIR(OptimalStrategy,VisitedStatesPre(iStatePre),iAgent,iPrice,1, &
                        numShockPeriodsPrint,PeriodsLengthPre,VisitedStatesPre(:PeriodsLengthPre), &
                        ShockStates,ShockPrices,ShockRealPrices,ShockProfits,AvgPostPrices,AvgPostProfits, &
                        ShockLength,SameCyclePrePost,PostLength)
                    !
                    ! Computing additional information
                    !
                    DO iPeriod = 1, numShockPeriodsPrint
                        !
                        p = RESHAPE(convertNumberBase(ShockStates(iPeriod)-1,numPrices,numAgents*DepthState),(/ DepthState,numAgents /))
                        ShockPrices(iPeriod,:) = p(1,:)
                        !
                        IF (iPeriod .EQ. 1) iPeriodState = VisitedStatesPre(iStatePre)
                        IF (iPeriod .GT. 1) iPeriodState = ShockStates(iPeriod-1)
                        !
                        DO jAgent = 1, numAgents
                            !
                            ! Find DynamicBR prices and Qs
                            CALL ComputeDynamicBestResponse(OptimalStrategy,iPeriodState,jAgent,delta, &
                                DynamicBRPrices(iPeriod,jAgent),DynamicBRQ(iPeriod,jAgent))
                            !
                            ! Find prices and Qs according to the strategy at convergence
                            CALL computeQCell(OptimalStrategy,iPeriodState,OptimalStrategy(iPeriodState,jAgent),jAgent,delta, &
                                OptStratQ(iPeriod,jAgent),VisitedStatesTMP,PreCycleLength,QCellCycleLength)                                !
                            !
                            ! Find StaticBR prices and PIs
                            CALL ComputeStaticBestResponse(OptimalStrategy,VisitedStatesPre(iStatePre),jAgent, &
                                StaticBRPrices(iPeriod,jAgent),PIStaticBR)
                            !
                        END DO
                        !
                    END DO
                    !
                    ! Printing results to output files
                    !
                    DO jAgent = 1, numAgents    ! Start of loop over observed agent
                        !
                        WRITE(100033,12) iGame, iPrice, &
                            NashProfits, CoopProfits, &
                            PeriodsLengthPre, &
                            AvgPrePrices, AvgPreProfits, ProfitGains, &
                            converged(iGame), timeToConvergence(iGame), iStatePre, &
                            flagEQAll,flagEQOnPath,flagEQOffPath, &
                            freqEQAll,freqEQOnPath,freqEQOffPath, &
                            flagBRAll,flagBROnPath,flagBROffPath, &
                            freqBRAll,freqBROnPath,freqBROffPath, &
                            QGapTotGame(0),QGapOnPathGame(0),QGapNotOnPathGame(0),QGapNotBRAllStatesGame(0), &
                                QGapNotBRonPathGame(0),QGapNotEqAllStatesGame(0),QGapNotEqonPathGame(0), &
                            QGapTotGame(1:numAgents),QGapOnPathGame(1:numAgents),QGapNotOnPathGame(1:numAgents), &
                                QGapNotBRAllStatesGame(1:numAgents), QGapNotBRonPathGame(1:numAgents), &
                                QGapNotEqAllStatesGame(1:numAgents),QGapNotEqonPathGame(1:numAgents), &
                            IndPrePrices(iStatePre,:), PreProfits(iStatePre,:), &
                            iAgent, jAgent, DeviationQ(jAgent), ShockLength, SameCyclePrePost, &
                            ShockPrices(:,jAgent), &
                            ShockProfits(:,jAgent), &
                            StaticBRPrices(:,jAgent), &
                            DynamicBRPrices(:,jAgent), &
                            OptStratQ(:,jAgent), &
                            DynamicBRQ(:,jAgent), &
                            PostLength, AvgPostPrices, AvgPostProfits
12                      FORMAT(I8, 1X, I11, 1X, &
                            <2*numAgents>(F12.5, 1X), &
                            I20, 1X, &
                            <numAgents>(F14.5, 1X), <numAgents>(F15.5, 1X), <numAgents>(F11.5, 1X), &
                            I9, 1X, F17.5, 1X, I19, 1X, &
                            I9, 1X, I12, 1X, I13, 1X, &
                            F9.5, 1X, F12.5, 1X, F13.5, 1X, &
                            <numAgents>(I10, 1X), <numAgents>(I13, 1X), <numAgents>(I14, 1X), &
                            <numAgents>(F10.5, 1X), <numAgents>(F13.5, 1X), <numAgents>(F14.5, 1X), &
                            F10.5, 1X, F10.5, 1X, F13.5, 1X, F18.5, 1X, F15.5, 1X, F18.5, 1X, F15.5, 1X, &
                            <numAgents>(F11.5, 1X), <numAgents>(F11.5, 1X), <numAgents>(F14.5, 1X), &
                                <numAgents>(F19.5, 1X), <numAgents>(F16.5, 1X), &
                                <numAgents>(F19.5, 1X), <numAgents>(F16.5, 1X), &
                            <numAgents>(I15, 1X), <numAgents>(F16.5, 1X), &
                            I11, 1X, I9, 1X, F12.5, 1X, I11, 1X, I16, 1X, &
                            <numShockPeriodsPrint>(I14, 1X), &
                            <numShockPeriodsPrint>(F15.5, 1X), &
                            <numShockPeriodsPrint>(I17, 1X), &
                            <numShockPeriodsPrint>(I18, 1X), &
                            <numShockPeriodsPrint>(F13.5, 1X), &
                                <numShockPeriodsPrint>(F14.5, 1X), &
                            I21, 1X, <numAgents>(F15.5, 1X), <numAgents>(F16.5, 1X))
                        !
                    END DO                  ! End of loop over observed agent
                    !
                END DO                      ! End of loop over deviating agent
                !
            END DO                          ! End of loop over pre-shock cycle states
            !
        END DO                              ! End of loop over deviation prices
        !
        ! Closing output file
        !
        CLOSE(UNIT = 100)
        !
    END DO                                  ! End of loop over games
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE ComputeDetailedAnalysis
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE DetailedAnalysis
