MODULE globals
!
USE generic_routines
!
! Declare global parameters and variables 
!
IMPLICIT NONE
!
! Parameters
!
INTEGER, PARAMETER :: numQLParameters = 3
INTEGER, PARAMETER :: numThresFrequencies = 5
INTEGER, PARAMETER :: thresFrequencies(numThresFrequencies) = (/ 1, 3, 5, 10, 20 /)
INTEGER, PARAMETER :: QTrajectoriesPeriodPrint = 120
REAL(8), PARAMETER :: twoOverPi = 0.636619772367581343075535053490d0    ! 2/Pi
REAL(8), PARAMETER :: piOverTwo = 1.57079632679489661923132169164d0     ! Pi/2
!
! Variables
!
INTEGER :: numModels, numCores, numGames, itersPerYear, maxNumYears, maxIters, &
    itersInPerfMeasPeriod, printExp, printQ, &
    PerfMeasPeriodTime, numPrices, lengthFormatActionPrint, numNashStrategies, &
    numPrintStrategies, typeExplorationMechanism, useNashStrategies, useOtherStrategies, &
    depthState, lengthStates, lengthStatesPrint, numStates, lengthStrategies, &
    lengthStrategiesPrint, typePayoffInput, &
    numAgents, numActions, numDemandParameters, numOtherStrategies, &
    numExplorationParameters, computeQLearningResults, computeConvergenceResults, computePreShockCycles, &
    computeImpulseResponseToBR, computeImpulseResponseToNash, computeImpulseResponseToAll, &
    computeEquilibriumCheck, computePIGapToMaximum, computeQGapToMaximum, computeRestart
REAL(8) :: PerfMeasPeriodLength, meanNashProfit, meanCoopProfit,gammaSinghVives
CHARACTER(len = 50) :: ModelNumber, FileNameIndexStrategies, FileNameIndexLastState
!
INTEGER, ALLOCATABLE :: converged(:), indexActions(:,:), indexLastState(:,:), indexStrategies(:,:), &
    cStates(:), cActions(:), &
    indexNashStrategies(:,:), indexPrintStrategies(:,:), &
    indexStates(:,:), indexEquivalentStates(:,:), indexNashPrices(:), indexCoopPrices(:)
REAL(8), ALLOCATABLE :: timeToConvergence(:), NashProfits(:), CoopProfits(:), &
    vecProfit(:,:), vecProfitQ(:,:), &
    vecAvgProfit(:), vecAvgProfitQ(:), freqStates(:,:), &
    maxValQ(:,:), meanFreqStates(:), NashPrices(:), CoopPrices(:), &
    avgFreqStatesMostObsStrategies(:,:), PI(:,:), PIQ(:,:), avgPI(:), avgPIQ(:), &
    freqMostObsStrategies(:), meanProfGainMostObsStrategies(:,:), meanAvgProfGainMostObsStrategies(:), &
    avgTTCMostObsStrategies(:), alpha(:), delta(:), & 
    meanProfit(:), seProfit(:), meanProfitGain(:), seProfitGain(:), DemandParameters(:), &
    NashMarketShares(:), CoopMarketShares(:), PricesGrids(:,:), MExpl(:), ExplorationParameters(:)
CHARACTER(len = :), ALLOCATABLE :: labelStates(:), labelStrategies(:)
LOGICAL, ALLOCATABLE :: maskConverged(:), IsTot(:,:), IsOnPath(:,:), IsBRAllStates(:,:), &
    IsBRonPath(:,:), IsEqAllStates(:,:), IsEqonPath(:,:) 
!
CONTAINS
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE readBatchVariables ( unitNumber )
    !
    ! Reads input variables
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: unitNumber
    !
    ! Declaring local variables
    !
    INTEGER :: i, iAgent, iPrices, jPrices, iState, iAction, jState
    INTEGER, ALLOCATABLE :: switchedState(:), indA1(:), indA2(:)
    !
    ! Beginning execution
    !
    READ(unitNumber,'(1X)') 
    READ(unitNumber,*) numModels
    READ(unitNumber,'(1X)') 
    READ(unitNumber,*) numCores
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) numGames
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) itersPerYear
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) maxNumYears
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) PerfMeasPeriodTime
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) PerfMeasPeriodLength
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) numAgents
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) depthState
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) numPrices
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) useNashStrategies
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) numNashStrategies
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) useOtherStrategies
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) numOtherStrategies
    !
    ! Global variables
    !
    maxIters = maxNumYears*itersPerYear
    itersInPerfMeasPeriod = INT(PerfMeasPeriodLength*itersPerYear)
    lengthStates = numAgents*depthState
    lengthStatesPrint = lengthStates*(1+FLOOR(LOG10(DBLE(numPrices))))+lengthStates-1
    numStates = numPrices**(numAgents*depthState)
    numActions = numPrices**numAgents   ! Actions contain combinations of prices;
                                        ! they coincide with states when depthState == 1
    lengthStrategies = numAgents*numStates
    lengthFormatActionPrint = FLOOR(LOG10(DBLE(numPrices)))+1
    lengthStrategiesPrint = lengthFormatActionPrint*lengthStrategies+lengthStrategies-1
    IF (useNashStrategies .EQ. 0) numPrintStrategies = numOtherStrategies
    IF (useNashStrategies .EQ. 1) numPrintStrategies = numNashStrategies+numOtherStrategies
    numPrintStrategies = MAX(1,numPrintStrategies)
    !
    ! Continue reading input settings
    !
    IF (useNashStrategies .EQ. 1) THEN
        !        
        ALLOCATE(indexNashStrategies(lengthStrategies,numNashStrategies))
        READ(unitNumber,'(1X)')
        DO i = 1, numNashStrategies
            !
            READ(unitNumber,*) indexNashStrategies(:,i)
            !
        END DO
        !
    ELSE IF (useNashStrategies .EQ. 0) THEN
        !
        READ(unitNumber,'(<numNashStrategies>(/))')
        !
    END IF
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) typeExplorationMechanism
    IF (typeExplorationMechanism .EQ. 1) numExplorationParameters = 1*numAgents           ! Constant 
    IF (typeExplorationMechanism .EQ. 2) numExplorationParameters = 1*numAgents           ! Exponentially decreasing
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) typePayoffInput
    IF (typePayoffInput .EQ. 0) numDemandParameters = numActions                ! Pi1 matrix
    IF (typePayoffInput .EQ. 1) numDemandParameters = 1                         ! sigma
    IF (typePayoffInput .EQ. 2) numDemandParameters = 2*numAgents+4             ! a0, ai, ci, sigma, extend
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) computeQLearningResults
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) computeConvergenceResults
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) computePreShockCycles
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) computeImpulseResponseToBR
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) computeImpulseResponseToNash
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) computeImpulseResponseToAll
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) computeEquilibriumCheck
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) computeQGapToMaximum
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) computePIGapToMaximum
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) computeRestart
    READ(unitNumber,'(1X)')
    !
    ! Allocating matrices and vectors
    !
    ALLOCATE(indexActions(numActions,numAgents), &
        indexStates(numStates,lengthStates), &
        indexEquivalentStates(numStates,numAgents), &
        timeToConvergence(numGames),vecProfit(numGames,numAgents),vecProfitQ(numGames,numAgents), &
        vecAvgProfit(numGames),vecAvgProfitQ(numGames), freqStates(numStates,numGames), &
        converged(numGames),maskConverged(numGames),cStates(lengthStates),cActions(numAgents), &
        maxValQ(numStates,numAgents),meanFreqStates(numStates), &
        avgFreqStatesMostObsStrategies(numStates,numPrintStrategies), &
        indexPrintStrategies(lengthStrategies,numPrintStrategies), &
        freqMostObsStrategies(numPrintStrategies),DemandParameters(numDemandParameters), &
        ExplorationParameters(numExplorationParameters), MExpl(numExplorationParameters), &
        meanProfGainMostObsStrategies(numPrintStrategies,numAgents), &
        meanAvgProfGainMostObsStrategies(numPrintStrategies), &
        avgTTCMostObsStrategies(numPrintStrategies), &
        alpha(numAgents),delta(numAgents),NashProfits(numAgents),CoopProfits(numAgents), &
        meanProfit(numAgents),seProfit(numAgents),meanProfitGain(numAgents),seProfitGain(numAgents), &
        PI(numActions,numAgents),PIQ(numActions,numAgents),avgPI(numActions),avgPIQ(numActions), &
        indexNashPrices(numAgents),indexCoopPrices(numAgents), &
        NashPrices(numAgents),CoopPrices(numAgents), &
        NashMarketShares(numAgents),CoopMarketShares(numAgents),PricesGrids(numPrices,numAgents), &
        IsTot(numStates,numAgents), IsOnPath(numStates,numAgents), IsBRAllStates(numStates,numAgents), &
        IsBRonPath(numStates,numAgents), IsEqAllStates(numStates,numAgents), IsEqonPath(numStates,numAgents))
    ALLOCATE(CHARACTER(len = 3+lengthStatesPrint) :: labelStates(numStates))
    ALLOCATE(CHARACTER(len = 9) :: labelStrategies(numPrintStrategies))
    !
    cStates = (/ (numPrices**i, i = lengthStates-1, 0, -1) /)
    cActions = (/ (numPrices**i, i = numAgents-1, 0, -1) /)
    !
    ! Actions contain the most recent prices of all agents. Actions and States
    ! coincide when depthState == 1
    !
    DO iAction = 1, numActions
        !
        indexActions(iAction,:) = convertNumberBase(iAction-1,numPrices,numAgents)
        !
    END DO
    !
    ! Find indexes of equivalent states, i.e. states that are the same when considered 
    ! from the point of view of the two agents
    ! E.g.: When numPrices = 3 and depthState = 2, state (1.2)-(1.3) for agent 1 is the same
    ! as state (2.1)-(3.1) for agent 2
    !
    IF (numAgents .EQ. 2) THEN
        !
        ALLOCATE(switchedState(lengthStates),indA1(depthState),indA2(depthState))
        DO iState = 1, numStates
            !
            indexStates(iState,:) = convertNumberBase(iState-1,numPrices,lengthStates)
            !
        END DO
        !
        indA1 = (/ (i, i = 1, lengthStates, numAgents) /)
        indA2 = indA1+1
        DO iState = 1, numStates
            !
            indexEquivalentStates(iState,1) = iState
            switchedState(indA1) = indexStates(iState,indA2)
            switchedState(indA2) = indexStates(iState,indA1)
            DO jState = 1, numStates
                !
                IF (ALL(switchedState .EQ. indexStates(jState,:))) THEN
                    !
                    indexEquivalentStates(iState,2) = jState
                    EXIT
                    !
                END IF
                !
            END DO
            !
        END DO
        DEALLOCATE(switchedState,indA1,indA2)        
    END IF
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE readBatchVariables
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE closeBatch ( )
    !
    ! Reads input variables
    !
    IMPLICIT NONE
    !
    ! Beginning execution
    !
    DEALLOCATE(freqStates,indexActions,timeToConvergence,vecProfit,vecProfitQ, &
        vecAvgProfit,vecAvgProfitQ,converged,maskConverged,cStates,cActions,maxValQ, &
        meanFreqStates,avgFreqStatesMostObsStrategies,labelStates,indexStates, &
        indexPrintStrategies,freqMostObsStrategies,meanProfGainMostObsStrategies, &
        meanAvgProfGainMostObsStrategies,avgTTCMostObsStrategies,NashProfits,CoopProfits, &
        labelStrategies,alpha,MExpl,ExplorationParameters,delta,indexEquivalentStates, &
        meanProfit,seProfit,meanProfitGain,seProfitGain,DemandParameters,PI,PIQ,avgPI,avgPIQ, &
        indexNashPrices,indexCoopPrices,NashMarketShares,CoopMarketShares,PricesGrids, &
        IsTot,IsOnPath,IsBRAllStates,IsBRonPath,IsEqAllStates,IsEqonPath)
    IF (useNashStrategies .EQ. 1) DEALLOCATE(indexNashStrategies)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE closeBatch
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE readModelVariables ( unitNumber )
    !
    ! Reads input variables
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: unitNumber
    !
    ! Declaring local variables
    !
    INTEGER :: i
    !
    ! Beginning execution
    !
    READ(unitNumber,*) i, printExp, printQ, alpha, MExpl, delta, &
        DemandParameters, NashPrices, CoopPrices
    ExplorationParameters = -DBLE(itersPerYear)/DBLE(numAgents+1)* &
        LOG(1.d0-(DBLE(numPrices-1)/DBLE(numPrices))**numAgents/(DBLE(numStates*numPrices)*MExpl))
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE readModelVariables
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE globals