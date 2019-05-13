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
REAL(8), PARAMETER :: twoOverPi = 0.636619772367581343075535053490d0    ! 2/Pi
REAL(8), PARAMETER :: piOverTwo = 1.57079632679489661923132169164d0     ! Pi/2
!
! Variables
!
INTEGER :: numModels, FirstModel, numCores, numGames, itersPerYear, maxNumYears, maxIters, &
    itersInPerfMeasPeriod, printQ, codModel, &
    PerfMeasPeriodTime, numPrices, lengthFormatActionPrint, &
    typeExplorationMechanism, &
    DepthState0, DepthState, LengthStates, lengthStatesPrint, numStates, lengthStrategies, &
    typePayoffInput, &
    numAgents, numActions, numDemandParameters, &
    numExplorationParameters, computeQLearningResults, computeConvergenceResults, computePreShockCycles, &
    computeImpulseResponseToBR, computeImpulseResponseToNash, computeImpulseResponseToAll, &
    computeEquilibriumCheck, computePIGapToMaximum, computeQGapToMaximum, computeRestart
REAL(8) :: PerfMeasPeriodLength, meanNashProfit, meanCoopProfit,gammaSinghVives
CHARACTER(len = 50) :: ModelNumber, FileNameIndexStrategies, FileNameIndexLastState, FileNamePriceCycles
!
INTEGER, ALLOCATABLE :: converged(:), indexActions(:,:), indexLastState(:,:), indexStrategies(:,:), &
    cStates(:), cActions(:), priceCycles(:,:), sampledIndexStrategies(:,:), sampledPriceCycles(:,:), &
    indexStates(:,:), indexEquivalentStates(:,:), indexNashPrices(:), indexCoopPrices(:), &
    computeMixedStrategies(:)
REAL(8), ALLOCATABLE :: timeToConvergence(:), NashProfits(:), CoopProfits(:), &
    vecProfit(:,:), vecProfitQ(:,:), &
    vecAvgProfit(:), vecAvgProfitQ(:), freqStates(:,:), &
    maxValQ(:,:), meanFreqStates(:), NashPrices(:), CoopPrices(:), &
    PI(:,:), PIQ(:,:), avgPI(:), avgPIQ(:), &
    alpha(:), delta(:), & 
    meanProfit(:), seProfit(:), meanProfitGain(:), seProfitGain(:), DemandParameters(:), &
    NashMarketShares(:), CoopMarketShares(:), PricesGrids(:,:), MExpl(:), ExplorationParameters(:)
CHARACTER(len = :), ALLOCATABLE :: labelStates(:)
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
    READ(unitNumber,*) FirstModel
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
    READ(unitNumber,*) DepthState0
    DepthState = MAX(1,DepthState0)          ! Accomodates the DepthState = 0 case
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) numPrices
    !
    ! Global variables
    !
    maxIters = maxNumYears*itersPerYear
    itersInPerfMeasPeriod = INT(PerfMeasPeriodLength*itersPerYear)
    LengthStates = MAX(1,numAgents*DepthState0)
    lengthStatesPrint = LengthStates*(1+FLOOR(LOG10(DBLE(numPrices))))+LengthStates-1
    numStates = numPrices**(numAgents*DepthState0)
    numActions = numPrices**numAgents   ! Actions contain combinations of prices;
                                        ! they coincide with states when DepthState == 1
    lengthStrategies = numAgents*numStates
    lengthFormatActionPrint = FLOOR(LOG10(DBLE(numPrices)))+1
    ALLOCATE(computeMixedStrategies(numAgents))
    !
    ! Continue reading input settings
    !
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) typeExplorationMechanism
    IF (typeExplorationMechanism .EQ. 1) numExplorationParameters = 1*numAgents           ! Constant 
    IF (typeExplorationMechanism .EQ. 2) numExplorationParameters = 1*numAgents           ! Exponentially decreasing (beta)
    IF (typeExplorationMechanism .EQ. 3) numExplorationParameters = 1*numAgents           ! Exponentially decreasing (m)
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) typePayoffInput
    IF (typePayoffInput .EQ. 0) numDemandParameters = numActions                ! Pi1 matrix
    IF (typePayoffInput .EQ. 1) numDemandParameters = 1+2                       ! gamma, extend
    IF (typePayoffInput .EQ. 2) numDemandParameters = 2*numAgents+4             ! a0, ai, ci, sigma, extend
    IF (typePayoffInput .EQ. 3) numDemandParameters = 2*numAgents+4             ! a0, ai, ci, sigma = 0, extend
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
    READ(unitNumber,*) computeMixedStrategies
    READ(unitNumber,'(1X)')
    !
    ! Allocating matrices and vectors
    !
    ALLOCATE(indexActions(numActions,numAgents), &
        indexStates(numStates,LengthStates), &
        indexEquivalentStates(numStates,numAgents), &
        timeToConvergence(numGames),vecProfit(numGames,numAgents),vecProfitQ(numGames,numAgents), &
        vecAvgProfit(numGames),vecAvgProfitQ(numGames), freqStates(numStates,numGames), &
        converged(numGames),maskConverged(numGames),cStates(LengthStates),cActions(numAgents), &
        maxValQ(numStates,numAgents),meanFreqStates(numStates), &
        DemandParameters(numDemandParameters), &
        ExplorationParameters(numExplorationParameters), MExpl(numExplorationParameters), &
        alpha(numAgents),delta(numAgents),NashProfits(numAgents),CoopProfits(numAgents), &
        meanProfit(numAgents),seProfit(numAgents),meanProfitGain(numAgents),seProfitGain(numAgents), &
        PI(numActions,numAgents),PIQ(numActions,numAgents),avgPI(numActions),avgPIQ(numActions), &
        indexNashPrices(numAgents),indexCoopPrices(numAgents), &
        NashPrices(numAgents),CoopPrices(numAgents), &
        NashMarketShares(numAgents),CoopMarketShares(numAgents),PricesGrids(numPrices,numAgents), &
        IsTot(numStates,numAgents),IsOnPath(numStates,numAgents),IsBRAllStates(numStates,numAgents), &
        IsBRonPath(numStates,numAgents),IsEqAllStates(numStates,numAgents),IsEqonPath(numStates,numAgents))
    ALLOCATE(CHARACTER(len = 3+lengthStatesPrint) :: labelStates(numStates))
    !
    cStates = (/ (numPrices**i, i = LengthStates-1, 0, -1) /)
    cActions = (/ (numPrices**i, i = numAgents-1, 0, -1) /)
    !
    ! Actions contain the most recent prices of all agents. Actions and States
    ! coincide when DepthState == 1
    !
    DO iAction = 1, numActions
        !
        indexActions(iAction,:) = convertNumberBase(iAction-1,numPrices,numAgents)
        !
    END DO
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
        meanFreqStates,labelStates,indexStates, &
        NashProfits,CoopProfits, &
        alpha,MExpl,ExplorationParameters,delta,indexEquivalentStates, &
        meanProfit,seProfit,meanProfitGain,seProfitGain,DemandParameters,PI,PIQ,avgPI,avgPIQ, &
        indexNashPrices,indexCoopPrices,NashMarketShares,CoopMarketShares,PricesGrids, &
        IsTot,IsOnPath,IsBRAllStates,IsBRonPath,IsEqAllStates,IsEqonPath,computeMixedStrategies)
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
    READ(unitNumber,*) codModel, printQ, alpha, MExpl, delta, &
        DemandParameters, NashPrices, CoopPrices
    IF (typeExplorationMechanism .EQ. 2) THEN
        !
        ExplorationParameters = MExpl
        !
    ELSE IF (typeExplorationMechanism .EQ. 3) THEN
        !
        ExplorationParameters = -DBLE(itersPerYear)/DBLE(numAgents+1)* &
            LOG(1.d0-(DBLE(numPrices-1)/DBLE(numPrices))**numAgents/(DBLE(numStates*numPrices)*MExpl))
        !
    END IF
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE readModelVariables
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE globals