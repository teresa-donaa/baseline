PROGRAM main
!
USE globals
USE LearningSimulation
USE LearningSimulationRestart
USE ConvergenceResults
USE PreShockCycles
USE ImpulseResponseToBR
USE ImpulseResponseToNash
USE ImpulseResponseToAll
USE EquilibriumCheck
USE QGapToMaximum
USE PIGapToMaximum
USE IO_routines
USE QL_routines
USE PI_routines
USE generic_routines
!
IMPLICIT NONE
!
! Declaring variables
!
INTEGER :: iModel, iter, i, j, h, iGames, numGamesConverged, iAgent
REAL(8) :: meanTimeToConvergence, seTimeToConvergence, medianTimeToConvergence
REAL(8) :: meanAvgProfit, seAvgProfit
REAL(8) :: meanAvgProfitGain, seAvgProfitGain
REAL(8) :: herfFreq, entropyFreq, giniFreq, freqSymmetricStrategies
INTEGER :: countThresFrequencies(numThresFrequencies)
CHARACTER(len = 50) :: ModelName, FileName
!
! Beginning execution
!
! Opening files
!
ModelName = "figure_3.txt"
FileName = "mod_" // ModelName
OPEN(UNIT = 10001,FILE = FileName)
CALL readBatchVariables(10001)
IF (computeQLearningResults .EQ. 1) THEN
    !
    FileName = "res_" // ModelName
    OPEN(UNIT = 10002,FILE = FileName)
    !
END IF
IF (computeConvergenceResults .EQ. 1) THEN
    !
    FileName = "ConvResults_" // ModelName
    OPEN(UNIT = 100022,FILE = FileName) 
    !
END IF
IF (computePreShockCycles .EQ. 1) THEN
    !
    FileName = "PreShockCycles_" // ModelName
    OPEN(UNIT = 100021,FILE = FileName) 
    !
END IF
IF (computeImpulseResponseToBR .EQ. 1) THEN
    !
    FileName = "irToBR_" // ModelName
    OPEN(UNIT = 10003,FILE = FileName) 
    !
END IF
IF (computeImpulseResponseToNash .NE. 0) THEN
    !
    FileName = "irToNash_" // ModelName
    OPEN(UNIT = 100031,FILE = FileName) 
    !
END IF
IF (computeImpulseResponseToAll .EQ. 1) THEN
    !
    FileName = "irToAll_" // ModelName
    OPEN(UNIT = 100032,FILE = FileName) 
    !
END IF
IF (computeEquilibriumCheck .EQ. 1) THEN
    !
    FileName = "ec_" // ModelName
    OPEN(UNIT = 10004,FILE = FileName) 
    !
END IF
IF (computeQGapToMaximum .EQ. 1) THEN
    !
    FileName = "qg_" // ModelName
    OPEN(UNIT = 10006,FILE = FileName) 
    !
END IF
IF (computePIGapToMaximum .EQ. 1) THEN
    !
    FileName = "pg_" // ModelName
    OPEN(UNIT = 10007,FILE = FileName) 
    !
END IF
labelStates = computeStatesCodePrint()
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Loop over models
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
DO iModel = 1, numModels
    !
    ! Read model parameters
    !
    CALL readModelVariables(10001)
    !
    ! Creating the PI matrix
    !
    ALLOCATE(indexStrategies(lengthStrategies,numGames),indexLastState(LengthStates,numGames))
    IF (typePayoffInput .EQ. 0) CALL computePIMatricesGiven(DemandParameters,NashPrices,CoopPrices,&
        PI,NashProfits,CoopProfits, &
        indexNashPrices,indexCoopPrices,NashMarketShares,CoopMarketShares,PricesGrids)
    IF (typePayoffInput .EQ. 1) CALL computePIMatricesSinghVives(DemandParameters,NashPrices,CoopPrices,&
        PI,NashProfits,CoopProfits, &
        indexNashPrices,indexCoopPrices,NashMarketShares,CoopMarketShares,PricesGrids)
    IF (typePayoffInput .EQ. 2) CALL computePIMatricesLogit(DemandParameters,NashPrices,CoopPrices,&
        PI,NashProfits,CoopProfits, &
        indexNashPrices,indexCoopPrices,NashMarketShares,CoopMarketShares,PricesGrids)
    IF (typePayoffInput .EQ. 3) CALL computePIMatricesLogitSigma0(DemandParameters,NashPrices,CoopPrices,&
        PI,NashProfits,CoopProfits, &
        indexNashPrices,indexCoopPrices,NashMarketShares,CoopMarketShares,PricesGrids)
    PIQ = PI**2
    avgPI = SUM(PI,DIM = 2)/numAgents
    avgPIQ = avgPI**2
    !
    ! Creating I/O filenames
    !
    i = 1+INT(LOG10(DBLE(numModels)))
    WRITE(ModelNumber, "(I0.<i>, A4)") iModel, ".txt"
    FileNameIndexStrategies = "indexStrategiesTransposed_" // ModelNumber
    FileNameIndexLastState = "indexLastState_" // ModelNumber
    !
    ! Print message
    !
    WRITE(*,11) iModel, numModels, numCores
11  FORMAT('model = ', I6, ' / numModels = ', I6, ' / numCores = ', I6)  
    !
    ! Compute QL strategy 
    !
    IF (computeQLearningResults .EQ. 1) THEN
        !
        IF (computeRestart .EQ. 0) THEN
            !
            CALL computeModel(iModel,alpha,ExplorationParameters,delta, &
                converged,indexLastState,timeToConvergence)
            !
        ELSE IF (computeRestart .EQ. 1) THEN
            !
            CALL computeModelRestart(iModel,alpha,ExplorationParameters,delta, &
                converged,indexLastState,timeToConvergence)
            !
        END IF
        !
        ! Output of results:
        ! Computing indicators
        !
        CALL computeIndicators(iModel,converged,timeToConvergence, &
            numGamesConverged,meanTimeToConvergence,seTimeToConvergence,medianTimeToConvergence) 
        !
        ! Printing standard output to file
        !
        CALL PrintResFile(iModel,numGamesConverged, &
            meanTimeToConvergence,seTimeToConvergence,medianTimeToConvergence)
        !
    END IF
    !
    ! Results at convergence
    ! 
    IF (computeConvergenceResults .EQ. 1) CALL ComputeConvResults(iModel)
    !
    ! Pre-shock cycles of prices and profits
    ! 
    IF (computePreShockCycles .EQ. 1) CALL computePSCycles(iModel)
    !
    ! Impulse Response analysis to one-period deviation to static best response
    ! 
    IF (computeImpulseResponseToBR .EQ. 1) CALL computeIRAnalysisToBR(iModel)
    !
    ! Impulse Response to a permanent or transitory deviation to Nash prices
    !
    IF (computeImpulseResponseToNash .NE. 0) CALL computeIRToNashAnalysis(iModel)
    !
    ! Impulse Response analysis to one-period deviation to all prices
    !
    IF (computeImpulseResponseToAll .EQ. 1) CALL computeIRToAllAnalysis(iModel)
    !
    ! Equilibrium Check
    !
    IF (computeEquilibriumCheck .EQ. 1) CALL computeEqCheck(iModel)
    !
    ! Q and Average PI Gap w.r.t. Maximum
    !
    IF (computeQGapToMaximum .EQ. 1) CALL computeQGapToMax(iModel)
    !
    ! Q and Average PI Gap w.r.t. Maximum
    !
    IF (computePIGapToMaximum .EQ. 1) CALL computeAvgPIGapToMax(iModel)
    !
    ! Deallocate arrays
    !
    DEALLOCATE(indexStrategies,indexLastState)
    !    
    ! End of loop over models
    !
END DO
!
! Deallocating arrays
!
CALL closeBatch()
!
! Closing output files
!
CLOSE(UNIT = 10001)
IF (computeQLearningResults .EQ. 1) CLOSE(UNIT = 10002)
IF (computeConvergenceResults .EQ. 1) CLOSE(UNIT = 100022)
IF (computePreShockCycles .EQ. 1) CLOSE(UNIT = 100021)
IF (computeImpulseResponseToBR .EQ. 1) CLOSE(UNIT = 10003)
IF (computeImpulseResponseToNash .NE. 0) CLOSE(UNIT = 100031)
IF (computeImpulseResponseToAll .EQ. 1) CLOSE(UNIT = 100032)
IF (computeEquilibriumCheck .EQ. 1) CLOSE(UNIT = 10004)
IF (computeQGapToMaximum .EQ. 1) CLOSE(UNIT = 10006)
IF (computePIGapToMaximum .EQ. 1) CLOSE(UNIT = 10007)
!
! End of execution
!
END PROGRAM main