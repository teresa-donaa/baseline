PROGRAM main
!
USE globals
USE LearningSimulation
USE LearningSimulationRestart
USE ConvergenceResults
USE PolicyImprovement
USE ImpulseResponse
USE EquilibriumCheck
USE QGapToMaximum
USE LearningTrajectory
USE DetailedAnalysis
USE QL_routines
USE PI_routines
USE generic_routines
!
IMPLICIT NONE
!
! Declaring variables and parameters
!
INTEGER :: iModel, i, iAgent
CHARACTER(len = 50) :: FileName
!
! Beginning execution
!
! Opening files
!
ModelName = "table_A5"
FileName = TRIM("A_mod_" // ModelName) // ".txt"
!
OPEN(UNIT = 10001,FILE = FileName)
CALL readBatchVariables(10001)
!
IF (SwitchQLearningResults .EQ. 1) THEN
    !
    FileName = TRIM("A_res_" // ModelName) // ".txt"
    OPEN(UNIT = 10002,FILE = FileName)
    !
END IF
IF (SwitchConvergenceResults .EQ. 1) THEN
    !
    FileName = TRIM("A_convResults_" // ModelName) // ".txt"
    OPEN(UNIT = 100022,FILE = FileName)
    !
END IF
IF (SwitchImpulseResponseToBR .EQ. 1) THEN
    !
    FileName = TRIM("A_irToBR_" // ModelName) // ".txt"
    OPEN(UNIT = 10003,FILE = FileName)
    !
END IF
IF (SwitchImpulseResponseToNash .GE. 1) THEN
    !
    FileName = TRIM("A_irToNash_" // ModelName) // ".txt"
    OPEN(UNIT = 100031,FILE = FileName)
    !
END IF
IF (SwitchImpulseResponseToAll .EQ. 1) THEn
    !
    FileName = TRIM("A_irToAll_" // ModelName) // ".txt"
    OPEN(UNIT = 100032,FILE = FileName)
    !
END IF
IF (SwitchEquilibriumCheck .EQ. 1) THEN
    !
    FileName = TRIM("A_ec_" // ModelName) // ".txt"
    OPEN(UNIT = 10004,FILE = FileName)
    !
END IF
IF (SwitchQGapToMaximum .EQ. 1) THEN
    !
    FileName = TRIM("A_qg_" // ModelName) // ".txt"
    OPEN(UNIT = 10006,FILE = FileName)
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
    ! When using trained Q matrices, store the name of the directory containing them
    !
    DO iAgent = 1, numAgents
        !
        IF (typeQInitialization(iAgent) .EQ. 'T') THEN
            !
            QFileFolderName(iAgent) = 'trained_Q/'
            !
        END IF
        !
    END DO
    !
    ! Creating the PI matrix
    !
    IF (typePayoffInput .EQ. 0) CALL computePIMatricesGiven(DemandParameters,NashPrices,CoopPrices,&
        PI,NashProfits,CoopProfits, &
        indexNashPrices,indexCoopPrices,NashMarketShares,CoopMarketShares,PricesGrids)
    IF (typePayoffInput .EQ. 1) CALL computePIMatricesSinghVives(DemandParameters,NashPrices,CoopPrices,&
        PI,NashProfits,CoopProfits, &
        indexNashPrices,indexCoopPrices,NashMarketShares,CoopMarketShares,PricesGrids)
    IF (typePayoffInput .EQ. 2) CALL computePIMatricesLogit(DemandParameters,NashPrices,CoopPrices,&
        PI,NashProfits,CoopProfits, &
        indexNashPrices,indexCoopPrices,NashMarketShares,CoopMarketShares,PricesGrids)
    IF (typePayoffInput .EQ. 3) CALL computePIMatricesLogitMu0(DemandParameters,NashPrices,CoopPrices,&
        PI,NashProfits,CoopProfits, &
        indexNashPrices,indexCoopPrices,NashMarketShares,CoopMarketShares,PricesGrids)
    PIQ = PI**2
    avgPI = SUM(PI,DIM = 2)/numAgents
    avgPIQ = avgPI**2
    !
    ! Computing profit gains
    !
    DO iAgent = 1, numAgents
        !
        PG(:,iAgent) = (PI(:,iAgent)-NashProfits(iAgent))/(CoopProfits(iAgent)-NashProfits(iAgent))
        !
    END DO
    PGQ = PG**2
    avgPG = SUM(PG,DIM = 2)/numAgents
    avgPGQ = avgPG**2
    !
    ! Creating I/O filenames
    !
    WRITE(ModelNumber, "(I0.<LengthFormatTotModelsPrint>, A4)") codModel, ".txt"
    FileNameInfoModel = "InfoModel_" // ModelNumber
    !
    ! Print message
    !
    WRITE(*,11) iModel, numModels, numCores
11  FORMAT('model = ', I6, ' / numModels = ', I6, ' / numCores = ', I6)  
    !
    ! Compute QL strategy 
    !
    IF (SwitchQLearningResults .EQ. 1) THEN
        !
        IF (SwitchRestart .EQ. 0) THEN
            !
            CALL computeModel(iModel,codModel,alpha,ExplorationParameters,delta)
            !
        ELSE IF (SwitchRestart .EQ. 1) THEN
            !
            CALL computeModelRestart(iModel,codModel,alpha,ExplorationParameters,delta)
            !
        END IF
        !
    END IF
    !
    ! Results at convergence
    ! 
    IF (SwitchConvergenceResults .EQ. 1) CALL ComputeConvResults(iModel)
!@SP
!CALL ImprovePolicy ( )
!@SP
    !
    ! Impulse Response analysis to one-period deviation to static best response
    ! NB: The last argument in computeIRAnalysis is "IRType", and it's crucial:
    ! IRType < 0 : One-period deviation to the price IRType
    ! IRType = 0 : One-period deviation to static BR
    ! IRType > 0 : IRType-period deviation to Nash 
    ! 
    IF (SwitchImpulseResponseToBR .EQ. 1) CALL computeIRAnalysis(iModel,10003,0)
    !
    ! Impulse Response to a permanent or transitory deviation to Nash prices
    !
    IF (SwitchImpulseResponseToNash .GE. 1) CALL computeIRAnalysis(iModel,100031,SwitchImpulseResponseToNash)
    !
    ! Impulse Response analysis to one-period deviation to all prices
    !
    IF (SwitchImpulseResponseToAll .EQ. 1) THEN
        !
        DO i = 1, numPrices
            !
            CALL computeIRAnalysis(iModel,100032,-i)
            !
        END DO
        !
    END IF
    !
    ! Equilibrium Check
    !
    IF (SwitchEquilibriumCheck .EQ. 1) CALL computeEqCheck(iModel)
    !
    ! Q Gap w.r.t. Maximum
    !
    IF (SwitchQGapToMaximum .EQ. 1) CALL computeQGapToMax(iModel)
    !
    ! Learning Trajectory analysis
    !
    IF (ParamsLearningTrajectory(1) .GT. 0) &
        CALL ComputeLearningTrajectory(iModel,codModel,alpha,ExplorationParameters,delta)
    !
    ! Detailed Impulse Response analysis to one-period deviation to all prices
    !
    IF (SwitchDetailedAnalysis .EQ. 1) CALL ComputeDetailedAnalysis(iModel)
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
IF (SwitchQLearningResults .EQ. 1) CLOSE(UNIT = 10002)
IF (SwitchConvergenceResults .EQ. 1) CLOSE(UNIT = 100022)
IF (SwitchImpulseResponseToBR .EQ. 1) CLOSE(UNIT = 10003)
IF (SwitchImpulseResponseToNash .GE. 1) CLOSE(UNIT = 100031)
IF (SwitchImpulseResponseToAll .EQ. 1) CLOSE(UNIT = 100032)
IF (SwitchEquilibriumCheck .EQ. 1) CLOSE(UNIT = 10004)
IF (SwitchQGapToMaximum .EQ. 1) CLOSE(UNIT = 10006)
!
! End of execution
!
END PROGRAM main