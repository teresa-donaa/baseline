MODULE LearningTrajectory
!
USE globals
USE generic_routines
USE QL_routines
USE LearningSimulation
USE ImpulseResponse
USE omp_lib
USE ifport
!
! Computes Monte Carlo Q-Learning Trajectories
!
IMPLICIT NONE
!
CONTAINS
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE computeLearningTrajectory ( iModel, codModel, alpha, ExplorationParameters, delta )
    !
    ! Computes statistics for one model
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: iModel, codModel
    REAL(8), DIMENSION(numAgents), INTENT(IN) :: alpha, delta
    REAL(8), DIMENSION(numExplorationParameters) :: ExplorationParameters
    !
    ! Declaring local variable
    !
    INTEGER :: idumIP, ivIP(32), iyIP, idum2IP, idum, iv(32), iy, idum2, idumQ, ivQ(32), iyQ, idum2Q
    INTEGER :: iIters, i, j, h, l, iGame, iAgent, jAgent, iCycle, iPrice, ss0, ss1
    INTEGER, DIMENSION(DepthState,numAgents) :: p, pLT1
    INTEGER, DIMENSION(numAgents) :: pPrime, pPrimeLT1, pPrimeLT2
    INTEGER :: state, statePrime, actionPrime
    INTEGER, DIMENSION(numStates,numAgents) :: strategy, strategyPrime
    INTEGER(8) :: numGames_I8
    INTEGER :: CycleLengthGame, CycleStatesGame(numPeriods)
    INTEGER, DIMENSION(numPeriods) :: NUVisitedStates
    INTEGER :: NUPreCycleLength, NUCycleLength
    !
    REAL(8), DIMENSION(numStates,numPrices,numAgents) :: Q
    REAL(8) :: uIniPrice(DepthState,numAgents,numGames), uExploration(2,numAgents)
    REAL(8) :: u(2), eps(numAgents)
    REAL(8) :: newq, oldq, profitgain(numAgents)
    REAL(8), DIMENSION(ParamsLearningTrajectory(1),numGames) :: PGmat, ICmat, IRmat
    REAL(8), DIMENSION(ParamsLearningTrajectory(1),9) :: PGss, ICss, IRss
    REAL(8) :: tmp_r(numGames), PGsum, ICsum, IRsum
    REAL(8) :: QRowValues(numPrices), avgPRatio
    REAL(8) :: NUPIStaticBR
    !
    CHARACTER(len = 200) :: LTrajectoryFileName
    !
    ! Beginning execution
    !
    ! Reading strategies and states at convergence from file
    !
    CALL ReadInfoModel()
    !
    ! Initializing various quantities
    !
    PGmat = 0.d0
    PGss = 0.d0
    !
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Loop over numGames
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !$ CALL OMP_SET_NUM_THREADS(numCores)
    !
    ! Generating uIniPrice
    !
    idumIP = -1
    idum2IP = 123456789
    ivIP = 0
    iyIP = 0
    CALL generate_uIniPrice(uIniPrice,idumIP,ivIP,iyIP,idum2IP)  
    !
    ! Starting loop over games
    !
    !$omp parallel do &
    !$omp private(idum,iv,iy,idum2,idumQ,ivQ,iyQ,idum2Q,Q,maxValQ, &
    !$omp   strategyPrime,pPrime,p,statePrime,actionPrime,iIters,iCycle,ss0,ss1, &
    !$omp   pLT1,pPrimeLT1,pPrimeLT2, &
    !$omp   state,strategy,eps,uExploration,u,oldq,newq,iAgent,jAgent,CycleLengthGame,CycleStatesGame, &
    !$omp   NUVisitedStates,NUPreCycleLength,NUCycleLength,NUPIStaticBR, &
    !$omp   PGsum,QRowValues,ICsum,avgPRatio,IRsum) &
    !$omp firstprivate(numGames,PI,PG,delta,uIniPrice,ExplorationParameters,itersPerYear,alpha, &
    !$omp   itersInPerfMeasPeriod,maxIters,profitgain)
    DO iGame = 1, numGames
        !
        PRINT*, 'Game = ', iGame, ' started'
        CycleLengthGame = CycleLength(iGame)
        CycleStatesGame = 0
        CycleStatesGame(1:CycleLengthGame) = CycleStates(:,iGame)
        !
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Learning phase
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !
        ! Initializing random number generators
        !
        idum = -iGame
        idum2 = 123456789
        iv = 0
        iy = 0
        !
        idumQ = -iGame
        idum2Q = 123456789
        ivQ = 0
        iyQ = 0
        !
        ! Initializing Q matrices
        !
        !$omp critical
        CALL initQMatrices(iGame,idumQ,ivQ,iyQ,idum2Q,PI,delta,Q,maxValQ,strategyPrime)
        !$omp end critical
        strategy = strategyPrime
        !
        ! Randomly initializing prices and state
        !
        CALL initState(uIniPrice(:,:,iGame),p,statePrime,actionPrime)
        state = statePrime
        !
        ! Loop
        !
        IF ((typeExplorationMechanism .EQ. 2) .OR. (typeExplorationMechanism .EQ. 3)) eps = ExplorationParameters(:numAgents)
        IF (typeExplorationMechanism .EQ. 4) eps = ExplorationParameters(:numAgents)
        PGsum = 0.d0
        ICsum = 0.d0
        IRsum = 0.d0
        DO iIters = 1, ParamsLearningTrajectory(1)*ParamsLearningTrajectory(2)
            !
            ! Generating exploration random numbers
            !            
            CALL generateUExploration(uExploration,idum,iv,iy,idum2)  
            !
            ! Compute pPrime by balancing exploration vs. exploitation
            !
            CALL computePPrime(ExplorationParameters,uExploration,strategyPrime,state,iIters,pPrime,Q,eps)
            !
            ! Defining the new state
            !
            IF (DepthState .GT. 1) p(2:DepthState,:) = p(1:DepthState-1,:)
            p(1,:) = pPrime
            statePrime = computeStateNumber(p)
            actionPrime = computeActionNumber(pPrime)
            !
            ! Store LTrajectories
            !
            ! Profit Gain
            profitgain = PG(actionPrime,:)
            PGsum = PGsum+SUM(profitgain)/DBLE(numAgents)
            !
            IF (MOD(iIters,ParamsLearningTrajectory(2)) .EQ. 0) THEN
                !
                ! Profit Gain
                PGmat(iIters/ParamsLearningTrajectory(2),iGame) = &
                    PGmat(iIters/ParamsLearningTrajectory(2),iGame)+PGsum/DBLE(ParamsLearningTrajectory(2))
                PGsum = 0.d0
                !
                ! Incentive Compatibility
                DO iCycle = 1, CycleLengthGame      ! Loop over states in path at convergence
                    !
                    ss1 = CycleStatesGame(iCycle)
                    DO iAgent = 1, numAgents        ! Loop over agents
                        !
                        DO iPrice = 1, numPrices    ! Loop over prices
                            !
                            ! Compute row of true Q matrix
                            CALL computeQCell(strategy,ss1,iPrice,iAgent,delta, &
                                QRowValues(iPrice),NUVisitedStates,NUPreCycleLength,NUCycleLength)
                            !
                        END DO
                        IF (AreEqualReals(MAXVAL(QRowValues),QrowValues(strategy(ss1,iAgent)))) &
                            ICsum = ICsum+1.d0/DBLE(numAgents*CycleLengthGame)
                        !
                    END DO
                    !
                END DO
                ICmat(iIters/ParamsLearningTrajectory(2),iGame) = ICsum
                ICsum = 0.d0
                !
                ! Punishment from nondeviating agent in period t+2
                DO iCycle = 1, CycleLengthGame      ! Loop over states in path at convergence
                    !
                    DO iAgent = 1, numAgents        ! Loop over agents
                        !
                        ! Period t: Initial state
                        ss0 = CycleStatesGame(iCycle)   
                        pLT1 = RESHAPE(convertNumberBase(ss0-1,numPrices,LengthStates),(/ DepthState,numAgents /))
                        !
                        ! Period t+1: Shock to static best response
                        pPrimeLT1 = strategy(ss0,:)
                        IF (DepthState .GT. 1) pLT1(2:DepthState,:) = pLT1(1:DepthState-1,:)
                        pLT1(1,:) = pPrimeLT1
                        CALL ComputeStaticBestResponse(strategy,ss0,iAgent,pLT1(1,iAgent),NUPIStaticBR)
                        ss1 = computeStateNumber(pLT1)
                        !
                        ! Period t+2
                        pPrimeLT2 = strategy(ss1,:)
                        avgPRatio = 0.d0
                        DO jAgent = 1, numAgents
                            !
                            IF (jAgent .EQ. iAgent) CYCLE
                            avgPRatio = avgPRatio+PricesGrids(pPrimeLT2(jAgent),jAgent)/PricesGrids(pPrimeLT1(jAgent),jAgent)
                            !
                        END DO
                        avgPRatio = avgPRatio/DBLE(numAgents-1)
                        IRsum = IRsum+avgPRatio/DBLE(numAgents*CycleLengthGame)                    
                        !
                    END DO
                    !
                END DO
                IRmat(iIters/ParamsLearningTrajectory(2),iGame) = IRsum
                IRsum = 0.d0
                !
            END IF
            !
            ! Each agent collects his payoff and updates
            !
            DO iAgent = 1, numAgents
                !
                ! Q matrices and strategies update
                !
                oldq = Q(state,pPrime(iAgent),iAgent)
                newq = oldq+alpha(iAgent)*(PI(actionPrime,iAgent)+delta(iAgent)*maxValQ(statePrime,iAgent)-oldq)
                Q(state,pPrime(iAgent),iAgent) = newq
                IF (newq .GT. maxValQ(state,iAgent)) THEN
                    !
                    maxValQ(state,iAgent) = newq
                    IF (strategyPrime(state,iAgent) .NE. pPrime(iAgent)) strategyPrime(state,iAgent) = pPrime(iAgent)
                    !
                END IF
                IF ((newq .LT. maxValQ(state,iAgent)) .AND. (strategyPrime(state,iAgent) .EQ. pPrime(iAgent))) THEN
                    !
                    CALL MaxLocBreakTies(numPrices,Q(state,:,iAgent),idumQ,ivQ,iyQ,idum2Q, &
                        maxValQ(state,iAgent),strategyPrime(state,iAgent))
                    !
                END IF
                !
            END DO
            !
            ! Update and iterate
            !
            strategy(state,:) = strategyPrime(state,:)
            state = statePrime
            !
            ! End of loop over iterations
            !
        END DO
        !
        ! End of loop over games
        !
    END DO
    !$omp end parallel do
    !
    ! Compute statistics on PG and IC trajectories
    !
    PGss = ComputeRowSummaryStatistics(ParamsLearningTrajectory(1),numGames,PGmat)
    ICss = ComputeRowSummaryStatistics(ParamsLearningTrajectory(1),numGames,ICmat)
    IRss = ComputeRowSummaryStatistics(ParamsLearningTrajectory(1),numGames,IRmat)
    !
    ! File name
    !
    LTrajectoryFileName = 'LTrajectories_' // ModelNumber 
    OPEN(UNIT = 9,FILE = LTrajectoryFileName,STATUS = "REPLACE")
    WRITE(9,8) 
8   FORMAT('         iter ', &
        '          AvgPG            SdPG           MinPG         q0025PG          q025PG           q05PG          q075PG         q0975PG           MaxPG ', &
        '          AvgIC            SdIC           MinIC         q0025IC          q025IC           q05IC          q075IC         q0975IC           MaxIC ', &
        '          AvgIR            SdIR           MinIR         q0025IR          q025IR           q05IR          q075IR         q0975IR           MaxIR ')
    DO i = 1, ParamsLearningTrajectory(1)
        !
        WRITE(9,9) i*ParamsLearningTrajectory(2), &
            PGss(i,:), ICss(i,:), IRss(i,:)
        !
    END DO
9   FORMAT(1X, I12, 1X, <9>(F15.5, 1X), <9>(F15.5, 1X), <9>(F15.5, 1X))
    CLOSE(UNIT = 9)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE computeLearningTrajectory
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE LearningTrajectory
