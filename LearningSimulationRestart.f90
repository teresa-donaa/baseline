MODULE LearningSimulationRestart
!
USE globals
USE QL_routines
USE LearningSimulation
USE omp_lib
!
! Computes Monte Carlo Q-Learning simulations
!
IMPLICIT NONE
!
CONTAINS
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE computeModelRestart ( iModel, codModel, alpha, ExplorationParameters, delta, &
        converged, indexLastState, timeToConvergence )
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
    INTEGER, DIMENSION(numGames), INTENT(OUT) :: converged
    INTEGER, INTENT(OUT) :: indexLastState(LengthStates,numGames)
    REAL(8), DIMENSION(numGames), INTENT(OUT) :: timeToConvergence
    !
    ! Declaring local variable
    !
    REAL(8), DIMENSION(numAgents) :: pricesGridsPrime
    INTEGER :: idumIP, ivIP(32), iyIP, idum2IP, idum, iv(32), iy, idum2, idumQ, ivQ(32), iyQ, idum2Q
    REAL(8), DIMENSION(numStates,numPrices,numAgents) :: Q
    REAL(8) :: uIniPrice(DepthState,numAgents,numGames), uExploration(2,numAgents)
    REAL(8) :: u(2), eps(numAgents)
    REAL(8) :: newq, oldq, profitgain
    INTEGER :: iIters, i, j, iGames, iItersInStrategy, convergedGame
    INTEGER :: state, statePrime, actionPrime
    INTEGER, DIMENSION(numStates,numAgents) :: strategy, strategyPrime
    INTEGER :: pPrime(numAgents), p(DepthState,numAgents)
    INTEGER :: iAgent, iState, iPrice, jAgent
    INTEGER :: minIndexStrategies, maxIndexStrategies
    CHARACTER(len = 25) :: QFileName
    CHARACTER(len = 5) :: iGamesChar
    CHARACTER(len = 5) :: codModelChar
    !
    ! Beginning execution
    !
    ! Initializing various quantities
    !
    converged = 0    
    indexLastState = 0
    timeToConvergence = 0.d0
    WRITE(codModelChar,'(I0.5)') codModel
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
    !$omp   strategyPrime,pPrime,p,statePrime,actionPrime,iIters,iItersInStrategy, &
    !$omp   convergedGame,pricesGridsPrime, &
    !$omp   state,strategy,eps,uExploration,u,oldq,newq,iAgent,iState,iPrice,jAgent, &
    !$omp   QFileName,iGamesChar) &
    !$omp firstprivate(numGames,PI,delta,uIniPrice,ExplorationParameters,itersPerYear,alpha, &
    !$omp   itersInPerfMeasPeriod,maxIters,printQ,printP,profitgain,codModelChar)
    DO iGames = 1, numGames
        !
        PRINT*, 'iGames = ', iGames
        !
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Start of first learning phase, random initialization of Q matrix
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !
        ! Initializing random number generator
        !
        idum = -iGames
        idum2 = 123456789
        iv = 0
        iy = 0
        !
        ! Initializing Q matrices
        !
        !$omp critical
        CALL initQMatrices(iGames,idumQ,ivQ,iyQ,idum2Q,PI,delta,Q,maxValQ,strategyPrime)
        !$omp end critical
        strategy = strategyPrime
        !
        ! Randomly initializing prices and state
        !
        CALL initState(uIniPrice(:,:,iGames),p,statePrime,actionPrime)
        state = statePrime
        !
        ! Loop
        !
        iIters = 0
        iItersInStrategy = 0
        DO 
            !
            ! Iterations counter
            !
            iIters = iIters+1
            !
            ! Generating exploration random numbers
            !            
            CALL generateUExploration(uExploration,idum,iv,iy,idum2)  
            !
            ! Compute pPrime by balancing exploration vs. exploitation
            !
            CALL computePPrime(ExplorationParameters,uExploration,strategyPrime,state,iIters,pPrime,Q)
            !
            ! Defining the new state
            !
            IF (DepthState .GT. 1) p(2:DepthState,:) = p(1:DepthState-1,:)
            p(1,:) = pPrime
            statePrime = computeStateNumber(p)
            actionPrime = computeActionNumber(pPrime)
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
!@SP                    maxValQ(state,iAgent) = MAXVAL(Q(state,:,iAgent))
!@SP                    strategyPrime(state,iAgent) = SUM(MAXLOC(Q(state,:,iAgent)))
                    !
                END IF
                !
            END DO
            !
            ! Measuring performance
            !
            IF ((PerfMeasPeriodTime .LE. 0) .AND. ((iIters-1)/itersPerYear .GE. ABS(PerfMeasPeriodTime))) THEN
                !
                ! PerfMeasPeriodTime < 0: at convergence after mandatory training of PerfMeasPeriodTime years
                !
                IF (ALL(strategyPrime(state,:) .EQ. strategy(state,:))) THEN
                    !
                    iItersInStrategy = iItersInStrategy+1
                    !
                ELSE
                    !
                    iItersInStrategy = 1
                    !
                END IF
                !
            ELSE 
                !
                ! PerfMeasPeriodTime > 1: Starts right after the completion of year PerfMeasPeriodTime - 1
                !
                IF ((PerfMeasPeriodTime .GE. 1) .AND. ((iIters-1)/itersPerYear .GE. (PerfMeasPeriodTime-1))) THEN
                    !
                    iItersInStrategy = iItersInStrategy+1
                    !
                END IF
                !
            END IF
            !
            ! Check for convergence in strategy
            !
            IF (iItersInStrategy .EQ. itersInPerfMeasPeriod) EXIT
            !
            ! Check for too many iterations
            !
            IF (iIters .GT. maxIters) THEN
                !
                EXIT
                !
            END IF
            !
            ! If no convergence yet, update and iterate
            !
            strategy(state,:) = strategyPrime(state,:)
            state = statePrime
            !
            ! End of loop over iterations
            !
        END DO          
        !
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! End of first learning phase, random initialization of Q matrix
        ! Start of second learning phase, expert initialization of Q matrix
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !
        ! Initializing random number generator
        !
        idum = -iGames
        idum2 = 123456789
        iv = 0
        iy = 0
        !
        ! Initializing Q matrices
        !
        !$omp critical
        CALL initQMatricesRestart(iGames,idumQ,ivQ,iyQ,idum2Q,PI,delta,Q,maxValQ,strategyPrime)
        !$omp end critical
        strategy = strategyPrime
        !
        ! Randomly initializing prices and state
        !
        CALL initState(uIniPrice(:,:,iGames),p,statePrime,actionPrime)
        state = statePrime
        !
        ! Loop
        !
        iIters = 0
        iItersInStrategy = 0
        pricesGridsPrime = 0.d0
        convergedGame = 1
        DO 
            !
            ! Iterations counter
            !
            iIters = iIters+1
            !
            ! Generating exploration random numbers
            !            
            CALL generateUExploration(uExploration,idum,iv,iy,idum2)  
            !
            ! Compute pPrime by balancing exploration vs. exploitation
            !
            CALL computePPrimeRestart(ExplorationParameters,uExploration,strategyPrime,state,iIters,pPrime,Q)
            !
            ! Defining the new state
            !
            DO iAgent = 1, numAgents
                !
                pricesGridsPrime(iAgent) = PricesGrids(pPrime(iAgent),iAgent)
                !
            END DO
            IF (DepthState .GT. 1) p(2:DepthState,:) = p(1:DepthState-1,:)
            p(1,:) = pPrime
            statePrime = computeStateNumber(p)
            actionPrime = computeActionNumber(pPrime)
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
!@SP                    maxValQ(state,iAgent) = MAXVAL(Q(state,:,iAgent))
!@SP                    strategyPrime(state,iAgent) = SUM(MAXLOC(Q(state,:,iAgent)))
                    !
                END IF
                !
            END DO
            !
            ! Measuring performance
            !
            IF ((PerfMeasPeriodTime .LE. 0) .AND. ((iIters-1)/itersPerYear .GE. ABS(PerfMeasPeriodTime))) THEN
                !
                ! PerfMeasPeriodTime < 0: at convergence after mandatory training of PerfMeasPeriodTime years
                !
                IF (ALL(strategyPrime(state,:) .EQ. strategy(state,:))) THEN
                    !
                    iItersInStrategy = iItersInStrategy+1
                    !
                ELSE
                    !
                    iItersInStrategy = 1
                    !
                END IF
                !
            ELSE 
                !
                ! PerfMeasPeriodTime > 1: Starts right after the completion of year PerfMeasPeriodTime - 1
                !
                IF ((PerfMeasPeriodTime .GE. 1) .AND. ((iIters-1)/itersPerYear .GE. (PerfMeasPeriodTime-1))) THEN
                    !
                    iItersInStrategy = iItersInStrategy+1
                    !
                END IF
                !
            END IF
            !
            ! Check for convergence in strategy
            !
            IF (iItersInStrategy .EQ. itersInPerfMeasPeriod) EXIT
            !
            ! Check for too many iterations
            !
            IF (iIters .GT. maxIters) THEN
                !
                convergedGame = 0
                EXIT
                !
            END IF
            !
            ! If no convergence yet, update and iterate
            !
            strategy(state,:) = strategyPrime(state,:)
            state = statePrime
            !
            ! End of loop over iterations
            !
        END DO          
        !
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! End of second learning phase, expert initialization of Q matrix
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !
        ! Write Q matrices to file
        !
        !$omp critical
        IF (printQ .EQ. 1) THEN
            !
            ! Open Q matrices output file
            !
            WRITE(iGamesChar,'(I0.5)') iGames
            QFileName = 'Q_' // codModelChar // '_' // iGamesChar // '.txt'
            !
            ! Write on Q matrices to file
            !
            OPEN(UNIT = iGames,FILE = QFileName,RECL = 10000)
            DO iAgent = 1, numAgents
                !
                DO iState = 1, numStates
                    !
                    WRITE(iGames,*) Q(iState,:,iAgent)
                    !
                END DO
                !
            END DO
            CLOSE(UNIT = iGames)
            !
        END IF
        !$omp end critical
        !
        ! Record results at convergence
        !
        converged(iGames) = convergedGame
        indexLastState(:,iGames) = convertNumberBase(state-1,numPrices,LengthStates)
        timeToConvergence(iGames) = DBLE(iIters-itersInPerfMeasPeriod)/itersPerYear
        !
        ! End of loop over games
        !
    END DO
    !$omp end parallel do
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE computeModelRestart
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE computePPrimeRestart ( ExplorationParameters, uExploration, strategyPrime, state, iIters, &
        pPrime, Q )
    !
    ! Computes pPrime by balancing exploration vs. exploitation
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: ExplorationParameters(numExplorationParameters)
    REAL(8), INTENT(IN) :: uExploration(2,numAgents)
    INTEGER, INTENT(IN) :: strategyPrime(numStates,numAgents)
    INTEGER, INTENT(IN) :: state, iIters
    INTEGER, INTENT(OUT) :: pPrime(numAgents)
    REAL(8), INTENT(IN) :: Q(numStates,numPrices,numAgents)
    !
    ! Declaring local variables
    !
    INTEGER :: iAgent, iPrice
    REAL(8) :: eps(numAgents), u(2)
    !
    ! Beginning execution
    !
    ! 1. Greedy with probability 1-epsilon, with constant epsilon
    !
    IF (typeExplorationMechanism .EQ. 1) THEN
        !
        DO iAgent = 1, numAgents-computeRestart
            !
            pPrime(iAgent) = strategyPrime(state,iAgent)        ! "Expert" agents do not experiment
            !
        END DO
        DO iAgent = numAgents-computeRestart+1, numAgents
            !
            eps(iAgent) = ExplorationParameters(iAgent)
            u = uExploration(:,iAgent)                          ! "Non-Expert" agents can experiment
            IF (u(1) .LE. eps(iAgent)) THEN
                !
                pPrime(iAgent) = 1+INT(numPrices*u(2))
                !
            ELSE
                !
                pPrime(iAgent) = strategyPrime(state,iAgent)
                !
            END IF
            !
        END DO
        !
    END IF
    !
    ! 2. Greedy with probability 1-epsilon, with exponentially decreasing epsilon
    !
    IF (typeExplorationMechanism .GE. 2) THEN
        !
        DO iAgent = 1, numAgents-computeRestart
            !
            pPrime(iAgent) = strategyPrime(state,iAgent)        ! "Expert" agents do not experiment
            !
        END DO
        DO iAgent = numAgents-computeRestart+1, numAgents
            !
            IF (ExplorationParameters(iAgent) .LT. 0.d0) THEN
                !
                eps(iAgent) = 0.d0
                !
            ELSE
                !
                eps(iAgent) = EXP(-ExplorationParameters(iAgent)*DBLE(iIters-1)/DBLE(itersPerYear))
                !
            END IF
            u = uExploration(:,iAgent)                          ! "Non-Expert" agents can experiment
            IF (u(1) .LE. eps(iAgent)) THEN
                !
                pPrime(iAgent) = 1+INT(numPrices*u(2))
                !
            ELSE
                !
                pPrime(iAgent) = strategyPrime(state,iAgent)
                !
            END IF
            !
        END DO
        !
    END IF
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE computePPrimeRestart
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE initQMatricesRestart ( iGames, idumQ, ivQ, iyQ, idum2Q, PI, delta, Q, maxValQ, maxLocQ )
    !
    ! Initializing Q matrices
    ! - at the Q matrix at convergence for the first numAgents-1 agents
    ! - randomly for the last agent
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: iGames
    INTEGER, INTENT(INOUT) :: idumQ, ivQ(32), iyQ, idum2Q
    REAL(8), DIMENSION(numActions,numAgents), INTENT(IN) :: PI
    REAL(8), DIMENSION(numAgents), INTENT(IN) :: delta
    REAL(8), DIMENSION(numStates,numPrices,numAgents), INTENT(OUT) :: Q
    INTEGER, DIMENSION(numStates,numAgents), INTENT(OUT) :: maxLocQ
    REAL(8), DIMENSION(numStates,numAgents), INTENT(OUT) :: maxValQ
    !
    ! Declaring local variables
    !
    INTEGER :: iAgent, iPrice, iState, i, h, status
    INTEGER :: tied(numPrices)
    REAL(8) :: den, u
    CHARACTER(len = 225) :: QFileName
    CHARACTER(len = 5) :: iChar
    CHARACTER(len = 5) :: codModelChar
    CHARACTER(len = 200) :: QFileFolderNameAgent
    !
    ! Beginning execution
    !
    ! Initializing the Q matrix, only for the last agents
    !
    DO iAgent = numAgents-computeRestart+1, numAgents
        !
        IF (typeQInitialization(iAgent) .EQ. 'O') THEN
            !
            ! Randomizing over the opponents decisions
            !
            DO iPrice = 1, numPrices
                !
                den = COUNT(indexActions(:,iAgent) .EQ. iPrice)*(1.d0-delta(iAgent))
                Q(:,iPrice,iAgent) = SUM(PI(:,iAgent),MASK = indexActions(:,iAgent) .EQ. iPrice)/den
                !
            END DO
            !
        ELSE IF (typeQInitialization(iAgent) .EQ. 'T') THEN
            !
            ! Start from a randomly drawn Q matrix at convergence 
            ! on model "QMatrixInitializationT(iAgent)"
            !
            WRITE(codModelChar,'(I0.5)') QMatrixInitializationT(iAgent)
            i = 1+INT(DBLE(numGames)*ran2(idumQ,ivQ,iyQ,idum2Q))
            WRITE(iChar,'(I0.5)') i
            QFileName = 'Q_' // codModelChar // '_' // iChar // '.txt'
            QFileFolderNameAgent = QFileFolderName(iAgent)
            QFileName = TRIM(QFileFolderNameAgent) // TRIM(QFileName)
            !
            ! Read Q matrices from file
            !
            OPEN(UNIT = iGames,FILE = QFileName,READONLY,RECL = 10000,IOSTAT = status)
            IF (iAgent .GT. 1) READ(iGames,100)
100         FORMAT(<(iAgent-1)*numStates-1>(/))
            DO iState = 1, numStates
                !
                READ(iGames,*) Q(iState,:,iAgent)
                !
            END DO
            CLOSE(UNIT = iGames)
            !
        ELSE IF (typeQInitialization(iAgent) .EQ. 'R') THEN
            !
            ! Randomly initialized Q matrix using a uniform distribution between 
            ! QMatrixInitializationR(2,iAgent) and QMatrixInitializationR(1,iAgent)
            !
            DO iState = 1, numStates
                !
                DO iPrice = 1, numPrices
                    !
                    Q(iState,iPrice,iAgent) = ran2(idumQ,ivQ,iyQ,idum2Q)
                    !
                END DO
                !
            END DO
            Q(:,:,iAgent) = QMatrixInitializationR(1,iAgent)+ &
                (QMatrixInitializationR(2,iAgent)-QMatrixInitializationR(1,iAgent))*Q(:,:,iAgent)
            !
        ELSE IF (typeQInitialization(iAgent) .EQ. 'U') THEN
            !
            ! Constant Q matrix with all elements set to QMatrixInitializationU(iAgent)
            !
            Q(:,:,iAgent) = QMatrixInitializationU(iAgent)
            !
        END IF
        !
    END DO
    !
    ! Find initial optimal strategy
    !
    DO iAgent = 1, numAgents
        !
        DO iState = 1, numStates
            !
            CALL MaxLocBreakTies(numPrices,Q(iState,:,iAgent),idumQ,ivQ,iyQ,idum2Q, &
                maxValQ(iState,iAgent),maxLocQ(iState,iAgent))
            !
        END DO
        !
    END DO
!@SP    maxValQ = MAXVAL(Q,DIM = 2)
!@SP    maxLocQ = MAXLOC(Q,DIM = 2)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE initQMatricesRestart
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE LearningSimulationRestart
