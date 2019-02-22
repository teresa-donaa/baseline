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
    SUBROUTINE computeModelRestart ( iModel, alpha, ExplorationParameters, delta, &
        converged, indexLastState, timeToConvergence )
    !
    ! Computes statistics for one model
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: iModel
    REAL(8), DIMENSION(numAgents), INTENT(IN) :: alpha, delta
    REAL(8), DIMENSION(numExplorationParameters) :: ExplorationParameters
    INTEGER, DIMENSION(numGames), INTENT(OUT) :: converged
    INTEGER, INTENT(OUT) :: indexLastState(lengthStates,numGames)
    REAL(8), DIMENSION(numGames), INTENT(OUT) :: timeToConvergence
    !
    ! Declaring local variable
    !
    REAL(8), DIMENSION(numAgents) :: pricesGridsPrime
    INTEGER :: idumIP, ivIP(32), iyIP, idum2IP, idum, iv(32), iy, idum2
    REAL(8), DIMENSION(numStates,numPrices,numAgents) :: Q
    REAL(8) :: uIniPrice(depthState,numAgents,numGames), uExploration(2,numAgents), u(2), eps(numAgents)
    REAL(8) :: newq, oldq
    INTEGER :: iIters, i, j, iGames, iItersInStrategy, convergedGame
    INTEGER :: state, statePrime, actionPrime
    INTEGER, DIMENSION(numStates,numAgents) :: strategy, strategyPrime
    INTEGER :: pPrime(numAgents), p(depthState,numAgents)
    INTEGER :: iAgent, iState, iPrice, jAgent
    INTEGER :: minIndexStrategies, maxIndexStrategies
    CHARACTER(len = 25) :: printQFileName
    CHARACTER(len = 3) :: iGamesChar
    REAL(8), ALLOCATABLE :: Qtrajectories(:,:)
    REAL(8), ALLOCATABLE :: Ptrajectories(:,:)
    REAL(8) :: RDeviation
    !
    ! Beginning execution
    !
    ! Initializing various quantities
    !
    converged = 0    
    indexLastState = 0
    timeToConvergence = 0.d0
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
    CALL generate_uIniPrice(numGames,uIniPrice,idumIP,ivIP,iyIP,idum2IP)  
    !
    ! Starting loop over games
    !
    !$omp parallel do &
    !$omp private(idum,iv,iy,idum2,Q,maxValQ, &
    !$omp   strategyPrime,pPrime,p,statePrime,actionPrime,iIters,iItersInStrategy, &
    !$omp   convergedGame,pricesGridsPrime, &
    !$omp   state,strategy,eps,uExploration,u,oldq,newq,iAgent,iState,iPrice,jAgent, &
    !$omp   printQFileName,iGamesChar) &
    !$omp firstprivate(numGames,PI,delta,uIniPrice,ExplorationParameters,itersPerYear,alpha, &
    !$omp   itersInPerfMeasPeriod,maxIters,printQ,Qtrajectories,Ptrajectories)
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
        CALL initQMatrices(PI,delta,Q,maxValQ,strategyPrime)
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
            IF (depthState .GT. 1) p(2:depthState,:) = p(1:depthState-1,:)
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
                    maxValQ(state,iAgent) = MAXVAL(Q(state,:,iAgent))
                    strategyPrime(state,iAgent) = SUM(MAXLOC(Q(state,:,iAgent)))
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
        CALL initQMatricesRestart(PI,delta,Q,maxValQ,strategyPrime)
        strategy = strategyPrime
        !
        ! Randomly initializing prices and state
        !
        CALL initState(uIniPrice(:,:,iGames),p,statePrime,actionPrime)
        state = statePrime
        !
        ! Allocate Qtrajectories and Ptrajectories matrices
        !
        IF (iGames .LE. printQ) THEN
            !
            ALLOCATE(Qtrajectories(itersInPerfMeasPeriod/QTrajectoriesPeriodPrint,lengthStrategies*numPrices), &
                     Ptrajectories(itersInPerfMeasPeriod,numAgents))
            !
        END IF
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
            IF (depthState .GT. 1) p(2:depthState,:) = p(1:depthState-1,:)
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
                    maxValQ(state,iAgent) = MAXVAL(Q(state,:,iAgent))
                    strategyPrime(state,iAgent) = SUM(MAXLOC(Q(state,:,iAgent)))
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
            ! Add a line to Qtrajectories and Ptrajectories output file
            !
            IF (iGames .LE. printQ) THEN
                !
                ! Qtrajectories
                !
                IF (MOD(iItersInStrategy,QTrajectoriesPeriodPrint) .EQ. 0) THEN
                    !
                    DO iAgent = 1, numAgents
                        !
                        Qtrajectories(iItersInStrategy/QTrajectoriesPeriodPrint,(iAgent-1)*numStates*numPrices+1:iAgent*numStates*numPrices) = &
                            RESHAPE(TRANSPOSE(Q(:,:,iAgent)),(/ numStates*numPrices /))
                        !
                    END DO
                    !
                END IF
                !
                ! Ptrajectories
                !
                Ptrajectories(iItersInStrategy,:) = pricesGridsPrime
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
        ! Write Qtrajectories and Ptrajectories output file
        !
        IF (iGames .LE. printQ) THEN
            !
            ! Opening Qtrajectories output file
            !
            WRITE(iGamesChar,'(I0.3)') iGames
            printQFileName = 'Qtrajectories_' // iGamesChar // '.txt'
            OPEN(UNIT = iGames,FILE = printQFileName)
            !
            ! Writing on Qtrajectories output file
            !
            WRITE(iGames,11) (((iAgent, iState, iPrice, iPrice = 1, numPrices), &
                iState = 1, numStates), iAgent = 1, numAgents)
11          FORMAT('     Iter ', &
                <6750>('  Q', I1, '(', I3, ',', I2, ') '))
!@SP USE NUMERIC VALUE <lengthStrategies*numPrices>
            WRITE(iGames,12) (iIters*QTrajectoriesPeriodPrint, Qtrajectories(iIters,:), iIters = 1, itersInPerfMeasPeriod/QTrajectoriesPeriodPrint)
12          FORMAT(<219>(I9, 1X, <6750>(F12.5, 1X), /))
!@SP USE NUMERIC VALUE <itersInPerfMeasPeriod/QTrajectoriesPeriodPrint>
            !
            ! Closing Qtrajectories output file
            !
            CLOSE(UNIT = iGames)
            !
            ! Opening output file
            !
            printQFileName = 'Ptrajectories_' // iGamesChar // '.txt'
            OPEN(UNIT = iGames,FILE = printQFileName)
            !
            ! Writing on output file
            !
            WRITE(iGames,13) (iAgent, iAgent = 1, numAgents)
13          FORMAT('     Iter    pPrime(', I1, ')    pPrime(', I1, ')')
            WRITE(iGames,14) (iIters, Ptrajectories(iIters,:), iIters = 1, itersInPerfMeasPeriod)
14          FORMAT(<26280>(I9, 1X, <2>(F12.5, 1X), /))
!@SP USE NUMERIC VALUE <itersInPerfMeasPeriod> AND <numAgents>
            !
            ! Closing output file
            !
            CLOSE(UNIT = iGames)
            !
            ! Deallocating Qtrajectories and Ptrajectories matrices
            !
            DEALLOCATE(Qtrajectories,Ptrajectories)
            !
        END IF
        !
        converged(iGames) = convergedGame
        indexLastState(:,iGames) = convertNumberBase(state-1,numPrices,lengthStates)
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
        eps = ExplorationParameters
        !
        DO iAgent = 1, numAgents-computeRestart
            !
            pPrime(iAgent) = strategyPrime(state,iAgent)        ! "Expert" agents do not experiment
            !
        END DO
        DO iAgent = numAgents-computeRestart+1, numAgents
            !
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
    IF (typeExplorationMechanism .EQ. 2) THEN
        !
        eps = EXP(-ExplorationParameters*DBLE(iIters-1)/DBLE(itersPerYear))
        !
        DO iAgent = 1, numAgents-computeRestart
            !
            pPrime(iAgent) = strategyPrime(state,iAgent)        ! "Expert" agents do not experiment
            !
        END DO
        DO iAgent = numAgents-computeRestart+1, numAgents
            !
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
    SUBROUTINE initQMatricesRestart ( PI, delta, Q, maxValQ, maxLocQ )
    !
    ! Initializing Q matrices
    ! - at the Q matrix at convergence for the first numAgents-1 agents
    ! - randomly for the last agent
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), DIMENSION(numActions,numAgents), INTENT(IN) :: PI
    REAL(8), DIMENSION(numAgents), INTENT(IN) :: delta
    REAL(8), DIMENSION(numStates,numPrices,numAgents), INTENT(INOUT) :: Q
    INTEGER, DIMENSION(numStates,numAgents), INTENT(OUT) :: maxLocQ
    REAL(8), DIMENSION(numStates,numAgents), INTENT(OUT) :: maxValQ
    !
    ! Declaring local variables
    !
    INTEGER :: iAgent, iPrice, p(depthState,numAgents), iAction, iNash, iDevNash
    INTEGER :: devPrices(numAgents)
    REAL(8) :: den
    !
    ! Beginning execution
    !
    ! Randomizing over the opponents decisions, only for the last agent
    !
    DO iAgent = numAgents-computeRestart+1, numAgents
        !
        DO iPrice = 1, numPrices
            !
            den = COUNT(indexActions(:,iAgent) .EQ. iPrice)*(1.d0-delta(iAgent))
            Q(:,iPrice,iAgent) = SUM(PI(:,iAgent),MASK = indexActions(:,iAgent) .EQ. iPrice)/den
            !
        END DO
        !
    END DO
    !
    maxValQ = MAXVAL(Q,DIM = 2)
    maxLocQ = MAXLOC(Q,DIM = 2)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE initQMatricesRestart
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE LearningSimulationRestart
