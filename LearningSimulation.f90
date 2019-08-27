MODULE LearningSimulation
!
USE globals
USE QL_routines
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
    SUBROUTINE computeModel ( iModel, codModel, alpha, ExplorationParameters, delta, &
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
    INTEGER :: idumIP, ivIP(32), iyIP, idum2IP, idum, iv(32), iy, idum2, idumQ, ivQ(32), iyQ, idum2Q
    REAL(8), DIMENSION(numStates,numPrices,numAgents) :: Q
    REAL(8) :: uIniPrice(DepthState,numAgents,numGames), uExploration(2,numAgents)
    REAL(8) :: u(2), eps(numAgents)
    REAL(8) :: newq, oldq, profitgain
    INTEGER :: iIters, i, j, h, iGame, iItersInStrategy, convergedGame
    INTEGER :: state, statePrime, actionPrime
    INTEGER, DIMENSION(numStates,numAgents) :: strategy, strategyPrime
    INTEGER :: pPrime(numAgents), p(DepthState,numAgents)
    INTEGER :: iAgent, iState, iPrice, jAgent
    INTEGER :: minIndexStrategies, maxIndexStrategies
    CHARACTER(len = 25) :: QFileName
    CHARACTER(len = 5) :: iGamesChar
    CHARACTER(len = 5) :: codModelChar
    REAL(8), ALLOCATABLE :: printPMat(:,:), printPMatQ(:,:)
    CHARACTER(len = 200) :: PTrajectoryFileName
    !
    ! Beginning execution
    !
    ! Initializing various quantities
    !
    converged = 0    
    indexStrategies = 0
    indexLastState = 0
    timeToConvergence = 0.d0
    WRITE(codModelChar,'(I0.5)') codModel
    !
    ! Allocate variables
    !
    IF (printP .GT. 0) THEN
        !
        ALLOCATE(printPMat(0:itersPerYear,2*numAgents))    ! [Iterations, (Prices and Profits)]
        ALLOCATE(printPMatQ(0:itersPerYear,2*numAgents))  
        printPMat = 0.d0
        printPMatQ = 0.d0
        !
    END IF
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
    !$omp   strategyPrime,pPrime,p,statePrime,actionPrime,iIters,iItersInStrategy,convergedGame, &
    !$omp   state,strategy,eps,uExploration,u,oldq,newq,iAgent,iState,iPrice,jAgent, &
    !$omp   QFileName,iGamesChar) &
    !$omp firstprivate(numGames,PI,delta,uIniPrice,ExplorationParameters,itersPerYear,alpha, &
    !$omp   itersInPerfMeasPeriod,maxIters,printQ,printP,profitgain,codModelChar) &
    !$omp reduction(+ : printPMat, printPMatQ)
    DO iGame = 1, numGames
        !
        PRINT*, 'Game = ', iGame, ' started'
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
        iIters = 0
        iItersInStrategy = 0
        convergedGame = 1
        IF (printP .GT. 0) THEN
            !
            actionPrime = computeActionNumber(p(1,:))
            DO iAgent = 1, numAgents
                !
                printPMat(iIters,iAgent) = printPMat(iIters,iAgent)+PricesGrids(p(1,iAgent),iAgent)
                profitgain = (PI(actionPrime,iAgent)-NashProfits(iAgent))/(CoopProfits(iAgent)-NashProfits(iAgent))
                printPMat(iIters,numAgents+iAgent) = printPMat(iIters,numAgents+iAgent)+profitgain
                printPMatQ(iIters,iAgent) = printPMatQ(iIters,iAgent)+PricesGrids(p(1,iAgent),iAgent)**2
                printPMatQ(iIters,numAgents+iAgent) = printPMatQ(iIters,numAgents+iAgent)+profitgain**2
                !
            END DO
            !
        END IF
        !
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
            ! Store PTrajectories
            !
            IF ((printP .GT. 0) .AND. (iIters .LE. itersPerYear*printP) .AND. (MOD(iIters,printP) .EQ. 0)) THEN
                !
                DO iAgent = 1, numAgents
                    !
                    printPMat(iIters/printP,iAgent) = printPMat(iIters/printP,iAgent)+ &    
                        PricesGrids(pPrime(iAgent),iAgent)
                    profitgain = (PI(actionPrime,iAgent)-NashProfits(iAgent))/(CoopProfits(iAgent)-NashProfits(iAgent))
                    printPMat(iIters/printP,numAgents+iAgent) = printPMat(iIters/printP,numAgents+iAgent)+profitgain
                    printPMatQ(iIters/printP,iAgent) = printPMatQ(iIters/printP,iAgent)+ &
                        PricesGrids(p(1,iAgent),iAgent)**2
                    printPMatQ(iIters/printP,numAgents+iAgent) = printPMatQ(iIters/printP,numAgents+iAgent)+profitgain**2
                    !
                END DO
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
            IF (iItersInStrategy .EQ. itersInPerfMeasPeriod) THEN
                !
                IF (PerfMeasPeriodTime .GE. 1) convergedGame = 0
                EXIT
                !
            END IF
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
        ! Write Q matrices to file
        !
        !$omp critical
        IF (printQ .EQ. 1) THEN
            !
            ! Open Q matrices output file
            !
            WRITE(iGamesChar,'(I0.5)') iGame
            QFileName = 'Q_' // codModelChar // '_' // iGamesChar // '.txt'
            !
            ! Write on Q matrices to file
            !
            OPEN(UNIT = iGame,FILE = QFileName,RECL = 10000)
            DO iAgent = 1, numAgents
                !
                DO iState = 1, numStates
                    !
                    WRITE(iGame,*) Q(iState,:,iAgent)
                    !
                END DO
                !
            END DO
            CLOSE(UNIT = iGame)
            !
        END IF
        !$omp end critical
        !
        ! Record results at convergence
        !
        converged(iGame) = convergedGame
        indexStrategies(:,iGame) = computeStrategyNumber(strategy)
        indexLastState(:,iGame) = convertNumberBase(state-1,numPrices,LengthStates)
        timeToConvergence(iGame) = DBLE(iIters-itersInPerfMeasPeriod)/itersPerYear
        !
        IF (convergedGame .EQ. 1) PRINT*, 'Game = ', iGame, ' converged'
        IF (convergedGame .EQ. 0) PRINT*, 'Game = ', iGame, ' did not converge'
        !
        ! End of loop over games
        !
    END DO
    !$omp end parallel do
    !
    ! Print P trajectories, if required
    !
    IF (printP .GT. 0) THEN
        !
        ! Statistics
        !
        printPMat = printPMat/DBLE(numGames)
        printPMatQ = printPMatQ/DBLE(numGames)
        printPmatQ = SQRT(ABS(printPMatQ-printPMat**2))
        !
        ! File name
        !
        WRITE(codModelChar,'(I0.5)') iModel
        PTrajectoryFileName = 'PTrajectories_' // codModelChar // '.txt'
        OPEN(UNIT = 991,FILE = PTrajectoryFileName,STATUS = "REPLACE")
        WRITE(991,990) (iAgent, iAgent = 1, numAgents), (iAgent, iAgent = 1, numAgents), & 
                       (iAgent, iAgent = 1, numAgents), (iAgent, iAgent = 1, numAgents)
990     FORMAT('      iter ', &
            <numAgents>('      Price', I1, ' '), <numAgents>('    sePrice', I1, ' '), &
            <numAgents>('     Profit', I1, ' '), <numAgents>('   seProfit', I1))
        DO i = 0, itersPerYear
            !
            WRITE(991,991) i*printP, printPMat(i,:numAgents), printPMatQ(i,:numAgents), &
                printPMat(i,numAgents+1:), printPMatQ(i,numAgents+1:)
            !
        END DO
991     FORMAT(I10, 1X, <4*numAgents>(F12.5,1X))
        CLOSE(UNIT = 991)
        !
    END IF
    !
    ! Print indexStrategies and indexLastState to file
    !
    OPEN(UNIT = 996,FILE = FileNameIndexStrategies,STATUS = "REPLACE")
    WRITE(996,992) converged
992 FORMAT(<numGames>(I1, 1X))
    WRITE(996,993) timeToConvergence
993 FORMAT(<numGames>(F9.2, 1X))        
    DO i = 1, lengthStrategies
        !
        WRITE(996,996) indexStrategies(i,:)
996     FORMAT(<numGames>(I<lengthFormatActionPrint>,1X))
        !
    END DO
    CLOSE(UNIT = 996)
    !    
    OPEN(UNIT = 999,FILE = FileNameIndexLastState,STATUS = "REPLACE")
    DO iGame = 1, numGames
        !
        WRITE(999,999) indexLastState(:,iGame), converged(iGame), timeToConvergence(iGame)
999     FORMAT(<LengthStates>(I<lengthFormatActionPrint>,1X), I3, 1X, ES12.5, 1X)
        !
    END DO
    CLOSE(UNIT = 999)
    !
    ! Deallocate variables
    !
    IF (printP .GT. 0) DEALLOCATE(printPMat)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE computeModel
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE computePPrime ( ExplorationParameters, uExploration, strategyPrime, state, iIters, pPrime, Q )
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
        DO iAgent = 1, numAgents
            !
            eps(iAgent) = ExplorationParameters(iAgent)
            u = uExploration(:,iAgent)
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
    ! 2. & 3. Greedy with probability 1-epsilon, with exponentially decreasing epsilon
    !
    IF (typeExplorationMechanism .GE. 2) THEN
        !
        DO iAgent = 1, numAgents
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
            u = uExploration(:,iAgent)
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
!@SP Sanity check    
!pPrime(1) = 14
!@SP Sanity check    
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE computePPrime
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE LearningSimulation
