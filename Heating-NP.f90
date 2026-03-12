!     Heating of the material in a less violent manner then 
!     just by velocities distribution (as in Vel.f [KFLAG=2])
!     Leonid Zhigilei, 2000
      SUBROUTINE HEATING_NP()
        INCLUDE 'common.h'
        INCLUDE 'commonTTM.h'
        INCLUDE 'mpif.h'

        tau_laser = DURAT ! ps
        t_zero = 2.5d0*tau_laser
        IF(TIME.GE.2*t_zero) GOTO 456
        
        E_incedent = ENRIN/(XL*YL)*XJOUTOEV !J/m2
        omegalas = 4.0d0*LOG(2.00)
        Refl = 0.96d0 ! Estimated from nTTM for the current pulse duration and the energy
        
        if(mypid.eq.0) print *,"The total absorbed fluence is:", &
             (1-Refl)*ENRIN/(XL*YL)*1.0d+19," mJ/cm2"
        
        Enmult = (1-Refl)*E_incedent/tau_laser*SQRT(omegalas/PI) &
             *EXP(-omegalas*(time - t_zero)**2/tau_laser**2)*delta
        
        Qnp = 0.0d0
        
        DO I = 1,NN1
           IF(KTYPE(I).EQ.1) THEN
              Qnp = Qnp + 0.5d0*XMASS(KTYPE(I))*(Q1D(1,I)*Q1D(1,I) + &
                   Q1D(2,I)*Q1D(2,I) + Q1D(3,I)*Q1D(3,I))/delta/delta
           ENDIF
        ENDDO

        CALL MPI_REDUCE(Qnp,QnpG,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
             MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(QnpG,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        
        DO I=1,NN1
           IF(KTYPE(I).EQ.1) THEN
              SC = SQRT(1.0d0 + Enmult*XL*YL/(QnpG*ENUNIT))
              Q1D(1,I)=Q1D(1,I)*SC
              Q1D(2,I)=Q1D(2,I)*SC
              Q1D(3,I)=Q1D(3,I)*SC
           ENDIF
        ENDDO

        F_total = F_total + Qnp*(SC*SC-1.0d0)*ENUNIT*EVTOJOU
        
456     CONTINUE
        
        RETURN
      END SUBROUTINE HEATING_NP
