!     Making the Energy and Force tables
!     for Jonhson's new EAM potential
!     Dewey Murdick and Leonid Zhigilei, 2000, 2002
!     Number of pair potential defined here should be 
!     Npots=Ntype*(Ntype+1)/2,  Ntype is number of particle types

      SUBROUTINE EF1_EAM(KTEF,KE1,KE2)
        INCLUDE 'common.h'
!       KTEF - type of the potential
!       KE1, KE2 - elements from the list Element(N)
!       Units: XT in Angs, UT in eV, FT in eV/Angs, mass in amu

        INTEGER I
        real*8, parameter :: rcc = 6.875d0   ! cutoff in angstroem
	real*8 xx, xc, rhomax

        REAL *8 XT, DXRho, DRho, Rhoin, Rhoout
        REAL *8 re, fe, Rhoe, Alpha, Beta, A, B, Kappa, Lambda
        REAL *8 Fi0,Fi1,Fi2,Fi3,Fm0,Fm1,Fm2,Fm3,niu,Fn
        
        CHARACTER*2, DIMENSION(8), PARAMETER:: Element =   &
             (/'Cu','Au','Ni','Fe','Zr','Al','ST','HO'/)
        
        REAL*8, DIMENSION(8), PARAMETER:: Tmass=   &
             (/63.6d0, 196.966569d0,  58.7d0,  55.9d0, 91.2d0, &
             26.981539d0, 183.49d0, 18.02d0/)

        
        IF(KE1.EQ.KE2) Then   
!       going from the type of the potential to the type of atom
           IF(KTEF.EQ.1) KF=1  
           IF(KTEF.EQ.2) KF=2

           IF(KTEF.EQ.1) THEN
              OPEN(UNIT = 111,FILE = 'gold.data')
              READ(111,*)(XT,UT(KTEF,I),FT(KTEF,I),DFT(KF,I),DDFT(KF,I), &
                   Rho,EFT(KF,I),DEFT(KF,I),I=1,NT)
              CLOSE(UNIT = 111)
           ENDIF

           IF(KTEF.EQ.2) THEN
              OPEN(UNIT = 111,FILE = 'water.data')
              READ(111,*)(XT,UT(KTEF,I),FT(KTEF,I),DFT(KF,I),DDFT(KF,I), &
                   Rho,EFT(KF,I),DEFT(KF,I),I=1,NT)
              CLOSE(UNIT = 111)
           ENDIF
           
           XMass(KF)=Tmass(KE1)      ! Mass of the particle in pr.unit [aem]
           RM(KTEF) = rcc		       ! Cutoff dist. for 2-body potential [A]
           IF(KTEF.EQ.1) alat = 4.065d0

        ELSEIF(KTEF.EQ.3.AND.KE1.NE.KE2) THEN

           DX = RM(2)/NT
           DRho=RhoM(2)/NT
           DO I = 1,NT

              UT(KTEF,I) = 0.5d0*(UT(1,I) + UT(2,I))
              FT(KTEF,I) = 0.5d0*(FT(1,I) + FT(2,I))

              XT = DX*I
              Rho = Drho*I

           ENDDO
           RM(KTEF) = rcc!RM(1)
           GOTO 789
        ELSE
           If(mypid.eq.0)Print *,"This combination of EAM potentials is not implemented yet. Program will stop."
           STOP
        ENDIF
        
        RhoM(KF) = 8.d0

        DX = RM(KTEF)/NT
        DRho=RhoM(KF)/NT

!          Transfer to program units
!          Electron density and its derivative is in eV/A, do not convert to pr.un.
        DO I = 1,NT
           XT = DX*I
           Rho = Drho*I

           UT(KTEF,I) = UT(KTEF,I)/ENUNIT
           FT(KTEF,I) = FT(KTEF,I)/ENUNIT
           EFT(KF,I)  =  EFT(KF,I)/ENUNIT
           DEFT(KF,I) = DEFT(KF,I)/ENUNIT
        END DO

789     CONTINUE
	
        RETURN
      END SUBROUTINE EF1_EAM
