!     Two-body part of Stillinger-Weber potential for Silicon
!     Avinash Dongare and Leonid Zhigilei, 2003
!     [F. H. Stillinger and T. A. Weber, Phys. Rev. B 31, 5262-5271 (1985)]
      SUBROUTINE EF1_SW(KTEF,KE1,KE2)
        INCLUDE 'common.h'

!       KTEF - type of the potential
!       KE1, KE2 - elements from the list Element(N)
!       Units: XT in Angs, UT in eV, FT in eV/Angs, mass in amu
!       Parameters for Stillinger-Weber potential for Silicon
!       [F. H. Stillinger and T. A. Weber, Phys. Rev. B 31, 5262-5271 (1985)]
!       There are two subsequent papers that give Si-F interactions
!       This is only coded for Si, but we are prepared for Si-F interactions
!       Paramerets for Stillinger-Weber potential for Ge:
!       [K. Ding and H. Andersen. Phys. Rev. B 34, 6987 (1986)]
!       Paramerets for Stillinger-Weber potential for Si-Ge:
!       [M. Laradji, D.P. Landau and B. Dunweg. Phys. Rev. B 51 4894 (1995)]
!       Stillinger and Weber use reduced units.  We use real ones.
 
        CHARACTER*2, DIMENSION(2), PARAMETER:: Element = (/'Si','Ge'/)

        REAL*8, DIMENSION(2), PARAMETER:: mass=   & 
             (/ 28.0855D+00, 72.64D+00/)
        
!       Eps is not used anywhere...
        REAL*8, DIMENSION(3), PARAMETER:: TSW_eps=   &    
             (/ 2.167222D+00, 1.93D+00, 2.0427D+00 /)
!       SW length units in A
        REAL*8, DIMENSION(3), PARAMETER:: TSW_sig=   &  
             (/2.0951D+00,2.181D+00,2.13805D+00/)  

!       Parameters for 2-body term
!       q and p are integers but we define them 
!       as reals to be ready for Si-F potential
        REAL*8, DIMENSION(3), PARAMETER:: SW_P= (/4.0D+00,4.0D+00,4.0D+00/)
        REAL*8, DIMENSION(3), PARAMETER:: SW_Q= (/0.0D+00,0.0D+00,0.0D+00/) 
!       A and B in non-reduced units
        REAL*8, DIMENSION(3), PARAMETER:: SW_A=   &
             (/15.2779D+00,13.6056D+00,14.40013D+00/)
        REAL*8, DIMENSION(3), PARAMETER:: SW_B=   &  
             (/0.60222D+00,0.60222D+00,0.60222D+00/)

!       Parameters for 3-body term - Lambda and Gamma in non-reduced units
        REAL*8, DIMENSION(2), PARAMETER:: TSW_LAM=   &  
             (/21.0D+00,31.0D+00/)
        REAL*8, DIMENSION(3), PARAMETER:: TSW_GAM=   &  
             (/1.2D+00,1.2D+00,1.2D+00/)
        REAL*8, DIMENSION(3), PARAMETER:: TSW_RM=   &  
             (/3.7712D+00,3.9258D+00,3.8437D+00 /)

        real*8 :: rcc = 6.875d0   ! cutoff in angstroem if Si and H2O are in the couple
        
!       For Si S-W give SW_eps = 3.4723D-19 J = 2.167222 eV
!       Thijsse (following Balamane) proposes to multiply the energy scale 
!       by 1.0676 so that cohesive energy would be 4.63 eV instead of 4.34 eV
!       Print *, "A = ", 7.049556277D+00*TSW_eps(1)*TSW_sig(1)**SW_q(1)
!       Print *, "B = ", 0.6022245584D+00*TSW_sig(1)**(SW_p(1)-SW_q(1))
!       Print *, "lamb = ", 21.0D+00*TSW_eps(1)
!       Print *, "gam = ", 1.20D+00*TSW_sig(1)
        !       Print *, "rm = ", 1.80D+00*TSW_sig(1)

        IF(mypid.eq.0) Write (17,*) 'Creating tables for ', Element(KE1),'-', &
             Element(KE2),' potential' 

        IF(KE1.EQ.KE2) Then
           RM(KTEF) = TSW_RM(KE1)
           SW_lam(KTEF) = TSW_lam(KE1)
           SW_gam(KTEF) = TSW_gam(KE1)
           SW_sig(KTEF) = TSW_sig(KE1)
           SW_eps(KTEF) = TSW_eps(KE1)
           A = SW_A(KE1)
           P = SW_P(KE1)
           B = SW_B(KE1)
           Q = SW_Q(KE1)
        ELSE
           RM(KTEF) = TSW_RM(3)
           SW_gam(KTEF) =TSW_gam(3)
           SW_sig(KTEF) =TSW_sig(3)
           SW_eps(KTEF) =TSW_eps(3)
           A = SW_A(3)
           P = SW_P(3)
           B = SW_B(3)
           Q = SW_Q(3)
        ENDIF
  
!       IF(mypid.eq.0) OPEN (UNIT = 444,FILE='SW-'//Element(KE1)//'-'//Element(KE2)//'.pot')
        IF(KTEF.EQ.1) alat = 5.437325d0! for 273K , but for 0K = 5.43096d0!5.437325d0

        XMass(1) = MASS(KE1)
!        XMass(2) = MASS(KE2)

        SIG = SW_SIG(KTEF)
        RMSW = RM(KTEF)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        DX=RM(KTEF)/NT
!        IF(KEYBS.EQ.5.AND.NTYPE.EQ.2) DX=rcc/NT ! In order to mix SW and EAM water at the same range
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        IF(mypid.EQ.0) OPEN(UNIT = 444, FILE = "SW-potential.out")

        UT(KTEF,1:NT) = 0.0d0
        FT(KTEF,1:NT) = 0.0d0
        
        XT=0.0d0  
        DO I=1,NT
           XT=XT+DX  
           
           IF(XT.LE.RMSW) Then  
              ee = exp(SIG/(XT-RMSW))
              rp = (XT/SIG)**(-P)
              rq = (XT/SIG)**(-Q)
              rpp = (XT/SIG)**(-P-1)
              rqq = (XT/SIG)**(-Q-1)
              
              UT(KTEF,I) = A*(B*rp-rq)*ee
              FT(KTEF,I) = -A*(-B*P*rpp+Q*rqq)*ee/SIG+ &
                   UT(KTEF,I)*SIG/(XT-RMSW)**2

!W            Uncomment this line if want to plot the potential
             IF(mypid.eq.0) WRITE (444,*) XT,UT(KTEF,I),FT(KTEF,I)

!             Transfer to program units
              UT(KTEF,I) = UT(KTEF,I)/ENUNIT
              FT(KTEF,I) = FT(KTEF,I)/ENUNIT
           ELSE 
              UT(KTEF,I) =0.0
              FT(KTEF,I) =0.0  
           ENDIF
        END DO
        
!W      Uncomment this line if want to plot the potential
       IF(mypid.eq.0) CLOSE (UNIT = 444)

       RETURN
      END SUBROUTINE EF1_SW
