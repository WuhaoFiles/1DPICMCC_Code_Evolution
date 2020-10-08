Module DiagnosticsEEPF
    use ModuleParticleBoundary
    Use ModuleGrid
    Use ModuleParticleBundle
    Implicit none

   Integer(4),Parameter,Private :: NeMax=1000
   Type ParticleEDF
         Integer(4) :: Ne=NeMax
         Real(8) :: EnergyInterval,Mass
         Real(8) :: EDF(NeMax),EDFNormalized(NeMax)
     EndType ParticleEDF
    contains
        subroutine DiagParticleEDFOne(GD,PB,PBDO,Mode)
            Implicit none
            Class(*),intent(inout)  :: GD
            Type(ParticleBundle),intent(in) :: PB
			   Type(ParticleBoundaryOne),intent(in) :: PBDO
            Integer(4),intent(in) ::  Mode
            Type(ParticleEDF) :: PEDF
				!Type(ParticleBundle),intent(in) :: PB(2)
            Integer(4) :: Shift,i

            If (PB%SO%SpecyIndex==0) Then
                PEDF%EnergyInterval = 0.1d0
            Else
                PEDF%EnergyInterval = 0.02d0
            End If
            PEDF%Mass=PB%Mass

            Select Type (GD)
                  Type is (Grid1D(*,*))
                     Select Case (Mode)
                        Case(-1)
                            !Call GridInitialization(GD,PB%Period,PEDF%EnergyInterval,PB%Dt)
                        case(0)
                            
                            PEDF%Ne=GD%Nx
                            Shift=1
									 
                            Call WeightingParticleEDF(PB,PEDF)
                            Call GD%Update(GD%Nx,PEDF%EDF,Shift)
                            Call GD%Update(GD%Nx,PEDF%EDFNormalized,Shift)
									 !AEDF
									 Call WeightingParticleEDF2(PBDO%PBLower,PEDF,1,PBDO%CountMinOne)
                            Call GD%Update(GD%Nx,PEDF%EDF,Shift)
                            Call GD%Update(GD%Nx,PEDF%EDFNormalized,Shift)
									 !Shift=Shift-2
									 Call WeightingParticleEDF2(PBDO%PBUpper,PEDF,1,PBDO%CountMaxOne)
                            Call GD%Update(GD%Nx,PEDF%EDF,Shift)
                            Call GD%Update(GD%Nx,PEDF%EDFNormalized,Shift)
									 !see EEDF
									 Call WeightingParticleEDF2(PBDO%PBLower,PEDF,PBDO%CountMinOne+1,PBDO%CountMinOne+PBDO%SEECountMinOne)
                            Call GD%Update(GD%Nx,PEDF%EDF,Shift)
                            Call GD%Update(GD%Nx,PEDF%EDFNormalized,Shift)
									 Shift=Shift-2
									 Call WeightingParticleEDF2(PBDO%PBUpper,PEDF,PBDO%CountMaxOne+1,PBDO%CountMaxOne+PBDO%SEECountMaxOne)
                            Call GD%Update(GD%Nx,PEDF%EDF,Shift)
                            Call GD%Update(GD%Nx,PEDF%EDFNormalized,Shift)
									 
									 
                            GD%Timer=GD%Timer+1
                         Case(1) 
                                 Call GD%Rescale
                                 Call GD%Dump(1)
                                 GD%Value=0.d0
                                 GD%Timer=0
                         Case(2)
                             Call GD%Dump(0)
                         case default
                         End Select
                    Type is (Grid2D(*,*,*))
                        Select Case (Mode)
                            Case(-1)
                                !Call GridInitialization(GD,PB%Period,PEDF%EnergyInterval,PB%Dt)
                            case(0)
                                Shift=1
                                Call WeightingParticleEDF(PB,PEDF)
                                Call GD%Update(GD%Nx,PEDF%EDF,Shift)
                                Call GD%Update(GD%Nx,PEDF%EDFNormalized,Shift)
                                GD%Timer=GD%Timer+1
                            Case(1) 
                                     Call GD%Rescale
                                     Call GD%Dump(1)
                                     GD%Value=0.d0
                                     GD%Timer=0
                             Case(2)
                                 Call GD%Dump(0)
                             case default
                         End Select     
                    End select
            return
             end subroutine DiagParticleEDFOne
             
        subroutine WeightingParticleEDF(PB,PEDF)
                implicit none
                Type(ParticleBundle),intent(in) :: PB
                Type(ParticleEDF),intent(inout) :: PEDF
                Real(8) :: Energy,Frac
                Integer(4) :: i,N
                
                PEDF%EDF=0.d0
                PEDF%EDFNormalized=0.d0
                
                Frac=1.d0/DBLE(PB%Npar)
                do i=1,PB%Npar
                   Energy=PB%PO(i)%Energy(PB%Mass,PB%VFactor)/JtoeV
                   N=Ceiling( Energy/PEDF%EnergyInterval)
                   If (N>=1.and.N<PEDF%Ne) Then
                           PEDF%EDF(N) = PEDF%EDF(N)+PB%Weight*PB%Dx
                           PEDF%EDFNormalized(N)=PEDF%EDFNormalized(N)+Frac
                   End IF
                end do
                return
		  end subroutine WeightingParticleEDF     
		  
		  subroutine WeightingParticleEDF2(PB,PEDF,NumLow,NumMax)
                implicit none
                Type(ParticleBundle),intent(in) :: PB
                Type(ParticleEDF),intent(inout) :: PEDF
					 integer,intent(in) :: NumLow
					 integer,intent(in) :: NumMax
                Real(8) :: Energy,Frac
                Integer(4) :: i,N
                
                PEDF%EDF=0.d0
                PEDF%EDFNormalized=0.d0
                Frac=1.d0/DBLE(PB%Npar)
                do i=NumLow,NumMax
                   Energy=PB%PO(i)%Energy(PB%Mass,PB%VFactor)/JtoeV
                   N=Ceiling( Energy/PEDF%EnergyInterval)
                   If (N>=1.and.N<PEDF%Ne) Then
                           PEDF%EDF(N) = PEDF%EDF(N)+PB%Weight*PB%Dx
                           PEDF%EDFNormalized(N)=PEDF%EDFNormalized(N)+Frac
                   End IF
                end do
                return
    end subroutine WeightingParticleEDF2     
		  
End Module DiagnosticsEEPF