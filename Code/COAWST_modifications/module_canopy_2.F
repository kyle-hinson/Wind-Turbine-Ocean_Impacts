!WRF:MEDIATION_LAYER:PHYSICS

MODULE module_canopy_2

USE module_model_constants
USE module_big_step_utilities_em, only: grid_config_rec_type

CONTAINS

!=======================================================================
! Sets the area density for canopy_opt = 2

  SUBROUTINE set_area_density_2( area_density, area_density_opt,   & 
                                 canopy_height, canopy_top, z,     &
                                 ids, ide, jds, jde, kds, kde,     &
                                 ims, ime, jms, jme, kms, kme,     &
                                 its, ite, jts, jte, kts, kte      )

!-----------------------------------------------------------------------
! Begin declarations.

  IMPLICIT NONE

  INTEGER, INTENT( IN )  &
  :: ids, ide, jds, jde, kds, kde,  &
     ims, ime, jms, jme, kms, kme,  &
     its, ite, jts, jte, kts, kte

  REAL, DIMENSION( kms:kme ), INTENT(INOUT)  &
  :: area_density

  INTEGER, INTENT( IN ) &
  :: area_density_opt

  INTEGER, INTENT(INOUT) &
  :: canopy_top

  REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN)  &
  :: z !at the midpoint
  
  REAL, INTENT( IN ) &
  :: canopy_height

! Local variables.
  REAL &
  :: area_density_decay

  INTEGER  &
  :: i, j, k

! End declarations.
!-----------------------------------------------------------------------

  !Define area_density function 
  ! IF ( area_density_opt .eq. 3 ) THEN
     ! !exp^2 to match canopy stress model (TEMPORARY!!!!)
     ! area_density_decay = 1.0 / ( 0.18 * canopy_height )
  ! ELSE
     area_density_decay = -log(1e-3) / canopy_height
  ! ENDIF

  !Set area_density array
   i = its
   j = jts
   DO k = kts, kte-1 !z(kte)=0
      !print*,"k=",k,"z=",z(i,k,j)
      IF ( z(i,k,j) .LE. canopy_height ) THEN
         if ( area_density_opt .eq. 0 ) then
           !exponential
           area_density(k) = exp( -area_density_decay * z(i,k,j) )
         else if ( area_density_opt .eq. 1 ) then
           !cos^2
           area_density(k) = cos( piconst * z(i,k,j) / 2.0 / canopy_height )**2
         else if ( area_density_opt .eq. 2 ) then
           !new exp
           area_density(k) = exp( -area_density_decay * z(i,k,j) / 2.0 ) - exp( -area_density_decay * canopy_height / 2.0 )
         else if ( area_density_opt .eq. 3 ) then
           !exp^2
           area_density(k) = exp( -2.0 * area_density_decay * z(i,k,j) )
         else if ( area_density_opt .eq. 4 ) then
           !linear
           area_density(k) = 1.0 - ( 1.0 / canopy_height ) * z(i,k,j)
         endif
      ELSE
         area_density(k) = 0.0
      ENDIF
   ENDDO

   !Find canopy_top
   DO k = kts, kte
     IF (area_density(k) .eq. 0.0) THEN
       canopy_top = k-1
       EXIT
     ENDIF
   ENDDO

   !RSAtest
   ! DO k = kts, kte
   !   print*,"k=",k," area_density=",area_density(k)
   ! ENDDO
   ! print*,"canopy_top=",canopy_top
   !RSAtest end
  
  END SUBROUTINE set_area_density_2

!=======================================================================
!=======================================================================
!Function to calculate drag coefficient for modified canopy model

  SUBROUTINE calculate_cd_canopy_2( Cd_canopy, ust_pseudo, ust,    &
                                    u, v, z, zw,                   &
                                    Vmag, z0_pseudo,               &
                                    canopy_height,                 &
                                    canopy_top,                    &
                                    area_density,                  &
                                    area_density_opt,              &
                                    freeslip,                      &
                                    ids, ide, jds, jde, kds, kde,  &
                                    ims, ime, jms, jme, kms, kme,  &
                                    its, ite, jts, jte, kts, kte   )

!-----------------------------------------------------------------------
! Begin declarations.

  IMPLICIT NONE

  INTEGER, INTENT( IN )  &
  :: ids, ide, jds, jde, kds, kde,  &
     ims, ime, jms, jme, kms, kme,  &
     its, ite, jts, jte, kts, kte

  REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) & 
  :: Cd_canopy, ust_pseudo

  REAL, DIMENSION( ims:ime, jms:jme ), INTENT( IN ) & 
  :: ust

  REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN ) &
  :: u, v, z, zw, Vmag

  REAL, INTENT( IN ) &
  :: z0_pseudo, canopy_height

  REAL, DIMENSION( kms:kme ), INTENT(INOUT)  &
  :: area_density

  INTEGER, INTENT( IN ) &
  :: area_density_opt, canopy_top, freeslip

! Local variables.

  INTEGER  &
  :: i, j, k, i_start, i_end, j_start, j_end

  REAL &
  :: Cd0, uc, vc, tau_13, tau_23, area_density_decay, &
     Vb, hPBL, zk, zkp1, wlo, whi, Cd_eff, area_density_int

  REAL, DIMENSION( kms:kme ) &
  :: Vhc

! End declarations.
!-----------------------------------------------------------------------
  i_start = its
  i_end   = MIN(ite,ide-1)
  j_start = jts
  j_end   = MIN(jte,jde-1)

  ! IF ( area_density_opt .eq. 3 ) THEN
     ! !exp^2 to match canopy stress model (TEMPORARY!!!!)
     ! area_density_decay = 1.0 / ( 0.18 * canopy_height )
  ! ELSE
     area_density_decay = -log(1e-3) / canopy_height
  ! ENDIF

  DO j = j_start, j_end
  DO i = i_start, i_end

     !Calculate effective surface drag
     IF (freeslip .eq. 0) THEN
        Cd_eff = ( karman / log( z(i,1,j) / z0_pseudo ) )**2 - ( karman / log( z(i,1,j) / 0.0001 ) )**2
     ELSEIF (freeslip .eq. 1) THEN
        Cd_eff = ( karman / log( z(i,1,j) / z0_pseudo ) )**2
     ENDIF

     !Numerically integrate area density
     area_density_int = 0.0
     DO k = kts, canopy_top
        area_density_int = area_density_int &
                         + ( area_density(k) * ( zw(i,k+1,j) - zw(i,k,j) ) )
     END DO

     !Calcule Cd_canopy (Cd_eff distributed over the canopy height)
     Cd_canopy(i,j) = Cd_eff / area_density_int

  END DO
  END DO

  END SUBROUTINE calculate_cd_canopy_2

!=======================================================================
!=======================================================================
! Adds canopy drag to the momentum equation based on Shaw & Patton 2003

  SUBROUTINE canopy_drag_pseudo( ru_tendf, rv_tendf, rw_tendf,         & 
                                 u, v, w, rho, mu,                     &
                                 area_density,                         &
                                 Vmag, Cd_canopy,                      & 
                                 canopy_top,                           &
                                 config_flags,                         &
                                 fnm, fnp,                             &
                                 ids, ide, jds, jde, kds, kde,         &
                                 ims, ime, jms, jme, kms, kme,         &
                                 its, ite, jts, jte, kts, kte          )

!-----------------------------------------------------------------------
! Begin declarations.

    IMPLICIT NONE

    TYPE( grid_config_rec_type ), INTENT( IN )  &
    :: config_flags

    INTEGER, INTENT( IN )  &
    :: ids, ide, jds, jde, kds, kde,  &
       ims, ime, jms, jme, kms, kme,  &
       its, ite, jts, jte, kts, kte

    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( INOUT )  &
    :: ru_tendf, rv_tendf, rw_tendf

    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN )  &
    :: u, v, w, rho, Vmag

    REAL , DIMENSION( ims:ime, jms:jme), INTENT( IN ) &
    :: mu, Cd_canopy

    REAL, DIMENSION( kms:kme ), INTENT( IN ) &
    :: area_density, fnm, fnp

    INTEGER, INTENT( IN ) &
    :: canopy_top

! Local variables.

    INTEGER  &
    :: i, j, k, i_start, i_end, j_start, j_end

    REAL  &
    :: tmpCd, tmpvmag, tmpmu

! End declarations.
!-----------------------------------------------------------------------

   !U TENDENCY
   i_start = its
   i_end   = ite
   j_start = jts
   j_end   = MIN(jte,jde-1)

   IF ( config_flags%open_xs .or. config_flags%specified .or. &
        config_flags%nested) i_start = MAX(ids+1,its)
   IF ( config_flags%open_xe .or. config_flags%specified .or. &
        config_flags%nested) i_end   = MIN(ide-1,ite)
   IF ( config_flags%open_ys .or. config_flags%specified .or. &
        config_flags%nested) j_start = MAX(jds+1,jts)
   IF ( config_flags%open_ye .or. config_flags%specified .or. &
        config_flags%nested) j_end   = MIN(jde-2,jte)
      IF ( config_flags%periodic_x ) i_start = its
      IF ( config_flags%periodic_x ) i_end = ite

    DO j = j_start, j_end
    DO k = kts, canopy_top
    DO i = i_start, i_end
       tmpCd = 0.5 * ( Cd_canopy(i-1,j) + Cd_canopy(i,j) ) 
       tmpvmag = 0.5 * ( Vmag(i-1,k,j) + Vmag(i,k,j) ) 
       tmpmu = 0.5 * ( mu(i-1,j) + mu(i,j) ) 
       ru_tendf(i,k,j) = ru_tendf(i,k,j) &
                         - tmpCd * area_density(k) * u(i,k,j) * tmpvmag * tmpmu 
    END DO
    END DO
    END DO

    !V TENDENCY
    i_start = its
    i_end   = MIN(ite,ide-1)
    j_start = jts
    j_end   = jte

    IF ( config_flags%open_xs .or. config_flags%specified .or. &
         config_flags%nested) i_start = MAX(ids+1,its)
    IF ( config_flags%open_xe .or. config_flags%specified .or. &
         config_flags%nested) i_end   = MIN(ide-2,ite)
    IF ( config_flags%open_ys .or. config_flags%specified .or. &
         config_flags%nested) j_start = MAX(jds+1,jts)
    IF ( config_flags%open_ye .or. config_flags%specified .or. &
         config_flags%nested) j_end   = MIN(jde-1,jte)
       IF ( config_flags%periodic_x ) i_start = its
       IF ( config_flags%periodic_x ) i_end = MIN(ite,ide-1)

    DO j = j_start, j_end
    DO k = kts, canopy_top
    DO i = i_start, i_end
       tmpCd = 0.5 * ( Cd_canopy(i,j-1) + Cd_canopy(i,j) )
       tmpvmag = 0.5 * ( Vmag(i,k,j-1) + Vmag(i,k,j) )
       tmpmu = 0.5 * ( mu(i,j-1) + mu(i,j) )
       rv_tendf(i,k,j) = rv_tendf(i,k,j) &
                         - tmpCd * area_density(k) * v(i,k,j) * tmpvmag * tmpmu
    END DO
    END DO
    END DO

    !W TENDENCY
    i_start = its
    i_end   = MIN(ite,ide-1)
    j_start = jts
    j_end   = MIN(jte,jde-1)
 
    IF ( config_flags%open_xs .or. config_flags%specified .or. &
         config_flags%nested) i_start = MAX(ids+1,its)
    IF ( config_flags%open_xe .or. config_flags%specified .or. &
         config_flags%nested) i_end   = MIN(ide-2,ite)
    IF ( config_flags%open_ys .or. config_flags%specified .or. &
         config_flags%nested) j_start = MAX(jds+1,jts)
    IF ( config_flags%open_ye .or. config_flags%specified .or. &
         config_flags%nested) j_end   = MIN(jde-2,jte)
       IF ( config_flags%periodic_x ) i_start = its
       IF ( config_flags%periodic_x ) i_end = MIN(ite,ide-1)

    DO j = j_start, j_end
    DO k = kts+1, canopy_top
    DO i = i_start, i_end
       tmpvmag = fnm(k) * Vmag(i,k,j) + fnp(k) * Vmag(i,k-1,j)
       rw_tendf(i,k,j) = rw_tendf(i,k,j) &
                         - Cd_canopy(i,j) * area_density(k) * w(i,k,j) * tmpvmag * mu(i,j)
    END DO
    END DO
    END DO

  END SUBROUTINE canopy_drag_pseudo

!=======================================================================

END MODULE module_canopy_2
