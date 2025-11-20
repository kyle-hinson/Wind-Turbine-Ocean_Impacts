!WRF:MEDIATION_LAYER:PHYSICS

MODULE module_canopy_1

USE module_model_constants
USE module_big_step_utilities_em, only: grid_config_rec_type, zero_tend, sixth_order_diffusion, set_tend, &
        horizontal_diffusion, vertical_diffusion, vertical_diffusion_mp
USE module_state_description, only: p_qv, tiedtkescheme, ntiedtkescheme, positivedef, monotonic, wenopd_scalar, weno_scalar, &
        gdscheme, g3scheme, gfscheme, kfetascheme, mskfscheme
USE module_advect_em, only: advect_scalar, advect_scalar_pd, advect_scalar_mono, advect_scalar_weno, advect_scalar_wenopd
USE module_diffusion_em, only: calc_l_scale, cal_dampkm, calculate_N2

IMPLICIT NONE

REAL, PARAMETER :: lf=0.1 !canopy element length scale

CONTAINS

!=======================================================================
! Sets the area density for canopy_opt = 1

  SUBROUTINE set_area_density_1( area_density,                     & 
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

  INTEGER, INTENT(INOUT) &
  :: canopy_top

  REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN)  &
  :: z !at the midpoint

  REAL, INTENT( IN ) &
  :: canopy_height

! Local variables.
  REAL, DIMENSION(1:11) &
  :: a_LAI_2, z_LAI_2_norm, z_LAI_2

  REAL &
  :: fac_z

  INTEGER  &
  :: i, j, k, k_canopy

! End declarations.
!-----------------------------------------------------------------------

  !Define area_density function from Dwyer et al 1997
  a_LAI_2 = (/ 0.0447, 0.0498, 0.0642, 0.0910, 0.1249, 0.1477, 0.1515, 0.1460, &
               0.1231, 0.0667, 0.0 /)

  z_LAI_2_norm = (/ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 /)

  z_LAI_2 = z_LAI_2_norm * canopy_height

  !Set area_density array
  !RSA this will only work for flat terrain and homogeneous canopy right now.
  !For terrain and/or heterogeneous canopy, need area_density and canopy_top to
  !be functions of x,y. Also will need to include HGT to get z AGL instead of
  !raw z.
!  DO j = jts, jte
!  DO i = its, ite
   i = its
   j = jts
   DO k = kts, kte-1 !z(kte)=0
      area_density(k) = 0.0
      !print*,"k=",k,"z=",z(i,k,j)
      DO k_canopy = 1, 10
        IF ( ( z(i,k,j) .GE. z_LAI_2(k_canopy) ) .AND. ( z(i,k,j) .LT. z_LAI_2(k_canopy+1) ) ) THEN 
          fac_z = ( z(i,k,j) - z_LAI_2(k_canopy) ) / ( z_LAI_2(k_canopy+1) - z_LAI_2(k_canopy) )
          area_density(k) = (1.0-fac_z) * a_LAI_2(k_canopy) + fac_z * a_LAI_2(k_canopy+1) 
          !print*,"k_canopy=",k_canopy,"z_LAI(k_canopy)=",z_LAI_2(k_canopy),"z_LAI(k_canopy+1)=",z_LAI_2(k_canopy+1)
          EXIT !move on to z(k+1)
       ENDIF
     ENDDO
   ENDDO
!  ENDDO
!  ENDDO

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
  
  END SUBROUTINE set_area_density_1

!=======================================================================
!=======================================================================
!Function to calculate velocity magnitude for canopy les model

  SUBROUTINE calculate_vmag( Vmag, u, v, w, canopy_top,            &
                             ids, ide, jds, jde, kds, kde,         &
                             ims, ime, jms, jme, kms, kme,         &
                             its, ite, jts, jte, kts, kte          )

!-----------------------------------------------------------------------
! Begin declarations.

  IMPLICIT NONE

    INTEGER, INTENT( IN )  &
    :: ids, ide, jds, jde, kds, kde,  &
       ims, ime, jms, jme, kms, kme,  &
       its, ite, jts, jte, kts, kte

    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(INOUT) &
    :: Vmag

    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN ) &
    :: u, v, w

    INTEGER, INTENT( IN ) &
    :: canopy_top

! Local variables.

    INTEGER &
    :: i, j, k, i_start, i_end, j_start, j_end

! End declarations.
!-----------------------------------------------------------------------
    i_start = its
    i_end   = MIN(ite,ide-1)
    j_start = jts
    j_end   = MIN(jte,jde-1)

    !Center points
    DO j = j_start, j_end
    DO k = kts, canopy_top
    DO i = i_start, i_end

       Vmag(i,k,j) = sqrt( ( 0.5 * ( u(i  ,k  ,j  ) + u(i+1,k  ,j  ) ) )**2 + &
                           ( 0.5 * ( v(i  ,k  ,j  ) + v(i  ,k  ,j+1) ) )**2 + & 
                           ( 0.5 * ( w(i  ,k  ,j  ) + w(i  ,k+1,j  ) ) )**2   )

    END DO
    END DO
    END DO

  END SUBROUTINE calculate_vmag

!=======================================================================
!=======================================================================
!Function to calculate drag coefficient for canopy les model

  SUBROUTINE calculate_cd_canopy_1( Cd_canopy,                     &
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

  REAL, INTENT(INOUT) & 
  :: Cd_canopy

! Local variables.

  INTEGER  &
  :: i, j, k

! End declarations.
!-----------------------------------------------------------------------

  Cd_canopy = 0.15

  END SUBROUTINE calculate_cd_canopy_1

!=======================================================================
!=======================================================================
!Function to calculate skin friction drag coefficient for canopy les model

  SUBROUTINE calculate_cd_skin_canopy( Cd_skin_canopy,                       &
                                       Vmag, canopy_top,                     &
                                       ids, ide, jds, jde, kds, kde,         &
                                       ims, ime, jms, jme, kms, kme,         &
                                       its, ite, jts, jte, kts, kte          )

!-----------------------------------------------------------------------
! Begin declarations.

    IMPLICIT NONE

    INTEGER, INTENT( IN )  &
    :: ids, ide, jds, jde, kds, kde,  &
       ims, ime, jms, jme, kms, kme,  &
       its, ite, jts, jte, kts, kte

    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(INOUT) &
    :: Cd_skin_canopy

    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN ) &
    :: Vmag

    INTEGER, INTENT(IN) &
    :: canopy_top
   
  ! Local variables.

    INTEGER  &
    :: i, j, k, i_start, i_end, j_start, j_end

    REAL &
    :: Rl, kinematic_viscosity

! End declarations.
!-----------------------------------------------------------------------

    i_start = its
    i_end   = MIN(ite,ide-1)
    j_start = jts
    j_end   = MIN(jte,jde-1)

    kinematic_viscosity = 1.0E-05

    DO j = j_start, j_end
    DO k = kts, canopy_top
    DO i = i_start, i_end

       Rl = Vmag(i,k,j) * lf / kinematic_viscosity
       Cd_skin_canopy(i,k,j) = 1.328 / sqrt(Rl)

    END DO
    END DO
    END DO
    
   Cd_skin_canopy = 0.0

  END SUBROUTINE calculate_cd_skin_canopy

!=======================================================================
!=======================================================================
! Adds canopy drag to the momentum equation based on Shaw & Patton 2003

  SUBROUTINE canopy_drag_sp2003( ru_tendf, rv_tendf, rw_tendf,         & 
                                 u, v, w, rho, mu,                     &
                                 area_density,                         &
                                 Vmag, Cd_canopy,                      & 
                                 Cd_skin_canopy,                       &  
                                 canopy_top,                           &
                                 config_flags,                         &
                                 c1h, c2h,                             & 
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
    :: u, v, w, rho, Vmag, Cd_skin_canopy

    REAL , DIMENSION( ims:ime, jms:jme), INTENT( IN ) &
    :: mu

    REAL, DIMENSION( kms:kme ), INTENT( IN ) &
    :: area_density, fnm, fnp

    REAL, INTENT( IN ) &
    :: Cd_canopy

    INTEGER, INTENT( IN ) &
    :: canopy_top

    REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: c1h, c2h
! Local variables.

    INTEGER  &
    :: i, j, k, i_start, i_end, j_start, j_end

    REAL  &
    :: tmpskin, tmpvmag, tmpmu

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
       tmpskin = 0.5 * ( Cd_skin_canopy(i-1,k,j) + Cd_skin_canopy(i,k,j) ) 
       tmpvmag = 0.5 * ( Vmag(i-1,k,j) + Vmag(i,k,j) ) 
       tmpmu = 0.5 * ( mu(i-1,j) + mu(i,j) ) * c1h(k) + c2h(k) 
       ru_tendf(i,k,j) = ru_tendf(i,k,j) &
                         - (Cd_canopy + tmpskin) * area_density(k) &
                           * u(i,k,j) * tmpvmag * tmpmu
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
       tmpskin = 0.5 * ( Cd_skin_canopy(i,k,j-1) + Cd_skin_canopy(i,k,j) )
       tmpvmag = 0.5 * ( Vmag(i,k,j-1) + Vmag(i,k,j) )
       tmpmu = 0.5 * ( mu(i,j-1) + mu(i,j) ) * c1h(k) + c2h(k)
       rv_tendf(i,k,j) = rv_tendf(i,k,j) &
                         - (Cd_canopy + tmpskin) * area_density(k) &
                           * v(i,k,j) * tmpvmag * tmpmu
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
       tmpskin = fnm(k) * Cd_skin_canopy(i,k,j) + fnp(k) * Cd_skin_canopy(i,k-1,j)
       tmpvmag = fnm(k) * Vmag(i,k,j) + fnp(k) * Vmag(i,k-1,j)
       rw_tendf(i,k,j) = rw_tendf(i,k,j) &
                         - (Cd_canopy + tmpskin) * area_density(k) &
                           * w(i,k,j) * tmpvmag * (mu(i,j)* c1h(k) + c2h(k))
    END DO
    END DO
    END DO

  END SUBROUTINE canopy_drag_sp2003



!=======================================================================

END MODULE module_canopy_1
