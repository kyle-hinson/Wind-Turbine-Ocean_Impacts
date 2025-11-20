!WRF:MEDIATION_LAYER:PHYSICS

MODULE module_canopy_driver

CONTAINS

!=============================================================================

  SUBROUTINE canopy_driver( grid, config_flags,                   &
                            ru_tendf, rv_tendf, rw_tendf,         &
                            Cd_canopy_1, canopy_top               )

  ! Driver layer modules
    USE module_domain
    USE module_model_constants
    USE module_configure
    USE module_tiles
    USE module_machine
    USE module_state_description
  ! Model layer modules
    USE module_bc
  ! Canopy modules
    USE module_canopy_1
    USE module_canopy_4  !, ONLY : set_area_density_4, calculate_cd_canopy_4

#ifdef DM_PARALLEL
    USE module_dm
    USE module_comm_dm, ONLY : HALO_EM_CANOPY_PARAMS_sub, PERIOD_BDY_EM_CANOPY_PARAMS_sub   
#endif

    IMPLICIT NONE

  ! Input data.

    TYPE(domain) , TARGET :: grid

    TYPE (grid_config_rec_type), INTENT(IN) :: config_flags

    REAL, DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33), INTENT(INOUT) &
    :: ru_tendf, rv_tendf, rw_tendf

    REAL, INTENT(INOUT) &
    :: Cd_canopy_1

    INTEGER, INTENT(INOUT) &
    :: canopy_top

  ! Local data

    INTEGER :: k_start , k_end, its, ite, jts, jte
    INTEGER :: ids , ide , jds , jde , kds , kde , &
               ims , ime , jms , jme , kms , kme , &
               ips , ipe , jps , jpe , kps , kpe

    INTEGER :: imsx, imex, jmsx, jmex, kmsx, kmex, &
               ipsx, ipex, jpsx, jpex, kpsx, kpex, &
               imsy, imey, jmsy, jmey, kmsy, kmey, &
               ipsy, ipey, jpsy, jpey, kpsy, kpey

    INTEGER :: ij, i, j, k

    CALL get_ijk_from_grid ( grid ,                              &
                             ids, ide, jds, jde, kds, kde,       &
                             ims, ime, jms, jme, kms, kme,       &
                             ips, ipe, jps, jpe, kps, kpe,       &
                             imsx, imex, jmsx, jmex, kmsx, kmex, &
                             ipsx, ipex, jpsx, jpex, kpsx, kpex, &
                             imsy, imey, jmsy, jmey, kmsy, kmey, &
                             ipsy, ipey, jpsy, jpey, kpsy, kpey  )

    k_start = kps
    k_end   = kpe

    !CHOOSE CANOPY DRAG OPTION
    SELECT CASE (config_flags%canopy_opt )

!-----------------------------------------------------------------------------------    
! Drag term from Shaw & Patton 2003
!-----------------------------------------------------------------------------------    
      CASE(1)

          !$OMP PARALLEL DO   &
          !$OMP PRIVATE ( ij )
          DO ij = 1 , grid%num_tiles

!               CALL set_area_density_1( grid%area_density,                      &
!                                        config_flags%canopy_height,       &
!                                        canopy_top, grid%z,                     &
!                                        ids, ide, jds, jde, kds, kde,           &
!                                        ims, ime, jms, jme, kms, kme,           &
!                                        grid%i_start(ij), grid%i_end(ij),       &
!                                        grid%j_start(ij), grid%j_end(ij),       &
!                                        k_start    , k_end                      )
               CALL set_area_density_4(  grid%area_density,                      &
                                         model_config_rec%a_lai_2,               &
                                         model_config_rec%z_lai_2_norm,          &
                                         config_flags%LAD_param_a,               &
                                         config_flags%LAD_param_l,               & 
                                         config_flags%LAI_canopy,                & 
                                         config_flags%canopy_height,             &
                                         canopy_top, grid%z,                     &
                                         ids, ide, jds, jde, kds, kde,           &
                                         ims, ime, jms, jme, kms, kme,           &
                                         grid%i_start(ij), grid%i_end(ij),       &
                                         grid%j_start(ij), grid%j_end(ij),       &
                                         k_start    , k_end                      )

               CALL calculate_vmag( grid%Vmag,                         &
                                    grid%u_2, grid%v_2, grid%w_2,      &
                                    canopy_top,                        &
                                    ids, ide, jds, jde, kds, kde,      &
                                    ims, ime, jms, jme, kms, kme,      &
                                    grid%i_start(ij), grid%i_end(ij),  &
                                    grid%j_start(ij), grid%j_end(ij),  &
                                    k_start    , k_end                 )

               CALL calculate_cd_canopy_1( Cd_canopy_1,                       &
                                           ids, ide, jds, jde, kds, kde,      &
                                           ims, ime, jms, jme, kms, kme,      &
                                           grid%i_start(ij), grid%i_end(ij),  &
                                           grid%j_start(ij), grid%j_end(ij),  &
                                           k_start    , k_end                 )
               
               CALL calculate_cd_skin_canopy( grid%Cd_skin_canopy,               & 
                                              grid%Vmag, canopy_top,             &
                                              ids, ide, jds, jde, kds, kde,      &
                                              ims, ime, jms, jme, kms, kme,      &
                                              grid%i_start(ij), grid%i_end(ij),  &
                                              grid%j_start(ij), grid%j_end(ij),  &
                                              k_start    , k_end                 )

          ENDDO
          !$OMP END PARALLEL DO

#ifdef DM_PARALLEL
#     include "HALO_EM_CANOPY_PARAMS.inc"
#     include "PERIOD_BDY_EM_CANOPY_PARAMS.inc"
#endif

          !$OMP PARALLEL DO   &
          !$OMP PRIVATE ( ij )
          DO ij = 1 , grid%num_tiles

               CALL set_physical_bc3d( grid%Vmag, 'p', config_flags,       &
                                       ids, ide, jds, jde, kds, kde,       &
                                       ims, ime, jms, jme, kms, kme,       &
                                       ips, ipe, jps, jpe, kps, kpe,       &
                                       grid%i_start(ij), grid%i_end(ij),   &
                                       grid%j_start(ij), grid%j_end(ij),   &
                                       k_start    , k_end                 )

               CALL set_physical_bc3d( grid%Cd_skin_canopy, 'p', config_flags,  &
                                       ids, ide, jds, jde, kds, kde,            &
                                       ims, ime, jms, jme, kms, kme,            &
                                       ips, ipe, jps, jpe, kps, kpe,            &
                                       grid%i_start(ij), grid%i_end(ij),        &
                                       grid%j_start(ij), grid%j_end(ij),        &
                                       k_start    , k_end                       )

          ENDDO
          !$OMP END PARALLEL DO
          
          !$OMP PARALLEL DO   &
          !$OMP PRIVATE ( ij )
          DO ij = 1 , grid%num_tiles

               CALL canopy_drag_sp2003( ru_tendf, rv_tendf, rw_tendf,           &
                                        grid%u_2, grid%v_2, grid%w_2,           &
                                        grid%rho, grid%mut,                     &
                                        grid%area_density,                      &
                                        grid%Vmag, Cd_canopy_1,                 & 
                                        grid%Cd_skin_canopy,                    &  
                                        canopy_top,                             &
                                        config_flags,                           &
                                        grid%c1h, grid%c2h,                     &
                                        grid%fnm, grid%fnp,                     &
                                        ids, ide, jds, jde, kds, kde,           &
                                        ims, ime, jms, jme, kms, kme,           &
                                        grid%i_start(ij), grid%i_end(ij),       &
                                        grid%j_start(ij), grid%j_end(ij),       &
                                        k_start    , k_end                      )
               
          ENDDO
          !$OMP END PARALLEL DO


!-----------------------------------------------------------------------------------    
! Drag term from Shaw & Patton 2003  LAI INFO input from namelist
!-----------------------------------------------------------------------------------    
      CASE(4)

           
          !$OMP PARALLEL DO   &
          !$OMP PRIVATE ( ij )
          DO ij = 1 , grid%num_tiles

               CALL set_area_density_4( grid%area_density,                      &
                                        model_config_rec%a_lai_2,               &
                                        model_config_rec%z_lai_2_norm,          &
                                        config_flags%lad_param_a,               &
                                        config_flags%lad_param_l,               &
                                        config_flags%lai_canopy,                & 
                                        config_flags%canopy_height,             &
                                        canopy_top, grid%z,                     &
                                        ids, ide, jds, jde, kds, kde,           &
                                        ims, ime, jms, jme, kms, kme,           &
                                        grid%i_start(ij), grid%i_end(ij),       &
                                        grid%j_start(ij), grid%j_end(ij),       &
                                        k_start    , k_end                      )
!PRINT*, 'Area density = ',grid%area_density(1:5)
               CALL calculate_vmag_4( grid%Vmag,                         &
                                    grid%u_2, grid%v_2, grid%w_2,      &
                                    canopy_top,                        &
                                    ids, ide, jds, jde, kds, kde,      &
                                    ims, ime, jms, jme, kms, kme,      &
                                    grid%i_start(ij), grid%i_end(ij),  &
                                    grid%j_start(ij), grid%j_end(ij),  &
                                    k_start    , k_end                 )

               CALL calculate_cd_canopy_4( grid%Cd_canopy,                       &
                                           ids, ide, jds, jde, kds, kde,      &
                                           ims, ime, jms, jme, kms, kme,      &
                                           grid%i_start(ij), grid%i_end(ij),  &
                                           grid%j_start(ij), grid%j_end(ij),  &
                                           k_start    , k_end                 )
               
               CALL calculate_cd_skin_canopy_4( grid%Cd_skin_canopy,               & 
                                              grid%Vmag, canopy_top,             &
                                              ids, ide, jds, jde, kds, kde,      &
                                              ims, ime, jms, jme, kms, kme,      &
                                              grid%i_start(ij), grid%i_end(ij),  &
                                              grid%j_start(ij), grid%j_end(ij),  &
                                              k_start    , k_end                 )

          ENDDO
          !$OMP END PARALLEL DO

#ifdef DM_PARALLEL
#     include "HALO_EM_CANOPY_PARAMS.inc"
#     include "PERIOD_BDY_EM_CANOPY_PARAMS.inc"
#endif

          !$OMP PARALLEL DO   &
          !$OMP PRIVATE ( ij )
          DO ij = 1 , grid%num_tiles

               CALL set_physical_bc3d( grid%Vmag, 'p', config_flags,       &
                                       ids, ide, jds, jde, kds, kde,       &
                                       ims, ime, jms, jme, kms, kme,       &
                                       ips, ipe, jps, jpe, kps, kpe,       &
                                       grid%i_start(ij), grid%i_end(ij),   &
                                       grid%j_start(ij), grid%j_end(ij),   &
                                       k_start    , k_end                 )

               CALL set_physical_bc3d( grid%Cd_skin_canopy, 'p', config_flags,  &
                                       ids, ide, jds, jde, kds, kde,            &
                                       ims, ime, jms, jme, kms, kme,            &
                                       ips, ipe, jps, jpe, kps, kpe,            &
                                       grid%i_start(ij), grid%i_end(ij),        &
                                       grid%j_start(ij), grid%j_end(ij),        &
                                       k_start    , k_end                       )

          ENDDO
          !$OMP END PARALLEL DO
          
          !$OMP PARALLEL DO   &
          !$OMP PRIVATE ( ij )
          DO ij = 1 , grid%num_tiles

               CALL canopy_drag_sp2003_4( ru_tendf, rv_tendf, rw_tendf,           &
                                        grid%u_2, grid%v_2, grid%w_2,           &
                                        grid%rho, grid%mut,                     &
                                        grid%area_density,                      &
                                        grid%Vmag, grid%Cd_canopy,              &  ! TL
                                        grid%Cd_skin_canopy,                    &  
                                        canopy_top,                             &
                                        config_flags,                           &
                                        grid%c1h, grid%c2h,                     &
                                        grid%fnm, grid%fnp,                     &
                                        ids, ide, jds, jde, kds, kde,           &
                                        ims, ime, jms, jme, kms, kme,           &
                                        grid%i_start(ij), grid%i_end(ij),       &
                                        grid%j_start(ij), grid%j_end(ij),       &
                                        k_start    , k_end                      )
               
          ENDDO
          !$OMP END PARALLEL DO

    END SELECT

!-----------------------------------------------------------------------------------    

  END SUBROUTINE canopy_driver

END MODULE module_canopy_driver
