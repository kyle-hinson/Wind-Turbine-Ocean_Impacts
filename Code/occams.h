/*
** svn $Id: estuary_test.h 429 2009-12-20 17:30:26Z arango $
*******************************************************************************
** Copyright (c) 2002-2010 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Estuary test.
**
** Application flag:   ESTUARY_TEST2
** Input script:       ocean_estuary_test2.in
*/

#define ROMS_MODEL

/* ROMS OPTIONS */
#define SALINITY
#define NONLIN_EOS
#define SOLVE3D
#define SPLINES_VVISC
#define SPLINES_VDIFF

#define DJ_GRADPS
#define CURVGRID

#define UV_ADV
#define UV_COR
#define UV_VIS2
#define MIX_S_UV
#define MIX_S_TS
#define TS_DIF2     /* Harmoic horizontal mixing of tracers */

/* Bottom Boundary Conditions */
#define UV_QDRAG

/* Turbulence Closure */
#define GLS_MIXING
#define LIMIT_VDIFF
#define LIMIT_VVISC
#if defined GLS_MIXING
# define KANTHA_CLAYSON
# define N2S2_HORAVG
# define RI_SPLINES
#endif

#define MASKING

/* Undefine if starting from a file */
#undef ANA_INITIAL

/* Surface Forcing */
#define BULK_FLUXES   /* Adapatation of COARE Fairall et al. 1996) */
#ifdef BULK_FLUXES
#  define EMINUSP       /* Evaporation and Precipitation */
#  define ATM_PRESS     /* Sea level pressure (inverse barometer) */
#  define LONGWAVE_OUT  /* Compute upward lwrad from SST (needs lwrad_down) */
#  undef LONGWAVE       /* Needs cloud */
#  define SOLAR_SOURCE /* Activate if the net shortwave radiation needs to penetrate layers */
#  undef COOL_SKIN     /* Correction if you input downward (not net) short wave radiation */
#  undef DIURNAL_SRFLUX /* use to impose shortwave radiation local diurnal cycle */
#else
#  define ANA_SMFLUX
#  define ANA_STFLUX
#  define ANA_SSFLUX
#  define ANA_SPFLUX
#  define ANA_SRFLUX   /* Analytical surface shortwave radiation flux */
#endif

/* Bottom Layer*/
#define ANA_BTFLUX
#define ANA_BSFLUX
#define ANA_BPFLUX

/* Open Boundary Conditions*/
#undef ANA_FSOBC   /* Free surface Open Boundary Condition */
#undef ANA_M2OBC   /* 2D Momentum Open Boundary Condition */
#undef ANA_TOBC    /* Analytical tracer open boundary conditions */

/* Tides */
#undef ADD_FSOBC  /* use to add tidal elevations to procssed OBC data */
#undef ADD_M2OBC  /* use to add tidal currents to processed OBC data */
#undef SSH_TIDES  /* use if imposing tidal elevation */
#undef UV_TIDES   /* use if imposing tidal currents */
#undef RAMP_TIDES /* use if ramping tidal forcing */

/* Output */
#define STATIONS  /* Extract station output */
