# =======================================================================================
# DAIS-fortran90 (# estimation by calling fortran routine)
# DAIS: Simple model for Antarctic ice-sheet volume [m sle] (Schaffer 2014)
# =======================================================================================
#
#  Requires (input variables):
#  - Ta        Antarctic mean surface temperature [degC]
#  - SL        Sea level [m]
#  - Toc       High latitude ocean subsurface temperatures [degC]
#
#  Simulates (output variables):
#  - Rad       Ice sheet radius [m]
#  - Vais      Volume of Antarctic ice sheet [m3]

#  Internal variables:
#  - b         Undisturbed bed profile [m]
#  - h         Ice sheet surface height [m]
#  - Btot      Total mass accumulation rate on the ice-sheet [ ??? ]
#  - F         Total ice flux across the grounding line [m3/yr]
#  - rc        Distance from the continent center to where the ice sheets enters the sea [m]
#  - hr        Height of runoff line above which precipitation accumulates as snow [m]
#  - P         Precipitation [m (ice equivalent)]
#  - beta      Mass balance gradient [ m-0.5 ]
#  - rR        Distance from the continent center to where the runoff line intersects the ice sheet surface [m]
#  - Hw        Water depth at grounding line
#
#  Coordinates:
#  - r         Radial coordinate
#
#  Parameters:
#  - b0        Undisturbed bed height at the continent center [m]
#  - slope     Slope of ice sheet bed before loading [-]
#  - mu        Profile parameter for parabolic ice sheet surface (related to ice stress) [m0.5]
#  - h0        hr(Ta=0): Height of runoff line at Ta = 0 [m]
#  - c         Sensitivity of Height of runoff line (hr) [m/degC]
#  - P0        P(Ta=0): Annual precipitation for Ta = 0 [m (ice equivalent)]
#  - kappa     Coefficient for the exponential dependency of precipitation on Ta [degC-1]
#  - nu        Proportionality constant relating runoff decrease with height to precipitation [m^(-1/2) yr^(-1/2)]
#  - f0        Proportionality constant for ice flow at grounding line [m/yr]
#  - gamma     Power for the relation of ice flow speed to water depth [-]
#  - alpha     Partition parameter for effect of ocean subsurface temperature on ice flux [-]
#  - Toc_0     Present-day, high latitude ocean subsurface temperature [degC]
#  - Rad0      Reference ice sheet radius [m]
#  - dSL0      Initial sea level rate
#  - tstep     time step
#
#  Constants:
#  - Tf        Freecing temperature sea water [degC]
#  - rho_w      (Sea) water density [kg/m3]
#  - rho_i      Ice density [kg/m3]
#  - rho_m      Rock density [kg/m3]
# =======================================================================================

# =======================================================================================
# load fortran subroutine (# to check if library is loaded is.loaded("run_dais") )
dyn.load("../fortran/dais.so")

daisF <- function(
  tstep = 1,
  b0    = 775,
  slope = 6 * 10^(-4),
  mu    = 8.7,
  h0    = 1471,
  c     = 95,
  P0    = 0.35,
  kappa = 4 * 10^(-2),
  nu    = 1.2 * 10^(-2),
  f0    = 1.2,
  gamma = 2.5,
  alpha = 0.5,
  Tf    = -1.8,
  rho_w  = 1030,
  rho_i  = 917,
  rho_m  = 4000,
  Toc_0 = 0.72,
  Rad0  = 1.864 * 10^6,
  dSL0  = 0,
  Ta,        # Antarctic mean surface temperature [degC]
  SL,        # (global mean) sea level [m]
  Toc        # High latitude ocean subsurface temperatures [degC]
) {

  # determine series length
  ns <- length(Ta)
  
  # call fortran
  f.output <- .Fortran("run_dais",
                  ns            = ns,
                  tstep         = as.double(tstep),
                  dais_b0       = as.double(b0),
                  dais_slope    = as.double(slope),
                  dais_mu       = as.double(mu),
                  dais_h0       = as.double(h0),
                  dais_c        = as.double(c),
                  dais_P0       = as.double(P0),
                  dais_kappa    = as.double(kappa),
                  dais_nu       = as.double(nu),
                  dais_f0       = as.double(f0),
                  dais_gamma    = as.double(gamma),
                  dais_alpha    = as.double(alpha),
                  dais_Toc_0    = as.double(Toc_0),
                  dais_Rad0     = as.double(Rad0),
                  dais_dSL0     = as.double(dSL0),
                  dais_Tf       = as.double(Tf),
                  dais_rho_w    = as.double(rho_w),
                  dais_rho_i    = as.double(rho_i),
                  dais_rho_m    = as.double(rho_m),
                  Ant_Temp      = as.double(Ta),
                  Ant_Sea_Level = as.double(SL),
                  Ant_Sur_Ocean_Temp = as.double(Toc),
                  AIS_Radius_out = as.double(rep(-999.99,ns)),
                  AIS_Volume_out = as.double(rep(-999.99,ns))
  )
  Vsle = 57*(1-f.output$AIS_Volume_out/Volo) #Takes steady state present day volume to correspond to 57m SLE
  
  
    return(Vsle)
  
}
# =================================================================================