# =======================================================================================
# DAIS: Simple model for Antarctic ice-sheet volume [m sle] (Schaffer 2014)
# =======================================================================================
# 
#  Requires (input variables):
#  - Ta        Antarctic mean surface temperature [degC]
#  - SL        (global mean) sea level [m]
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
#
#  Coordinates:
#  - r         Radial coordinate
#
#  Parameters:
#  - b0        Undisturbed bed height at the continent center [m]
#  - slope     Slope of ice sheet bed before loading [-]
#  - mu        Profile parameter for parabolic ice sheet surface (related to ice stress) [m0.5]
#  - ho        hr(Ta=0): Height of runoff line at Ta = 0 [m]
#  - c         Sensitivity of Height of runoff line (hr) [m/degC]
#  - P0        P(Ta=0): Annual precipitation for Ta = 0 [m (ice equivalent)]
#  - kappa     Coefficient for the exponential dependency of precipitation on Ta [degC-1]
#  - nu        Proportionality constant relating runoff decrease with height to precipitation [m^(-1/2) yr^(-1/2)]
#  - f0        Proportionality constant for ice flow at grounding line [m/yr]
#  - gamma     Power for the relation of ice flow speed to water depth [-]
#  - alpha     Partition parameter for effect of ocean subsurface temperature on ice flux [-]
#  - tstep     time step
#
#  Constants:
#  - Tf        Freecing temperature sea water [degC]
#  - ro_w      (Sea) water density [kg/m3]
#  - ro_i      Ice density [kg/m3]
#  - ro_m      Rock density [kg/m3]
#  - Toc_0     Present-day, high latitude ocean subsurface temperature [degC]
#  - Rad0      Reference ice sheet radius [m]
#
# =======================================================================================

dais <- function(
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
  ro_w  = 1030,
  ro_i  = 917,
  ro_m  = 4000,
  Toc_0 = 0.72,
  Rad0  = 1.864 * 10^6,
  tstep = 1) {
  
  # Initialize intermediate parameters
  del  <- ro_w/ro_i                # Ratio sea water and ice density [-]
  eps1 <- ro_i/(ro_m - ro_i)       # Ratio ice density and density difference between rock and ice [-]
  eps2 <- ro_i/(ro_m - ro_i)       # Ratio sea water density and density difference between rock and ice [-]
  
  # Define vectors with state variables
  np     <- length(Ta)
  Rad    <- rep(NA,np)               # Radius of ice sheet
  Vais   <- rep(NA,np)               # Ice volume
#  SLE    <- rep(NA,np)               # Sea-level equivalent [m]
  
  # Initial conditions
  R       <- Rad0                  # gets updated at end of loop
  rc      <- (b0 - SL[1])/slope    # application of equation 1 (paragraph after eq3)
  Vais[1] <- pi * (1+eps1) * ( (8/15) *  mu^0.5 * R^2.5 - 1/3*slope*R^3 ) -
             pi*eps2 * ( (2/3) * slope*(R^3-rc^3)-b0*(R^2-rc^2) )

  # Run model
  for(i in 2:np) { 
  
    hr   <- ho + c * Ta[i-1]        # equation 5
    rc   <- (b0 - SL[i-1])/slope    # application of equation 1 (paragraph after eq3)
    P    <- P0 * exp(kappa*Ta[i-1]) # equation 6
    beta <- nu * P^(0.5)            # equation 7 (corrected with respect to text)
    
    # Total mass accumulation on ice sheet (equation 8)
    if(hr > 0) {
      rR   <- R - ((hr - b0 + slope*R)^2) / mu

      Btot <- P * pi * R^2 - 
        pi * beta * (hr - b0 + slope*R) * (R*R - rR*rR) -
        (4 * pi * beta * mu^0.5 * (R-rR)^2.5) / 5  +
        (4 * pi * beta * mu^0.5 * R*(R-rR)^1.5) / 3 
    } else {
      Btot <- P * pi*R^2
    }
    
    # Ice flux at grounding line (F), ISO and fac
    if (R <= rc) {
      # no grounding line / ice shelves
      mit <- 0 # marine ice term (if mit=0 marine ice sheet terms are ignored)
      F   <- 0                                                     # no ice flux
      ISO <- 0                                                     # (third term equation 14) NAME?
    } else {
      # marine ice sheet / grounding line
      mit <- 1
      
      Hw <- slope*R - b0 + SL[i-1] # (equation 10) 
      
      # Ice speed at grounding line (equation 11)
      # KELSEY: last term is different than in manuscript! (slope*Rad0 - b0) rather than (b0 - slope*Rad0)
      # KELSEY: this seems to corrected for in equation 9 that misses a minus sign (wrt ms)
      Speed <- f0 * 
        ((1-alpha) + alpha* ((Toc[i-1] - Tf)/(Toc_0 - Tf))^2) *
        (Hw^gamma) / ( (slope*Rad0 - b0)^(gamma-1) )

      F     <- 2*pi*R * del * Hw * Speed # (equation 9)
      
      ISO   <- 2*pi*eps2* (slope*rc^2 - (b0/slope)*rc) * (SL[i] - SL[i-1]) # third term equation 14 !! NAME?
    }

    # dV/dR (equation 14)
    fac   <- pi * (1+eps1) * (4/3 * mu^0.5 * R^1.5 - slope*R^2) -
             mit * ( 2*pi*eps2 * (slope*R^2 - b0*R) )               # in case of marine ice sheet (mit=1)
    
    # Ice sheet volume (equation 13)
    # KELSEY: it seems some parantheses were missing in your code
    R       <- R + (Btot-F+ISO)/fac
    Rad[i]  <- R
    Vais[i] <- pi * (1+eps1) * ( (8/15) *  mu^0.5 * R^2.5 - 1/3*slope*R^3 ) -
               mit * pi*eps2 * ( (2/3) * slope*(R^3-rc^3)-b0*(R^2-rc^2) )
  }
  Vsle = 57*(1-Vais/Volo) #Takes steady state present day volume to correspond to 57m SLE

#return(Vsle)
  
return(Vais)

}

