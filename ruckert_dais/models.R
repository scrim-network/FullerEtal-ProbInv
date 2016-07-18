useCmodel <- F
useFmodel <- T


if (!useCmodel && !useFmodel) {
    source("Scripts/DAIS_IceFlux_model.R")
}


if (useFmodel) {
    source("daisF.R")

    iceflux <- function(iceflux, forcings, standards)
    {
        Volume_F <- daisF(
          tstep = 1,
          b0    = iceflux[10],
          slope = iceflux[11],
          mu    = iceflux[3],
          h0    = iceflux[8],
          c     = iceflux[9],
          P0    = iceflux[5],
          kappa = iceflux[6],
          nu    = iceflux[4],
          f0    = iceflux[7],
          gamma = iceflux[1],
          alpha = iceflux[2],
          Tf    = -1.8,             #Freezing temperature of sea water
          rho_w = 1030,             #Density of sea water [g/cm^3]
          rho_i = 917,              #Density of ice water [g/cm^3]
          rho_m = 4000,             #Density of rock [g/cm^3]
          Toc_0 = 0.72,             #Present day high latitude ocean subsurface temperature [K]
          Rad0  = 1.8636e6,         #Steady state AIS radius for present day Ta and SL [m]
          dSL0   = 0,
          Ta     = forcings[, 1], 
          SL     = forcings[, 4],
          Toc    = forcings[, 2])

        return (Volume_F)
    }
}


if (useCmodel) {
    source("roblib.R")
    dynReload("dais", srcname=c("dais.c", "r.c"), extrasrc="r.h")

    iceflux <- function(iceflux, forcings, standards)
    {
        mp <- c(
          b0    = iceflux[10],
          slope = iceflux[11],
          mu    = iceflux[3],
          h0    = iceflux[8],
          c     = iceflux[9],
          P0    = iceflux[5],
          kappa = iceflux[6],
          nu    = iceflux[4],
          f0    = iceflux[7],
          gamma = iceflux[1],
          alpha = iceflux[2],
          Tf    = -1.8,             #Freezing temperature of sea water
          rho_w = 1030,             #Density of sea water [g/cm^3]
          rho_i = 917,              #Density of ice water [g/cm^3]
          rho_m = 4000,             #Density of rock [g/cm^3]
          Toc_0 = 0.72,             #Present day high latitude ocean subsurface temperature [K]
          Rad0  = 1.8636e6          #Steady state AIS radius for present day Ta and SL [m]
        )

        np     <- nrow(forcings)
        Rad    <- numeric(length=np)               # Radius of ice sheet
        Vais   <- numeric(length=np)               # Ice volume
        SLE    <- numeric(length=np)               # Sea-level equivalent [m]

        .Call("daisOdeC", list(mp=mp, frc=forcings, out=list(Rad, Vais, SLE)), PACKAGE = "dais")

        return(SLE)
    }
}
