useCmodel <- T
useFmodel <- F


if (!useCmodel && !useFmodel) {
    source("Scripts/DAIS_IceFlux_model.R")
    source("Scripts/iceflux.mult_func_outRHF.R")
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
          Ta     = forcings[, 1], 
          SL     = forcings[, 4],
          Toc    = forcings[, 2],
          dSL    = forcings[, 3])

        return (Volume_F)
    }
}


if (useCmodel) {
    source("roblib.R")
    dynReload("dais", srcname=c("dais.c", "r.c"), extrasrc="r.h")


    iceflux_RHF <- function(iceflux, forcings, standards)
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
        SLE    <- numeric(length=np)               # Sea-level equivalent [m]
        Vais   <- numeric(length=np)               # Ice volume
        Rad    <- numeric(length=np)               # Radius of ice sheet
        Flow   <- numeric(length=np)               # Ice flow
        Depth  <- numeric(length=np)               # Water depth

        .Call("daisOdeC", list(mp=mp, frc=forcings, out=list(SLE, Vais, Rad, Flow, Depth)), PACKAGE = "dais")

        out <- list(SLE=SLE, Vol=Vais, Rad=Rad, Flow=Flow, WatDepth=Depth)

        return(out)
    }


    iceflux <- function(iceflux, forcings, standards)
    {
        out <- iceflux_RHF(iceflux, forcings, standards)
        return (out$SLE)
    }
}
