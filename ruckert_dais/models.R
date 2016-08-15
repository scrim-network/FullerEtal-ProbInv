cModel <- "rob"         # use Rob's C model
#cModel <- "kelsey"     # use Kelsey's C model
#cModel <- NULL         # use Alex's Fortran model
#cModel <- NA           # use Kelsey's R model


if (exists("daisassimctx") && exists("cModel", env=daisassimctx)) {
    cModel <- daisassimctx$cModel
    print(paste("cModel is", daisassimctx$cModel))
} else {
    if (!is.null(cModel) && is.na(cModel)) {
        source("Scripts/DAIS_IceFlux_model.R")
        source("Scripts/iceflux.mult_func_outRHF.R")
    }
}


# for dynReload() function
source("roblib.R")


daisLoadModel <- function(cModel="rob")
{
    if (is.null(cModel)) {
        dynReload("../fortran/dais", makevars='PKG_FCFLAGS="-I../fortran -J../fortran"',
            srcname=paste("../fortran/src/", c("dais.f90", "run_dais.f90", "global.f90"), sep=""))
    } else {
        daisLib <- paste(cModel, "_dais", sep="")
        dynReload(daisLib, srcname=c(paste(daisLib, ".c", sep=""), "r.c"), extrasrc="r.h")
    }
}


if (is.null(cModel) || !is.na(cModel)) {

    # need an iceflux_RHF() for Fortran.  Rob's C model is close enough.
    model <- ifelse(is.null(cModel), "rob", cModel)
    print(paste("model is", model))

    daisLoadModel(model)
    daisCmodel <- paste("dais", toupper(substring(model, 1, 1)), substring(model, 2), "OdeC", sep="")

    iceflux_RHF <- function(iceflux, forcings, standards)
    {
        # avoid names like "b0.b0"
        iceflux <- unname(iceflux)
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

        .Call(daisCmodel, list(mp=mp, frc=forcings, out=list(SLE, Vais, Rad, Flow, Depth)))

        out <- list(SLE=SLE, Vol=Vais, Rad=Rad, IceFlux=Flow, WatDepth=Depth)

        return(out)
    }

    iceflux <- function(iceflux, forcings, standards)
    {
        out <- iceflux_RHF(iceflux, forcings, standards)
        return (out$SLE)
    }
}


if (is.null(cModel)) {
    print("loading Fortran model")
    daisLoadModel(NULL)
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
