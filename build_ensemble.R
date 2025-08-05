library(nleqslv)  
library(terra)
library(stringi)
library(nnls)

osmWeight <- 0.95  
GHSdiv <- 3
MSBFdiv <- 3  
cell_area  <- 100   

shp_path <- "DATA/Covasna_county.shp"
polygons <- vect(shp_path)

raster_files <- list(
  OSM  = "DATA/OSM_CV.tif",
  GHS  = "DATA/GHS_CV.tif",
  MSBF = "DATA/MSBF_CV.tif"
)
rasters <- lapply(raster_files, rast)

stat_df <- data.frame(
  natCode    = character(),
  name       = character(),
  wOSM       = numeric(),
  wGHS       = numeric(),
  wMSBF      = numeric(),
  ensArea    = numeric(),
  OSMarea    = numeric(),
  GHSarea    = numeric(),
  MSBFarea   = numeric()
)

for (i in seq_len(nrow(polygons))) {
  poly    <- polygons[i, ]
  poly_nm <- stri_trans_general(poly$name, "Latin-ASCII")
  cat(">>> Working on:", poly_nm, "\n")
  
  masked <- lapply(rasters, function(r) mask(crop(r, poly), poly))
  osm  <- masked$OSM
  ghs  <- masked$GHS
  msbf <- masked$MSBF
  
  miss <- values(ghs)>0 & values(msbf)>0 & values(osm)==0
  values(osm)[miss] <- NA
  if (sum(values(osm), na.rm=TRUE)==0) next
  
  stacked <- c(osm, ghs, msbf)
  vals    <- values(stacked)
  valid   <- complete.cases(vals)
  raw_x   <- vals[valid,1]  # OSM m²
  raw_y   <- vals[valid,2]  # GHS m²
  raw_z   <- vals[valid,3]  # MSBF m²
  
  x0 <- raw_x - mean(raw_x)
  y0 <- raw_y - mean(raw_y)
  z0 <- raw_z - mean(raw_z)
  
  etc <- function(x,y,z) {
    C12 <- cov(x,y); C13 <- cov(x,z); C23 <- cov(y,z)
    list(
      a       = c(a1=C12/C23, a2=C13/C23, a3=1),
      err.var = c(C12*C13/C23, C12*C23/C13, C13*C23/C12)
    )
  }
  et     <- etc(x0,y0,z0)
  a      <- et$a
  err.var<- et$err.var
  
  T0 <- (x0/a[1] + y0/a[2] + z0/a[3]) / 3
  b  <- c(
    mean(raw_x - a[1]*T0),
    mean(raw_y - a[2]*T0),
    mean(raw_z - a[3]*T0)
  )
  
  x_c <- (raw_x - b[1]) / a[1]
  y_c <- (raw_y - b[2]) / a[2]
  z_c <- (raw_z - b[3]) / a[3]
  
  ctc <- function(x,y,z) {
    C   <- cov(cbind(x,y,z))
    V1  <- C[1,1]; V2 <- C[2,2]; V3 <- C[3,3]
    C12 <- C[1,2];  C23 <- C[2,3]
    f   <- function(p) c(
      p[5] - C12,
      p[5] + p[4] - C23,
      p[5] + p[1] - V1,
      p[5] + p[2] - V2,
      p[5] + p[3] - V3
    )
    init <- c(0.1*V1,0.1*V2,0.1*V3,0, mean(c(V1,V2,V3)))
    sol  <- nleqslv(init, f)$x
    list(err.var=sol[1:3],
         err.cov=matrix(c(0,0,0,0,0,sol[4],0,sol[4],0),3,3))
  }
  ct  <- ctc(x_c,y_c,z_c)
  Sig <- ct$err.cov; diag(Sig) <- ct$err.var
  
  w    <- nnls(Sig, rep(1,3))$x; w    <- w/sum(w)
  w12  <- {v<-solve(Sig[1:2,1:2],rep(1,2)); v/sum(v)}
  w13  <- {v<-solve(Sig[c(1,3),c(1,3)],rep(1,2)); v/sum(v)}
  w23  <- {v<-solve(Sig[2:3,2:3],rep(1,2)); v/sum(v)}
  
  p1_r <- (!is.na(osm)) & (osm > 0)
  p2_r <- (!is.na(ghs)) & (ghs > 0)
  p3_r <- (!is.na(msbf)) & (msbf > 0)
  
  T_prop_r <- lapp(
    c(osm, ghs, msbf, p1_r, p2_r, p3_r),
    fun = function(o, g, m, p1, p2, p3) {
      o[is.na(o)] <- 0
      g[is.na(g)] <- 0
      m[is.na(m)] <- 0
      
      n <- p1 + p2 + p3
      
      w1 <- numeric(length(n))
      w2 <- numeric(length(n))
      w3 <- numeric(length(n))
      
      idx3 <- (n == 3)
      w1[idx3] <- w[1]
      w2[idx3] <- w[2]
      w3[idx3] <- w[3]
      
      idx12 <- (n == 2 & p1 & p2)
      w1[idx12] <- w12[1]
      w2[idx12] <- w12[2]
      
      w1[idx12] <- osmWeight
      w2[idx12] <- 1-osmWeight
      
      idx13 <- (n == 2 & p1 & p3)
      w1[idx13] <- w13[1]
      w3[idx13] <- w13[2]
      
      w1[idx13] <- osmWeight
      w3[idx13] <- 1-osmWeight
      
      idx23 <- (n == 2 & p2 & p3)
      w2[idx23] <- w23[1]
      w3[idx23] <- w23[2]
      
      idx1o <- (n == 1 & p1)
      idx1g <- (n == 1 & !p1 & p2)
      idx1m <- (n == 1 & !p1 & !p2 & p3)
      w1[idx1o] <- osmWeight
      w2[idx1g] <- w[2] / GHSdiv
      w3[idx1m] <- w[3] / MSBFdiv
      
      w1*o + w2*g + w3*m
    }
  )
  
  var_rast <- lapp(
    c(p1_r, p2_r, p3_r),
    fun = function(p1, p2, p3) {
      n <- p1 + p2 + p3
      w1 <- w2 <- w3 <- numeric(length(n))
      
      idx3 <- (n == 3)
      w1[idx3] <- w[1]
      w2[idx3] <- w[2]
      w3[idx3] <- w[3]
      
      idx12 <- (n == 2 & p1 & p2); w1[idx12] <- osmWeight; w2[idx12] <- 1-osmWeight
      idx13 <- (n == 2 & p1 & p3); w1[idx13] <- osmWeight; w3[idx13] <- 1-osmWeight
      
      idx23 <- (n == 2 & p2 & p3); w2[idx23] <- w23[1]; w3[idx23] <- w23[2]
      
      idx1o <- (n == 1 & p1);      w1[idx1o] <- osmWeight
      idx1g <- (n == 1 & !p1 & p2); w2[idx1g] <- w[2] / GHSdiv
      idx1m <- (n == 1 & !p1 & !p2 & p3); w3[idx1m] <- w[3] / MSBFdiv
      
      w1^2 * Sig[1,1] +
        w2^2 * Sig[2,2] +
        w3^2 * Sig[3,3] +
        2 * w2 * w3 * Sig[2,3]
    }
  )
  sd_rast <- sqrt(var_rast)
  T_lo_r  <- T_prop_r - 1.96 * sd_rast
  T_hi_r  <- T_prop_r + 1.96 * sd_rast
  
  zero_mask <- (values(osm)==0 & values(ghs)==0 & values(msbf)==0)
  values(T_prop_r)[zero_mask] <- 0
  values(T_lo_r)[zero_mask]   <- 0
  values(T_hi_r)[zero_mask]   <- 0
  
  base <- paste0("ENS/", poly$natCode, "_", poly_nm)
  writeRaster(T_prop_r, paste0(base, "_est.tif"),   overwrite=TRUE)
  writeRaster(T_lo_r,   paste0(base, "_ci_lo.tif"), overwrite=TRUE)
  writeRaster(T_hi_r,   paste0(base, "_ci_hi.tif"), overwrite=TRUE)
  
  ci_width <- T_hi_r - T_lo_r
  values(ci_width)[zero_mask] <- 0
  writeRaster(ci_width,
              paste0("ENS/", poly$natCode, "_", poly_nm, "_ci_width.tif"),
              overwrite = TRUE)
  
  v_osm  <- values(osm);  v_ghs <- values(ghs);  v_msbf <- values(msbf)
  v_ens  <- values(T_prop_r)
  
  stat_df <- rbind(
    stat_df,
    data.frame(
      natCode    = poly$natCode,
      name       = poly$name,
      wOSM       = w[1],
      wGHS       = w[2],
      wMSBF      = w[3],
      ensArea    = sum(v_ens, na.rm=TRUE),
      OSMarea    = sum(v_osm,  na.rm=TRUE),
      GHSarea    = sum(v_ghs,  na.rm=TRUE),
      MSBFarea   = sum(v_msbf, na.rm=TRUE)
    )
  )
}

write.csv(stat_df, "stat.csv", row.names=FALSE)
