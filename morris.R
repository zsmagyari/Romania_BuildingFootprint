library(terra)
library(nnls)
library(nleqslv)
library(sensitivity)

shp_path <- "DATA/Covasna_county.shp"
raster_files <- list(
  OSM  = "DATA/OSM_CV.tif",
  GHS  = "DATA/GHS_CV.tif",
  MSBF = "DATA/MSBF_CV.tif"
)

polygons <- vect(shp_path)
rasters <- lapply(raster_files, rast)

classify_nonlinearity <- function(r) {
  if (is.na(r)) return(NA_character_)
  if (r < 0.1) return("linear")
  if (r < 0.5) return("moderate")
  return("nonlinear")
}

results <- vector("list", length = nrow(polygons))

morris_r <- 30
factor_names <- c("osmWeight", "GHSdiv", "MSBFdiv")

for (i in seq_len(nrow(polygons))) {
  poly <- polygons[i, ]
  natCode <- poly$natCode
  name <- poly$name
  cat(sprintf("[%d/%d] Processing: %s (%s)\n", i, nrow(polygons), natCode, name))
  
  out_row <- list(
    natCode = natCode,
    name = name,
    mu_osmWeight = NA_real_,
    mu_star_osmWeight = NA_real_,
    sigma_osmWeight = NA_real_,
    mu_GHSdiv = NA_real_,
    mu_star_GHSdiv = NA_real_,
    sigma_GHSdiv = NA_real_,
    mu_MSBFdiv = NA_real_,
    mu_star_MSBFdiv = NA_real_,
    sigma_MSBFdiv = NA_real_,
    dominant_factor = NA_character_,
    dominance_strength = NA_real_,
    dominance_relative = NA_real_,
    dominance_rank_osmWeight = NA_real_,
    dominance_rank_GHSdiv = NA_real_,
    dominance_rank_MSBFdiv = NA_real_,
    linearity_osmWeight = NA_character_,
    linearity_GHSdiv = NA_character_,
    linearity_MSBFdiv = NA_character_,
    max_interaction_ratio = NA_real_,
    insufficient_data = FALSE
  )
  
  tryCatch({
    masked <- lapply(rasters, function(r) mask(crop(r, poly), poly))
    osm  <- masked$OSM
    ghs  <- masked$GHS
    msbf <- masked$MSBF
    
    miss <- values(ghs) > 0 & values(msbf) > 0 & values(osm) == 0
    values(osm)[miss] <- NA
    
    if ((all(is.na(values(osm))) || sum(values(osm), na.rm=TRUE)==0) &&
        (all(is.na(values(ghs))) || sum(values(ghs), na.rm=TRUE)==0) &&
        (all(is.na(values(msbf))) || sum(values(msbf), na.rm=TRUE)==0)) {
      out_row$insufficient_data <- TRUE
      results[[i]] <- out_row
      next
    }
    
    stacked <- c(osm, ghs, msbf)
    vals    <- values(stacked)
    valid   <- complete.cases(vals)
    raw_x   <- vals[valid,1]
    raw_y   <- vals[valid,2]
    raw_z   <- vals[valid,3]
    
    if (length(raw_x) < 5) {
      out_row$insufficient_data <- TRUE
      results[[i]] <- out_row
      next
    }
    
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
    et      <- etc(x0,y0,z0)
    a       <- et$a
    err.var <- et$err.var
    
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
    
    o_vals <- values(osm)
    g_vals <- values(ghs)
    m_vals <- values(msbf)
    p1 <- (!is.na(o_vals)) & (o_vals > 0)
    p2 <- (!is.na(g_vals)) & (g_vals > 0)
    p3 <- (!is.na(m_vals)) & (m_vals > 0)
    
    ensemble_area <- function(osmWeight, GHSsingleSourceDivisor, MSBFsingleSourceDivisor) {
      o <- o_vals; g <- g_vals; m <- m_vals
      o[is.na(o)] <- 0; g[is.na(g)] <- 0; m[is.na(m)] <- 0
      
      n <- as.numeric(p1) + as.numeric(p2) + as.numeric(p3)
      w1 <- numeric(length(n))
      w2 <- numeric(length(n))
      w3 <- numeric(length(n))
      
      idx3 <- (n == 3)
      w1[idx3] <- w[1]
      w2[idx3] <- w[2]
      w3[idx3] <- w[3]
      
      idx12 <- (n == 2 & p1 & p2)
      w1[idx12] <- osmWeight
      w2[idx12] <- 1 - osmWeight
      
      idx13 <- (n == 2 & p1 & p3)
      w1[idx13] <- osmWeight
      w3[idx13] <- 1 - osmWeight
      
      idx23 <- (n == 2 & p2 & p3)
      w2[idx23] <- w23[1]
      w3[idx23] <- w23[2]
      
      # egy forrás
      idx1o <- (n == 1 & p1)
      idx1g <- (n == 1 & !p1 & p2)
      idx1m <- (n == 1 & !p1 & !p2 & p3)
      w1[idx1o] <- osmWeight
      w2[idx1g] <- w[2] / GHSsingleSourceDivisor
      w3[idx1m] <- w[3] / MSBFsingleSourceDivisor
      
      cell_vals <- w1 * o + w2 * g + w3 * m
      zero_mask <- (o == 0 & g == 0 & m == 0)
      cell_vals[zero_mask] <- 0
      sum(cell_vals, na.rm = TRUE)
    }
    
    model_wrapper <- function(U) {
      apply(U, 1, function(u) {
        osmW <- u[1]
        GHSsingleDiv <- 1 + 4 * u[2]
        MSBFsingleDiv <- 1 + 4 * u[3]
        ensemble_area(osmW, GHSsingleDiv, MSBFsingleDiv)
      })
    }
    
    set.seed(10000 + i)
    morris_obj <- morris(
      model = NULL,
      factors = 3,
      r = morris_r,
      design = list(type = "oat", levels = 10, grid.jump = 1),
      factor.names = factor_names
    )
    Y_morris <- model_wrapper(morris_obj$X)
    Y_morris <- as.numeric(Y_morris)
    tell(morris_obj, Y_morris)
    
    ee <- morris_obj$ee  
    
    mu      <- apply(ee, 2, mean)
    mu_star <- apply(abs(ee), 2, mean)
    sigma   <- apply(ee, 2, sd)
    
    out_row$mu_osmWeight      <- mu[1]
    out_row$mu_star_osmWeight <- mu_star[1]
    out_row$sigma_osmWeight   <- sigma[1]
    
    out_row$mu_GHSdiv         <- mu[2]
    out_row$mu_star_GHSdiv    <- mu_star[2]
    out_row$sigma_GHSdiv      <- sigma[2]
    
    out_row$mu_MSBFdiv        <- mu[3]
    out_row$mu_star_MSBFdiv   <- mu_star[3]
    out_row$sigma_MSBFdiv     <- sigma[3]
    
    mu_star_vec <- c(mu_star[1], mu_star[2], mu_star[3])
    names(mu_star_vec) <- factor_names
    total_mu_star <- sum(mu_star_vec, na.rm = TRUE)
    if (total_mu_star > 0) {
      dominant <- names(mu_star_vec)[which.max(mu_star_vec)]
      out_row$dominant_factor <- dominant
      out_row$dominance_strength <- mu_star_vec[dominant]
      out_row$dominance_relative <- mu_star_vec[dominant] / total_mu_star
    }
    ranks <- rank(-mu_star_vec, ties.method = "min")
    out_row$dominance_rank_osmWeight <- ranks["osmWeight"]
    out_row$dominance_rank_GHSdiv    <- ranks["GHSdiv"]
    out_row$dominance_rank_MSBFdiv   <- ranks["MSBFdiv"]
    
    ratio <- rep(NA_real_, 3)
    for (j in seq_len(3)) {
      if (mu_star[j] != 0) ratio[j] <- sigma[j] / abs(mu_star[j]) else ratio[j] <- NA
    }
    out_row$linearity_osmWeight <- classify_nonlinearity(ratio[1])
    out_row$linearity_GHSdiv    <- classify_nonlinearity(ratio[2])
    out_row$linearity_MSBFdiv   <- classify_nonlinearity(ratio[3])
    out_row$max_interaction_ratio <- max(ratio, na.rm = TRUE)
    
    results[[i]] <- out_row
    
  }, error = function(e) {
    warning(sprintf("Error on %s (%s): %s", natCode, name, e$message))
    out_row$insufficient_data <- TRUE
    results[[i]] <- out_row
  })
}

df_full <- do.call(rbind, lapply(results, as.data.frame))
write.csv(df_full, "morris_stat.csv", row.names = FALSE)

