

PPM_OSP <- function(bbp, chl, irr, k490, MLD,cafe, doy,year){
  
  # Set bounding box for OSP region
  
  
  osp <- extent(-150.25, -149.75,49.75,50.25)

  # Remove negative fill values
  
  bbp[bbp<0]   <- NA  
  chl[chl<0]   <- NA
  irr[irr<0]   <- NA
  k490[k490<0] <- NA
  MLD[MLD<0] <- NA
  cafe[cafe<0] <- NA
  MLD[MLD<10] = 10
  
  # Set bounding box for OSP region
  
  bbp <- bbp * (443/700)

  # Clip input parameters 
  
  chl <- raster(t(chl), xmn=-180, xmx=180, ymn=-90, ymx=90, crs="+proj=longlat")
  chl  <- crop(chl, osp)
  chl <- as.matrix(chl)
  
  k490 <- raster(t(k490), xmn=-180, xmx=180, ymn=-90, ymx=90, crs="+proj=longlat")
  k490  <- crop(k490, osp)
  k490 <- as.matrix(k490)
  
  irr <- raster(t(irr), xmn=-180, xmx=180, ymn=-90, ymx=90, crs="+proj=longlat")
  irr  <- crop(irr, osp)
  irr <- as.matrix(irr)
  
  MLD <- raster(t(MLD), xmn=-180, xmx=180, ymn=-90, ymx=90, crs="+proj=longlat")
  MLD  <- crop(MLD, osp)
  MLD <- as.matrix(MLD)
  
  bbp <- raster(t(bbp), xmn=-180, xmx=180, ymn=-90, ymx=90, crs="+proj=longlat")
  bbp  <- crop(bbp, osp)
  bbp <- as.matrix(bbp)
  
  cafe <- raster(t(cafe), xmn=-180, xmx=180, ymn=-90, ymx=90, crs="+proj=longlat")
  cafe  <- crop(cafe, osp)
  cafe <- as.matrix(cafe)
  
  ocean <- which(is.finite(bbp) & is.finite(k490) & is.finite(chl) 
                 &  is.finite(irr) & is.finite(MLD) & is.finite(cafe)
                 ,arr.ind = T)
  
  bbp   <- bbp[ocean]  
  cafe   <- cafe[ocean]  
  chl   <- chl[ocean]
  irr   <- irr[ocean]
  k490  <- k490[ocean]
  MLD   <- MLD[ocean]
  
  # Estimate kd 
  
  kdpar <- ifelse(MLD <= 1/k490, 0.0864 + 0.884*k490 - 0.00137/k490,
                  0.0665 + 0.874*k490 - 0.00121/k490)
  
  # Create depth resolved vector, calculate PCD and Ez 
  
  z    = c(0:1000)
  Carbon.z     <- array(NA, c(length(chl), 1001))
  
  r <- 0.427
  
  decay_depth  =   round_any((-log(r /(irr))/kdpar),1)
  decay_depth = ifelse(decay_depth < MLD, MLD,
                       decay_depth)
  decay_depth = ifelse(decay_depth > 500, MLD,
                       decay_depth)
  
  Ez_0.1 <-   (irr)/100*0.1
  Ez_0.1 <- -log(Ez_0.1/(irr))/kdpar
  Ez <- -log(0.01/(irr))/kdpar
  Ez[Ez>300] <- NA
  zeu <- 4.6/kdpar
  
  # Calculate POC for discrete depths 
  
  surface_poc <- bbp*(43255.75)
  decay_depth_c <-  bbp * (12100+(43255.75-12100)*exp(-0.00657*(decay_depth-MLD)))
  decay_depth_c[decay_depth_c > bbp*(43255.75)] = bbp*(43255.75)
  
  ez_depth_norm <- Ez-decay_depth
  ez_depth_norm[ez_depth_norm < 0] <- 0
  ez_depth_norm_100 <- ez_depth_norm + 100
  
  tz_depth_norm <- 1000-decay_depth
  tz_depth_norm[tz_depth_norm < 0] <- 0
  
  A  =  decay_depth_c
  B = 450.29093 * exp(-(0.0708*A))+19
  C = 0.325 # for non blank corrected profiles
  
  tz_depth_c_martin <- A * ((1000)/decay_depth)^-0.858
  tz_depth_c <- A * exp(-((tz_depth_norm) / B)^C)
  ez_depth_c <- A * exp(-((ez_depth_norm) / B)^C)
  ez_depth_c_100 <- A * exp(-((ez_depth_norm_100) / B)^C)
  decay_depth_c_100 <- A * exp(-((100) / B)^C)
  decay_depth_c_100_martin <- A * ((decay_depth+100)/decay_depth)^-0.858
  Ez_0.1 = ifelse(Ez_0.1 < MLD, MLD,
                  Ez_0.1)
  
 
  # Calculate POC through depth 
  
  Carbon     <- array(NA, c(length(MLD), 1001))
  
  for (k in 1:1001){Carbon[,k] <-  ifelse(k > MLD, 
                                          bbp * (12100+(43255.75-12100)*exp(-0.00657*(k-MLD))) ,
                                          bbp * (43255.75*(exp(-0.00657*(0)))))
  } 
  
  for (k in 1:1001){
    
    deep <- which(decay_depth<=z[k], arr.ind = T)
    
    A[deep] = decay_depth_c[deep]
    B[deep]  = 19  + 450.29093 * exp(-(0.0708*A[deep]))
    C = 0.325
    
    Carbon[deep,k] =   A[deep] * exp(-((z[k]- decay_depth[deep]) / B[deep])^C)
    
      print(k)
    
  }
  
  Carbon_martin     <- array(NA, c(length(MLD), 1001))
  for (k in 1:1001){Carbon_martin[,k] <-  ifelse(k > MLD, 
                                                 bbp * (12100+(43255.75-12100)*exp(-0.00657*(k-MLD))) ,
                                                 bbp * (43255.75*(exp(-0.00657*(0)))))
  } 
  
  for (k in 1:1001){
    
    deep <- which(Ez<=z[k], arr.ind = T)
    
    

    Carbon_martin[deep,k] =    A[deep] * (z[k]/decay_depth[deep])^-0.858
    
    print(k)
    
  }
  
  
  # Create DF and save file 
  
  Carbon_df <- as.data.frame(t(Carbon))
  new_colnames <- gsub("V", "Pix", colnames(Carbon_df))
  Carbon_martin_df <- as.data.frame(t(Carbon_martin))
  new_colnames_martin <- gsub("V", "Pix_Martin", colnames(Carbon_martin_df))
  colnames(Carbon_df) <- new_colnames
  colnames(Carbon_martin_df) <- new_colnames_martin
  Carbon_df <- cbind(Carbon_df,Carbon_martin_df)
  Carbon_df$doy <- doy
  Carbon_df$date <-  as.Date(Carbon_df$doy, format = "%j", origin = paste("1.1.",yr,sep = ""))
  day <- mday(Carbon_df$date) 
  Carbon_df$month <- month(Carbon_df$date) 
  Carbon_df$date <- paste(day,"01",yr, sep = "/")
  Carbon_df$year <- yr
  Carbon_df$site <- "osp"
  Carbon_df$lat <- 50
  Carbon_df$lon <- -150
  Carbon_df$mld <- mean(MLD)
  Carbon_df$surface_poc <- mean(surface_poc)
  Carbon_df$Ez <- mean(Ez)
  Carbon_df$Ez_C <- mean(ez_depth_c)
  Carbon_df$DP <- mean(decay_depth)
  Carbon_df$DP_C <- mean(decay_depth_c)
  Carbon_df$DP_C_100 <- mean(decay_depth_c_100)
  Carbon_df$Ez_C_100 <- mean(ez_depth_c_100)
  Carbon_df$tz_depth_c <- mean(tz_depth_c)
  Carbon_df$cafe <- mean(cafe)
  
  write.csv(Carbon_df,paste("c:/Users/foxja/Documents/OSU/Active projects/GLIDER_PP_SH/DATA/processed_files/sat_poc_estimates/OSP_iv/8day/AQUA_MODIS.",year,doy,".8day.poc.9km.csv", sep=""),row.names = F )
  
}


PPM_PAP <- function(bbp, chl, irr, k490, MLD,cafe, doy,year){
  
  #Remove negative fill values
  pap <- extent(-16.75,-16.25,48.75,49.25)
  
  bbp[bbp<0]   <- NA  
  chl[chl<0]   <- NA
  irr[irr<0]   <- NA
  k490[k490<0] <- NA
  MLD[MLD<0] <- NA
  MLD[MLD<10] = 10
  
  bbp <- bbp * (443/700)
  #Stretch MLD and ZNOx
  
  # this can be used if the dimensions of the MLD data grid do not match. Check to see whether
  # the MLD data are compatible with 4km or 9km sat data
  POC_surf <- bbp*40800
  
  chl <- raster(t(chl), xmn=-180, xmx=180, ymn=-90, ymx=90, crs="+proj=longlat")
  chl  <- crop(chl, pap)
  chl <- as.matrix(chl)
  
  k490 <- raster(t(k490), xmn=-180, xmx=180, ymn=-90, ymx=90, crs="+proj=longlat")
  k490  <- crop(k490, pap)
  k490 <- as.matrix(k490)
  
  irr <- raster(t(irr), xmn=-180, xmx=180, ymn=-90, ymx=90, crs="+proj=longlat")
  irr  <- crop(irr, pap)
  irr <- as.matrix(irr)
  
  MLD <- raster(t(MLD), xmn=-180, xmx=180, ymn=-90, ymx=90, crs="+proj=longlat")
  MLD  <- crop(MLD, pap)
  MLD <- as.matrix(MLD)
  
  bbp <- raster(t(bbp), xmn=-180, xmx=180, ymn=-90, ymx=90, crs="+proj=longlat")
  bbp  <- crop(bbp, pap)
  bbp <- as.matrix(bbp)
  
  cafe <- raster(t(cafe), xmn=-180, xmx=180, ymn=-90, ymx=90, crs="+proj=longlat")
  cafe  <- crop(cafe, pap)
  cafe <- as.matrix(cafe)
  
  ocean <- which(is.finite(bbp) & is.finite(k490) & is.finite(chl) 
                 &  is.finite(irr) & is.finite(MLD) & is.finite(cafe)
                 ,arr.ind = T)
  
  bbp   <- bbp[ocean]  
  cafe   <- cafe[ocean]  
  chl   <- chl[ocean]
  irr   <- irr[ocean]
  k490  <- k490[ocean]
  MLD   <- MLD[ocean]
  
  # Estimate kd 
  
  kdpar <- ifelse(MLD <= 1/k490, 0.0864 + 0.884*k490 - 0.00137/k490,
                  0.0665 + 0.874*k490 - 0.00121/k490)
  
  # Create depth resolved vector, calculate PCD and Ez 
  
  z    = c(0:1000)
  Carbon.z     <- array(NA, c(length(chl), 1001))
  
  r <- 0.427
  
  decay_depth  =   round_any((-log(r /(irr))/kdpar),1)
  decay_depth = ifelse(decay_depth < MLD, MLD,
                       decay_depth)
  decay_depth = ifelse(decay_depth > 500, MLD,
                       decay_depth)
  
  Ez_0.1 <-   (irr)/100*0.1
  Ez_0.1 <- -log(Ez_0.1/(irr))/kdpar
  Ez <- -log(0.01/(irr))/kdpar
  Ez[Ez>300] <- NA
  zeu <- 4.6/kdpar
  
  # Calculate POC for discrete depths 
  
  surface_poc <- bbp*(43255.75)
  decay_depth_c <-  bbp * (12100+(43255.75-12100)*exp(-0.00657*(decay_depth-MLD)))
  decay_depth_c[decay_depth_c > bbp*(43255.75)] = bbp*(43255.75)
  
  ez_depth_norm <- Ez-decay_depth
  ez_depth_norm[ez_depth_norm < 0] <- 0
  ez_depth_norm_100 <- ez_depth_norm + 100
  
  tz_depth_norm <- 1000-decay_depth
  tz_depth_norm[tz_depth_norm < 0] <- 0
  
  A  =  decay_depth_c
  B = 450.29093 * exp(-(0.0708*A))+19
  C = 0.325 # for non blank corrected profiles
  
  tz_depth_c_martin <- A * ((1000)/decay_depth)^-0.858
  tz_depth_c <- A * exp(-((tz_depth_norm) / B)^C)
  ez_depth_c <- A * exp(-((ez_depth_norm) / B)^C)
  ez_depth_c_100 <- A * exp(-((ez_depth_norm_100) / B)^C)
  decay_depth_c_100 <- A * exp(-((100) / B)^C)
  decay_depth_c_100_martin <- A * ((decay_depth+100)/decay_depth)^-0.858
  Ez_0.1 = ifelse(Ez_0.1 < MLD, MLD,
                  Ez_0.1)
  
  
  # Calculate POC through depth 
  
  Carbon     <- array(NA, c(length(MLD), 1001))
  
  for (k in 1:1001){Carbon[,k] <-  ifelse(k > MLD, 
                                          bbp * (12100+(43255.75-12100)*exp(-0.00657*(k-MLD))) ,
                                          bbp * (43255.75*(exp(-0.00657*(0)))))
  } 
  
  for (k in 1:1001){
    
    deep <- which(decay_depth<=z[k], arr.ind = T)
    
    A[deep] = decay_depth_c[deep]
    B[deep]  = 19  + 450.29093 * exp(-(0.0708*A[deep]))
    C = 0.325
    
    Carbon[deep,k] =   A[deep] * exp(-((z[k]- decay_depth[deep]) / B[deep])^C)
    
    print(k)
    
  }
  
  Carbon_martin     <- array(NA, c(length(MLD), 1001))
  for (k in 1:1001){Carbon_martin[,k] <-  ifelse(k > MLD, 
                                                 bbp * (12100+(43255.75-12100)*exp(-0.00657*(k-MLD))) ,
                                                 bbp * (43255.75*(exp(-0.00657*(0)))))
  } 
  
  for (k in 1:1001){
    
    deep <- which(Ez<=z[k], arr.ind = T)
    
        Carbon_martin[deep,k] =    A[deep] * (z[k]/decay_depth[deep])^-0.858
   
    print(k)
    
  }
  
  
  # Create DF and save file 
  
  Carbon_df <- as.data.frame(t(Carbon))
  new_colnames <- gsub("V", "Pix", colnames(Carbon_df))
  Carbon_martin_df <- as.data.frame(t(Carbon_martin))
  new_colnames_martin <- gsub("V", "mart_p", colnames(Carbon_martin_df))
  colnames(Carbon_df) <- new_colnames
  colnames(Carbon_martin_df) <- new_colnames_martin
  Carbon_df <- cbind(Carbon_df,Carbon_martin_df)
  Carbon_df$doy <- doy
  Carbon_df$date <-  as.Date(Carbon_df$doy, format = "%j", origin = paste("1.1.",yr,sep = ""))
  day <- mday(Carbon_df$date) 
  Carbon_df$month <- month(Carbon_df$date) 
  Carbon_df$date <- paste(day,"01",yr, sep = "/")
  Carbon_df$year <- yr
  Carbon_df$site <- "pap"
  Carbon_df$lat <- 49
  Carbon_df$lon <- -16
  Carbon_df$mld <- mean(MLD)
  Carbon_df$surface_poc <- mean(surface_poc)
  Carbon_df$Ez <- mean(Ez)
  Carbon_df$Ez_C <- mean(ez_depth_c)
  Carbon_df$DP <- mean(decay_depth)
  Carbon_df$DP_C <- mean(decay_depth_c)
  Carbon_df$DP_C_100 <- mean(decay_depth_c_100)
  Carbon_df$Ez_C_100 <- mean(ez_depth_c_100)
  Carbon_df$DP_C_MART <- mean(tz_depth_c_martin)
  Carbon_df$DP_C_100_MART <- mean(decay_depth_c_100_martin)
  Carbon_df$tz_depth_c <- mean(tz_depth_c)
  Carbon_df$cafe <- mean(cafe)

  write.csv(Carbon_df,paste("c:/Users/foxja/Documents/OSU/Active projects/GLIDER_PP_SH/DATA/processed_files/sat_poc_estimates/PAP_viii/AQUA_MODIS.", year,doy,".8day.poc.9km.csv", sep=""),row.names = F )

}





