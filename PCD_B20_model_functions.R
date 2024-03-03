
pcd_model <- function(x){
  
  # Load data 
  
  bbp   = x$bbp_med 
  chl   = x$chl_cor
  irr   = x$par_8day
  mld= x$mld
  kd490 <- 0.0166+0.07242*(chl^0.68955)
  kdpar <- ifelse(mld <= 1/kd490, 0.0864 + 0.884*kd490 - 0.00137/kd490,
                  0.0665 + 0.874*kd490 - 0.00121/kd490)
  carbon_od = bbp*43255.75  
  
  # Define PCD and Ez 
  
  PCD <- 0.427 

  Ez_0.1pc <-   irr/100*0.1
  
  decay_depth  =   round_any((-log(PCD /irr)/kdpar),1)
  
  if (decay_depth < mld){
    decay_depth = mld}
  
  if (decay_depth > 500){
    decay_depth = mld}
  
  # Create vector for depth-resolved carbon
  
  z    = c(0:1000)
  Carbon.z    <- rep(carbon_od, 1001)
  
  # Calculate POC through depth

  for (k in 1:1001){Carbon.z[k] <-  ifelse(k > mld, 
                                           bbp * (12100+(43255.75-12100)*exp(-0.00657*(k-mld))) ,
                                           bbp * (43255.75*(exp(-0.00657*(0)))))
  } 
  
  
  # Apply algorithm to overwrite data below the PCD
  
    for (k in decay_depth:1001){
    
    decay_depth_c <-  bbp * (12100+(43255.75-12100)*exp(-0.00657*(decay_depth-mld)))
    if (decay_depth_c > bbp*(43255.75)){
      decay_depth_c = bbp*(43255.75)}

    depth_norm = k-decay_depth
    
    A  =  decay_depth_c
    B = 450.29093 * exp(-(0.0708*A))+19
    C = 0.325 
    
    carbon.a =    A * exp(-((depth_norm) / B)^C)
    
    Carbon.z[k] <- round(carbon.a,1)

    print(k)
    
  }
  
  return(Carbon.z)

}



B20_ez <- function(x){
  
  # Load data 
  
  bbp   = x$bbp_med 
  chl   = x$chl_cor
  irr   = x$par_8day
  mld= x$mld
  kd490 <- 0.0166+0.07242*(chl^0.68955)
  kdpar <- ifelse(mld <= 1/kd490, 0.0864 + 0.884*kd490 - 0.00137/kd490,
                  0.0665 + 0.874*kd490 - 0.00121/kd490)
  kdpar_deep  = 0.0665 + 0.874*kd490 - 0.00121/kd490 # equn 9'

  # Define Ez 
  
  Ez_0.1pc <-   irr/100*0.1
 
   decay_depth  =   round_any((-log(Ez_0.1pc /irr)/kdpar),1)
  
  if (decay_depth < mld){
    decay_depth = mld}
  if (decay_depth > 500){
    decay_depth = mld}
  
   # Create vector for depth-resolved carbon
   
  z    = c(0:1000)
  
  Carbon.z    <- rep(carbon_od, 1001)
  
  # Calculate POC through depth
  
  
  for (k in 1:1001){Carbon.z[k] <-  ifelse(k > mld, 
                                           bbp * (12100+(43255.75-12100)*exp(-0.00657*(k-mld))) ,
                                           bbp * (43255.75*(exp(-0.00657*(0)))))
  } 

   # Apply algorithm to overwrite data below the PCD
  
  for (k in decay_depth:1001){
    
    decay_depth_c <-  bbp * (12100+(43255.75-12100)*exp(-0.00657*(decay_depth-mld)))
    if (decay_depth_c > bbp*(43255.75)){
      decay_depth_c = bbp*(43255.75)}
    
    carbon.a =    decay_depth_c * (k/decay_depth)^(-0.858)
    
     Carbon.z[k] <- round(carbon.a,1)
    
    print(k)
    
  }
  
  return(Carbon.z)

}

