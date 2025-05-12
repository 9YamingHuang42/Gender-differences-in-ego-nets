library( dagitty )
library( dplyr )
library( rethinking )
library( tibble )
library( gt )
library( gto )
library( officer )
library( ggplot2 )
library( ggdist )
options( digits = 8 )

# Load data ----
load( "Data_gender-difference_in_ego_nets.RData")

#1 DAG ----
Dag <- 
  dagitty( 'dag {
bb="-0.5,-0.5,0.5,0.5"
Age [pos="0.000,-0.350"]
Density [outcome,pos="0.000,0.320"]
Gender [pos="0.160,-0.250"]
Market [pos="0.210,0.050"]
Residence [pos="-0.160,-0.250"]
U_Individual [latent,pos="-0.210,0.150"]
U_Village [latent,pos="-0.280,-0.050"]
Age -> Density
Age -> Market
Age -> Residence
Gender -> Density
Gender -> Market
Market -> Density
Residence -> Density
Residence -> Market
U_Individual -> Density
U_Village -> Density
}
')

drawdag( Dag , radius = 8.8 )

#2 Data preparation ----
# Standardise age 
# Natural logarithmic transformation and standardisation of market
Ego.data <- Ego.data %>% 
  mutate( Age_std = ( Age - mean( Age ) )/sd( Age ) ,
          Market_log = log( Market + 1 ) , 
          Market_log_std = ( Market_log - mean( Market_log ) )/sd( Market_log ) )

#3 Models ----
##3.1 Market participation ----
# Other variable should be included into models
adjustmentSets( Dag , 
                exposure = "Market" , 
                outcome = "Density" ,
                effect = "total" )

###3.1.1 Models of all density ----
Model_GM_all_list <- with( Ego.data , list(
  Possible_tie =  as.integer( Possible_tie ) , 
  Tie_all = as.integer( Tie_all ) , 
  Age = as.numeric( sprintf( "%0.4f" , Age_std ) ) ,
  Gender = as.integer( Gender ) , 
  PMR = as.integer( Residence ) , 
  Market = as.numeric( sprintf("%0.4f" , Market_log_std ) ) ,
  ID = as.integer( ID ) ,
  VID = as.integer( VID ) ) )

{set.seed(123)
  Model_GM_all <- ulam(
    alist(
      Tie_all ~ binomial( Possible_tie , p ) ,
      logit( p ) <- bGA[Gender] * Age + 
        bGR[Gender,PMR] + 
        bM[Gender] * Market +
        V[VID] + # village intercepts 
        I[ID] , # individual intercepts 
      
      # define effects using other parameters
      transpars> vector[VID]: V <<- V_bar + z_V*sigma_V,
      transpars> vector[ID]: I <<- z_I*sigma_I,
      
      bGA[Gender] ~ normal( 0 , 0.5 ) ,
      
      matrix[Gender,PMR]: bGR ~ normal( 0 , 0.5 ) ,
      
      bM[Gender] ~ normal( 0 , 0.5 ) ,
      
      V_bar ~ normal( 0 , 0.5 ) ,
      z_V[VID] ~ normal( 0 , 1 ) ,
      sigma_V ~ exponential( 1 ) ,
      z_I[ID] ~ normal( 0 , 4 ) ,
      sigma_I ~ exponential( 1 )  
    ) , data = Model_GM_all_list , 
    iter = 1500 , warmup = 500 , chains = 6 , cores = 6 , 
    log_lik = TRUE , control = list(adapt_delta = 0.99)
  )}

Pre_model_GM_all <- precis( Model_GM_all , depth = 3 , prob = 0.90 , 
                            pars = c( "bGA" , "bGR" , "bM" ) )
Pre_model_GM_all

###3.1.2 Models of biological kin density ----
Model_GM_genetic_list <- with( Ego.data , list(
  Possible_tie =  as.integer( Possible_tie ) , 
  Tie_genetic = as.integer( Tie_genetic ) ,
  Age = as.numeric( sprintf( "%0.4f" , Age_std ) ) ,
  Gender = as.integer( Gender ) , 
  PMR = as.integer( Residence ) , 
  Market = as.numeric( sprintf("%0.4f" , Market_log_std ) ) ,
  ID = as.integer( ID ) ,
  VID = as.integer( VID ) ) )

{set.seed(123)
  Model_GM_genetic <- ulam(
    alist(
      Tie_genetic ~ binomial( Possible_tie , p ) ,
      logit( p ) <- bGA[Gender] * Age + 
        bGR[Gender,PMR] + 
        bM[Gender] * Market +
        V[VID] + # village intercepts 
        I[ID] , # individual intercepts 
      
      # define effects using other parameters
      transpars> vector[VID]: V <<- V_bar + z_V*sigma_V,
      transpars> vector[ID]: I <<- z_I*sigma_I,
      
      bGA[Gender] ~ normal( 0 , 0.5 ) ,
      
      matrix[Gender,PMR]: bGR ~ normal( 0 , 0.5 ) ,
      
      bM[Gender] ~ normal( 0 , 0.5 ) ,
      
      V_bar ~ normal( 0 , 0.5 ) ,
      z_V[VID] ~ normal( 0 , 1 ) ,
      sigma_V ~ exponential( 1 ) ,
      z_I[ID] ~ normal( 0 , 4 ) ,
      sigma_I ~ exponential( 1 ) 
    ) , data = Model_GM_genetic_list , 
    iter = 1500 , warmup = 500 , chains = 6 , cores = 6 , 
    log_lik = TRUE , control = list(adapt_delta = 0.99)
  )}
Pre_model_GM_genetic <- precis( Model_GM_genetic , depth = 3 , prob = 0.90 , 
                                pars = c( "bGA" , "bGR" , "bM" ) )
Pre_model_GM_genetic

###3.1.3 Models of affinal kin density ----
Model_GM_affine_list <- with( Ego.data , list(
  Possible_tie =  as.integer( Possible_tie ) , 
  Tie_affine = as.integer( Tie_affine ) ,
  Age = as.numeric( sprintf( "%0.4f" , Age_std ) ) ,
  Gender = as.integer( Gender ) , 
  PMR = as.integer( Residence ) , 
  Market = as.numeric( sprintf("%0.4f" , Market_log_std ) ) ,
  ID = as.integer( ID ) ,
  VID = as.integer( VID ) ) )

{set.seed(123)
  Model_GM_affine <- ulam(
    alist(
      Tie_affine ~ binomial( Possible_tie , p ) ,
      logit( p ) <- bGA[Gender] * Age + 
        bGR[Gender,PMR] + 
        bM[Gender] * Market +
        V[VID] + # village intercepts 
        I[ID] , # individual intercepts 
      
      # define effects using other parameters
      transpars> vector[VID]: V <<- V_bar + z_V*sigma_V,
      transpars> vector[ID]: I <<- z_I*sigma_I,
      
      bGA[Gender] ~ normal( 0 , 0.5 ) ,
      
      matrix[Gender,PMR]: bGR ~ normal( 0 , 0.5 ) ,
      
      bM[Gender] ~ normal( 0 , 0.5 ) ,
      
      V_bar ~ normal( 0 , 0.5 ) ,
      z_V[VID] ~ normal( 0 , 1 ) ,
      sigma_V ~ exponential( 1 ) ,
      z_I[ID] ~ normal( 0 , 4 ) ,
      sigma_I ~ exponential( 1 )
    ) , data = Model_GM_affine_list , 
    iter = 1500 , warmup = 500 , chains = 6 , cores = 6 , 
    log_lik = TRUE , control = list(adapt_delta = 0.99)
  )}
Pre_model_GM_affine <- precis( Model_GM_affine , depth = 3 , prob = 0.90 , 
                               pars = c( "bGA" , "bGR" , "bM" ) )
Pre_model_GM_affine

###3.1.4 Models of friend density ----
Model_GM_fri_list <- with( Ego.data , list(
  Possible_tie =  as.integer( Possible_tie ) , 
  Tie_fri = as.integer( Tie_fri ) ,
  Age = as.numeric( sprintf( "%0.4f" , Age_std ) ) ,
  Gender = as.integer( Gender ) , 
  PMR = as.integer( Residence ) , 
  Market = as.numeric( sprintf("%0.4f" , Market_log_std ) ) ,
  ID = as.integer( ID ) ,
  VID = as.integer( VID ) ) )

{set.seed(123)
  Model_GM_fri <- ulam(
    alist(
      Tie_fri ~ binomial( Possible_tie , p ) ,
      logit( p ) <- bGA[Gender] * Age + 
        bGR[Gender,PMR] + 
        bM[Gender] * Market +
        V[VID] + # village intercepts 
        I[ID] , # individual intercepts 
      
      # define effects using other parameters
      transpars> vector[VID]: V <<- V_bar + z_V*sigma_V,
      transpars> vector[ID]: I <<- z_I*sigma_I,
      
      bGA[Gender] ~ normal( 0 , 0.5 ) ,
      
      matrix[Gender,PMR]: bGR ~ normal( 0 , 0.5 ) ,
      
      bM[Gender] ~ normal( 0 , 0.5 ) ,
      
      V_bar ~ normal( 0 , 0.5 ) ,
      z_V[VID] ~ normal( 0 , 1 ) ,
      sigma_V ~ exponential( 1 ) ,
      z_I[ID] ~ normal( 0 , 4 ) ,
      sigma_I ~ exponential( 1 )
    ) , data = Model_GM_fri_list , 
    iter = 1500 , warmup = 500 , chains = 6 , cores = 6 , 
    log_lik = TRUE , control = list(adapt_delta = 0.99)
  )}
Pre_model_GM_fri <- precis( Model_GM_fri , depth = 3 , prob = 0.90 , 
                            pars = c( "bGA" , "bGR" , "bM" ) )
Pre_model_GM_fri

###3.1.5 Figures ----
# Posterior samples of estimates
Post_model_GM_all <- extract.samples( Model_GM_all )
Post_model_GM_genetic <- extract.samples( Model_GM_genetic )
Post_model_GM_affine <- extract.samples( Model_GM_affine )
Post_model_GM_fri <- extract.samples( Model_GM_fri )

####3.1.5.1 Posterior estimates for gender-specific effects of market participation ----
jpeg( "Estimates of market and gender.jpeg" , 
      width = 200 , height = 160 , units = "mm" , res = 300 )
{
  par( mfrow = c( 2 , 2 ) , mar = c( 5 , 4 , 1 , 1 ) )
  dens( Post_model_GM_all$bM[,1] ,
        xlim = c( -0.25 , 0.2 ) ,
        ylim = c( 0 , 14 ) ,
        xlab = "" , 
        ylab = "" ,
        lwd = 3 , col = 2 , cex.axis = 1.2 )
  dens( Post_model_GM_all$bM[,2] ,
        lwd = 3 , col = 4 , cex.axis = 1.2 , 
        add = T )
  mtext( "All" , side = 1 , line = 2.5 , cex = 1.2 )
  mtext( "Density" , side = 2 , line = 2.5 , cex = 1.2 )
  
  dens( Post_model_GM_genetic$bM[,1] ,
        xlim = c( -0.4 , 0.5 ) ,
        ylim = c( 0 , 7 ) ,
        xlab = "" , 
        ylab = "" ,
        lwd = 3 , col = 2 , cex.axis = 1.2 ,
        cex.lab = 1.5 )
  dens( Post_model_GM_genetic$bM[,2] ,
        lwd = 3 , col = 4 , cex.axis = 1.2 ,
        cex.lab = 1.5 , 
        add = T )
  mtext( "Biological kin" , side = 1 , line = 2.5 , cex = 1.2 )
  
  dens( Post_model_GM_affine$bM[,1] ,
        xlim = c( -0.4 , 0.4 ) ,
        ylim = c( 0 , 9 ) ,
        xlab = "" , 
        ylab = "" ,
        lwd = 3 , col = 2 , cex.axis = 1.2 ,
        cex.lab = 1.5 )
  dens( Post_model_GM_affine$bM[,2] ,
        lwd = 3 , col = 4 , cex.axis = 1.2 ,
        cex.lab = 1.5 , 
        add = T )
  mtext( "Affinal kin" , side = 1 , line = 2.5 , cex = 1.2 )
  mtext( "Density" , side = 2 , line = 2.5 , cex = 1.2 )
  
  dens( Post_model_GM_fri$bM[,1] ,
        xlim = c( -0.6 , 0.4 ) ,
        ylim = c( 0 , 9 ) ,
        xlab = "" , 
        ylab = "" ,
        lwd = 3 , col = 2 , cex.axis = 1.2 ,
        cex.lab = 1.5 )
  dens( Post_model_GM_fri$bM[,2] ,
        lwd = 3 , col = 4 , cex.axis = 1.2 ,
        cex.lab = 1.5 , 
        add = T )
  mtext( "Friend" , side = 1 , line = 2.5 , cex = 1.2 )
}

dev.off()

####3.1.5.2 Predicted effects of market participation by gender ----
# Function to predict p-matrix based on posterior samples of estimates and new data
predict_p_matrix_GM <- function( post_samples , new_data , n_samples = 6000 ) {
  n_market <- length( unique( new_data$N_market ) )  # Number of markets
  p_matrix <- matrix( NA , nrow = n_samples , ncol = n_market )  # Initialize matrix
  
  # Loop over posterior samples
  for (s in 1:n_samples) {
    
    # Extract posterior sample for current iteration
    bGA <- post_samples$bGA[s, ]
    bGR <- post_samples$bGR[s, , ]
    bM <- post_samples$bM[s, ]
    V_bar <- post_samples$V_bar[s, 1]
    
    # Loop over each market
    for (m in 1:n_market) {
      market_subset <- new_data[new_data$N_market == m, ]  # Subset data for market
      logit_p_vals <- with( market_subset,
                            bGA[Gender] * Age + 
                              bGR[cbind(Gender, PMR)] + 
                              bM[Gender] * Market + 
                              V_bar)
      
      # Convert to probabilities and average over Age and PMR
      p_matrix[s, m] <- mean( rethinking::inv_logit( logit_p_vals ) , na.rm = TRUE)
    }
  }
  
  return(p_matrix)
}

# New data for women and men, respectively 
New_data_GM_F <- data.frame( Age_act = as.integer( rep( 20:86 , 222 ) ) , 
                             Age = rep( ( as.integer( 20:86 ) - mean( Ego.data$Age ) )/ sd( Ego.data$Age ) , 222) , 
                             Gender = as.integer( rep( 1 , 14874 ) ) ,
                             PMR = as.integer( rep( 1:3 , 4958 ) ) ,
                             N_market = rep( 1:74 , 201 ) ,
                             Market_act = rep( seq( 0 , 730 , 10 ) , 201 ) ,
                             Market = rep( c( ( log( c( seq( 0 , 730 , 10 ) + 1 ) ) - mean( Ego.data$Market_log ) ) / sd( Ego.data$Market_log ) ) ) , 201 ) 

New_data_GM_M <- data.frame( Age_act = as.integer( rep( 20:86 , 222 ) ), 
                             Age = rep( ( as.integer( 20:86 ) - mean( Ego.data$Age ) )/ sd( Ego.data$Age ) , 222) , 
                             Gender = as.integer( rep( 2 , 14874 ) ) ,
                             PMR = as.integer( rep( 1:3 , 4958 ) ) ,
                             N_market = rep( 1:74 , 201 ) ,
                             Market_act = rep( seq( 0 , 730 , 10 ) , 201 ) ,
                             Market = rep( c( ( log( c( seq( 0 , 730 , 10 ) + 1 ) ) - mean( Ego.data$Market_log ) ) / sd( Ego.data$Market_log ) ) ) , 201 ) 

# Predicted data based on new data
p_matrix_GM_all_F <- predict_p_matrix_GM( post_samples = Post_model_GM_all , 
                                          new_data = New_data_GM_F, 
                                          n_samples = 6000 )
p_matrix_GM_all_M <- predict_p_matrix_GM( post_samples = Post_model_GM_all , 
                                          new_data = New_data_GM_M, 
                                          n_samples = 6000 )

p_matrix_GM_genetic_F <- predict_p_matrix_GM( post_samples = Post_model_GM_genetic , 
                                              new_data = New_data_GM_F, 
                                              n_samples = 6000 )
p_matrix_GM_genetic_M <- predict_p_matrix_GM( post_samples = Post_model_GM_genetic , 
                                              new_data = New_data_GM_M, 
                                              n_samples = 6000 )

p_matrix_GM_affine_F <- predict_p_matrix_GM( post_samples = Post_model_GM_affine , 
                                             new_data = New_data_GM_F, 
                                             n_samples = 6000 )
p_matrix_GM_affine_M <- predict_p_matrix_GM( post_samples = Post_model_GM_affine , 
                                             new_data = New_data_GM_M, 
                                             n_samples = 6000 )

p_matrix_GM_fri_F <- predict_p_matrix_GM( post_samples = Post_model_GM_fri , 
                                          new_data = New_data_GM_F, 
                                          n_samples = 6000 )
p_matrix_GM_fri_M <- predict_p_matrix_GM( post_samples = Post_model_GM_fri , 
                                          new_data = New_data_GM_M, 
                                          n_samples = 6000 )

#####3.1.5.2.1 Gender-specific effects of market participation ----
Market_seq <- c( ( log( c( seq( 0 , 730 , 10 ) + 1 ) ) - mean( Ego.data$Market_log ) ) / sd( Ego.data$Market_log ) )

jpeg( "Market_stratified_by_gender.jpg" , 
      width = 200 , height = 180 , units = "mm" , res = 300 )

{par(mfrow=c(2,2), mar = c(4.5, 4.5, 1.8, 1))
  # All
  plot( NULL,
        xlim = c(-2.5, 2.2),
        ylim = c(0,1),
        xlab = "" ,
        ylab = "Predicted density" ,
        main = "All" ,
        cex.axis = 1.2 ,
        cex.lab = 1.5 ,
        cex.main = 1.5 )
  # women
  lines( Market_seq , 
         apply( p_matrix_GM_all_F , 2 , mean ) ,
         lwd = 2 , col = "#CB2313" )
  shade( apply( p_matrix_GM_all_F , 2 , PI , prob = 0.90 ) , 
         Market_seq , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  shade( apply( p_matrix_GM_all_F , 2 , PI , prob = 0.60 ) , 
         Market_seq , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  shade( apply( p_matrix_GM_all_F , 2 , PI , prob = 0.30 ) , 
         Market_seq , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  # men
  lines( Market_seq , 
         apply( p_matrix_GM_all_M , 2 , mean ) ,
         lwd = 2 , col = "#046C9A" )
  shade( apply( p_matrix_GM_all_M , 2 , PI , prob = 0.90 ) , 
         Market_seq , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  shade( apply( p_matrix_GM_all_M , 2 , PI , prob = 0.60 ) , 
         Market_seq , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  shade( apply( p_matrix_GM_all_M , 2 , PI , prob = 0.30 ) , 
         Market_seq , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  
  # Biological kin
  plot( NULL,
        xlim = c(-2.5, 2.2),
        ylim = c(0,1),
        xlab = "" ,
        ylab = "" ,
        main = "Biological kin" ,
        cex.axis = 1.2 ,
        cex.lab = 1.5 ,
        cex.main = 1.5 )
  legend( x = 0.8 , y = 1.0 ,  
          box.col = "white",
          legend = c( "Female" , "Male" ) , 
          lty = c( 1 , 1 ) , 
          col = c( "#CB2313" , "#046C9A" ) , 
          lwd = 2 ,
          cex = 1.2 , 
          bty = "n" ,
          y.intersp = 1.2 ,
          x.intersp = 0.3 ,
          seg.len = 0.8  )
  
  # women
  lines( Market_seq , 
         apply( p_matrix_GM_genetic_F , 2 , mean ) ,
         lwd = 2 , col = "#CB2313" )
  shade( apply( p_matrix_GM_genetic_F , 2 , PI , prob = 0.90 ) , 
         Market_seq , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  shade( apply( p_matrix_GM_genetic_F , 2 , PI , prob = 0.60 ) , 
         Market_seq , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  shade( apply( p_matrix_GM_genetic_F , 2 , PI , prob = 0.30 ) , 
         Market_seq , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  # men
  lines( Market_seq , 
         apply( p_matrix_GM_genetic_M , 2 , mean ) ,
         lwd = 2 , col = "#046C9A" )
  shade( apply( p_matrix_GM_genetic_M , 2 , PI , prob = 0.90 ) , 
         Market_seq , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  shade( apply( p_matrix_GM_genetic_M , 2 , PI , prob = 0.60 ) , 
         Market_seq , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  shade( apply( p_matrix_GM_genetic_M , 2 , PI , prob = 0.30 ) , 
         Market_seq , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  
  # Affinal kin
  plot( NULL,
        xlim = c(-2.5, 2.2),
        ylim = c(0,1),
        xlab = "Market participation (std)" ,
        ylab = "Predicted density" ,
        main = "Affinal kin" ,
        cex.axis = 1.2 ,
        cex.lab = 1.5 ,
        cex.main = 1.5 )
  # women
  lines( Market_seq , 
         apply( p_matrix_GM_affine_F , 2 , mean ) ,
         lwd = 2 , col = "#CB2313" )
  shade( apply( p_matrix_GM_affine_F , 2 , PI , prob = 0.90 ) , 
         Market_seq , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  shade( apply( p_matrix_GM_affine_F , 2 , PI , prob = 0.60 ) , 
         Market_seq , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  shade( apply( p_matrix_GM_affine_F , 2 , PI , prob = 0.30 ) , 
         Market_seq , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  # men
  lines( Market_seq , 
         apply( p_matrix_GM_affine_M , 2 , mean ) ,
         lwd = 2 , col = "#046C9A" )
  shade( apply( p_matrix_GM_affine_M , 2 , PI , prob = 0.90 ) , 
         Market_seq , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  shade( apply( p_matrix_GM_affine_M , 2 , PI , prob = 0.60 ) , 
         Market_seq , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  shade( apply( p_matrix_GM_affine_M , 2 , PI , prob = 0.30 ) , 
         Market_seq , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  
  # Friend kin
  plot( NULL,
        xlim = c(-2.5,2.2),
        ylim = c(0,1),
        xlab = "Market participation (std)" ,
        ylab = "" ,
        main = "Friend" ,
        cex.axis = 1.2 ,
        cex.lab = 1.5 ,
        cex.main = 1.5 )
  # women
  lines( Market_seq , 
         apply( p_matrix_GM_fri_F , 2 , mean ) ,
         lwd = 2 , col = "#CB2313" )
  shade( apply( p_matrix_GM_fri_F , 2 , PI , prob = 0.90 ) , 
         Market_seq , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  shade( apply( p_matrix_GM_fri_F , 2 , PI , prob = 0.60 ) , 
         Market_seq , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  shade( apply( p_matrix_GM_fri_F , 2 , PI , prob = 0.30 ) , 
         Market_seq , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  # men
  lines( Market_seq , 
         apply( p_matrix_GM_fri_M , 2 , mean ) ,
         lwd = 2 , col = "#046C9A" )
  shade( apply( p_matrix_GM_fri_M , 2 , PI , prob = 0.90 ) , 
         Market_seq , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  shade( apply( p_matrix_GM_fri_M , 2 , PI , prob = 0.60 ) , 
         Market_seq , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  shade( apply( p_matrix_GM_fri_M , 2 , PI , prob = 0.30 ) , 
         Market_seq , 
         col = col.alpha( "#046C9A" , 0.15 ) )
}

dev.off()

#####3.1.5.2.2 Gender difference in effects of market participation ----
jpeg( "Market_gender_diff.jpg" , 
      width = 200 , height = 160 , units = "mm" , res = 300 )

{par(mfrow=c(2,2), mar = c(4.5, 4.5, 1.8, 1))
  # All
  plot( NULL,
        xlim = c(-2.5, 2.2),
        ylim = c( -0.4 , 0.4 ),
        xlab = "" ,
        ylab = "Predicted density" ,
        main = "All" ,
        cex.axis = 1.2 ,
        cex.lab = 1.5 ,
        cex.main = 1.5 )
  abline( h = 0 , lty = 2 , lwd = 2, col = 1 )
  lines( Market_seq , 
         apply( p_matrix_GM_all_F - p_matrix_GM_all_M , 2 , mean ) ,
         lwd = 2 , col = "black" )
  shade( apply( p_matrix_GM_all_F - p_matrix_GM_all_M , 2 , PI , prob = 0.90 ) , 
         Market_seq , 
         col = col.alpha( "dimgray" , 0.15 ) )
  shade( apply( p_matrix_GM_all_F - p_matrix_GM_all_M , 2 , PI , prob = 0.60 ) , 
         Market_seq , 
         col = col.alpha( "dimgray" , 0.15 ) )
  shade( apply( p_matrix_GM_all_F - p_matrix_GM_all_M , 2 , PI , prob = 0.30 ) , 
         Market_seq , 
         col = col.alpha( "dimgray" , 0.15 ) )
  
  # Genetic kin
  plot( NULL,
        xlim = c(-2.5, 2.2),
        ylim = c( -0.4 , 0.4 ),
        xlab = "" ,
        ylab = "" ,
        main = "Biological kin" ,
        cex.axis = 1.2 ,
        cex.lab = 1.5 ,
        cex.main = 1.5 )
  abline( h = 0 , lty = 2 , lwd = 2, col = 1 )
  lines( Market_seq , 
         apply( p_matrix_GM_genetic_F - p_matrix_GM_genetic_M , 2 , mean ) ,
         lwd = 2 , col = "black" )
  shade( apply( p_matrix_GM_genetic_F - p_matrix_GM_genetic_M , 2 , PI , prob = 0.90 ) , 
         Market_seq , 
         col = col.alpha( "dimgray" , 0.15 ) )
  shade( apply( p_matrix_GM_genetic_F - p_matrix_GM_genetic_M , 2 , PI , prob = 0.60 ) , 
         Market_seq , 
         col = col.alpha( "dimgray" , 0.15 ) )
  shade( apply( p_matrix_GM_genetic_F - p_matrix_GM_genetic_M , 2 , PI , prob = 0.30 ) , 
         Market_seq , 
         col = col.alpha( "dimgray" , 0.15 ) )
  
  # Affinal kin
  plot( NULL,
        xlim = c(-2.5, 2.2),
        ylim = c( -0.4 , 0.4 ),
        xlab = "Market participation (std)" ,
        ylab = "Predicted density" ,
        main = "Affinal kin" ,
        cex.axis = 1.2 ,
        cex.lab = 1.5 ,
        cex.main = 1.5 )
  abline( h = 0 , lty = 2 , lwd = 2, col = 1 )
  lines( Market_seq , 
         apply( p_matrix_GM_affine_F - p_matrix_GM_affine_M , 2 , mean ) ,
         lwd = 2 , col = "black" )
  shade( apply( p_matrix_GM_affine_F - p_matrix_GM_affine_M , 2 , PI , prob = 0.90 ) , 
         Market_seq , 
         col = col.alpha( "dimgray" , 0.15 ) )
  shade( apply( p_matrix_GM_affine_F - p_matrix_GM_affine_M , 2 , PI , prob = 0.60 ) , 
         Market_seq , 
         col = col.alpha( "dimgray" , 0.15 ) )
  shade( apply( p_matrix_GM_affine_F - p_matrix_GM_affine_M , 2 , PI , prob = 0.30 ) , 
         Market_seq , 
         col = col.alpha( "dimgray" , 0.15 ) )
  
  # Friend kin
  plot( NULL,
        xlim = c(-2.5,2.2),
        ylim = c( -0.4 , 0.4 ),
        xlab = "Market participation (std)" ,
        ylab = "" ,
        main = "Friend" ,
        cex.axis = 1.2 ,
        cex.lab = 1.5 ,
        cex.main = 1.5 )
  abline( h = 0 , lty = 2 , lwd = 2, col = 1 )
  lines( Market_seq , 
         apply( p_matrix_GM_fri_F - p_matrix_GM_fri_M , 2 , mean ) ,
         lwd = 2 , col = "black" )
  shade( apply( p_matrix_GM_fri_F - p_matrix_GM_fri_M , 2 , PI , prob = 0.90 ) , 
         Market_seq , 
         col = col.alpha( "dimgray" , 0.15 ) )
  shade( apply( p_matrix_GM_fri_F - p_matrix_GM_fri_M , 2 , PI , prob = 0.60 ) , 
         Market_seq , 
         col = col.alpha( "dimgray" , 0.15 ) )
  shade( apply( p_matrix_GM_fri_F - p_matrix_GM_fri_M , 2 , PI , prob = 0.30 ) , 
         Market_seq , 
         col = col.alpha( "dimgray" , 0.15 ) )
}

dev.off()

###3.1.6 Outputs ----
Output_model_GM <- data.frame( Mean = c( Pre_model_GM_all$mean ,
                                         Pre_model_GM_genetic$mean ,
                                         Pre_model_GM_affine$mean ,
                                         Pre_model_GM_fri$mean ) ,
                               CI5 = c( Pre_model_GM_all$`5%` ,
                                        Pre_model_GM_genetic$`5%` ,
                                        Pre_model_GM_affine$`5%` ,
                                        Pre_model_GM_fri$`5%` ) ,
                               CI95 = c( Pre_model_GM_all$`95%` ,
                                         Pre_model_GM_genetic$`95%` ,
                                         Pre_model_GM_affine$`95%` ,
                                         Pre_model_GM_fri$`95%` ) ,
                               Type = rep( c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) , each = 10 ) ,
                               Variable = rep( c( "Age[female]" , "Age[male]" , 
                                                  "Female: matrilocal" , 
                                                  "Male: matrilocal" , 
                                                  "Female: bilocal" , 
                                                  "Male: bilocal" , 
                                                  "Female: patrilocal" , 
                                                  "Male: patrilocal" , 
                                                  "Market[female]" , "Market[male]" ) , 4 ) ) %>% 
  mutate( Type = factor( Type , 
                         levels = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) ,
                         labels = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) ) ,
          Variable = factor( Variable , 
                             levels = c( "Age[female]" , "Age[male]" , 
                                         "Female: matrilocal" , 
                                         "Female: bilocal" , 
                                         "Female: patrilocal" , 
                                         "Male: patrilocal" , 
                                         "Male: bilocal" , 
                                         "Male: matrilocal" , 
                                         "Market[female]" , "Market[male]" ) ,
                             labels = c( "Age[female]" , "Age[male]" , 
                                         "Female: matrilocal" , 
                                         "Female: bilocal" , 
                                         "Female: patrilocal" , 
                                         "Male: patrilocal" , 
                                         "Male: bilocal" , 
                                         "Male: matrilocal" , 
                                         "Market[female]" , "Market[male]" ) ) ,
          Model = "Market participation" , 
          "Mean [CI]" = paste( format( round( Mean , 2 ), nsmall = 2 ) , 
                               " [" , 
                               format( round( CI5 , 2 ), nsmall = 2 ) , 
                               ", " , 
                               format( round( CI95 , 2 ), nsmall = 2 ) , 
                               "]", 
                               sep = "") )

GM.output <- tibble( Variable =  c( "Age[female]" , "Age[male]" , 
                                    "Female: matrilocal" , 
                                    "Female: bilocal" , 
                                    "Female: patrilocal" , 
                                    "Male: matrilocal" , 
                                    "Male: bilocal" , 
                                    "Male: patrilocal" ,
                                    "Market[female]" , "Market[male]" ) ,
                     Type = "All" ) %>% 
  left_join( Output_model_GM[ , c( 4 , 5 , 7 ) ] ,
             by = c( "Variable" , "Type" ) ) %>% 
  dplyr::rename( `All` = "Mean [CI]" ) %>% 
  mutate( Type = "Biological kin" ) %>% 
  left_join( Output_model_GM[ , c( 4 , 5 , 7 ) ] ,
             by = c( "Variable" , "Type" ) ) %>% 
  dplyr::rename( `Biological kin` = "Mean [CI]" ) %>% 
  mutate( Type = "Affinal kin" ) %>% 
  left_join( Output_model_GM[ , c( 4 , 5 , 7 ) ] ,
             by = c( "Variable" , "Type" ) ) %>% 
  dplyr::rename( `Affinal kin` = "Mean [CI]" ) %>% 
  mutate( Type = "Friend" ) %>% 
  left_join( Output_model_GM[ , c( 4 , 5 , 7 ) ] ,
             by = c( "Variable" , "Type" ) ) %>% 
  dplyr::rename( `Friend` = "Mean [CI]" ) %>% 
  select( - Type ) %>% 
  gt() %>% 
  fmt_number( ) %>% 
  sub_missing( columns = everything() , 
               missing_text = " " ) %>% 
  tab_style( style = list( cell_text( weight = "bold" ) ) ,
             locations = cells_column_labels( everything() ) ) %>% 
  tab_style( style = list( cell_text( font = "Times New Roman" ) ) ,
             locations = cells_body() ) %>% 
  cols_align( align = "left",
              columns = "Variable" ) %>%
  cols_align( align = "center",
              columns = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) ) %>% 
  opt_table_lines( extent = "default" ) %>%
  tab_options( column_labels.border.top.color = "black" ,
               column_labels.border.top.width = px( 3 ) ,
               column_labels.border.bottom.color = "black" ,
               table_body.hlines.color = "lightgrey" ,
               table.border.bottom.color = "black" ,
               table.border.bottom.width = px( 3 ) ) %>% 
  cols_width( "Variable" ~ px( 170 ) , 
              everything() ~ px( 150 ) )
GM.output

# Proportions of mean difference in βM values comparing women to men
tibble( Type = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) ,
        Woman = c( sprintf("%0.4f" , 
                           length( which( Post_model_GM_all$bM[,1] > 0 ) ) /6000 ) , 
                   sprintf("%0.4f" , 
                           length( which( Post_model_GM_genetic$bM[,1] > 0 ) ) /6000 ) , 
                   sprintf("%0.4f" , 
                           length( which( Post_model_GM_affine$bM[,1] > 0 ) ) /6000 ) , 
                   sprintf("%0.4f" , 
                           length( which( Post_model_GM_fri$bM[,1] > 0 ) ) /6000 ) ) , 
        Man = c( sprintf("%0.4f" , 
                         length( which( Post_model_GM_all$bM[,2] > 0 ) ) /6000 ) , 
                 sprintf("%0.4f" , 
                         length( which( Post_model_GM_genetic$bM[,2] > 0 ) ) /6000 ) , 
                 sprintf("%0.4f" , 
                         length( which( Post_model_GM_affine$bM[,2] > 0 ) ) /6000 ) , 
                 sprintf("%0.4f" , 
                         length( which( Post_model_GM_fri$bM[,2] > 0 ) ) /6000 ) ) , 
        `Gender difference` = c( sprintf("%0.4f" , 
                                         length( which( c( Post_model_GM_all$bM[,1] - Post_model_GM_all$bM[,2] ) > 0 ) ) /6000 ) , 
                                 sprintf("%0.4f" , 
                                         length( which( c( Post_model_GM_genetic$bM[,1] - Post_model_GM_genetic$bM[,2] ) > 0 ) ) /6000 ) , 
                                 sprintf("%0.4f" , 
                                         length( which( c( Post_model_GM_affine$bM[,1] - Post_model_GM_affine$bM[,2] ) > 0 ) ) /6000 ) , 
                                 sprintf("%0.4f" , 
                                         length( which( c( Post_model_GM_fri$bM[,1] - Post_model_GM_fri$bM[,2] ) > 0 ) ) /6000 ) ) ) %>% 
  mutate( 
    Woman = paste( round( as.numeric( Woman ) * 100 , 2 ) , "%" ) , 
    Man = paste( round( as.numeric( Man ) * 100 , 2 ) , "%" ) , 
    `Gender difference` = paste( round( as.numeric( `Gender difference` ) * 100 , 2 ) , "%" ) ) %>% 
  gt() %>% 
  fmt_number( ) %>% 
  sub_missing( columns = everything() , 
               missing_text = " " ) %>% 
  tab_style( style = list( cell_text( weight = "bold" ) ) ,
             locations = cells_column_labels(everything() ) ) %>% 
  tab_style( style = list( cell_text(font = "Times New Roman") ) ,
             locations = cells_body() ) %>% 
  cols_align( align = "center",
              columns = c( "Type" , Woman , Man , `Gender difference` ) ) %>% 
  opt_table_lines( extent = "default" ) %>%
  tab_options( column_labels.border.top.color = "black" ,
               column_labels.border.top.width = px( 3 ) ,
               column_labels.border.bottom.color = "black" ,
               table_body.hlines.color = "lightgrey" ,
               table.border.bottom.color = "black" ,
               table.border.bottom.width = px( 3 ) ) %>% 
  cols_width( everything() ~ px( 150 ) )

##3.2 Residence ----
# Other variable should be included into models
adjustmentSets( Dag , 
                exposure = "Residence" , 
                outcome = "Density" ,
                effect = "total" )

###3.2.1 Models of all density ----
Model_GR_all_list <- with( Ego.data , list(
  Possible_tie =  as.integer( Possible_tie ) , 
  Tie_all = as.integer( Tie_all ) , 
  Age = as.numeric( sprintf( "%0.4f" , Age_std ) ) ,
  Gender = as.integer( Gender ) , 
  PMR = as.integer( Residence ) , 
  ID = as.integer( ID ) ,
  VID = as.integer( VID ) ) )

{set.seed(123)
  Model_GR_all <- ulam(
    alist(
      Tie_all ~ binomial( Possible_tie , p ) ,
      logit( p ) <- bGA[Gender] * Age + 
        bGR[Gender,PMR] + 
        V[VID] + # village intercepts 
        I[ID] , # individual intercepts 
      
      # define effects using other parameters
      transpars> vector[VID]: V <<- V_bar + z_V*sigma_V,
      transpars> vector[ID]: I <<- z_I*sigma_I,
      
      bGA[Gender] ~ normal( 0 , 0.5 ) ,
      
      matrix[Gender,PMR]: bGR ~ normal( 0 , 0.5 ) ,
      
      V_bar ~ normal( 0 , 0.5 ) ,
      z_V[VID] ~ normal( 0 , 1 ) ,
      sigma_V ~ exponential( 1 ) ,
      z_I[ID] ~ normal( 0 , 4 ) ,
      sigma_I ~ exponential( 1 ) 
    ) , data = Model_GR_all_list , 
    iter = 1500 , warmup = 500 , chains = 6 , cores = 6 , 
    log_lik = TRUE , control = list(adapt_delta = 0.99)
  )}

Pre_model_GR_all <- precis( Model_GR_all , depth = 3 , prob = 0.90 , 
                            pars = c( "bGA" , "bGR" ) )
Pre_model_GR_all

###3.2.2 Models of biological kin density ----
Model_GR_genetic_list <- with( Ego.data , list(
  Possible_tie =  as.integer( Possible_tie ) , 
  Tie_genetic = as.integer( Tie_genetic ) ,
  Age = as.numeric( sprintf( "%0.4f" , Age_std ) ) ,
  Gender = as.integer( Gender ) , 
  PMR = as.integer( Residence ) , 
  ID = as.integer( ID ) ,
  VID = as.integer( VID ) ) )

{set.seed(123)
  Model_GR_genetic <- ulam(
    alist(
      Tie_genetic ~ binomial( Possible_tie , p ) ,
      logit( p ) <- bGA[Gender] * Age + 
        bGR[Gender,PMR] + 
        V[VID] + # village intercepts 
        I[ID] , # individual intercepts 
      
      # define effects using other parameters
      transpars> vector[VID]: V <<- V_bar + z_V*sigma_V,
      transpars> vector[ID]: I <<- z_I*sigma_I,
      
      bGA[Gender] ~ normal( 0 , 0.5 ) ,
      
      matrix[Gender,PMR]: bGR ~ normal( 0 , 0.5 ) ,
      
      V_bar ~ normal( 0 , 0.5 ) ,
      z_V[VID] ~ normal( 0 , 1 ) ,
      sigma_V ~ exponential( 1 ) ,
      z_I[ID] ~ normal( 0 , 4 ) ,
      sigma_I ~ exponential( 1 ) 
    ) , data = Model_GR_genetic_list , 
    iter = 1500 , warmup = 500 , chains = 6 , cores = 6 , 
    log_lik = TRUE , control = list(adapt_delta = 0.99)
  )}
Pre_model_GR_genetic <- precis( Model_GR_genetic , depth = 3 , prob = 0.90 , 
                                pars = c( "bGA" , "bGR" ) )
Pre_model_GR_genetic

###3.2.3 Models of affinal kin density ----
Model_GR_affine_list <- with( Ego.data , list(
  Possible_tie =  as.integer( Possible_tie ) , 
  Tie_affine = as.integer( Tie_affine ) ,
  Age = as.numeric( sprintf( "%0.4f" , Age_std ) ) ,
  Gender = as.integer( Gender ) , 
  PMR = as.integer( Residence ) , 
  ID = as.integer( ID ) ,
  VID = as.integer( VID ) ) )

{set.seed(123)
  Model_GR_affine <- ulam(
    alist(
      Tie_affine ~ binomial( Possible_tie , p ) ,
      logit( p ) <- bGA[Gender] * Age + 
        bGR[Gender,PMR] + 
        V[VID] + # village intercepts 
        I[ID] , # individual intercepts 
      
      # define effects using other parameters
      transpars> vector[VID]: V <<- V_bar + z_V*sigma_V,
      transpars> vector[ID]: I <<- z_I*sigma_I,
      
      bGA[Gender] ~ normal( 0 , 0.5 ) ,
      
      matrix[Gender,PMR]: bGR ~ normal( 0 , 0.5 ) ,
      
      V_bar ~ normal( 0 , 0.5 ) ,
      z_V[VID] ~ normal( 0 , 1 ) ,
      sigma_V ~ exponential( 1 ) ,
      z_I[ID] ~ normal( 0 , 4 ) ,
      sigma_I ~ exponential( 1 ) 
    ) , data = Model_GR_affine_list , 
    iter = 1500 , warmup = 500 , chains = 6 , cores = 6 , 
    log_lik = TRUE , control = list(adapt_delta = 0.99)
  )}
Pre_model_GR_affine <- precis( Model_GR_affine , depth = 3 , prob = 0.90 , 
                               pars = c( "bGA" , "bGR" ) )
Pre_model_GR_affine

###3.2.4 Models of friend density ----
Model_GR_fri_list <- with( Ego.data , list(
  Possible_tie =  as.integer( Possible_tie ) , 
  Tie_fri = as.integer( Tie_fri ) ,
  Age = as.numeric( sprintf( "%0.4f" , Age_std ) ) ,
  Gender = as.integer( Gender ) , 
  PMR = as.integer( Residence ) , 
  ID = as.integer( ID ) ,
  VID = as.integer( VID ) ) )

{set.seed(123)
  Model_GR_fri <- ulam(
    alist(
      Tie_fri ~ binomial( Possible_tie , p ) ,
      logit( p ) <- bGA[Gender] * Age + 
        bGR[Gender,PMR] + 
        V[VID] + # village intercepts 
        I[ID] , # individual intercepts 
      
      # define effects using other parameters
      transpars> vector[VID]: V <<- V_bar + z_V*sigma_V,
      transpars> vector[ID]: I <<- z_I*sigma_I,
      
      bGA[Gender] ~ normal( 0 , 0.5 ) ,
      
      matrix[Gender,PMR]: bGR ~ normal( 0 , 0.5 ) ,
      
      V_bar ~ normal( 0 , 0.5 ) ,
      z_V[VID] ~ normal( 0 , 1 ) ,
      sigma_V ~ exponential( 1 ) ,
      z_I[ID] ~ normal( 0 , 4 ) ,
      sigma_I ~ exponential( 1 ) 
    ) , data = Model_GR_fri_list , 
    iter = 1500 , warmup = 500 , chains = 6 , cores = 6 , 
    log_lik = TRUE , control = list(adapt_delta = 0.99)
  )}
Pre_model_GR_fri <- precis( Model_GR_fri , depth = 3 , prob = 0.90 , 
                            pars = c( "bGA" , "bGR" ) )
Pre_model_GR_fri

###3.2.5 Figures ----
# Posterior samples of estimates 
Post_model_GR_all <- extract.samples( Model_GR_all )
Post_model_GR_genetic <- extract.samples( Model_GR_genetic )
Post_model_GR_affine <- extract.samples( Model_GR_affine )
Post_model_GR_fri <- extract.samples( Model_GR_fri )

####3.2.5.1 Posterior estimates for gender-specific effects of residence ----
jpeg( "Estimates of residence and gender.jpeg" , 
      width = 180 , height = 180 , units = "mm" , res = 300 )

{
  par( mfrow = c( 4 , 3 ) , mar = c( 2 , 3.2 , 2 , 1 ) )
  # All
  dens( Post_model_GR_all$bGR[,1,1] ,
        xlim = c( -1.6 , 1.6 ) ,
        ylim = c( 0 , 2.5 ) ,
        xlab = "" , 
        ylab = "" ,
        lwd = 3 , col = 2 , cex.axis = 1.2 )
  dens( Post_model_GR_all$bGR[,2,1] ,
        lwd = 3 , col = 4 , 
        add = T )
  mtext( "All" , side = 2 , line = 1.9 , cex = 0.9 )
  mtext( "Matrilocal" , side = 3 , line = 0.2 , cex= 0.9 )
  
  dens( Post_model_GR_all$bGR[,1,2] ,
        xlim = c( -1.6 , 1.6 ) ,
        ylim = c( 0 , 2.5 ) ,
        xlab = "" , 
        ylab = "" ,
        lwd = 3 , col = 2 , cex.axis = 1.2 )
  dens( Post_model_GR_all$bGR[,2,2] ,
        lwd = 3 , col = 4 , 
        add = T )
  mtext( "Bilocal" , side = 3 , line = 0.2 , cex= 0.9 )
  
  dens( Post_model_GR_all$bGR[,1,3] ,
        xlim = c( -1.6 , 1.6 ) ,
        ylim = c( 0 , 2.5 ) ,
        xlab = "" , 
        ylab = "" ,
        lwd = 3 , col = 2 , cex.axis = 1.2 )
  dens( Post_model_GR_all$bGR[,2,3] ,
        lwd = 3 , col = 4 , 
        add = T )
  mtext( "Patrilocal" , side = 3 , line = 0.2 , cex= 0.9 )
  
  # Biological kin
  dens( Post_model_GR_genetic$bGR[,1,1] ,
        xlim = c( -1.6 , 1.6 ) ,
        ylim = c( 0 , 2.5 ) ,
        xlab = "" , 
        ylab = "" ,
        lwd = 3 , col = 2 , cex.axis = 1.2 )
  dens( Post_model_GR_genetic$bGR[,2,1] ,
        lwd = 3 , col = 4 , 
        add = T )
  mtext( "Biological kin" , side = 2 , line = 1.9 , cex = 0.9 )
  
  dens( Post_model_GR_genetic$bGR[,1,2] ,
        xlim = c( -1.6 , 1.6 ) ,
        ylim = c( 0 , 2.5 ) ,
        xlab = "" , 
        ylab = "" ,
        lwd = 3 , col = 2 , cex.axis = 1.2 )
  dens( Post_model_GR_genetic$bGR[,2,2] ,
        lwd = 3 , col = 4 , 
        add = T )
  
  dens( Post_model_GR_genetic$bGR[,1,3] ,
        xlim = c( -1.6 , 1.6 ) ,
        ylim = c( 0 , 2.5 ) ,
        xlab = "" , 
        ylab = "" ,
        lwd = 3 , col = 2 , cex.axis = 1.2 )
  dens( Post_model_GR_genetic$bGR[,2,3] ,
        lwd = 3 , col = 4 , 
        add = T )
  
  # Affinal kin
  dens( Post_model_GR_affine$bGR[,1,1] ,
        xlim = c( -1.6 , 1.6 ) ,
        ylim = c( 0 , 2.5 ) ,
        xlab = "" , 
        ylab = "" ,
        lwd = 3 , col = 2 , cex.axis = 1.2 )
  dens( Post_model_GR_affine$bGR[,2,1] ,
        lwd = 3 , col = 4 , 
        add = T )
  mtext( "Affinal kin" , side = 2 , line = 1.9 , cex = 0.9 )
  
  dens( Post_model_GR_affine$bGR[,1,2] ,
        xlim = c( -1.6 , 1.6 ) ,
        ylim = c( 0 , 2.5 ) ,
        xlab = "" , 
        ylab = "" ,
        lwd = 3 , col = 2 , cex.axis = 1.2 )
  dens( Post_model_GR_affine$bGR[,2,2] ,
        lwd = 3 , col = 4 , 
        add = T )
  
  dens( Post_model_GR_affine$bGR[,1,3] ,
        xlim = c( -1.6 , 1.6 ) ,
        ylim = c( 0 , 2.5 ) ,
        xlab = "" , 
        ylab = "" ,
        lwd = 3 , col = 2 , cex.axis = 1.2 )
  dens( Post_model_GR_affine$bGR[,2,3] ,
        lwd = 3 , col = 4 , 
        add = T )
  
  # Friend
  dens( Post_model_GR_fri$bGR[,1,1] ,
        xlim = c( -1.6 , 1.6 ) ,
        ylim = c( 0 , 2.5 ) ,
        xlab = "" , 
        ylab = "" ,
        lwd = 3 , col = 2 , cex.axis = 1.2 )
  dens( Post_model_GR_fri$bGR[,2,1] ,
        lwd = 3 , col = 4 , 
        add = T )
  mtext( "Friend" , side = 2 , line = 1.9 , cex = 0.9 )
  
  dens( Post_model_GR_fri$bGR[,1,2] ,
        xlim = c( -1.6 , 1.6 ) ,
        ylim = c( 0 , 2.5 ) ,
        xlab = "" , 
        ylab = "" ,
        lwd = 3 , col = 2 , cex.axis = 1.2 )
  dens( Post_model_GR_fri$bGR[,2,2] ,
        lwd = 3 , col = 4 , 
        add = T )
  
  dens( Post_model_GR_fri$bGR[,1,3] ,
        xlim = c( -1.6 , 1.6 ) ,
        ylim = c( 0 , 2.5 ) ,
        xlab = "" , 
        ylab = "" ,
        lwd = 3 , col = 2 , cex.axis = 1.2 )
  dens( Post_model_GR_fri$bGR[,2,3] ,
        lwd = 3 , col = 4 , 
        add = T )
}

dev.off()

####3.2.5.2 Gender difference across residence ----
# F - M
tibble( Value = c( c( Post_model_GR_all$bGR[,1,1] - Post_model_GR_all$bGR[,2,1] ) , 
                   c( Post_model_GR_all$bGR[,1,2] - Post_model_GR_all$bGR[,2,2] ) , 
                   c( Post_model_GR_all$bGR[,1,3] - Post_model_GR_all$bGR[,2,3] ) , 
                   
                   c( Post_model_GR_genetic$bGR[,1,1] - Post_model_GR_genetic$bGR[,2,1] ) , 
                   c( Post_model_GR_genetic$bGR[,1,2] - Post_model_GR_genetic$bGR[,2,2] ) , 
                   c( Post_model_GR_genetic$bGR[,1,3] - Post_model_GR_genetic$bGR[,2,3] ) , 
                   
                   c( Post_model_GR_affine$bGR[,1,1] - Post_model_GR_affine$bGR[,2,1] ) , 
                   c( Post_model_GR_affine$bGR[,1,2] - Post_model_GR_affine$bGR[,2,2] ) , 
                   c( Post_model_GR_affine$bGR[,1,3] - Post_model_GR_affine$bGR[,2,3] ) , 
                   
                   c( Post_model_GR_fri$bGR[,1,1] - Post_model_GR_fri$bGR[,2,1] ) , 
                   c( Post_model_GR_fri$bGR[,1,2] - Post_model_GR_fri$bGR[,2,2] ) , 
                   c( Post_model_GR_fri$bGR[,1,3] - Post_model_GR_fri$bGR[,2,3] ) ) ,
        Residence = c(rep( c( rep( "Matrilocal" , 6000 ) , 
                              rep( "Bilocal" , 6000 ) , 
                              rep( "Patrilocal" , 6000 ) ) , 4 ) ) ,
        Type = c( rep( "All" , 18000 ) ,
                  rep( "Biological kin" , 18000 ) ,
                  rep( "Affinal kin" , 18000 ) ,
                  rep( "Friend" , 18000 ) ) ) %>% 
  mutate( Residence = factor( Residence , 
                              levels = c( "Matrilocal" , "Bilocal",  "Patrilocal" ) , 
                              labels = c( "Matrilocal" , "Bilocal",  "Patrilocal" ) ) , 
          Type = factor( Type , 
                         levels = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) , 
                         labels = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) ) ) %>% 
  ggplot( aes( x = Value , y = Residence , 
               fill = Residence , color = Residence) ) +
  facet_wrap( ~ Type ) +
  geom_vline( xintercept = 0, linetype = 2 , linewidth = 1 ,
              color = "dimgray" ) +
  stat_halfeye( .width = .90, height = 1.0 ) + 
  labs( title = "Gender difference" ) +
  scale_fill_manual( values = alpha( c( "#c5272d" , "#037f77" , "#0001a1" ) , 0.3 ) ) +
  scale_color_manual( values = c( "#c5272d" , "#037f77" , "#0001a1" ) ) +
  scale_y_discrete( NULL, labels = ggplot2:::parse_safe ) +
  xlab( "Posterior estimates" ) +
  coord_cartesian( ylim = c( 1.5, 3.6 ) ) +
  scale_x_continuous( limits = c( -1 , 1.5 ) , breaks = c( -1 , 0 , 1 ) ) +
  theme(plot.margin = margin(10, 10, 10, 20),
        strip.background = element_rect(color = "black", fill = "white") ,
        strip.text.x = element_text(size = 16),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(colour = "Gainsboro"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black",fill=NA),
        panel.background = element_blank(),
        plot.title = element_text(size=16),
        legend.position = "none",
        axis.title.x = element_text(size = 16,
                                    margin = margin(t = 0.3, r = 0, b = 0, l = 0,unit = "cm")),
        axis.text.x = element_text(colour="black",size=14,
                                   margin = margin(t = 0.1, r = 0, b = 0, l = 0,unit = "cm")),
        axis.title.y = element_text(size = 16,
                                    margin = margin(t = 0, r = 0.5, b = 0, l = 0,unit = "cm")),
        axis.text.y = element_text(colour="black",size=14,
                                   margin = margin(t = 0, r = 0.2, b = 0, l = 0,unit = "cm")),
        legend.title=element_text(size=14,face = "bold"),
        legend.text=element_text(size=12,face = "bold"))

ggsave( filename = "Gender_diff_by_residence.jpeg" , 
        width = 150 , height = 120 , units = "mm" , dpi = 300 )

###3.2.6 Outputs ----
Output_model_GR <- data.frame( Mean = c( Pre_model_GR_all$mean ,
                                         Pre_model_GR_genetic$mean ,
                                         Pre_model_GR_affine$mean ,
                                         Pre_model_GR_fri$mean ) ,
                               CI5 = c( Pre_model_GR_all$`5%` ,
                                        Pre_model_GR_genetic$`5%` ,
                                        Pre_model_GR_affine$`5%` ,
                                        Pre_model_GR_fri$`5%` ) ,
                               CI95 = c( Pre_model_GR_all$`95%` ,
                                         Pre_model_GR_genetic$`95%` ,
                                         Pre_model_GR_affine$`95%` ,
                                         Pre_model_GR_fri$`95%` ) ,
                               Type = rep( c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) , each = 8 ) ,
                               Variable = rep( c( "Age[female]" , "Age[male]" , 
                                                  "Female: matrilocal" , 
                                                  "Male: matrilocal" , 
                                                  "Female: bilocal" , 
                                                  "Male: bilocal" , 
                                                  "Female: patrilocal" , 
                                                  "Male: patrilocal" ) , 4 ) ) %>% 
  mutate( Type = factor( Type , 
                         levels = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) ,
                         labels = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) ) ,
          Variable = factor( Variable , 
                             levels = c( "Age[female]" , "Age[male]" , 
                                         "Female: matrilocal" , 
                                         "Female: bilocal" , 
                                         "Female: patrilocal" , 
                                         "Male: matrilocal" , 
                                         "Male: bilocal" , 
                                         "Male: patrilocal" ) ,
                             labels = c( "Age[female]" , "Age[male]" , 
                                         "Female: matrilocal" , 
                                         "Female: bilocal" , 
                                         "Female: patrilocal" , 
                                         "Male: matrilocal" , 
                                         "Male: bilocal" , 
                                         "Male: patrilocal" ) ) ,
          Model = "Residence" , 
          "Mean [CI]" = paste( format( round( Mean , 2 ), nsmall = 2 ) , 
                               " [" , 
                               format( round( CI5 , 2 ), nsmall = 2 ) , 
                               ", " , 
                               format( round( CI95 , 2 ), nsmall = 2 ) , 
                               "]", 
                               sep = "") )

Res.gender.output <- tibble( Variable = c( "Age[female]" , "Age[male]" , 
                                           "Female: matrilocal" , 
                                           "Female: bilocal" , 
                                           "Female: patrilocal" , 
                                           "Male: matrilocal" , 
                                           "Male: bilocal" , 
                                           "Male: patrilocal" ) ,
                             Type = "All" ) %>% 
  left_join( Output_model_GR[ , c( 4 , 5 , 7 ) ] ,
             by = c( "Variable" , "Type" ) ) %>% 
  dplyr::rename( `All` = "Mean [CI]" ) %>% 
  mutate( Type = "Biological kin" ) %>% 
  left_join( Output_model_GR[ , c( 4 , 5 , 7 ) ] ,
             by = c( "Variable" , "Type" ) ) %>% 
  dplyr::rename( `Biological kin` = "Mean [CI]" ) %>% 
  mutate( Type = "Affinal kin" ) %>% 
  left_join( Output_model_GR[ , c( 4 , 5 , 7 ) ] ,
             by = c( "Variable" , "Type" ) ) %>% 
  dplyr::rename( `Affinal kin` = "Mean [CI]" ) %>% 
  mutate( Type = "Friend" ) %>% 
  left_join( Output_model_GR[ , c( 4 , 5 , 7 ) ] ,
             by = c( "Variable" , "Type" ) ) %>% 
  dplyr::rename( `Friend` = "Mean [CI]" ) %>% 
  select( - Type ) %>% 
  gt() %>% 
  fmt_number( ) %>% 
  sub_missing( columns = everything() , 
               missing_text = " " ) %>% 
  tab_style( style = list( cell_text( weight = "bold" ) ) ,
             locations = cells_column_labels( everything() ) ) %>% 
  tab_style( style = list( cell_text( font = "Times New Roman" ) ) ,
             locations = cells_body() ) %>% 
  cols_align( align = "left",
              columns = "Variable" ) %>%
  cols_align( align = "center",
              columns = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) ) %>% 
  opt_table_lines( extent = "default" ) %>%
  tab_options( column_labels.border.top.color = "black" ,
               column_labels.border.top.width = px( 3 ) ,
               column_labels.border.bottom.color = "black" ,
               table_body.hlines.color = "lightgrey" ,
               table.border.bottom.color = "black" ,
               table.border.bottom.width = px( 3 ) ) %>% 
  cols_width( "Variable" ~ px( 170 ) , 
              everything() ~ px( 150 ) )
Res.gender.output

# Proportions of mean difference in βR values comparing women to men
tibble( Type = c( rep( "All" , 3 ) ,
                  rep( "Biological kin" , 3 ) ,
                  rep( "Affinal kin" , 3 ) ,
                  rep( "Friend" , 3 ) ) ,
        Residence = c(rep( c( "Matrilocal" , 
                              "Bilocal" , 
                              "Patrilocal" ) , 4 ) ) ,
        Woman = c( sprintf("%0.4f" , 
                           length( which( Post_model_GR_all$bGR[,1,1] > 0 ) ) /6000 ) , 
                   sprintf("%0.4f" , 
                           length( which( Post_model_GR_all$bGR[,1,2] > 0 ) ) /6000 ) , 
                   sprintf("%0.4f" , 
                           length( which( Post_model_GR_all$bGR[,1,3] > 0 ) ) /6000 ) , 
                   
                   sprintf("%0.4f" , 
                           length( which( Post_model_GR_genetic$bGR[,1,1] > 0 ) ) /6000 ) , 
                   sprintf("%0.4f" , 
                           length( which( Post_model_GR_genetic$bGR[,1,2] > 0 ) ) /6000 ) , 
                   sprintf("%0.4f" , 
                           length( which( Post_model_GR_genetic$bGR[,1,3] > 0 ) ) /6000 ) , 
                   
                   sprintf("%0.4f" , 
                           length( which( Post_model_GR_affine$bGR[,1,1] > 0 ) ) /6000 ) , 
                   sprintf("%0.4f" , 
                           length( which( Post_model_GR_affine$bGR[,1,2] > 0 ) ) /6000 ) , 
                   sprintf("%0.4f" , 
                           length( which( Post_model_GR_affine$bGR[,1,3] > 0 ) ) /6000 ) , 
                   
                   sprintf("%0.4f" , 
                           length( which( Post_model_GR_fri$bGR[,1,1] > 0 ) ) /6000 ) , 
                   sprintf("%0.4f" , 
                           length( which( Post_model_GR_fri$bGR[,1,2] > 0 ) ) /6000 ) , 
                   sprintf("%0.4f" , 
                           length( which( Post_model_GR_fri$bGR[,1,3] > 0 ) ) /6000 ) ) , 
        Man = c( sprintf("%0.4f" , 
                         length( which( Post_model_GR_all$bGR[,2,1] > 0 ) ) /6000 ) , 
                 sprintf("%0.4f" , 
                         length( which( Post_model_GR_all$bGR[,2,2] > 0 ) ) /6000 ) , 
                 sprintf("%0.4f" , 
                         length( which( Post_model_GR_all$bGR[,2,3] > 0 ) ) /6000 ) , 
                 
                 sprintf("%0.4f" , 
                         length( which( Post_model_GR_genetic$bGR[,2,1] > 0 ) ) /6000 ) , 
                 sprintf("%0.4f" , 
                         length( which( Post_model_GR_genetic$bGR[,2,2] > 0 ) ) /6000 ) , 
                 sprintf("%0.4f" , 
                         length( which( Post_model_GR_genetic$bGR[,2,3] > 0 ) ) /6000 ) , 
                 
                 sprintf("%0.4f" , 
                         length( which( Post_model_GR_affine$bGR[,2,1] > 0 ) ) /6000 ) , 
                 sprintf("%0.4f" , 
                         length( which( Post_model_GR_affine$bGR[,2,2] > 0 ) ) /6000 ) , 
                 sprintf("%0.4f" , 
                         length( which( Post_model_GR_affine$bGR[,2,3] > 0 ) ) /6000 ) , 
                 
                 sprintf("%0.4f" , 
                         length( which( Post_model_GR_fri$bGR[,2,1] > 0 ) ) /6000 ) , 
                 sprintf("%0.4f" , 
                         length( which( Post_model_GR_fri$bGR[,2,2] > 0 ) ) /6000 ) , 
                 sprintf("%0.4f" , 
                         length( which( Post_model_GR_fri$bGR[,2,3] > 0 ) ) /6000 ) ) ,
        `Gender difference` = c( sprintf("%0.4f" , 
                                         length( which( c( Post_model_GR_all$bGR[,1,1] - Post_model_GR_all$bGR[,2,1] ) > 0 ) ) /6000 ) , 
                                 sprintf("%0.4f" , 
                                         length( which( c( Post_model_GR_all$bGR[,1,2] - Post_model_GR_all$bGR[,2,2] ) > 0 ) ) /6000 ) , 
                                 sprintf("%0.4f" , 
                                         length( which( c( Post_model_GR_all$bGR[,1,3] - Post_model_GR_all$bGR[,2,3] ) > 0 ) ) /6000 ) , 
                                 
                                 sprintf("%0.4f" , 
                                         length( which( c( Post_model_GR_genetic$bGR[,1,1] - Post_model_GR_genetic$bGR[,2,1] ) > 0 ) ) /6000 ) , 
                                 sprintf("%0.4f" , 
                                         length( which( c( Post_model_GR_genetic$bGR[,1,2] - Post_model_GR_genetic$bGR[,2,2] ) > 0 ) ) /6000 ) , 
                                 sprintf("%0.4f" , 
                                         length( which( c( Post_model_GR_genetic$bGR[,1,3] - Post_model_GR_genetic$bGR[,2,3] ) > 0 ) ) /6000 ) , 
                                 
                                 sprintf("%0.4f" , 
                                         length( which( c( Post_model_GR_affine$bGR[,1,1] - Post_model_GR_affine$bGR[,2,1] ) > 0 ) ) /6000 ) , 
                                 sprintf("%0.4f" , 
                                         length( which( c( Post_model_GR_affine$bGR[,1,2] - Post_model_GR_affine$bGR[,2,2] ) > 0 ) ) /6000 ) , 
                                 sprintf("%0.4f" , 
                                         length( which( c( Post_model_GR_affine$bGR[,1,3] - Post_model_GR_affine$bGR[,2,3] ) > 0 ) ) /6000 ) , 
                                 
                                 sprintf("%0.4f" , 
                                         length( which( c( Post_model_GR_fri$bGR[,1,1] - Post_model_GR_fri$bGR[,2,1] ) > 0 ) ) /6000 ) , 
                                 sprintf("%0.4f" , 
                                         length( which( c( Post_model_GR_fri$bGR[,1,2] - Post_model_GR_fri$bGR[,2,2] ) > 0 ) ) /6000 ) , 
                                 sprintf("%0.4f" , 
                                         length( which( c( Post_model_GR_fri$bGR[,1,3] - Post_model_GR_fri$bGR[,2,3] ) > 0 ) ) /6000 ) ) ) %>% 
  mutate( 
    Woman = paste( round( as.numeric( Woman ) * 100 , 2 ) , "%" ) , 
    Man = paste( round( as.numeric( Man ) * 100 , 2 ) , "%" ) ,  
    `Gender difference` = paste( round( as.numeric( `Gender difference` ) * 100 , 2 ) , "%" ) ) %>% 
  gt() %>% 
  fmt_number( ) %>% 
  sub_missing( columns = everything() , 
               missing_text = " " ) %>% 
  tab_style( style = list( cell_text( weight = "bold" ) ) ,
             locations = cells_column_labels(everything() ) ) %>% 
  tab_style( style = list( cell_text(font = "Times New Roman") ) ,
             locations = cells_body() ) %>% 
  cols_align( align = "center",
              columns = c( "Type" , "Residence" , Woman , Man , `Gender difference` ) ) %>% 
  opt_table_lines( extent = "default" ) %>%
  tab_options( column_labels.border.top.color = "black" ,
               column_labels.border.top.width = px( 3 ) ,
               column_labels.border.bottom.color = "black" ,
               table_body.hlines.color = "lightgrey" ,
               table.border.bottom.color = "black" ,
               table.border.bottom.width = px( 3 ) ) %>% 
  cols_width( everything() ~ px( 150 ) )
