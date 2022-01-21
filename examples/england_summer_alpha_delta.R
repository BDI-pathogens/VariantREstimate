######################################################################
# Example: Estimation of relative transmission strengths of the summer
#         (B.1.177), Alpha (B.1.1.7) and Delta (B.1.617.2) variants using 
#         data in England
#
# Case data is taken from the UK Dashboard and variant from COG-UK
#
# Model uses the renewal equation with a common R(t) which follows up
# diffusion process with jumps at start/end of lockdowns
# 
# Each strain has its own multiplier which multiplies R(t) in the renewal
# equation
#
# Seeding times are pre-determined with seeding rates estimated
######################################################################

library( VariantREstimate )
library( rstan )
library( matrixStats)
library( plotly)
library( data.table)

# define length of the modelled period and get the data
start_date <- as.Date( "2020-09-05")
end_date   <- as.Date( "2021-07-16")
cases      <- example.data.case_cog_summer_alpha_delta()

# define the dates when R is casn jump (start and end of Lockdowns)
jump_dates <- c(
  as.Date( "2020-11-05"),
  as.Date( "2020-12-03"),
  as.Date( "2021-01-04")
)

# define the seeding period for the st
seed_start_dates <- c(
  as.Date( "2020-09-06"),
  as.Date( "2020-09-06"),
  as.Date( "2020-09-16"),
  as.Date( "2021-03-16")
) 
seed_end_dates <- c(
  as.Date( "2020-09-06"),
  as.Date( "2020-09-24"),
  as.Date( "2020-10-16"),
  as.Date( "2021-04-16")
)
prior_daily_seed_max <- c( 0.001, 20, 20, 20 )

# get the Stan model 
model = model.strain_multiplier()

# build the input data for the model
data = list(
  # define the fitted data dimensions
  n_dates   = cases[,.N],
  n_strains = 4L,
  data_freq = 7L,
  cases     = as.matrix( cases[, .(cases_WT, cases_summer, cases_alpha, cases_delta)] ),
  
  # jumps in R due to lockdowns
  n_jumps    = length( jump_dates ),
  jump_times = as.integer( jump_dates - cases[, min(date)]),
  
  # seeding windows
  seed_start_time = as.integer( seed_start_dates - cases[, min(date)] ),
  seed_end_time   = as.integer( seed_end_dates - cases[, min(date)]) ,
  
  # fixed parameters
  generation_mean = 5.5,
  generation_sd   = 2.14,
  generation_max  = 20,
  test_mean       = 6.5,
  test_sd         = 2.7,
  ascertainment_factor = 0.35,
  
  # priors
  prior_R0_min    = 0.5,
  prior_R0_max    = 2,
  prior_dR_sd_min = 0.0001,
  prior_dR_sd_max = 0.02,
  prior_phi_od_max = 1,
  prior_jump_up_max = 1.5,
  prior_jump_down_max = 0.5,
  prior_daily_seed_max = prior_daily_seed_max,
  prior_strain_multipliers_min = c( 0.9, 1.25, 1.25 ),
  prior_strain_multipliers_max = c( 1.2, 1.75, 1.75 ),
  
  # options on underlying fitted 
  multiply_strain_multipliers = 1L,
  mask_cases_less_than = 100L
)

# MCMC parameters
n_chains <- 3
n_iter   <- 1000

# initialize with some reasonable parameters for a few variables
init <- lapply( 1:n_chains, function(x) list( 
  R0    = runif( 1, 1.2,1.4),
  jumps = c( runif(1,0.7,0.9), runif(1,1.3,1.5), runif(1,0.7,0.9)),
  dR    = rep( -0.03, ( cases[,.N] - 1 )  )
) )

# sample from model and collect key parameters
raw <- sampling( 
  model, 
  data = data, 
  chains = n_chains, 
  iter = n_iter, 
  pars = c("R", "R_comb", "jumps", "strain_multipliers", "daily_seed"), 
  cores = n_chains, 
  init = init
)
res <- extract(raw)

# generate plots of results
p1 <- plot_ly( 
  x = res$strain_multipliers[,1], 
  type = "histogram", 
  nbinsx = 50,
  showlegend = FALSE,
  histnorm = "probability",
  marker = list(color = 'rgb(158,202,225)',
                line = list(color = 'rgb(8,48,107)', width = 1.5))
) %>%
  layout( 
    xaxis = list( title = list( text = "B.1.177 transmissibility relative to WT", standoff = 1) ),
    yaxis = list( title = list( text = "posterior probability" ) )
  )

p2 <- plot_ly( 
  x = res$strain_multipliers[,2] / res$strain_multipliers[,1], 
  type = "histogram", 
  nbinsx = 50,
  histnorm = "probability",
  showlegend = FALSE,
  marker = list(color = 'rgb(158,202,225)',
                line = list(color = 'rgb(8,48,107)', width = 1.5))
) %>%
  layout( 
    xaxis = list( title = list( text = "Alpha transmissibility realtive to B.1.177" ), standoff = 1 ),
    yaxis = list( title = list( text = "posterior probability"))
  )

p3 <- plot_ly( 
  x = res$strain_multipliers[,3] / res$strain_multipliers[,2], 
  type = "histogram", 
  nbinsx = 50,
  histnorm = "probability",
  showlegend = FALSE,
  marker = list(color = 'rgb(158,202,225)',
                line = list(color = 'rgb(8,48,107)', width = 1.5))
) %>%
  layout( 
    xaxis = list( title = list( text = "Delta transmissibility realtive Alpha" ), standoff = 1 ),
    yaxis = list( title = list( text = "posterior probability"))
  )

p4 <- plot_ly( 
  x = start_date + (1:ncol(res$R)), 
  y = colQuantiles(res$R, probs = 0.05), 
  type = "scatter", 
  mode = "lines",
  line = list( color = rgb(0,0,0.5), width= 0 ),
  showlegend = FALSE
) %>%
  add_trace( y = ~colQuantiles(res$R, probs = 0.95), fill = "tonexty", fillcolor = "rgba(0,0,0.5,0.3)", showlegend = TRUE, name = "CI 5%-95%") %>%
  add_trace( y = colMedians(res$R ) ,line = list( width = 5), name = "median", showlegend = TRUE ) %>%
  layout( 
    yaxis  = list( title = list( text = "wild-type R(t)" ), range = c( 0.25,1.75) ),
    xaxis  = list( title = list( text = "date" ) ),
    legend = list( x = 0.9, y = 1 )
  )

p5 <- plot_ly( 
  x = start_date + (1:ncol(res$R_comb)), 
  y = colQuantiles(res$R_comb, probs = 0.05), 
  type = "scatter", 
  mode = "lines",
  line = list( color = rgb(0,0,0.5), width= 0 ),
  showlegend = FALSE
) %>%
  add_trace( y = ~colQuantiles(res$R_comb, probs = 0.95), fill = "tonexty", fillcolor = "rgba(0,0,0.5,0.3)", name = "CI 5%-95%") %>%
  add_trace( y = colMedians(res$R_comb ) ,line = list( width = 5), name = "median" ) %>%
  layout( 
    yaxis  = list( title = list( text = "composite R(t)"), range = c( 0.25,1.75)),
    xaxis  = list( title = list( text = "date" )),
    legend = list( x = 0.9, y = 1)
)

s1 = subplot( list( p1, p2, p3), nrows = 3, shareX = FALSE, titleX = TRUE, titleY = TRUE, margin = 0.05)
s2 = subplot( list( p4, p5), nrows = 2, shareX = TRUE, shareY = FALSE, titleY = TRUE, margin = 0.05 )
subplot( list( s1, s2 ), shareX = FALSE, titleX = TRUE, titleY = TRUE, margin = 0.05)
