data {
  // data to model
  int<lower=0> n_dates;         // number case dates
  int<lower=1> n_strains;       // number of strains
  int<lower=0> cases[n_dates,n_strains];  // number of cases
  int<lower=0> data_freq;       // freqency of data points (i.e. weekly data then data_freq = 7)
  
  // seed times
  int<lower=0> seed_start_time[n_strains]; // time when a strain is started to be seeded
  int<lower=0> seed_end_time[n_strains];   // time when a strain is stopped being seeded
  
  // jump days (i.e. lockdowns)
  int<lower=0> n_jumps;              // number of points with jumps
  int<lower=0> jump_times[n_jumps];  // days of jump
  
  // fixed model parameters
  real generation_mean;         // mean generation time
  real generation_sd;           // sd of generation time
  int generation_max;           // maximum number of generation time    
  real test_mean;               // mean time to test from infection
  real test_sd;                 // sd time to test from infection
  real ascertainment_factor;    // ascertainment of cases from underlying infections
  
  // priors on other parameters
  real prior_R0_min;            // minimum of prior of R0
  real prior_R0_max;            // maximum of prior of R0
  real prior_dR_sd_min;         // minimum of prior on sd of dR (log-normal process)    
  real prior_dR_sd_max;         // minimum of prior on sd of dR (log-normal process)    
  real prior_phi_od_max;        // maximum of obersvation over dispersion parameter
  real prior_jump_up_max;       // maximum size of upward jump on jump days (e.g. 2)
  real prior_jump_down_max;     // maximum size of upward jump on jump days (e.g. 0.5)
  real prior_daily_seed_max[n_strains]; // maximum prior on daily seed
  real prior_strain_multipliers_min[n_strains-1]; // minimum prior on strain multipliers
  real prior_strain_multipliers_max[n_strains-1]; // minimum prior on strain multipliers
  
  // options
  int<lower=0,upper=1> multiply_strain_multipliers; // if true then multiply strain factors (improved paramertisation)
  int<lower=0> mask_cases_less_than;                // do not fit if the number of cases is very low
}

// transformed generation data
transformed data{
  int<lower=0> n_inf;
  int<lower=0> n_R;
  int<lower=0> n_strain_multipliers;
  int<lower=0> test_data_max;
  real<lower=0> inf0_data[n_strains];
  int<lower=1> idx_dp[n_dates];
  real<lower=0> generation_time[generation_max];  
  real<lower=0> test_time[generation_max];
  real<lower=0> test_time_data[generation_max+data_freq];
  real g_alpha = generation_mean * generation_mean / generation_sd/ generation_sd;
  real g_beta = generation_mean / generation_sd / generation_sd;
  real t_alpha = test_mean * test_mean / test_sd/ test_sd;
  real t_beta = test_mean / test_sd / test_sd;
  vector[ n_dates ] ones_dates;    
  int<lower=0> seed_start_time_adj[n_strains]; 
  int<lower=0> seed_end_time_adj[n_strains];
  real prior_daily_seed_max_max;
  int<lower=0> n_used_data[n_strains];
  int used_data[n_dates,n_strains];
  
  // set up generation and test kernel
  for( t in 1:generation_max ) {
    generation_time[t] = exp( gamma_lpdf( t | g_alpha, g_beta ) );
    test_time[t]       = exp( gamma_lpdf( t | t_alpha, t_beta ) );
  }
  
  // calculate the test kernel for the aggregated observed data points
  test_data_max = generation_max + data_freq;
  for( t in 1:test_data_max ) 
    test_time_data[t] = 0;
  for( t in 1:generation_max )
    for( td in 1:data_freq )
      test_time_data[t+td-1] += test_time[ t ];
  test_data_max = generation_max + data_freq - 1;
    
  // number of R abd infection days modelled and corresponding to data points
  n_R = ( n_dates - 1 ) * data_freq;
  n_inf = n_R + generation_max;
  for( t in 1:n_dates )
    idx_dp[t] = ( t- 1 ) * data_freq + generation_max;
  
  // initial cases
  for( sdx in 1:n_strains ) 
    inf0_data[sdx] = cases[1,sdx] / ascertainment_factor / data_freq;
  
  // useful array of ones
  for( t in 1:n_dates )
    ones_dates[t] = 1;
    
  // base strain has a multipler of 1
  n_strain_multipliers = n_strains - 1;
  
  // offset the seed times by the generation time (since infections are)
  for( sdx in 1:n_strains ) {
    seed_start_time_adj[sdx] = seed_start_time[sdx] + generation_max;
    seed_end_time_adj[sdx] = seed_end_time[sdx] + generation_max;
  }
  
  prior_daily_seed_max_max = max( prior_daily_seed_max);
  
  // work out which data is used in the fit (don't fit very small data)
  for( sdx in 1:n_strains ) {
    n_used_data[sdx] = 0;
    for( t in 1:n_dates ) {
      if( cases[ t, sdx ] >= mask_cases_less_than ) {
        n_used_data[sdx] += 1;
        used_data[ n_used_data[sdx], sdx ] = t;
      }
    }
  }
}


// fitted paramters
parameters {
  real<lower=prior_R0_min,upper=prior_R0_max> R0;
  real<lower=prior_dR_sd_min,upper=prior_dR_sd_max> dR_sd; 
  real dR[ n_dates - 1];
  real<lower=prior_jump_down_max,upper=prior_jump_up_max> jumps[ n_jumps];
  real<lower=0,upper=prior_phi_od_max> phi_od;
  real<lower=0,upper=1> strain_multipliers_raw[n_strain_multipliers];
  real<lower=0,upper=prior_daily_seed_max_max> daily_seed[n_strains];
}

// transformed parameters
transformed parameters{
  real<lower=0> R[ n_R ];
  real<lower=0> infections[n_inf,n_strains];
  real<lower=0> cases_expected_dp[n_dates,n_strains];
  real<lower=0> strain_multipliers[n_strain_multipliers];

  // inflate the raw strain multipliers base on the given priors
  for( sdx in 1:n_strain_multipliers )
    strain_multipliers[sdx] = prior_strain_multipliers_min[sdx] + 
      (prior_strain_multipliers_max[sdx]- prior_strain_multipliers_min[sdx]) * strain_multipliers_raw[sdx];
  
  if( multiply_strain_multipliers == 1 ) {
    for( sdx in 2:n_strain_multipliers )
      strain_multipliers[sdx]  = strain_multipliers[sdx] *strain_multipliers[sdx-1];
  }
  
  // calculate R as a log-process with jumps
  for( idx in 1:(n_dates-1))
    for( jdx in 1:data_freq )
      R[ ( idx - 1 ) * data_freq + jdx ] = dR[ idx ];
  R[1] = log(R0);
  for( idx in 1:n_jumps )
    R[jump_times[idx]] = log(jumps[idx]);
  R = exp( cumulative_sum( R) );
  
  // calculate infection series from R series
  for(sdx in 1:n_strains )
  {
    for( t_inf in 1:generation_max)
      infections[t_inf,sdx] = inf0_data[sdx];
    
    if( seed_start_time_adj[sdx] > generation_max )
    {
      for( t_inf in (generation_max+1):(seed_start_time_adj[sdx]-1))
      {
        infections[t_inf,sdx] = 0;
        for( t_g in 1:generation_max )
          infections[t_inf,sdx] += infections[t_inf-t_g,sdx] * generation_time[t_g];
        infections[t_inf,sdx] *= R[ t_inf - generation_max ];
        if( sdx > 1 )
          infections[t_inf,sdx] *= strain_multipliers[sdx-1];
      }
  
      for( t_inf in seed_start_time_adj[sdx]:seed_end_time_adj[sdx]) 
      {
        infections[t_inf,sdx] = daily_seed[sdx];
        for( t_g in 1:generation_max )
          infections[t_inf,sdx] += infections[t_inf-t_g,sdx] * generation_time[t_g];
        infections[t_inf,sdx] *= R[ t_inf - generation_max ];
        if( sdx > 1 )
          infections[t_inf,sdx] *= strain_multipliers[sdx-1];
      }
    }
    
    for( t_inf in (seed_end_time_adj[sdx]+1):n_inf) 
    {
      infections[t_inf,sdx] = 0;
      for( t_g in 1:generation_max )
        infections[t_inf,sdx] += infections[t_inf-t_g,sdx] * generation_time[t_g];
      infections[t_inf,sdx] *= R[ t_inf - generation_max ];
      if( sdx > 1 )
       infections[t_inf,sdx] *= strain_multipliers[sdx-1];
    }

    // calculate tests from infections series 
    cases_expected_dp[1,sdx] = fmax( cases[1,sdx], 0.0001);
    for( t in 2:n_dates)
    {
      cases_expected_dp[t,sdx] = 0;
      for( t_g in 1:test_data_max )
        cases_expected_dp[t,sdx] += infections[idx_dp[t]-t_g,sdx] * test_time_data[t_g];
        
      cases_expected_dp[t,sdx] *= ascertainment_factor;
      cases_expected_dp[t,sdx] += 0.1;
    }
  }
}

model {
  // priors
  // dR has a range prior
  // R0 has a range prior
  // jumps have a range prior
  // strain multipliers have range priors
  
  // infection process model (log normal)
  dR ~ normal( -0.5*dR_sd*dR_sd, dR_sd);
  
  // observation model
  for( sdx in 1:n_strains ) {
    cases[ used_data[ 1:n_used_data[ sdx ],sdx],sdx] ~ neg_binomial_2( cases_expected_dp[used_data[ 1:n_used_data[ sdx ],sdx],sdx], 
                                                                       ones_dates[ 1:n_used_data[ sdx ]] / phi_od );
  }
}

generated quantities{
  real<lower=0> R_comb[ n_R ];
  real<lower=0> R_denom;
  real<lower=0> R_numer;

  // calculate the composite R
  for( t_inf in (generation_max+1):n_inf )
  {
    R_denom = infections[ t_inf,1];
    R_numer = infections[ t_inf,1]; 
    for( sdx in 2:n_strains ){
      R_denom += infections[ t_inf,sdx] / strain_multipliers[ sdx - 1 ];
      R_numer += infections[ t_inf,sdx]; 
    }
    R_comb[ t_inf - generation_max ] = R[ t_inf - generation_max ] * R_numer / R_denom; 
  }
}

