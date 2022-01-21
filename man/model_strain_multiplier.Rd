\name{model.strain_multiplier}
\alias{model.strain_multiplier}
\title{
Estimate R and relative transmissibility of multiple variants of a virus
}
\description{
Estimate R(t) and the relative transmissibility of multiple variants of a virus from case numbers. The model is a SIS model with: a gamma distirubtion for the infectiousness; R(t) following a jump-diffusion process; and variants having a constant relative transmissibility. The models fit is fit to time-series of cases for each variant, with the posteriror disitribution of the parameters being sampled with Stan.
}
\usage{
\code{model.strain_multiplier() }
}
\value{ A Stan model object}
\details{
Input data is provided to the Stan model as a list containing the fields:
\itemize{
\item \code{n_dates} - the number of dates on which cases data is provided
\item \code{n_strains} - the number of variants modelled
\item \code{data_freq} - the frequency of the data (e.g. 7L for weekly)
\item \code{cases} - matrix of cases counts ( \code{n_dates} rows and \code{n_strains} columns)    
\item \code{n_jumps} - number of jumps in the R(t) process
\item \code{jump_times} - times at which jumps are allowed to take place (in days)
\item \code{seed_start_time} - start of seeding periods (integer array of length  \code{n_strains})
\item \code{seed_end_time} - end of seeding periods (integer array of length  \code{n_strains})  
\item \code{generation_mean} - mean generation time of the infection
\item \code{generation_sd} - standard deviation of generation time of the infection
\item \code{generation_max} - maximum length of infections period
\item \code{test_mean} - mean time after infection somebody is tested
\item \code{test_sd} - standard deviation of time after infection somebody is tested
\item \code{ascertainment_factor} - case ascertainment factor
\item \code{prior_R0_min} - prior minimum R(0)
\item \code{prior_R0_max} - prior maximum R(0)
\item \code{prior_dR_sd_min} - prior minimum (absolute) daily change in R(t)
\item \code{prior_dR_sd_max} - prior maximum (absolute) daily change in R(t)
\item \code{prior_phi_od_max} - prior maximum over-dispersion of negative-binmial in observation model
\item \code{prior_jump_up_max} - prior maximum factor of increase in R(t) at a jump
\item \code{prior_jump_up_min} - prior minimum factor of decrease in R(t) at a jump
\item \code{prior_daily_seed_max} - prior on maximum number of daily seeds during the seeding window (array of length \code{n_strains})
\item \code{prior_strain_multipliers_min} - prior minimum of strain multiplier (array of length \code{n_strains-1})
\item \code{prior_strain_multipliers_max} - prior maximum of strain multiplier (array of length \code{n_strains-1})
\item \code{multiply_strain_multipliers} - strain multipliers are relative to the previous strain or the first strain
\item \code{mask_cases_less_than} - mask cases in fit if less than this threshold
}
The sampled parameters are :
\itemize{
\item \code{strain_multipliers} - the relative transmissibility of the strains to the wild-type strains (output is relative to wild-type regardless of \code{multiply_strain_multipliers} - which controls the underlying parameter which is sampled and improves samplin)
\item \code{R} - array containing the daily estimated value of R(t) for the wild-type
\item \code{R_comb} - array containing the daily estimate of R(t) weighted by the number of infections across variants
\item \code{jumps} - array of estimated size of jumps (length \code{n_jumps})
\item \code{daily_seed} - array of estimated daily seeds(length \code{n_strains})
}
}
\examples{
Look at https://github.com/BDI-pathogens/VariantREstimate/blob/main/examples/england_summer_alpha.R.

\code{> library( rstan )}
\code{> model <- model.strain_multiplier() }
\code{> samples <- sampling( model, data = data, pars = pars )}
}
\references{
Hinch R., Panovska-Griffiths J., Probert W., Ferretti L., Wymant C., Di Lauro F., Baya N., Ghafari M., Abeler-Dörner, COG-UK, Fraser C. (2021), Estimating SARS-CoV-2 Variant Fitness and the Impact of Interventions in England using Statistical and Geo-Spatial Agent-Based Models. Phil. Trans Roy Soc A (under review)
}