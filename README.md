# VariantREstimate
This is a R package to estimate R(t) and the relative transmissibility of multiple variants of a virus from case numbers,
The model is a SIS model with: a gamma distirubtion for the infectiousness; R(t) following a jump-diffusion process; and variants having ta constant relative transmissibility [[1]](#1).
The models fit is fit to time-series of cases for each variant, with the posteriror disitribution of the parameters being sampled with Stan.

## Example
An example of a model fit is contained in ```examples\england_summer_alpha.R```, which takes data for SARS-Cov2 infections in England from the UK Dashboard and COG-UK between September 2020 and Januaray 2021.
The model contains 3 variants of the virus  (wild-type, B.1.177, Alpha) and estimates the relative transmissibility of the variant to one another.


## References
<a id="1">[1]</a> 
Hinch R., Panovska-Griffiths J., Probert W., Ferretti L., Wymant C., Di Lauro F., Baya N., Ghafari M., Abeler-Dörner, COG-UK, Fraser C. (2021), Estimating SARS-CoV-2 Variant Fitness and the Impact of Interventions in England using Statistical and Geo-Spatial Agent-Based Models. Phil. Trans Roy Soc A (under review)

