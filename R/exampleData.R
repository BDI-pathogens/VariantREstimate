###################################################################################/
# example.data.case_cog_summer_alpha
#
###################################################################################/
example.data.case_cog_summer_alpha = function()
{
  file <- "data/cog_dashboard_combined_1.csv"
  file <- system.file( file, package = "VariantREstimate" )
  return( fread( file ) )
}