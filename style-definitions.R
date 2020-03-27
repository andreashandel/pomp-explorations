# The following colors are useful for state variables, etc.
# They are all a bit transparent, so useful for spaghetti plots, etc.
# A very trasnparent version of each color can be used for confidence interval fill color
# Try to stick to the Red for Infectious, and the blue-green for exposed.
# The other colors may need to shuffle a bit depending on your needs.
# Black, Red and blue-purple stand out the most.

# Depends: Eric Marty's datacolor package
# install.packages("devtools")
# devtools::install_github("allopole/datacolor")
#library(datacolor)

# Base R
## black-ish for cases, fatalities, etc (given data)
col.cases <- rgb(0, 0, 0, 191, "black", 255)
col.cases.ci <- rgb(0, 0, 0, 63,"lt-black", 255) ## confidence intervals, etc.
## Red for Infectious class or similar
col.I <- rgb(230, 7, 7, 191, "red", 255)            
col.I.ci <- rgb(230, 7, 7, 63, "lt-red", 255) ## confidence intervals, etc.
## Bluish Green for Exposed class or similar
col.E <- rgb(7, 164, 181, 191, "green", 255)        
col.E.ci <- rgb(7, 164, 181, 63, "lt-green", 255) ## confidence intervals, etc.
## Blue-Purple for nowcast (E + I) or similar.
## Could be used for hospitalized class
col.nowcast <- rgb(7, 7, 230, 191, "blue", 255)
col.nowcast.ci <- rgb(7, 7, 230, 63, "lt-blue", 255)  ## confidence intervals, etc.
## Magenta for other
col.other <- rgb(164, 0, 181, 191, "magenta", 255)        
col.other.ci <- rgb(164, 0, 181, 63, "lt-magenta", 255)
palette <- c(col.cases, col.cases.ci, col.I, col.I.ci, col.E, col.E.ci, col.nowcast, col.nowcast.ci, col.other, col.other.ci)

#datacolor::colorbar(palette)

# Plotly
## black-ish for cases, fatalities, etc (given data)
# col.cases <- 'rgb(0, 0, 0, .75)'
# col.cases.ci <- 'rgb(0, 0, 0, .25)' ## confidence intervals, etc.
# ## Red for Infectious class or similar
# col.I <- 'rgb(230, 7, 7, .75)'            
# col.I.ci <- 'rgb(230, 7, 7, .25)' ## confidence intervals, etc.
# ## Bluish Green for Exposed class or similar
# col.E <- 'rgb(7, 164, 181, .75)'        
# col.E.ci <- 'rgb(7, 164, 181, .25)' ## confidence intervals, etc.
# ## Blue-Purple for nowcast (E + I) or similar.
# ## Could be used for hospitalized class
# col.nowcast <- 'rgb(7, 7, 230, .75)'
# col.nowcast.ci <- 'rgb(7, 7, 230, .25)'  ## confidence intervals, etc.
# ## Magenta for other
# col.other <- 'rgb(164, 0, 181, .75)'        
# col.other.ci <- 'rgb(164, 0, 181, .25)'
