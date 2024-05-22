### This script runs a simulation of directional selection for different strengths of selection and calculates
### rate of change in genetic mean expressed in genetic standard deviations.

## Load packages
library(simsimqt)
library(ggplot2)
library(purrr)
library(tibble)


# Read in population and simulation parameters from previous simulations
population <- readRDS("/Users/Admin/Documents/Examensarbete/Stabilizing_populationhistory/replicates_s0.1/01_population.Rda")
simparam <- readRDS( "/Users/Admin/Documents/Examensarbete/Stabilizing_populationhistory/replicates_s0.1/01_simparam.Rda")

## Mutation rate
mu <- 1e-5
## Keep the effective population size constant
n_crosses = length(population$id)



## Directional selection

# Fitness function directional selection
fitness_directional <- function(trait, s) {
  
  exp(s * trait)
  
}

# Strength of selection parameters directional selection
s_prim <- list(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)

n_gen = 12
generations_dir <- vector(mode = "list",
                          length = n_gen)

generations_dir[[1]] <- population

# Results of calculated rate of change
num_sd <- c()

# Selection mutation loop, each s is run 10 times trough the simulation and the mean rate of change is calculated for each s.
for (s in s_prim) {
  for (i in 1:10) {
    for (gen in 2:n_gen) {
  
      generations_dir[[gen]] <-
        select_cross_fitness(generations_dir[[gen - 1]],
                             n_crosses = n_crosses,
                             n_progeny = 1,
                             fitness_function = fitness_directional,
                             s = s,
                             simparam = simparam)
  
      generations_dir[[gen]] <-
       mutate_per_locus(generations_dir[[gen]],
                        mutation_rate = mu)
    }
    
    # Calculate mean rate of change per generation
    change <- (meanG(generations_dir[[n_gen]])-meanG(generations_dir[[1]]))/n_gen
    # Express in terms of number of genetic standard deviations
    num_sd[i] <- change/sqrt(varG(population))[1]
    
  }
  
  # Print mean rate of change for each s
  print(paste('s=', s, sep=' ')) 
  print(mean(num_sd))
}
