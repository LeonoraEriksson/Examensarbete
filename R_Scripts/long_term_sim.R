### This script runs long-term simulation of a population under stabilizing selection
### and a simulation of a population under stabilizing selection with adaptation
### to repeated shifts in optimum. Simulation is run for 5000 generations and
### population size is constant at 120 individuals

## Load packages
library(simsimqt)
library(ggplot2)
library(purrr)
library(tibble)

#Read in population and simulation parameters from previous simulation
population <- readRDS("/home/leonora/pophistory_simulation/Replicate_0/01_population.Rda")
simparam <- readRDS("/home/leonora/pophistory_simulation/Replicate_0/01_simparam.Rda")

## Fitness function stabilizing selection
fitness_stabilizing <- function(trait, s, trait_optimum) {
  exp(-0.5 *  s * (trait - trait_optimum)^2)
}


## Fitness function directional selection
fitness_directional <- function(trait, s) {
  exp(s * trait)
}


# Mutation rate
mu <- 1e-5
# Strength of selection parameters stabilizing selection
s <- 0.1
# Strength of selection parameters directional selection
s_prim = 1

# Keep the effective population size constant throughout simulation
n_crosses = length(population$id)

## Run 4 replicates of simulation
for (rep_num in 1:4) {
  
  ## Stabilizing selection with adaptation to repeated shifts in optimum
  
  n_gen_adapt <- 5000
  
  # Data frame to gather results
  generations_adapt <- vector(length = n_gen_adapt,
                            mode = "list")

  generations_adapt[[1]] <- population

  ## Set an initial optimum shift one genetic standard deviation away
  ## from the current genetic mean
  optimum <- meanG(population) +
    1 * sqrt(varG(population))[1]

  ## Create data frame with the new optimum for each shift
  optimum_shifts <- vector(length = n_gen_adapt/100,
                           mode = "list")
  optimum_shifts[[1]] <- optimum
  for (n in 2:(n_gen_adapt/100)){
    optimum_shifts[[n]] <- optimum_shifts[[n-1]] + 1 * sqrt(varG(population))[1]
  }
  
  # Create a data frame where the optimum in each generation is specified
  # Optimum changes every 100 generations
  optimi <- c()
  for(i in optimum_shifts) {
    optimi <- c(optimi, rep(i, 100))
  }


  ## Selection, mating and mutation loop, adaptation to repeated shifts in optimum
  for (gen in 2:n_gen_adapt) {
  
    if (gen %% 10 == 0) { print(gen) }
  
    generations_adapt[[gen]] <-
      select_cross_fitness(generations_adapt[[gen - 1]],
                           n_crosses = n_crosses,
                           n_progeny = 1,
                           fitness_function = fitness_stabilizing,
                           trait_optimum = optimi[gen],
                           s = s,
                           simparam = simparam)
  
    generations_adapt[[gen]] <-
      mutate_per_locus(generations_adapt[[gen]],
                       mutation_rate = mu)
  
  }

  
  
  ## Stabilizing selection with adaptation to a single shift in optimum
  
  n_gen_adapt1 <- 5000

  generations_adapt1 <- vector(length = n_gen_adapt1,
                               mode = "list")
  generations_adapt1[[1]] <- population

  ## Set the optimum shift at one genetic standard deviation away
  ## from the current genetic mean
  optimum <- meanG(population) +
    1 * sqrt(varG(population))[1]


  ## Selection, mating and mutation loop

  print("Adapting to a new optimum, stabilizing selection:")
  for (gen in 2:n_gen_adapt1) {
  
    if (gen %% 10 == 0) { print(gen) }
  
    generations_adapt1[[gen]] <-
      select_cross_fitness(generations_adapt1[[gen - 1]],
                           n_crosses = n_crosses,
                           n_progeny = 1,
                           fitness_function = fitness_stabilizing,
                           trait_optimum = optimum,
                           s = s,
                           simparam = simparam)
  
    generations_adapt1[[gen]] <-
      mutate_per_locus(generations_adapt1[[gen]],
                       mutation_rate = mu)
  
  }



  ## Simulation of a population under directional selection

  n_gen_dir = 5000

  generations_dir <- vector(mode = "list",
                             length = n_gen_dir)
  generations_dir[[1]] <- population

  ## Selection, mating, mutation loop
  for (gen in 2:n_gen_dir) {
  
    if (gen %% 10 == 0) { print(gen) }
  
    generations_dir[[gen]] <-
      select_cross_fitness(generations_dir[[gen - 1]],
                           n_crosses = n_crosses,
                           n_progeny = 1,
                           fitness_function = fitness_directional,                                                                                                                                                   
                           s = s_prim,
                           simparam = simparam)
  
    generations_dir[[gen]] <-
      mutate_per_locus(generations_dir[[gen]],
                       mutation_rate = mu)
  }
  
  
  ## Save results in files
  
  saveRDS(generations_adapt, paste("/home/leonora/long_term_sim/repeated_shifts/", "4", rep_num, "_longstab_optimumshifts_pop.Rda", sep=""))
  saveRDS(generations_adapt1, paste("/home/leonora/long_term_sim/one_shift/", "4", rep_num, "_longstab_pop.Rda", sep=""))
  saveRDS(generations_dir, paste("/home/leonora/long_term_sim/directional/", "4", rep_num, "_longdir_pop.Rda", sep=""))
  saveRDS(optimi, paste("/home/leonora/long_term_sim/repeated_shifts/", "4", rep_num, "_optimi.Rda", sep=""))
  
}
