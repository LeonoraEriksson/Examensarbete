### This script runs a simulation of a population under stabilizing selection
### and a simulation of a population under stabilizing selection with adaptation
### to a shift in optimum. Simulation is run for 1000 generations and
### population size is constant at 120 individuals

## Load packages
library(simsimqt)
library(ggplot2)
library(purrr)
library(tibble)

# Read in population and simulation parameters from previous simulation
population <- readRDS("/home/leonora/pophistory_simulation/Replicate_0/01_population.Rda")
simparam <- readRDS("/home/leonora/pophistory_simulation/Replicate_0/01_simparam.Rda")

## Run 4 replicates of simulation
for (rep_num in 1:4) {
  
  # Fitness function for stabilizing selection
  fitness_stabilizing <- function(trait, s, trait_optimum) {
    exp(-0.5 *  s * (trait - trait_optimum)^2)
  }
  
  
  # Mutation rate
  mu <- 1e-5
  # Strength of selection stabilizing selection
  s <- 0.1
  # Keep the effective population size constant
  n_crosses = length(population$id)
  
  
  ## Adaptation towards a new optimum, stabilizing selection
  
  n_gen_adapt <- 1000
  
  # Data frame to gather results
  generations_adapt <- vector(length = n_gen_adapt,
                              mode = "list")
  
  generations_adapt[[1]] <- population
  
  # Set the optimum shift at two genetic standard deviations away
  # from the current genetic mean
  new_optimum <- meanG(population) +
    2 * sqrt(varG(population))[1]
  
  
  ## Loop for selection, mating and mutation, the
  ## whole population is saved in each each generation.

  for (gen in 2:n_gen_adapt) {
    
    if (gen %% 10 == 0) { print(gen) }
    
    generations_adapt[[gen]] <-
      select_cross_fitness(generations_adapt[[gen - 1]],
                           n_crosses = n_crosses,
                           n_progeny = 1,
                           fitness_function = fitness_stabilizing,
                           trait_optimum = new_optimum,
                           s = s,
                           simparam = simparam)
    
    # Add mutations, and save population for next iteration
    generations_adapt[[gen]] <-
      mutate_per_locus(generations_adapt[[gen]],
                       mutation_rate = mu)
    
  }
  
  
  
  ## Directional selection
  
  # Fitness function directional selection
  fitness_directional <- function(trait, s) {
    
    exp(s * trait)
    
  }
  
  ## Strength of selection directional selection
  s_prim = 1
  
  n_gen_dir = 1000
  
  
  generations_dir <- vector(mode = "list",
                             length = n_gen_dir)
  
  generations_dir[[1]] <- population
  
  ## Same selection loop, but different fitness function
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
  
  saveRDS(generations_adapt, paste("/home/leonora/stabilizing_directional_sim/optimum_shift_simulation/", "4", rep_num, "_stabilizing_pop.Rda", sep=""))
  saveRDS(generations_dir, paste("/home/leonora/stabilizing_directional_sim/directional_simulation/", "4", rep_num, "_directional_pop.Rda", sep=""))
  saveRDS(new_optimum, paste("/home/leonora/stabilizing_directional_sim/optimum_shift_simulation/", "4", rep_num, "_optimum.Rda", sep=""))
  
}
