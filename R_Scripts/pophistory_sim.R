### This script runs a simulation of a population under stabilizing selection
### The population size changes according to the population history of Holstein cattle.

## Load packages
library(simsimqt)
library(ggplot2)
library(purrr)
library(tibble)

## Make population history

# Read population history file
Holstein_pophistory <- read.delim("/home/leonora/pophistory_simulation/Holstein_pophistory.txt")
pop_history <- c()

# Generate population history
for(i in 1:nrow(Holstein_pophistory)) {
  Ne <- Holstein_pophistory[i,1]
  n_gen <- Holstein_pophistory[i,2]
  pop_history <- c(pop_history, rep(Ne, n_gen))
}

### Set up of simulation
# Create 4 replicates of each simulation
for (rep_num in 1:4) {
  
  # Gaussian fitness function for stabilizing selection
  fitness_stabilizing <- function(trait, s, trait_optimum) {
    exp(-0.5 *  s * (trait - trait_optimum)^2)
  }

  # Mutation rate
  mu <- 1e-5
  # Strength of selection parameters
  s <- 0.1
  # Number of loci
  n_loci <- 200

  ## Simulate founders with genotype frequencies drawn from a beta
  ## distribution (to get a U-shape)
  founder_genotypes <- draw_founder_genotypes(n_ind = 10000,
                                              n_loci = 1000,
                                              distribution = "beta",
                                              parameters = list(shape1 = 1/10,
                                                                shape2 = 1/10))
  
  # Draw random causative variants that segregate in founders
  founder_p <- colSums(founder_genotypes)/2/nrow(founder_genotypes)
  candidates <- which(founder_p > 0 & founder_p < 1)
  locus_ix <- sample(candidates, n_loci)

  ## Set up a trait with additive effects drawn from an exponential distribution
  trait <- make_trait_manual(founder_genotypes,
                             loci_ix = locus_ix,
                             a = rexp(n = n_loci, rate = 1/0.1),
                             d = rep(0, n_loci), 
                             Vg = 1,
                             trait_mean = 0)

  ## Set up simulation parameters and population object
  simparam <- SimParam$new(traits = list(trait), Ve = 1)
  population <- new_population(founder_genotypes,
                               simparam)

  

  # Number of generations
  n_gen <- length(pop_history)

  # Data frame to gather results
  results <- tibble(gen = 1:n_gen, mean = numeric(n_gen), var = numeric(n_gen), varPheno = numeric(n_gen), freq = vector(length = n_gen, mode = "list"))

  # Add starting population to results
  results$freq[[1]] <- pull_qtl_freq(population,
                                     simparam = simparam)
  results$mean[1] <- meanG(population)
  results$var[1] <- varG(population)
  results$varPheno[1] <- var(population$pheno)

  
  ## Selection, mating and mutation loop
  for (gen in 2:n_gen) {
  
    if (gen %% 10 == 0) { print(gen) }
  
    offspring <- select_cross_fitness(population,
                                      n_crosses = pop_history[gen],
                                      n_progeny = 1,
                                      fitness_function = fitness_stabilizing,
                                      trait_optimum = 0,
                                      s = s,
                                      simparam = simparam)
  
    ## Get results
    results$freq[[gen]] <- pull_qtl_freq(offspring,
                                       simparam = simparam)
    results$mean[gen] <- meanG(offspring)
    results$var[gen] <- varG(offspring)
    results$varPheno[gen] <- var(offspring$pheno)

  
    ## Add mutations, and save population for next iteration
    population <- mutate_per_locus(offspring,
                                   mutation_rate = mu)
  }


  ## Save results in files
  saveRDS(pop_history, paste("/home/leonora/pophistory_simulation/Replicate_5/", "5", rep_num, "_pop_history.Rda", sep=""))
  saveRDS(population, paste("/home/leonora/pophistory_simulation/Replicate_5/", "5", rep_num, "_population.Rda", sep=""))
  saveRDS(results, paste("/home/leonora/pophistory_simulation/Replicate_5/", "5", rep_num, "_results.Rda", sep=""))
  saveRDS(simparam, paste("/home/leonora/pophistory_simulation/Replicate_5/", "5", rep_num, "_simparam.Rda", sep=""))

}
