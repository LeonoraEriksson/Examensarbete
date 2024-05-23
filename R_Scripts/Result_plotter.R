### Plot results from all simulations

## Load packages
library(simsimqt)
library(ggplot2)
library(purrr)
library(tibble)

## Create list of indexes for result files
indexes <- list('01', '02', '03', '04', '11', '12', '13', '14', '21', '22', '23', '24', '31', '32', '33', '34', '41', '42', '43', '44')



## Create plots for results of simulation of stabilizing selection through population history
for (i in 1:length(indexes)) {
  #Read result file
  results <- readRDS(paste("./Examensarbete/Stabilizing_populationhistory/replicates_s0.1/", indexes[i], "_results.Rda", sep=""))
  # Plot and save results as png
  x <- 1:length(results$mean)
  png(file = paste("./Examensarbete/Stabilizing_populationhistory/Plots/Pophistory_Gmean", indexes[i], ".png", sep=""))
  plot(x, results$mean, type = 'l', xlab = 'Generation', ylab = 'Genetic mean', main = paste("Genetic mean through population history", indexes[i], sep = " "), cex.lab=1.5, cex.axis=1.5)
  dev.off()
  png(file = paste("./Examensarbete/Stabilizing_populationhistory/Plots/Pophistory_Gvar", indexes[i], ".png", sep=""))
  plot(x, results$var, type = 'l', xlab = 'Generation', ylab = 'Genetic variance', main = paste("Genetic variance through population history", indexes[i], sep = " "), cex.lab=1.5, cex.axis=1.5)
  dev.off()
}



## Create plots for results of population history with burn-in period
results <- readRDS("./Examensarbete/Stabilizing_populationhistory/Burnin/Burnin1_results.Rda")
x <- 1:length(results$mean)
png(file= "./Examensarbete/Stabilizing_populationhistory/Burnin/Burnin_Gmean.png")
plot(x, results$mean, type = 'l', xlab = 'Generation', ylab = 'Genetic mean', main = "Genetic mean through population history")
dev.off()
png(file = "./Examensarbete/Stabilizing_populationhistory/Burnin/Burnin_Gvar.png")
plot(x, results$var, type = 'l', xlab = 'Generation', ylab = 'Genetic variance', main = "Genetic variance through population history")
dev.off()



## Create plots for simulation of directional selection and stabilizing selection with a shift in optimum

# Directional selection
for (i in 1:length(indexes)) {
  results <- readRDS(paste("./Examensarbete/Stabilizing_directional_sim/Directional_sim/", indexes[i], "_directional_pop.Rda", sep=""))
  png(file = paste("./Examensarbete/Stabilizing_directional_sim/Plots/Directional/Directional_Gmean", indexes[i], ".png", sep=""))
  plot(map_dbl(results, meanG), type = 'l', xlab = 'Generation', ylab = 'Genetic mean', main = paste("Genetic mean directional selection", indexes[i], sep = " "), cex.lab=1.5, cex.axis=1.5)
  dev.off()
  png(file = paste("./Examensarbete/Stabilizing_directional_sim/Plots/Directional/Directional_Gvar", indexes[i], ".png", sep=""))
  plot(map_dbl(results, varG), type = 'l', xlab = 'Generation', ylab = 'Genetic variance', main = paste("Genetic variance directional selection", indexes[i], sep = " "), cex.lab=1.5, cex.axis=1.5)
  dev.off()
}

# Stabilizing selection with a shift in optimum
for (i in 1:length(indexes)) {
  results <- readRDS(paste("./Examensarbete/Stabilizing_directional_sim/Optimum_shift_sim/", indexes[i], "_stabilizing_pop.Rda", sep=""))
  optimum <- readRDS(paste("./Examensarbete/Stabilizing_directional_sim/Optimum_shift_sim/", indexes[i], "_optimum.Rda", sep=""))
  png(file = paste("./Examensarbete/Stabilizing_directional_sim/Plots/Optimum_shift/Stabilizing_Gmean", indexes[i], ".png", sep=""))
  plot(map_dbl(results, meanG), type = 'l', xlab = 'Generation', ylab = 'Genetic mean', main = paste("Genetic mean stabilizing selection with a shift in optimum", indexes[i], sep = " "), cex.lab=1.5, cex.axis=1.5)
  abline(h = optimum, col = 'red')
  dev.off()
  png(file = paste("./Examensarbete/Stabilizing_directional_sim/Plots/Optimum_shift/Stabilizing_Gvar", indexes[i], ".png", sep=""))
  plot(map_dbl(results, varG), type = 'l', xlab = 'Generation', ylab = 'Genetic variance', main = paste("Genetic variance stabilizing selection with a shift in optimum", indexes[i], sep = " "), cex.lab=1.5, cex.axis=1.5)
  dev.off()
}



## Create plots for long-term simulation of directional and stabilizing selection

# Directional selection
for (i in 1:length(indexes)) {
  results <- readRDS(paste("./Examensarbete/Long_term_sim/Directional/", indexes[i], "_longdir_pop.Rda", sep=""))
  png(file = paste("./Examensarbete/Long_term_sim/Plots/Directional/Long_term_directional_Gmean", indexes[i], ".png", sep=""))
  plot(map_dbl(results, meanG), type = 'l', xlab = 'Generation', ylab = 'Genetic mean', main = paste("Genetic mean directional selection", indexes[i], sep = " "), cex.lab=1.5, cex.axis=1.5)
  dev.off()
  png(file = paste("./Examensarbete/Long_term_sim/Plots/Directional/Directional_Gvar", indexes[i], ".png", sep=""))
  plot(map_dbl(results, varG), type = 'l', xlab = 'Generation', ylab = 'Genetic variance', main = paste("Genetic variance directional selection", indexes[i], sep = " "), cex.lab=1.5, cex.axis=1.5)
  dev.off()

  # Retrieve allele frequencies and create plot
  allele_freqs <- data.frame()
  for (population in results) {
    allele_freqs <- rbind(allele_freqs, pull_qtl_freq(population, 1, simparam))
  }
  png(file = paste("./Examensarbete/Long_term_sim/Plots/Directional/Long_term_directional_Afreq", indexes[i], ".png", sep=""), width = 697, height = 345)
  matplot(1:length(results), allele_freqs, lty = 1, type = 'l', col = 'black', cex.lab=1.5, cex.axis=1.5, xlab = 'Generation', ylab = 'Allele frequency', main = paste('Allele frequencies directional selection', indexes[i], sep = ' '))
  dev.off()

  # Identify locus with small and large effect sizes, and plot
  a <- simparam$traits[[1]]$a
  small <- which(a<0.1)
  large <- which(a>0.1)
  png(file = paste("./Examensarbete/Long_term_sim/Plots/Directional/Long_term_directional_Asmall", indexes[i], ".png", sep=""), width = 697, height = 345)
  matplot(1:length(results), allele_freqs[ , small], type = 'l', lty = 1, lwd = 0.5, col = 'black', cex.lab=1.5, cex.axis=1.5, xlab = 'Generation', ylab = 'Allele frequency', main = paste('Allele frequencies directional selection small', indexes[i], sep = ' '))
  dev.off()
  png(file = paste("./Examensarbete/Long_term_sim/Plots/Directional/Long_term_directional_Alarge", indexes[i], ".png", sep=""), width = 697, height = 345)
  matplot(1:length(results), allele_freqs[ , large], type = 'l', lty = 1, lwd = 0.5, col = 'black', cex.lab=1.5, cex.axis=1.5, xlab = 'Generation', ylab = 'Allele frequency', main = paste('Allele frequencies directional selection large', indexes[i], sep = ' '))
  dev.off()
}

 # Stabilizing selection with repeated shifts in optimum
for (i in 1:length(indexes)) {
  results <- readRDS(paste("./Examensarbete/Long_term_sim/Repeated_shifts/", indexes[i], "_longstab_optimumshifts_pop.Rda", sep=""))
  optimi <- readRDS(paste("./Examensarbete/Long_term_sim/Repeated_shifts/", indexes[i], "_optimi.Rda", sep=""))
  png(file = paste("./Examensarbete/Long_term_sim/Plots/Repeated_shifts/Repeated_shifts_Gmean", indexes[i], ".png", sep=""))
  plot(map_dbl(results, meanG), type = 'l', xlab = 'Generation', ylab = 'Genetic mean', main = paste("Genetic mean stabilizing selection repeated shifts", indexes[i], sep = " "), cex.lab=1.5, cex.axis=1.5)
  lines(1:(length(optimi)), optimi, type = 'l', col = 'red')
  dev.off()
  png(file = paste("./Examensarbete/Long_term_sim/Plots/Repeated_shifts/Repeated_shifts_Gvar", indexes[i], ".png", sep=""))
  plot(map_dbl(results, varG), type = 'l', xlab = 'Generation', ylab = 'Genetic variance', main = paste("Genetic variance stabilizing selection repeated shifts", indexes[i], sep = " "), cex.lab=1.5, cex.axis=1.5)
  dev.off()
  
  allele_freqs <- data.frame()
  for (population in results) {
    allele_freqs <- rbind(allele_freqs, pull_qtl_freq(population, 1, simparam))
  }
  png(file = paste("./Examensarbete/Long_term_sim/Plots/Repeated_shifts/Repeated_shifts_A_freqs", indexes[i], ".png", sep=""), width = 697, height = 345)
  matplot(1:length(results), allele_freqs, lty = 1, type = 'l', col = 'black', cex.lab=1.5, cex.axis=1.5, xlab = 'Generation', ylab = 'Allele frequency', main = paste('Allele frequencies repeated shifts', indexes[i], sep = ' '))
  dev.off()

  a <- simparam$traits[[1]]$a
  small <- which(a<0.1)
  large <- which(a>0.1)
  png(file = paste("./Examensarbete/Long_term_sim/Plots/Repeated_shifts/Repeated_shifts_A_small", indexes[i], ".png", sep=""), width = 697, height = 345)
  matplot(1:length(results), allele_freqs[ , small], type = 'l', lty = 1, lwd = 0.5, col = 'black', cex.lab=1.5, cex.axis=1.5, xlab = 'Generation', ylab = 'Allele frequency', main = paste('Allele frequencies repeated shifts, small', indexes[i], sep = ' '))
  dev.off()
  png(file = paste("./Examensarbete/Long_term_sim/Plots/Repeated_shifts/Repeated_shifts_A_large", indexes[i], ".png", sep=""), width = 697, height = 345)
  matplot(1:length(results), allele_freqs[ , large], type = 'l', lty = 1, lwd = 0.5, col = 'black', cex.lab=1.5, cex.axis=1.5, xlab = 'Generation', ylab = 'Allele frequency', main = paste('Allele frequencies repeated shifts, large', indexes[i], sep = ' '))
  dev.off()
}

# Stabilizing selection with a single shift in optimum
for (i in 1:length(indexes)) {
  results <- readRDS(paste("./Examensarbete/Long_term_sim/Stabilizing/", indexes[i], "_longstab_pop.Rda", sep=""))
  optimi <- readRDS(paste("./Examensarbete/Long_term_sim/Repeated_shifts/", indexes[i], "_optimi.Rda", sep=""))
  png(file = paste("./Examensarbete/Long_term_sim/Plots/One_shift/Stabilizing_longterm_Gmean", indexes[i], ".png", sep=""))
  plot(map_dbl(results, meanG), type = 'l', xlab = 'Generation', ylab = 'Genetic mean', main = paste("Genetic mean stabilizing selection", indexes[i], sep = " "), cex.lab=1.5, cex.axis=1.5)
  abline(h=optimi[1], col = 'red')
  dev.off()
  png(file = paste("./Examensarbete/Long_term_sim/Plots/One_shift/Stabilizing_longterm_Gvar", indexes[i], ".png", sep=""))
  plot(map_dbl(results, varG), type = 'l', xlab = 'Generation', ylab = 'Genetic variance', main = paste("Genetic variance stabilizing selection", indexes[i], sep = " "), cex.lab=1.5, cex.axis=1.5)
  dev.off()
  
  allele_freqs <- data.frame()
  for (population in results) {
    allele_freqs <- rbind(allele_freqs, pull_qtl_freq(population, 1, simparam))
  }
  png(file = paste("./Examensarbete/Long_term_sim/Plots/One_shift/Stabilizing_longterm_Afreqs", indexes[i], ".png", sep=""), width = 697, height = 345)
  matplot(1:length(results), allele_freqs, lty = 1, type = 'l', col = 'black', cex.lab=1.5, cex.axis=1.5, xlab = 'Generation', ylab = 'Allele frequency', main = paste('Allele frequencies one shift', indexes[i], sep = ' '))
  dev.off()

  a <- simparam$traits[[1]]$a
  small <- which(a<0.1)
  large <- which(a>0.1)
  png(file = paste("./Examensarbete/Long_term_sim/Plots/One_shift/Stabilizing_longterm_Asmall", indexes[i], ".png", sep=""), width = 697, height = 345)
  matplot(1:length(results), allele_freqs[ , small], type = 'l', lty = 1, lwd = 0.5, col = 'black', cex.lab=1.5, cex.axis=1.5, xlab = 'Generation', ylab = 'Allele frequency', main = paste('Allele frequencies one shift, small', indexes[i], sep = ' '))
  dev.off()
  png(file = paste("./Examensarbete/Long_term_sim/Plots/One_shift/Stabilizing_longterm_Asmalarge", indexes[i], ".png", sep=""), width = 697, height = 345)
  matplot(1:length(results), allele_freqs[ , large], type = 'l', lty = 1, lwd = 0.5, col = 'black', cex.lab=1.5, cex.axis=1.5, xlab = 'Generation', ylab = 'Allele frequency', main = paste('Allele frequencies one shifts, large', indexes[i], sep = ' '))
  dev.off()
}



## Create plots for results of stabilizing selection throughout population history simulations, with heredity
indexes = list(51, 52, 53, 54)
for (i in 1:length(indexes)) {
  results <- readRDS(paste("./Examensarbete/Stabilizing_populationhistory/Heredity/", indexes[i], "_results.Rda", sep=""))
  x <- 1:length(results$mean)
  png(file = paste("./Examensarbete/Stabilizing_populationhistory/Heredity/Plots/Pophistory_Gmean", indexes[i], ".png", sep=""))
  plot(x, results$mean, type = 'l', xlab = 'Generation', ylab = 'Genetic mean', main = paste("Genetic mean through population history", indexes[i], sep = " "), cex.lab=1.5, cex.axis=1.5)
  dev.off()
  png(file = paste("./Examensarbete/Stabilizing_populationhistory/Heredity/Plots/Pophistory_Gvar", indexes[i], ".png", sep=""))
  plot(x, results$var, type = 'l', xlab = 'Generation', ylab = 'Genetic variance', main = paste("Genetic variance through population history", indexes[i], sep = " "), cex.lab=1.5, cex.axis=1.5)
  dev.off()
  png(file = paste("./Examensarbete/Stabilizing_populationhistory/Heredity/Plots/Pophistory_her", indexes[i], ".png", sep=""))
  plot(x, (results$var/results$varPheno), type = 'l', xlab = 'Generation', ylab = 'Heredity', main = paste("Heredity through population history", indexes[i], sep = " "), cex.lab=1.5, cex.axis=1.5)
  dev.off()
}
