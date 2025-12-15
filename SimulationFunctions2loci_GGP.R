# List of functions to simulate single locus bi-allelic population genetics model
# Sylvain Glémin 2023
# sylvain.glemin@univ-rennes.fr

library(ggplot2)

# Notation
# Gametes 1:A0B0 2:A0B1 3:A1B0 4:A1B1
# Genotypes noted G11 = A0A0;B0B0 etc...
# Vector of genotypes geno=c(G11,G12,G13,G14,G22,G23,G24,G33,G34,G44)

# G11 A0A0;B0B0   1
# G12 A0A0;B1B0   1 + hb*sb
# G13 A0A1;B0B0   1 + ha*sa   
# G14 A0A1;B0B1   1 + ha*sa + hb*sb + ehh
# G22 A0A0;B1B1   1 + sb
# G23 A0A1;B0B1   1 + ha*sa + hb*sb + ehh
# G24 A0A1;B1B1   1 + ha*sa + sb + ehb
# G33 A1A1;B0B0   1 + sa
# G34 A1A1;B0B1   1 + sa + hb*sb + eah
# G44 A1A1;B1B1   1 + sa + sb + eab



########################### #
# Life cycle functions ######
########################### #



##################### MUTATION  #########################################
#' @title mutation
#' @description Function that return the genotype frequencies after mutation
#' Multiple mutation are not allowed
#' @param geno: a vector of the three genotype frequencies
#' @param mu11: mutation rate from A0 to A1
#' @param mu12: mutation rate from A1 to A0
#' @param mu21: mutation rate from B0 to B1
#' @param mu22: mutation rate from B1 to B0
#' @return the genotype frequencies after mutation
mutation <- function(geno,mu11,mu12,mu21,mu22){
  # Before mutation
  G11 <- geno[1]
  G12 <- geno[2]
  G13 <- geno[3]
  G14 <- geno[4]
  G22 <- geno[5]
  G23 <- geno[6]
  G24 <- geno[7]
  G33 <- geno[8]
  G34 <- geno[9]
  G44 <- geno[10]
  # After mutation
  g11 <- G11*(1-2*mu11-2*mu21) + G12*mu22 + G13*mu12
  g12 <- G12*(1-2*mu11-mu21-mu22) + G11*2*mu21 + G14*mu12 + G22*2*mu22 + G23*mu12
  g13 <- G13*(1-mu11-mu12-2*mu21) + G11*2*mu11 + G14*mu22 + G23*mu22 + G33*2*mu12
  g14 <- G14*(1-mu11-mu12-mu21-mu22) + G12*mu11 + G13*mu21 + G24*mu22 + G34*mu12
  g22 <- G22*(1-2*mu11-2*mu22) + G12*mu21 + G24*mu12
  g23 <- G23*(1-mu11-mu12-mu21-mu22) + G12*mu11 + G13*mu21 + G24*mu22 + G34*mu12
  g24 <- G24*(1-mu11-mu12-2*mu22)  + G14*mu21 + G22*2*mu11 + G23*mu21 + G44*2*mu12
  g33 <- G33*(1-2*mu12-2*mu21) + G13*mu11 + G34*mu22
  g34 <- G34*(1-2*mu12-mu21-mu22) + G14*mu11 + G23*mu11 + G33*2*mu21 + G44*2*mu22
  g44 <- G44*(1-2*mu12-2*mu22) + G24*mu11 + G34*mu21
  newgeno <- c(g11,g12,g13,g14,g22,g23,g24,g33,g34,g44)
  return(newgeno)
}


###################### REPRODUCTION  #######################################
#' @title reproduction
#' @description Function that return the genotype frequencies after reproduction where selfing is allowed
#' @param geno: vector of genotype frequencies
#' @param self: selfing rate
#' @param rec: recombination rate
#' @return the genotype frequencies after reproduction

# Equation from Hedrick 1980 (Genetics 94: 791-808)
reproduction <- function(geno,self,rec) {  
  #Before reproduction
  G11 <- geno[1]
  G12 <- geno[2]
  G13 <- geno[3]
  G14 <- geno[4]
  G22 <- geno[5]
  G23 <- geno[6]
  G24 <- geno[7]
  G33 <- geno[8]
  G34 <- geno[9]
  G44 <- geno[10]
  X1 <- G11 + (G12+G13+G14)/2 - rec*(G14-G23)/2
  X2 <- G22 + (G12+G23+G24)/2 + rec*(G14-G23)/2
  X3 <- G33 + (G13+G23+G34)/2 + rec*(G14-G23)/2
  X4 <- G44 + (G14+G24+G34)/2 - rec*(G14-G23)/2
  # After reproduction
  g11 <- self*(G11+(G12+G13+G23*rec^2+G14*(1-rec)^2)/4) + (1-self)*X1^2
  g22 <- self*(G22+(G12+G24+G14*rec^2+G23*(1-rec)^2)/4) + (1-self)*X2^2
  g33 <- self*(G33+(G13+G34+G14*rec^2+G23*(1-rec)^2)/4) + (1-self)*X3^2
  g44 <- self*(G44+(G24+G34+G23*rec^2+G14*(1-rec)^2)/4) + (1-self)*X4^2
  g12 <- self*(G12+rec*(1-rec)*(G14+G23))/2 + (1-self)*2*X1*X2
  g13 <- self*(G13+rec*(1-rec)*(G14+G23))/2 + (1-self)*2*X1*X3
  g24 <- self*(G24+rec*(1-rec)*(G14+G23))/2 + (1-self)*2*X2*X4
  g34 <- self*(G34+rec*(1-rec)*(G14+G23))/2 + (1-self)*2*X3*X4
  g14 <- self*(G14*(1-rec)^2+G23*rec^2)/2 + (1-self)*2*X1*X4
  g23 <- self*(G23*(1-rec)^2+G14*rec^2)/2 + (1-self)*2*X2*X3
  newgeno <- c(g11,g12,g13,g14,g22,g23,g24,g33,g34,g44)
  return(newgeno)
}



##################### SELECTION   #########################################
#' @title selection
#' @description Function that return the genotype frequencies after selection
#' @param geno: vector of the three genotype frequencies
#' @param fit: vector of fitness
#' @return the genotype frequencies after selection

selection <- function(geno,fit) {
  geno.sel <- geno * fit
  return(geno.sel/sum(geno.sel))
}

# Function to set fitness matrix from selecion,dominance and epistasis coefficients#
# The fitness matrix can also be set manually #
fitness <- function(ha,sa,hb,sb,ehh,eah,ehb,eab){
  c(1,
    1 + hb*sb,
    1 + ha*sa,   
    1 + ha*sa + hb*sb + ehh,
    1 + sb,
    1 + ha*sa + hb*sb + ehh,
    1 + ha*sa + sb + ehb,
    1 + sa,
    1 + sa + hb*sb + eah,
    1 + sa + sb + eab  
  )
}



##################### DRIFT   #########################################
#' @title drift
#' @description Function that return the genotype frequencies after drift
#' We assume that Ne = N.
#' @param geno: vector of the three genotype frequencies
#' @param Npop: poopulation size
#' @return the genotype frequencies after drift
drift <- function(geno,Npop){
  newgeno <- rmultinom(n=1,size=Npop,prob = geno)/Npop
  return(newgeno)
}




############################### #
# Examples of simulations ######
############################# #


######################## ONE GENERATION ##################
#' @title generation
#' @description Function that return the genotype frequencies after one generation
#' @param geno: vector of the three genotype frequencies
#' @param mu11: mutation rate from A0 to A1
#' @param mu12: mutation rate from A1 to A0
#' @param mu21: mutation rate from B0 to B1
#' @param mu22: mutation rate from B1 to B0
#' @param self: selfing rate
#' @param rec: recombination rate
#' @param fit: vector of fitness
#' @param Npop: population size
#' @return the genotype frequencies after one generation
generation <- function(geno,Npop,self,rec,mu11,mu12,mu21,mu22,fit){
  g.mut <- mutation(geno,mu11,mu12,mu21,mu22)
  g.repro <- reproduction(g.mut,self,rec)
  g.sel <- selection(g.repro,fit)
  newgeno <- drift(g.sel,Npop)
  return(newgeno)
}



######################## ONE SIMULATION ##################
# Exemple 1: simulation of one population through time

#' @title simulation
#' @description Function that return vectors of allelic and genotypic frequencies thought time
#' @param geno: vector of the three genotype frequencies
#' @param mu11: mutation rate from A0 to A1
#' @param mu12: mutation rate from A1 to A0
#' @param mu21: mutation rate from B0 to B1
#' @param mu22: mutation rate from B1 to B0
#' @param self: selfing rate
#' @param rec: recombination rate
#' @param fit: vector of fitness
#' @param Npop: population size
#' @param Tmax: the total number of generations to simulate
#' @return list of: frequencies of allele a th
simulation <- function(geno0,Npop,self,rec,mu11,mu12,mu21,mu22,fit,Tmax){
  geno <- geno0
  tab.a <- c()
  tab.b <- c()
  tab.LD <- c()
  for(t in 1:Tmax){
    newgeno <- generation(geno,Npop,self,rec,mu11,mu12,mu21,mu22,fit)
    # Genotype frequencies
    G11 <- newgeno[1]
    G12 <- newgeno[2]
    G13 <- newgeno[3]
    G14 <- newgeno[4]
    G22 <- newgeno[5]
    G23 <- newgeno[6]
    G24 <- newgeno[7]
    G33 <- newgeno[8]
    G34 <- newgeno[9]
    G44 <- newgeno[10]
    # Haplotype frequencies
    X1 <- G11 + (G12+G13+G14)/2 - rec*(G14-G23)/2
    X2 <- G22 + (G12+G23+G24)/2 + rec*(G14-G23)/2
    X3 <- G33 + (G13+G23+G34)/2 + rec*(G14-G23)/2
    X4 <- G44 + (G14+G24+G34)/2 - rec*(G14-G23)/2
    # Allele frequencies
    fa <- X3 + X4
    fb <- X2 + X4
    LD <- X1*X4 - X2*X3
    tab.a <- c(tab.a,fa)
    tab.b <- c(tab.b,fb)
    tab.LD <- c(tab.LD,LD)
    geno <- newgeno
  }
  return(list("a"=tab.a,"b"=tab.b,"LD"=tab.LD))
}


# Vector of genotypes geno=c(G11,G12,G13,G14,G22,G23,G24,G33,G34,G44)

# G11 A0A0;B0B0   1
# G12 A0A0;B1B0   1 + hb*sb
# G13 A0A1;B0B0   1 + ha*sa   
# G14 A0A1;B0B1   1 + ha*sa + hb*sb + ehh
# G22 A0A0;B1B1   1 + sb
# G23 A0A1;B0B1   1 + ha*sa + hb*sb + ehh
# G24 A0A1;B1B1   1 + ha*sa + sb + ehb
# G33 A1A1;B0B0   1 + sa
# G34 A1A1;B0B1   1 + sa + hb*sb + eah
# G44 A1A1;B1B1   1 + sa + sb + eab

# Example
# set.seed(123)

Npop = 1000   # population size
nb = 40      # number of s values considered
rec = 0.01   # recombination rate

# 1 example of temporal evolution of neutral allele frequency

fit.matrix <- fitness(0.5,0.1,0.5,0,0,0,0,0) # Two loci with additive selection
sim <- simulation(geno0 = c(0.25*0.995,0,0.5*0.995,0.5*0.005,0,0,0,0.25*0.995,0.25*0.005,0),Npop = Npop,self = 0,rec=rec,mu11 = 0,mu12 = 0,mu21 = 0,mu22=0,fit = fit.matrix,Tmax = 400)
plot(NULL,xlim = c(0,200),ylim=c(-0.1,1), xlab = "Generations", ylab = "Allele frequencies")
lines(sim$a,col="darkblue")
lines(sim$b,col="orange")
lines(sim$LD,col="black")
abline(0,0,lty=2)

# Fixation probability graph for different values of s

freq <- function(nb,rec) {
  s = c(1:nb)/(5*nb)
  final = numeric(length(s))
  for (i in (1:length(s))) {
    results =logical(10)
    for (repetition in (1:10)) {
      fit.matrix <- fitness(0.5,s[i],0.5,0,0,0,0,0) # Two loci with additive selection
      sim <- simulation(geno0 = c(0.25*0.995,0,0.5*0.995,1*0.005,0,0,0,0.25*0.995,0,0),Npop = Npop,self = 0,rec=rec,mu11 = 0,mu12 = 0,mu21 = 0,mu22=0,fit = fit.matrix,Tmax = 400)
      results[repetition] <- (sim$a[400] >= 0.95)
    }
    final[i] <- mean(results)
  }
  return(final)
}

# rapport = s/rec
s = c(1:nb)/(5*nb)
final = freq(nb = nb,rec = rec)
smoothingSpline = smooth.spline(x = s, y = final, spar=0.35)
plot(s,final,  xlab="Selection coefficient", ylab="Fixation probability",xlim = c(0,0.2))
lines(smoothingSpline)
title(main="Neutral allele fixation depending on intensity of selection")


