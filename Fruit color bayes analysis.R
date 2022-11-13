#Bayesian hierarchical modelling of visits in function of fruit color and fruit size
library(rstan)
library(ape)

dados <- read.table("Full interaction data.txt", h=T)
dados <- subset(dados, Community == "Argentina") #Selecting only Argentina 
dados$Fruit_color <- as.factor(dados$Fruit_color)


#Loading and prunning the phylogenetic trees
##Plants
pltphy <- read.nexus("planttreeBIN.nexus")
plt_sp <- unique(dados$Plantsp) 
#Prunning the tree
pltphy <- drop.tip(pltphy, pltphy$tip.label[-match(plt_sp, pltphy$tip.label)])

##Birds
brdphy <- read.nexus("birdtreeBIN.nexus")
brd_sp <- unique(dados$Bird_phylo)
#Prunning the tree
brdphy <- drop.tip(brdphy, brdphy$tip.label[-match(brd_sp, brdphy$tip.label)])

#Correlation matrix of the phylogenetic distances (cophenetic matrix) - plants

coP <- cophenetic(pltphy)

tmp <- dimnames(coP)
idsp <- as.numeric(as.factor(tmp[[1]]))

DistP <- matrix(NA, ncol(coP), ncol(coP))
for(i in 1:ncol(coP)){
  for(j in 1:ncol(coP)){
    DistP[idsp[i],idsp[j]] <- coP[i,j]
  }
}

#Corelation matrix for the birds phylogenetic distances (varcov matrix)
CC = vcv(brdphy, corr = TRUE)
tmp <- dimnames(CC)
ids <- as.numeric(as.factor(tmp[[1]]))

C <- matrix(NA, ncol(CC), ncol(CC))
for(i in 1:ncol(CC)){
  for(j in 1:ncol(CC)){
    C[ids[i], ids[j]] <- CC[i,j]
  }
}

plants <- dimnames(coP)[[1]]
birds <- dimnames(CC)[[1]]
N = length(plants)  # number of plant spp
N_J = length(birds) # number of bird spp
visits = matrix(0, N_J, N)


for(i in 1:N_J){
  for(j in 1:N){
    tmp = which(dados$Plantsp == plants[j] & dados$Birdsp == birds[i])
    if(length(tmp) > 0){
      visits[i,j] = dados$Visits[tmp]
    }
  }
}

colores = unique(dados$Fruit_color)

K = length(colores) + 1       # number of predictive variables

X = matrix(0, N, K)

for(i in 1:N){
  tmp = dados$Fruit_color[which(dados$Plantsp == plants[i])]
  tmp1 = dados$Fruit_diam_mm[which(dados$Plantsp == plants[i])]
  X[i, which(colores == tmp[1])] = 1
  X[i,K] = tmp1[1]
}

X[,K] = scale(X[,K])

# stan_dat <- list(
#   K = K,
#   N = N,
#   D = D,
#   y = visits,
#   x = X
# )
# 
# library(cmdstanr)
# 
# mod <- cmdstan_model('multinomial_reg.stan')
# 
# fit <- mod$sample(
#   data = stan_dat, 
#   seed = 123, 
#   chains = 4, 
#   parallel_chains = 4,
#   iter_warmup = 1000,
#   iter_sampling = 1000,
#   thin = 1,
#   refresh = 200 # print update every 500 iters
# )
# 
# fit_summary = fit$summary()

# betas are in covariates by bird spp

#------------------------------------------------------------------------------

Xvar <- model.matrix(~ Fruit_color + scale(Fruit_diam_mm), data = dados) #catergorical predictors as dummy variables + fruit diam
brdspp <- dados$Bird_phylo #bird species
brdsp <- sort(unique(brdspp)) 
pltspp <- dados$Plantsp #plant species
pltsp <- unique(pltspp)

np <- dim(Xvar)[2] #number of parameters (fixed effects)
nt <- 2 #number of bird traits to model the intercept
npp <- 1 #(?)
ntp <- 1 #(?)

nbrdsp <- length(brdsp) #number of bird species
npltsp <- length(pltsp) #number of plant species

n_obs <- nrow(dados) #number of observations

#Matrix of bird traits
TT <- matrix(1, nrow = nbrdsp, ncol = nt+1)
for (i in 1:nbrdsp) {
  TT[i, 2] <- dados$Fruit_diet[brdspp == brdsp[i]][1]
  TT[i, 3] <- log(dados$Body_mass[brdspp == brdsp[i]][1])
}

fd <- TT[,2]
bm <- TT[,3]

TT[,2] <- scale(TT[,2])
TT[,3] <- scale(TT[,3])


jj <- as.numeric(as.factor(dados$Bird_phylo)) #turning the brd spp into factor
J <- max(jj)

J_1 <- as.numeric(as.factor(dados$Plantsp)) #turning the plt spp into factor
N_1 <- max(J_1)

L <- ncol(TT) #number of intercepts based on the n of brd spp

ones <- numeric(J) + 1 #vector of ones for the diagonal of phylo distance

mod <- cmdstan_model('multinomial_reg_JSM.stan')

stan_dat <- list(
  N = N,
  N_J = N_J,
  L_J = L,
  K = K,
  y = visits,
  x = X,
  TT = TT,
  C = C,
  ones = ones
)

fit <- mod$sample(
  data = stan_dat, 
  #seed = 123, 
  chains = 4, 
  parallel_chains = 4,
  iter_warmup = 50000,
  iter_sampling = 50000,
  thin = 50,
  refresh = 1000 # print update every 500 iters
)

fit_summary = fit$summary()

betas <- fit_summary[grepl("b_m", fit_summary$variable),]
zetas <- fit_summary[grepl("Z", fit_summary$variable),]
plot(betas$mean)
