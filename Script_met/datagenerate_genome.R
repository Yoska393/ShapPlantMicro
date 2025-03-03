library(MASS)
library(rrBLUP)

vec<- c()
l <- 200 # the number of SNPs
n <- 100 # the number of samples
h2 <- 0.8 # heritability
X <- matrix(sample(c(-1, 0, 1), size = n * l, replace = T), 
            nrow = n, ncol = l)

# Genetic relationship matrix
amat <- A.mat(X = X)
# image(amat)

u <- mvrnorm(n = 1, mu = rep(0, n), Sigma = amat)
varE <- (var(u) / h2) - var(u)
# heritability
print(var(u) / (var(u) + varE))
e <- rnorm(n = n, mean = 0, sd = sqrt(varE))
e <- mvrnorm(n = 1, mu = rep(0, n), 
             Sigma = diag(sqrt(varE), nrow = n))
var(e)

y <- u + e
model <- mixed.solve(y, K = amat)
Gh2<-model$Vu / (model$Vu + model$Ve)
print(Gh2)
vec<-c(vec,Gh2)

#######################################
# Using real genome data
amatPath <- "../Desktop/amat.csv"
amat <- read.csv(file = amatPath, row.names = 1)
amat <- as.matrix(amat)
# image(amat)

n <- nrow(amat)
u <- mvrnorm(n = 1, mu = rep(0, n), Sigma = amat)
varE <- (var(u) / h2) - var(u)
# heritability
# print(var(u) / (var(u) + varE))
e <- rnorm(n = n, mean = 0, sd = sqrt(varE))

y <- u + e
model <- mixed.solve(y, K = amat)
model$Vu / (model$Vu + model$Ve)


###################

library(MASS)
library(rrBLUP)


# Generate X

n <- 100 #Samples
a <- 20
a2<- 20
b <- 20
c <- 20
d <- 40

A  <-  matrix(rnorm(n = n*a, mean = 0, sd = 1),n,a)
#A2 <- matrix(rnorm(n = n*a2, mean = 0, sd = 1 ),n,a)
D  <-  matrix(rnorm(n = n*d, mean = 0, sd = 1),n,d)


amat.a  <-  A.mat(X = A)
#amat.a2 <-  A.mat(X = A2)
amat.d  <-  A.mat(X = D)
#image(amat.a)

# herit
h2a <- 0.3
#h2a2 <- 0.3
h2d <- 0.3
h2e <- 0.4

u_a <- mvrnorm(n = b, mu = rep(0, n), Sigma = h2a * amat.a)
u_d <- mvrnorm(n = b, mu = rep(0, n), Sigma = h2d * amat.d)

varA <- apply(u_a,1,var)
varE <- varA * (h2e / h2a)

e <- mvrnorm(n = b, mu = rep(0, n), Sigma = diag(sqrt(varE), nrow = n))

B <- scale(t(u_a + u_d + e))
#print(varA / apply(B,2,var))

amat.b <-  A.mat(X = B)

#check 
model <- mixed.solve(B[,1], K = amat.a)
Gh2<-model$Vu / (model$Vu + model$Ve)
print(Gh2)


# generate C

h2ac <- 0.2
#h2a2c <- 0.2
h2bc <- 0.3
h2dc <- 0.2
h2ec <- 0.3

#vec<-c()
#for (i in 1:100){

u_a <- mvrnorm(n = b, mu = rep(0, n), Sigma = h2ac * amat.a)
u_d <- mvrnorm(n = b, mu = rep(0, n), Sigma = h2dc * amat.d)
u_b <- mvrnorm(n = b, mu = rep(0, n), Sigma = h2bc * amat.b)

varA <- apply(u_a,1,var)
varE <- varA * (h2ec / h2ac)


e <- mvrnorm(n = b, mu = rep(0, n), Sigma = diag(sqrt(varE), nrow = n))

C <- scale(t(u_a + u_b + u_d + e))

print(varA / (varA + varB + varD + varE))

model <- mixed.solve(C[,1], K = amat.a)
Gh2_C<-model$Vu / (model$Vu + model$Ve)
print(Gh2_C)
vec<-c(vec,Gh2_C)
#}
#boxplot(vec)

fol<-paste0("out/TwostepSimulation_b","n",n)

if (!dir.exists(fol)) {
  dir.create(fol)
}





