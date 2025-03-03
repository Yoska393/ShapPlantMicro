library(rrBLUP)
library(MASS)

n <- 100 #Samples
a <- 200
b <- 200
d <- 200

A <- matrix(rnorm(n = n*a, mean = 0), nrow = n, ncol = a)
D <- matrix(rnorm(n = n*d, mean = 0), nrow = n, ncol = d)

amat.a  <-  tcrossprod(A)
amat.d  <-  tcrossprod(D)
#image(amat.a)

# herit
h2a <- 0.60
h2d <- 0.30
h2e <- 0.1

u_a <- mvrnorm(n = 1, mu = rep(0, n), Sigma = amat.a)
varA <- var(u_a)
varD <- (varA / h2a) - varA - varA*(h2e/h2a)
varE <- (varA / h2a) - varA - varD

varA / (varA + varD + varE)
varD / (varA + varD + varE)
varE / (varA + varD + varE)


t1 <- sum(diag(amat.d)) / n
t2 <- (sum(diag(amat.d)) + 2 * sum(amat.d[upper.tri(amat.d)])) / n^2

sigmaD <- varD / (t1 - t2)

u_d <- mvrnorm(n = 1, mu = rep(0, n), Sigma = sigmaD * amat.d)
var(u_d)
varD

e <- mvrnorm(n = 1, mu = rep(0, n), Sigma = diag(varE, nrow = n))
var(e)
var(u_a) / (var(u_a) + var(u_d) + var(e))
var(u_d) / (var(u_a) + var(u_d) + var(e))
var(e) / (var(u_a) + var(u_d) + var(e))

b <- u_a + u_d + e

#check 
model <- mixed.solve(b, K = amat.a)
Gh2 <- model$Vu / (model$Vu + model$Ve)
print(Gh2)


```{r ingherit model,include=F,eval=F}
#heritability model
n <- 5 #Samples
a <- 2
a2<- 2
b <- 2
c <- 2
d <- 2

#A  <-  matrix(rnorm(n = n*a, mean = 0, sd = 1),n,a)
#A2 <- matrix(rnorm(n = n*a2, mean = 0, sd = 1 ),n,a)
#D  <-  matrix(rnorm(n = n*d, mean = 0, sd = 1),n,d)

A <- matrix(rnorm(n = n*a, mean = 0, sd = 1 / (a**(1/2))),n,a)
#A2 <- matrix(rnorm(n = n*a, mean = 0, sd = 1 / (a2**(1/2))),n,a)
D <- matrix(rnorm(n = n*d, mean = 0, sd = 1 / (d**(1/2))),n,d)

amat.a  <-  A.mat(X = A)
#amat.a2 <-  A.mat(X = A2)
amat.d  <-  A.mat(X = D)
#image(amat.a)

# herit
h2a <- 0.61
#h2a2 <- 0.3
h2d <- 0.29
h2e <- 0.1

u_a <- mvrnorm(n = b, mu = rep(0, n), Sigma = h2a * amat.a)
u_d <- mvrnorm(n = b, mu = rep(0, n), Sigma = h2d * amat.d)

varA <- apply(u_a,1,var)
varD <- apply(u_d,1,var)
varE <- varA * (h2e / h2a)
#
print(varA / (varA + varD + varE))

e <- mvrnorm(n = b, mu = rep(0, n), Sigma = diag(sqrt(varE), nrow = n))
B <- t(u_a + u_d + e)
#B <- scale(B)
B <- scale(B) / (b**(1/2))
#print(varA / apply(B,2,var))

amat.b <-  A.mat(X = B)

#check 
model <- mixed.solve(B[,1], K = amat.a)
Gh2<-model$Vu / (model$Vu + model$Ve)
print(Gh2)


# generate C

h2ac <- 0.2
#h2a2c <- 0.2
h2bc <- 0.44
h2dc <- 0.26
h2ec <- 0.1

#vec<-c()
#for (i in 1:100){

u_a <- mvrnorm(n = b, mu = rep(0, n), Sigma = h2ac * amat.a)
u_d <- mvrnorm(n = b, mu = rep(0, n), Sigma = h2dc * amat.d)
u_b <- mvrnorm(n = b, mu = rep(0, n), Sigma = h2bc * amat.b)

varA <- apply(u_a,1,var)
varB <- apply(u_b,1,var)
varD <- apply(u_d,1,var)
varE <- varA * (h2ec / h2ac)


e <- mvrnorm(n = b, mu = rep(0, n), Sigma = diag(sqrt(varE), nrow = n))

C <- scale(t(u_a + u_b + u_d + e))

print(varA / (varA + varB + varD + varE))

model <- mixed.solve(C[,1], K = amat.a)
Gh2_C<-model$Vu / (model$Vu + model$Ve)
print(Gh2_C)
#}
#boxplot(vec)

fol<-paste0("out/TSSblup","nabcd",n,a,b,c,d)

if (!dir.exists(fol)) {
  dir.create(fol)
}
if (!dir.exists(paste0(fol,"/pdf"))) {
  dir.create(paste0(fol,"/pdf"))
}
```

