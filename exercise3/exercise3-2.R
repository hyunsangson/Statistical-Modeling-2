titanic <- read.csv('titanic.csv')

#rename the first column
summary(titanic)
names(titanic)[1] <- "Name"

titanic <- titanic[!is.na(titanic$Age), ] ##titanic dataset has 557 missing in age. Remove the missing age rows
summary(titanic) ##no missing in age

X <- matrix(titanic$Age, nrow(titanic), 1)
y <- as.numeric(titanic$Survive == "Yes")

##MAP
log_posterior <- function(beta, X, y){
  y %*% log(1/(1+exp(-X %*% beta))) + (1 - y) %*% log(1 - (1/(1+exp(-X %*% beta)))) - 0.5 * t(beta) %*% beta
}


MAP <- optim(0, function(beta) -log_posterior(beta, X, y), method = "Brent", lower = -1, upper = 1)
MAP


## Exercise 3-3. Plotting
betas <- seq(MAP$par -0.25, map$par + 0.25, 0.001)
values <- rep(NA, length(betas))
for (i in 1:length(betas)) values[i] <- exp(log_posterior(betas[i], X, y))

plot(betas, values, type = "l", main = "Posterior for betas")
