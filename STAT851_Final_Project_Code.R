############################################################################################
############################################################################################
################################## Jason's Code ############################################
############################################################################################
############################################################################################

# Section 2.1 Case I: Fixed Grand Total with Melanoma example
melanoma=read.table("melanoma.txt",header=T)

#Make Type and Site factors:
melanoma$Type=factor(melanoma$Type)
melanoma$Site=factor(melanoma$Site)

xtabs(Count ~ Type + Site, melanoma)

#Fit the log-linear model to the melanoma data:
fit=glm(Count~Type*Site,poisson,melanoma)

#Find the p-value associated with the interaction effect:
anova(fit, test="Chisq")


# Section 2.2 Case II: Fixed Row Totals with Diet example
library(ggplot2)

diet <- read.csv(file = "Fiber.csv")
diet$fiber<-factor(x = diet$fiber, levels = c("none", "bran", "gum", "both"))

xtabs(formula = count ~ fiber + bloat, data = diet)

fit_diet=glm(count~fiber*bloat,poisson,diet)
anova(fit_diet,test="Chisq")
summary(fit_diet)



# Section 2.3  Use multinom() function in the nnet package

# We build the multinomial logit model in three different ways, and to demonstrate that
# the multinom() function is unstable in this case as it produces different coefficient
# estimations if we specify the model differently.

fit_multinom_typesite <- multinom(formula = Type ~ Site, weights=Count, data = melanoma)
coef(fit_multinom_typesite)
Anova(fit_multinom_typesite)

fit_multinom_sitetype <- multinom(formula = Site ~ Type, weights=Count, data = melanoma)
coef(fit_multinom_sitetype)
Anova(fit_multinom_sitetype)

fit.glm <- glm(formula = Count ~ Type * Site, poisson, data = melanoma)
coef(fit.glm)
anova(fit.glm)

# Next, we apply the multinom function to the Wheat Kernel data set:

wheat <- read.csv(file  = "wheat.csv")
head(wheat)

levels(wheat$type)
mod.fit<-multinom(formula = type ~ class + density + hardness + size + weight + moisture, data=wheat)
coef(mod.fit)
Anova(mod.fit)
predict(mod.fit, type='probs')[1:5,]


############################################################################################
############################################################################################
################################## Michael's Code ##########################################
############################################################################################
############################################################################################

##Proportional model using VGAM package
library(VGAM)
library(nnet)
library(MASS)
library(car)


# Figure 1: Graph for proportional odds illustration
beta<-c(0, 2, 4, 2,3,4) #beta10, beta20, beta30, beta1
x.range<-c(-5, 3)

x11(width = 10, height = 6, pointsize = 12)

## CDF
par(mfrow = c(1, 2))
curve(expr = plogis(q = beta[1] + beta[4]*x), xlim = x.range, ylab = expression(P(Y<=j)),
      xlab = expression(x[1]), main = "Cumulative probabilities for Y", lwd = 2)
curve(expr = plogis(q = beta[2] + beta[4]*x), add = TRUE, lty = "dashed", col = "red", lwd = 2)
curve(expr = plogis(q = beta[3] + beta[4]*x), add = TRUE, lty = "dotted", , col = "blue", lwd = 2)
legend(x = -5.5, y = 0.9, legend = c(expression(P(Y<=1)), expression(P(Y<=2)), expression(P(Y<=3))),
       lty = c("solid", "dashed", "dotted", "dotdash"), col = c("black", "red", "blue"),
       bty = "n", lwd = 2)
## pdf
curve(expr = plogis(q = beta[1] + beta[4]*x), xlim = x.range, ylab = expression(pi[j]),
      xlab = expression(x[1]), main = "Probabilities for Y", lwd = 2)
curve(expr = plogis(q = beta[2] + beta[4]*x) - plogis(q = beta[1] + beta[4]*x), add = TRUE,
      lty = "dashed", col = "red", lwd = 2)
curve(expr = plogis(q = beta[3] + beta[4]*x) - plogis(q = beta[2] + beta[4]*x), add = TRUE,
      lty = "dotted", , col = "blue", lwd = 2)
curve(expr = 1 - plogis(q = beta[3] + beta[4]*x), add = TRUE,
      lty = "dotdash", col = "green", lwd = 2)
legend(x = -5.5, y = 0.9, legend = c(expression(pi[1]), expression(pi[2]), expression(pi[3]), expression(pi[4])),
       lty = c("solid", "dashed", "dotted", "dotdash"), col = c("black", "red", "blue", "green"),
       bty = "n", lwd = 2)



## Figure 3: Non-proportional odds graph illustration

## CDF
par(mfrow = c(1, 2))
curve(expr = plogis(q = beta[1] + beta[4]*x), xlim = x.range, ylab = expression(P(Y<=j)),
      xlab = expression(x[1]), main = "Cumulative probabilities for Y", lwd = 2)
curve(expr = plogis(q = beta[2] + beta[5]*x), add = TRUE, lty = "dashed", col = "red", lwd = 2)
curve(expr = plogis(q = beta[3] + beta[6]*x), add = TRUE, lty = "dotted", , col = "blue", lwd = 2)
legend(x = -5.5, y = 0.9, legend = c(expression(P(Y<=1)), expression(P(Y<=2)), expression(P(Y<=3))),
       lty = c("solid", "dashed", "dotted", "dotdash"), col = c("black", "red", "blue"),
       bty = "n", lwd = 2, cex=0.8)

## pdf
curve(expr = plogis(q = beta[1] + beta[4]*x), xlim = x.range, ylab = expression(pi[j]),
      xlab = expression(x[1]), main = "Probabilities for Y", lwd = 2)
curve(expr = plogis(q = beta[2] + beta[5]*x) - plogis(q = beta[1] + beta[4]*x), add = TRUE,
      lty = "dashed", col = "red", lwd = 2)
curve(expr = plogis(q = beta[3] + beta[6]*x) - plogis(q = beta[2] + beta[5]*x), add = TRUE,
      lty = "dotted", , col = "blue", lwd = 2)
curve(expr = 1 - plogis(q = beta[3] + beta[4]*x), add = TRUE,
      lty = "dotdash", col = "green", lwd = 2)
legend(x = -7.2, y = 0.9, legend = c(expression(pi[1]), expression(pi[2]), expression(pi[3]), expression(pi[4])),
       lty = c("solid", "dashed", "dotted", "dotdash"), col = c("black", "red", "blue", "green"),
       bty = "n", lwd = 2)


#Cumulative logit model, proportional odds, Wheat Kernel data

# Need to reorder factor so that Scab < Sprout < Healthy
wheat$type.order<-factor(wheat$type, levels = c("Scab",  "Sprout", "Healthy"))

mod.fit.ord<-polr(formula = type.order ~ class + density + hardness + size + weight + moisture, data = wheat, method = "logistic")
summary(mod.fit.ord)

Anova(mod.fit.ord)
# Estimate probability of being in a particular category

pi.hat.ord<-predict(object = mod.fit.ord, type = "probs")
head(pi.hat.ord)
wheat[1,]




# Confidence intervals for pi_j

# Obtain observation values, ";" means end of expression
x1<-0;          x2<-wheat[1,2]; x3<-wheat[1,3]
x4<-wheat[1,4]; x5<-wheat[1,5]; x6<-wheat[1,6]


## Rewriete the deltaMethod to correct negative sign issue

x1<-0;          x2<-wheat[1,2]; x3<-wheat[1,3]
x4<-wheat[1,4]; x5<-wheat[1,5]; x6<-wheat[1,6]


deltaMethod.polr2<-function(object, g)  {
  # All beta^'s where the slope estimates are adjusted
  beta.hat<-c(-object$coefficients, object$zeta)
  
  # Count the number of slopes and intercepts
  numb.slope<-length(object$coefficients)
  numb.int<-length(object$zeta)
  
  # Name corresponding parameters
  names(beta.hat)<-c(paste("b", 1:numb.slope, sep=""), paste("b", 1:numb.int, "0", sep=""))
  
  # Fix covariance matrix - All covariances between slopes and intercepts
  #  need to be multiplied by -1
  cov.mat<-vcov(object)
  # Above diagonal
  cov.mat[1:numb.slope, (numb.slope + 1):(numb.slope + numb.int)]<-
    -cov.mat[1:numb.slope, (numb.slope + 1):(numb.slope + numb.int)]
  # Below diagonal
  cov.mat[(numb.slope + 1):(numb.slope + numb.int), 1:numb.slope]<-
    -cov.mat[(numb.slope + 1):(numb.slope + numb.int), 1:numb.slope]
  
  # deltaMethod.default() method function completes calculations
  deltaMethod(object = beta.hat, g = g, vcov. = cov.mat)
}

alpha<-0.05

## Parts of character string

## Specify the model for each response level
scab<-"exp(b10 + b1*x1 + b2*x2 + b3*x3 + b4*x4 + b5*x5 + b6*x6)"
sprout<-"exp(b20 + b1*x1 + b2*x2 + b3*x3 + b4*x4 + b5*x5 + b6*x6)"

# pi^_Scab
g.scab<-paste(scab, "/ (1 + ", scab, ")")
g.scab
calc.scab<-deltaMethod.polr2(object = mod.fit.ord, g = g.scab)
calc.scab$Estimate + qnorm(p = c(alpha/2, 1-alpha/2))*calc.scab$SE

# pi^_Sprout
g.sprout<-paste(sprout, "/ (1 + ", sprout, ")", " - ", scab, "/ (1 + ", scab, ")")
g.sprout
calc.sprout<-deltaMethod.polr2(object = mod.fit.ord, g = g.sprout)
calc.sprout$Estimate + qnorm(p = c(alpha/2, 1-alpha/2))*calc.sprout$SE

# pi^_Healthy
g.healthy<-paste("1 - ", sprout, "/ (1 + ", sprout, ")")
g.healthy
# Alternatively, directly enter the string:
g.healthy<-"1 - exp(b20 + b1*x1 + b2*x2 + b3*x3 + b4*x4 + b5*x5 + b6*x6) /
(1 + exp(b20 + b1*x1 + b2*x2 + b3*x3 + b4*x4 + b5*x5 + b6*x6))"
calc.healthy<-deltaMethod.polr2(object = mod.fit.ord, g = g.healthy)
calc.healthy$Estimate  # pi^_Healthy
calc.healthy$SE  # sqrt(Var^(pi^_Healthy))
calc.healthy$Estimate + qnorm(p = c(alpha/2, 1-alpha/2))*calc.healthy$SE



## Comparison graph for cumulative logit model and nominal model
## Only has 'density' as predictor, treat all other variables as constant

# Estimate model with density only
mod.fit.dens<-polr(formula = type.order ~ density, data = wheat, method = "logistic")
summary(mod.fit.dens)


min(wheat$density)
max(wheat$density)


mod.fit.nom.density<-multinom(formula = type ~ density, data = wheat)
summary(mod.fit.nom.density)
beta.hat<-coefficients(mod.fit.nom.density)

# Figure 2: Estimated proportional odds (thicker line) and multinomial (thinner line) 
# regression models for the wheat data where density is the only explanatory variable 
# included in the model.
curve(expr = 1/(1 + exp(beta.hat[1,1] + beta.hat[1,2]*x) + exp(beta.hat[2,1] + beta.hat[2,2]*x)), ylab = expression(hat(pi)), xlab = "Density",
      xlim = c(min(wheat$density), max(wheat$density)), col = "black", lty = "solid", lwd = 2, n = 1000, type = "n",
      panel.first = grid(col = "gray", lty = "dotted"))
# Plot each pi_j
curve(expr = 1/(1 + exp(beta.hat[1,1] + beta.hat[1,2]*x) + exp(beta.hat[2,1] + beta.hat[2,2]*x)),
      col = "black", lty = "solid", lwd = 2, n = 1000, add = TRUE,
      xlim = c(min(wheat$density[wheat$type == "Healthy"]), max(wheat$density[wheat$type == "Healthy"])))  # Healthy
curve(expr = exp(beta.hat[1,1] + beta.hat[1,2]*x)/(1 + exp(beta.hat[1,1] + beta.hat[1,2]*x) + exp(beta.hat[2,1] + beta.hat[2,2]*x)),
      col = "green", lty = "dotdash", lwd = 2, n = 1000, add = TRUE,
      xlim = c(min(wheat$density[wheat$type == "Scab"]), max(wheat$density[wheat$type == "Scab"])))  # Scab
curve(expr = exp(beta.hat[2,1] + beta.hat[2,2]*x)/(1 + exp(beta.hat[1,1] + beta.hat[1,2]*x) + exp(beta.hat[2,1] + beta.hat[2,2]*x)),
      col = "red", lty = "longdash", lwd = 2, n = 1000, add = TRUE,
      xlim = c(min(wheat$density[wheat$type == "Sprout"]), max(wheat$density[wheat$type == "Sprout"])))  # Sprout
legend(x = 1.35, y = 0.8, legend=c("Healthy", "Sprout", "Scab"), lty=c("solid","longdash","dotdash"),
       col=c("black","red","green"), bty="n", lwd = c(2,2,2), seg.len = 4)


lwd.po<-4
# Plot each pi_j for proportional odds model
curve(expr = plogis(q = mod.fit.dens$zeta[1] - mod.fit.dens$coefficients*x), col = "green",
      type = "l", xlim = c(min(wheat$density[wheat$type.order == "Scab"]), max(wheat$density[wheat$type.order == "Scab"])),
      add = TRUE, lty = "dotdash", lwd = lwd.po, n = 1000)  # Scab
curve(expr = plogis(q = mod.fit.dens$zeta[2] - mod.fit.dens$coefficients*x) - plogis(q =mod.fit.dens$zeta[1] - mod.fit.dens$coefficients*x), col = "red",
      type = "l", xlim = c(min(wheat$density[wheat$type.order == "Sprout"]), max(wheat$density[wheat$type.order == "Sprout"])),
      add = TRUE, lty = "longdash", lwd = lwd.po, n = 1000)  # Sprout
curve(expr = 1 - plogis(q = mod.fit.dens$zeta[2] - mod.fit.dens$coefficients*x), col = "black",
      type = "l", xlim = c(min(wheat$density[wheat$type.order == "Healthy"]), max(wheat$density[wheat$type.order == "Healthy"])),
      add = TRUE, lty = "solid", lwd = lwd.po, n = 1000)  # Healthy
legend(x = 1.35, y = 0.8, legend=c("Healthy", "Sprout", "Scab"), lty=c("solid","longdash","dotdash"),
       col=c("black","red","green"), bty="n", lwd = c(2,2,2), seg.len = 4)



## Fiber Example

# Match order given at DASL
diet$fiber<-factor(x = diet$fiber, levels = c("none", "bran", "gum", "both"))
diet$bloat<-factor(x = diet$bloat, levels = c("none", "low", "medium", "high"))

diet.table<-xtabs(formula = count ~ fiber + bloat, data = diet)
diet.table


#Fit proportional odds model using VGAM package
mod.fit.po<-vglm(formula = bloat ~ fiber, family = cumulative(parallel = TRUE),
                 weights = count, data = diet[diet$count != 0,])
summary(mod.fit.po)
slotNames(mod.fit.po)  # Like names( ) in S3
mod.fit.po@coefficients
options(width = 80)
mod.fit.po@df.residual
#showMethods(class = "vglm") #Like method(class = " ") in S3

#Fit proportional odds model using polr
mod.fit.ord<-polr(formula = bloat ~ fiber, weights = count, data=diet, method = "logistic")
summary(mod.fit.ord)
Anova(mod.fit.ord)

#non-proportional odds model
mod.fit.npo<-vglm(formula = bloat ~ fiber, family = cumulative(parallel = FALSE),
                  weights = count, data = diet[diet$count != 0,])
summary(mod.fit.npo)



#anova(mod.fit.po, mod.fit.npo)
#LRT to compare two models to test proportional odds assumption
tran.LR<-deviance(mod.fit.po) - deviance(mod.fit.npo)
df<-mod.fit.po@df.residual - mod.fit.npo@df.residual
p.value<-1 - pchisq(q = tran.LR, df = df)
data.frame(tran.LR, df, p.value)


#adjacet-categories models, non-proportional
mod.fit.full <- vglm(formula = bloat ~ fiber , weights=count , family = acat ( parallel = FALSE ), 
                     data =diet[diet$count != 0,])
summary(mod.fit.full)



############################################################################################
############################################################################################
################################## John's Code #############################################
############################################################################################
############################################################################################

tonsils <-data.frame(carrier = c(1, 0),
                     y1  = c(19, 497),
                     y2  = c(29, 560), 
                     y3  = c(24, 269))
tonsils
# CR logit model 
library(VGAM)
fit <- vglm(cbind(y1, y2, y3) ~ carrier, 
            family = cratio(reverse = F, parallel = T),
            data = tonsils) # proportional odds 
summary(fit)
fitted(fit) # probabilities
constraints(fit)
coef(fit)

fit1 <- vglm(cbind(y1, y2, y3) ~ carrier, 
             family = cratio(reverse = F, parallel = F),
             data = tonsils) # non-proportional odds 
summary(fit1)