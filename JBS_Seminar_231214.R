library(plyr)
library(dplyr)
library(rstan)
library(cmdstanr)
library(bayesplot)
library(ggplot2)
library(RBesT)
library(SAMprior)

set.seed(954)

n_h <- c(107, 44, 51, 39, 139, 20, 78, 35)
x_h <- c( 23, 12, 19,  9,  39,  6,  9, 10)
n_CC <- 6
x_CC <- 1
n_CT <- 23
x_CT <- 14

Historical <- data.frame(
  Study = 1:length(n_h),
  N = n_h, 
  X = x_h
)

Current <- data.frame(
  TRT = c(0,1),
  N = c(n_CC, n_CT),
  X = c(x_CC, x_CT)
)

ANA <- list(H = length(Historical$N),
            x_h = Historical$X,
            n_h = Historical$N,
            x_CC = Current$X[1],
            n_CC = Current$N[1],
            x_CT = Current$X[2],
            n_CT = Current$N[2]
            )

prior <- c(0.5, 0.5)

SAMPLE <- 16000
BURN <- 4000
THIN <- 10

## Current Data analysis
post.cd <- rbeta(SAMPLE, prior[1] + ANA$x_CC, prior[2] + ANA$n_CC - ANA$x_CC)
post.treat <- rbeta(SAMPLE, prior[1] + ANA$x_CT, prior[2] + ANA$n_CT - ANA$x_CT)
post.g.pi.cd <- post.treat - post.cd

## Pooled Data analysis
post.pd <- rbeta(SAMPLE, prior[1] + sum(ANA$x_h) + ANA$x_CC, 
                 prior[2] + sum(ANA$n_h - ANA$x_h) + ANA$n_CC - ANA$x_CC)
post.g.pi.pd <- post.treat - post.pd


## MAP prior
map_mcmc <- gMAP(cbind(X, N - X) ~ 1 | Study,
                 family = binomial,
                 data = Historical,
                 tau.dist = "HalfNormal",
                 tau.prior = 1,
                 beta.prior = cbind(0, 10),
                 warmup = BURN, iter = BURN + SAMPLE, thin = THIN)
map <- mixfit(map_mcmc, Nc = 3, constrain_gt1 = TRUE)
plot(map)$mixdist
fit_map <- postmix(map, r = Current$X[1], n = Current$N[1])
treat_prior <- mixbeta(c(1, 0.5, 0.5))
post.pi.e <- postmix(treat_prior, r = ANA$x_CT, n = ANA$n_CT)
post.g.pi.map <- rmixdiff(post.pi.e, fit_map, SAMPLE)

## R-MAP prior
map_robust <- robustify(map, weight = 0.1, n = 0)
plot(map_robust)
fit_rmap <- postmix(map_robust, r = Current$X[1], n = Current$N[1])
post.g.pi.rmap <- rmixdiff(post.pi.e, fit_rmap, SAMPLE)

## SAM prior
wSAM <- SAM_weight(if.prior = map,
                   delta = 0.15,        ## Clinically significant difference
                   data = c(rep(1, ANA$x_CC), rep(0, ANA$n_CC - ANA$x_CC))   ## Control arm data
                   )
### Assume beta(0.5, 0.5) as the non-informative prior used for mixture
nf.prior  <- mixbeta(nf.prior = c(1, 0.5, 0.5))
### Generate the SAM prior
map_sam <- SAM_prior(if.prior = map, ## Informative prior
                     nf.prior = nf.prior,         ## Non-informative prior
                     weight = wSAM                ## Mixture weight of the SAM prior
                     )
plot(map_sam)
fit_sammap <- postmix(map_sam, r = ANA$x_CC, n = ANA$n_CC)
post.g.pi.SAMmap <- rmixdiff(post.pi.e, fit_sammap, SAMPLE)


DMPP_stan         <- cmdstan_model(paste0(getwd(),'/Stan/DMPP.stan'))
Commensurate_stan <- cmdstan_model(paste0(getwd(),'/Stan/Commensurate.stan'))
HS_stan           <- cmdstan_model(paste0(getwd(),'/Stan/HS.stan'))
UIPDir_stan       <- cmdstan_model(paste0(getwd(),'/Stan/UIPDir.stan'))

## dependent MPP by Banbeta et al (2019)
fit_DMPP   <- DMPP_stan$sample(data=ANA, seed = 999, chains = 4, 
                               parallel_chains = 4, iter_warmup = 4000, 
                               iter_sampling = 20000, thin = 5,
                               adapt_delta = 0.99, refresh = 0)

## Commensurate prior with spike-and-slab
ANA_CMP <- list(x_h = sum(Historical$X),
               n_h = sum(Historical$N),
               x_CC = Current$X[1],
               n_CC = Current$N[1],
               x_CT = Current$X[2],
               n_CT = Current$N[2], 
               w = 0.9, lower_slab = 0.0001, upper_slab = 10, spike = 10
               )
fit_CMP     <- Commensurate_stan$sample(data=ANA_CMP, seed = 999, chains = 4, 
                                       parallel_chains = 4, iter_warmup = 4000, 
                                       iter_sampling = 20000, thin = 5,
                                       max_treedepth = 12, 
                                       adapt_delta = 0.999, refresh = 0)

## Horseshoe prior
ANA_HS <- ANA
ANA_HS[c("betascale", "nu")] <- list(1, 1)
fit_HS     <- HS_stan$sample(data=ANA_HS, seed = 999, chains = 4, 
                             parallel_chains = 4, iter_warmup = 4000, 
                             iter_sampling = 20000, thin = 5,
                             adapt_delta = 0.999, refresh = 0)


## UIP
gamma_h = rep(NA, ANA$H)
for(h in 1:ANA$H){
  gamma_h[h] = min(1, ANA$n_h[h] / ANA$n_CC)
}

theta_h <- rep(NA, ANA$H)
I_U <- rep(NA, ANA$H)
for(h in 1:ANA$H){
  theta_h[h] = ANA$x_h[h] * 1.0 / ANA$n_h[h];
  I_U[h] = 1 / (theta_h[h] * (1 - theta_h[h]));
}  

ANA_UIPDir <- ANA
ANA_UIPDir[c("gamma_h", "theta_h", "I_U")] <- list(gamma_h, theta_h, I_U)
fit_UIPDir  <- UIPDir_stan$sample(data=ANA_UIPDir, seed = 999, chains = 4, 
                                  parallel_chains = 4, iter_warmup = 4000, 
                                  iter_sampling = 8000, adapt_delta = 0.999, 
                                  max_treedepth = 12, refresh = 0)


# Method, Mean, SD, LCI, UCI
Result_ALL <- matrix(0, nrow = 9, ncol = 5)
## Current data analysis
Result_ALL[1,]  <- c(1, mean(post.g.pi.cd), sd(post.g.pi.cd), quantile(post.g.pi.cd, 0.025, names = F), quantile(post.g.pi.cd, 0.975, names = F))
## Pooled data analysis
Result_ALL[2,]  <- c(2, mean(post.g.pi.pd), sd(post.g.pi.pd), quantile(post.g.pi.pd, 0.025, names = F), quantile(post.g.pi.pd, 0.975, names = F))
## MAP
Result_ALL[3,]  <- c(3, mean(post.g.pi.map), sd(post.g.pi.map), quantile(post.g.pi.map, 0.025, names = F), quantile(post.g.pi.map, 0.975, names = F))
## R-MAP
Result_ALL[4,]  <- c(4, mean(post.g.pi.rmap), sd(post.g.pi.rmap), quantile(post.g.pi.rmap, 0.025, names = F), quantile(post.g.pi.rmap, 0.975, names = F))
## SAM-MAP
Result_ALL[5,]  <- c(5, mean(post.g.pi.SAMmap), sd(post.g.pi.SAMmap), quantile(post.g.pi.SAMmap, 0.025, names = F), quantile(post.g.pi.SAMmap, 0.975, names = F))
## Power prior (DMPP)
Result_ALL[6,]  <- c(6, mean(fit_DMPP$draws("g_pi")), sd(fit_DMPP$draws("g_pi")), quantile(fit_DMPP$draws("g_pi"),0.025,names=F),quantile(fit_DMPP$draws("g_pi"),0.975,names=F))
## Commensurate prior
Result_ALL[7,]  <- c(7, mean(fit_CMP$draws("g_pi")), sd(fit_CMP$draws("g_pi")), quantile(fit_CMP$draws("g_pi"),0.025,names=F),quantile(fit_CMP$draws("g_pi"),0.975,names=F))
## Potential bias model with horseshoe prior
Result_ALL[8,]  <- c(8, mean(fit_HS$draws("g_pi")), sd(fit_HS$draws("g_pi")), quantile(fit_HS$draws("g_pi"),0.025,names=F),quantile(fit_HS$draws("g_pi"),0.975,names=F))
## UIP
Result_ALL[9,]  <- c(9, mean(fit_UIPDir$draws("g_pi")), sd(fit_UIPDir$draws("g_pi")), quantile(fit_UIPDir$draws("g_pi"),0.025,names=F),quantile(fit_UIPDir$draws("g_pi"),0.975,names=F))

Result_ALL1 <- as.data.frame(Result_ALL,stringsAsFactors = F) 
colnames(Result_ALL1) <- c("MethodN", "Mean", "SD", "LCI", "UCI")
Result_ALL1 <- Result_ALL1 %>% dplyr::mutate(Range = (UCI-LCI)) %>%
  mutate(Method = case_when(
    MethodN == 1 ~ "Current data",
    MethodN == 2 ~ "Pooled data",
    MethodN == 3 ~ "MAP",
    MethodN == 4 ~ "R-MAP",
    MethodN == 5 ~ "SAM-MAP",
    MethodN == 6 ~ "DMPP",
    MethodN == 7 ~ "Commensurate",
    MethodN == 8 ~ "Horseshoe",
    MethodN == 9 ~ "UIP",
    TRUE ~ ""
  ))

options(digits = 3)
Result_ALL10 <- Result_ALL1 %>% 
  mutate(Mean100 = Mean * 100, SD100 = SD * 100,
         LCI100 = LCI * 100, UCI100 = UCI * 100, Range100 = Range * 100)
Result_ALL10[,c(7, 8, 9, 10, 11, 12)]

Result_ALL10$Method <- factor(Result_ALL10$Method, 
                             levels = rev(c("Current data", "Pooled data", 
                                            "MAP", "R-MAP", "SAM-MAP", "DMPP", 
                                            "Commensurate", "Horseshoe", "UIP")))

Res_forest <- ggplot(Result_ALL10, aes(x = Mean100, y = Method)) +
  xlab("Treatment effect in terms of the response rate (%)") + 
  geom_point(size = 4) +
  geom_linerange(aes(xmin = LCI100, xmax = UCI100), size = 1.2) +
  coord_cartesian(xlim = c(0, 100), expand = T) +
  theme(
    plot.title = element_text(colour = "black", size = 24),
    axis.text.x = element_text(colour = "black", size = 24),
    axis.text.y = element_text(colour = "black", size = 24),
    axis.ticks.x = element_line(colour = "black", size=1),
    axis.ticks.y = element_line(colour = "black", size=1),
    axis.title.x = element_text(colour = "black", size = 24),
    axis.title.y = element_blank(),
    axis.line.x = element_line(colour = "#000000", size = 1, linetype = "solid", lineend = "round"),
    axis.line.y = element_line(colour = "#000000", size = 1, linetype = "solid", lineend = "round"),
    panel.background = element_rect(fill = "transparent", colour = NA),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 22),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.key.width = unit(0.8,"cm"),
    legend.position = "bottom"
  )
Res_forest
ggsave("Res_forestplot.png", Res_forest, width = 12, height = 7,bg = "white")

