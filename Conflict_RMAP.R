library(RBesT)
library(dplyr)
library(ggplot2)
set.seed(954)

n_h <- rep(60, 4)
x_h <- c(10, 12, 14, 16)

Historical <- data.frame(
  Study = 1:length(n_h),
  N = n_h, 
  X = x_h
)

SAMPLE <- 16000
BURN <- 4000
THIN <- 10

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

map_robust <- robustify(map, weight = 0.1, n = 0)
plot(map_robust)


prior_plot <- ggplot(data.frame(theta=c(0, 1)), aes(theta)) +
  stat_function(fun=dmix, args=list(mix=map), aes(linetype="MAP"), linewidth = 2) +
  stat_function(fun=dmix, args=list(mix=map_robust), aes(linetype="R-MAP"), linewidth = 2) +
  xlab(expression(theta[CC])) +
  ylab("Density") +
  theme(
    plot.title = element_text(size = 24),
    axis.text = element_text(colour = "black", size = 24),
    axis.ticks=element_line(colour = "black", size=1),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.line = element_line(colour = "#000000", size = 1, linetype = "solid", lineend = "round"),
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
    legend.text = element_text(size = 24),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.key.width = unit(2,"cm"),
    legend.position = c(0.6,0.8),
    legend.key.size = unit(2, 'lines'),
    legend.spacing.x = unit(0.5, 'cm')
  )
ggsave("Prior_plot.png", prior_plot, width = 8, height = 7,bg = "white")

# No conflict
n_CC <- 20
x_CC <- 4

fit_map1 <- postmix(map, r = x_CC, n = n_CC)
fit_rmap1 <- postmix(map_robust, r = x_CC, n = n_CC)

noconflict_plot <- ggplot(data.frame(theta=c(0, 1)), aes(theta)) +
  stat_function(fun=dmix, args=list(mix=fit_map1), aes(linetype="MAP"), linewidth = 2) +
  stat_function(fun=dmix, args=list(mix=fit_rmap1), aes(linetype="R-MAP"), linewidth = 2) +
  xlab(expression(theta[CC])) +
  ylab("Density") +
  theme(
    plot.title = element_text(size = 24),
    axis.text = element_text(colour = "black", size = 24),
    axis.ticks=element_line(colour = "black", size=1),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.line = element_line(colour = "#000000", size = 1, linetype = "solid", lineend = "round"),
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
    legend.text = element_text(size = 24),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.key.width = unit(2,"cm"),
    legend.position = c(0.6,0.8),
    legend.key.size = unit(2, 'lines'),
    legend.spacing.x = unit(0.5, 'cm')
  )
ggsave("No_conflict_plot.png", noconflict_plot, width = 8, height = 7,bg = "white")

# Conflict case
n_CC <- 20
x_CC <- 10

fit_map2 <- postmix(map, r = x_CC, n = n_CC)
fit_rmap2 <- postmix(map_robust, r = x_CC, n = n_CC)


conflict_plot <- ggplot(data.frame(theta=c(0, 1)), aes(theta)) +
  stat_function(fun=dmix, args=list(mix=fit_map2), aes(linetype="MAP"), linewidth = 2) +
  stat_function(fun=dmix, args=list(mix=fit_rmap2), aes(linetype="R-MAP"), linewidth = 2) +
  xlab(expression(theta[CC])) +
  ylab("Density") +
  theme(
    plot.title = element_text(size = 24),
    axis.text = element_text(colour = "black", size = 24),
    axis.ticks=element_line(colour = "black", size=1),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.line = element_line(colour = "#000000", size = 1, linetype = "solid", lineend = "round"),
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
    legend.text = element_text(size = 24),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.key.width = unit(2,"cm"),
    legend.position = c(0.6,0.8),
    legend.key.size = unit(2, 'lines'),
    legend.spacing.x = unit(0.5, 'cm')
  )
ggsave("Conflict_plot.png", conflict_plot, width = 8, height = 7,bg = "white")
