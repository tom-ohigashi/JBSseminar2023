
library(ggplot2)
library(dplyr)
library(fdrtool)
set.seed(954)


data_beta <- data.frame(
  theta = seq(0, 1, by = 0.01)
) %>% mutate(y = dbeta(theta, 4, 16))

fig_beta <- ggplot(data_beta, aes(x = theta, y = y)) +
  geom_line(linewidth = 2) +
  xlab(expression(theta)) +
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
    legend.position = c(0.8,0.5),
    legend.key.size = unit(2, 'lines'),
    legend.spacing.x = unit(0.5, 'cm')
  )
fig_beta
ggsave("beta1.png", fig_beta, width = 9, height = 7,bg = "white")


data_beta <- data.frame(
  theta = seq(0, 1, by = 0.01)
) %>% mutate(y = dbeta(theta, 0.5, 0.5))

fig_beta <- ggplot(data_beta, aes(x = theta, y = y)) +
  geom_line(linewidth = 2) +
  xlab(expression(theta)) +
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
    legend.position = c(0.8,0.5),
    legend.key.size = unit(2, 'lines'),
    legend.spacing.x = unit(0.5, 'cm')
  )
fig_beta
ggsave("beta2.png", fig_beta, width = 9, height = 7,bg = "white")



data_beta <- data.frame(
  theta = seq(0, 1, by = 0.01)
) %>% mutate(y = dbeta(theta, 4.5, 16.5))

fig_beta <- ggplot(data_beta, aes(x = theta, y = y)) +
  geom_line(linewidth = 2) +
  xlab(expression(theta)) +
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
    legend.position = c(0.8,0.5),
    legend.key.size = unit(2, 'lines'),
    legend.spacing.x = unit(0.5, 'cm')
  )
fig_beta
ggsave("beta3.png", fig_beta, width = 9, height = 7,bg = "white")


data_pool1 <- rbeta(50000, 10.5, 20.5)
data_pool2 <- rbeta(50000, 5.5, 5.5)

data_pool_CT <- rbeta(50000, 5.5, 5.5)

data_diff1 <- data_pool_CT - data_pool1
data_diff2 <- data_pool_CT - data_pool2

data_diff <- data.frame(
  diff1 = data_diff1,
  diff2 = data_diff2
)

fig_diff1 <- ggplot(data_diff, aes(x=diff1)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#FF0099", linewidth = 2) + 
  geom_histogram() +
  xlab(expression(theta[CT] - theta[CC])) +
  ylab("Density") +
  scale_x_continuous(breaks = seq(-0.5,0.5,by=0.5),limits = c(-0.8,0.8)) +
  scale_y_continuous(breaks = c(),limits = c(0, 8000)) +
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
    legend.position = c(0.8,0.5),
    legend.key.size = unit(2, 'lines'),
    legend.spacing.x = unit(0.5, 'cm')
  )
fig_diff1
ggsave("diff1.png", fig_diff1, width = 8, height = 6,bg = "white")

fig_diff2 <- ggplot(data_diff, aes(x=diff2)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#FF0099", linewidth = 2) + 
  geom_histogram() +
  xlab(expression(theta[CT] - theta[CC])) +
  ylab("Density") +
  scale_x_continuous(breaks = seq(-0.5,0.5,by=0.5),limits = c(-0.8,0.8)) +
  scale_y_continuous(breaks = c(),limits = c(0, 8000)) +
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
    legend.position = c(0.8,0.5),
    legend.key.size = unit(2, 'lines'),
    legend.spacing.x = unit(0.5, 'cm')
  )
fig_diff2
ggsave("diff2.png", fig_diff2, width = 8, height = 6,bg = "white")

mean(data_diff$diff1>0)
mean(data_diff$diff2>0)


data_hn <- data.frame(
  theta = rep(seq(0.001, 5, by = 0.01), 3),
  sigma = c(rep(0.5, length(seq(0.001, 5, by = 0.01))), 
            rep(1.0, length(seq(0.001, 5, by = 0.01))), 
            rep(2.0, length(seq(0.001, 5, by = 0.01)))),
  sigmap = c(rep("0.5", length(seq(0.001, 5, by = 0.01))), 
            rep("1.0", length(seq(0.001, 5, by = 0.01))), 
            rep("2.0", length(seq(0.001, 5, by = 0.01))))
) %>% mutate(y = dhalfnorm(theta, sqrt(pi/2)/sigma))

fig_hn <- ggplot(data_hn, aes(x = theta, y = y, linetype = sigmap)) +
  geom_line(linewidth = 2) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash"))+
  xlab(expression(tau)) +
  ylab("Density") +
  labs(linetype = "Scale") +
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
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 24),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.key.width = unit(2,"cm"),
    legend.position = c(0.8,0.5),
    legend.key.size = unit(2, 'lines'),
    legend.spacing.x = unit(0.5, 'cm')
  )
fig_hn
ggsave("hn.png", fig_hn, width = 6, height = 4,bg = "white")



data_ss <- data.frame(
  theta = seq(0.001, 5, by = 0.01)
) %>% mutate(y = (0.1 * dunif(theta, 0, 5) + 0.9 * dnorm(theta, 5, 0.01)))

fig_ss <- ggplot(data_ss, aes(x = theta, y = y)) +
  geom_line(linewidth = 2) +
  xlab(expression(1/tau[CP]^2)) +
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
    legend.position = c(0.8,0.5),
    legend.key.size = unit(2, 'lines'),
    legend.spacing.x = unit(0.5, 'cm')
  )
fig_ss
ggsave("ss.png", fig_ss, width = 6, height = 4,bg = "white")
