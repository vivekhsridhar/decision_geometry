rm(list = ls())

## Load packages
library(ggplot2)
library(viridis)

## Set working directory
setwd("/Users/viveksridhar/Documents/Code/multi-choice_decision_geometry/multi-choice_decision_geometry/output/")
dir <- getwd()

theme_cus <- function(base_size = 12, base_family = "Helvetica"){
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(axis.title.x = element_text(size = 15, margin = margin(15,0,0,0)), 
          axis.title.y = element_text(size = 15, margin = margin(0,15,0,0)), 
          axis.text = element_text(size = 12),
          axis.ticks.length = unit(0.3, "lines"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          legend.position="right",
          legend.background = element_rect(colour = "black"),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          axis.line.x = element_line(color="black", size = 0.5), 
          axis.line.y = element_line(color="black", size = 0.5)) 
}

## Load data
df <- read.csv("decision_geometry_spin1d.csv")

## Input parameters
symmetric = FALSE
n_cues = 2

targets_x <- rep(0, n_cues)
targets_y <- rep(0, n_cues)
if (symmetric) {
  max_angle = 2*pi
  theta = max_angle / n_cues
  
  targets_x <- 500 + 1000 * cos((i-2) * theta)
  targets_x <- 500 + 1000 * sin((i-2) * theta)
} else {
  for (i in 1:n_cues) {
    max_angle = pi/2
    theta = max_angle / (n_cues - 1)
    
    targets_x[i] <- 1000 * cos((i-1) * theta - max_angle / 2)
    targets_y[i] <- 500 + 1000 * sin((i-1) * theta - max_angle / 2)
  }
}

# for 1D spin model only
tmp <- aggregate(df$time, list(df$temp, df$replicate), max)
df_out <- merge(df, tmp[,c("Group.1", "Group.2", "x")], by.x=c("temp", "replicate", "time"), by.y=c("Group.1", "Group.2", "x"))
head(df_out)

ggplot(df_out, aes(x = temp, y = V)) + geom_point(size = 0.5, alpha = 0.5) + theme_cus()

# for 2D spin model
df_low <- df[which(df$temp == 0.0),]
df_mid <- df[which(df$temp == 0.05),]
df_high <- df[which(df$temp == 0.1),]

# for 1D spin model only
ggplot(df_low, aes(x = time, y = x)) + geom_point(size = 0.5, alpha = 0.05) + theme_cus()
ggplot(df_mid, aes(x = time, y = x)) + geom_point(size = 0.5, alpha = 0.05) + theme_cus()
ggplot(df_high, aes(x = time, y = x)) + geom_point(size = 0.5, alpha = 0.05) + theme_cus()

ggplot(df_low, aes(x = V)) + geom_histogram(binwidth = 0.02) + theme_cus()
ggplot(df_mid, aes(x = V)) + geom_histogram(binwidth = 0.02) + theme_cus()
ggplot(df_high, aes(x = V)) + geom_histogram(binwidth = 0.02) + theme_cus()

# for 2D spin model
ggplot(df_low, aes(x = x, y = y)) + geom_point(size = 0.5, alpha = 0.2) + theme_cus() + xlim(-500, 1500) + ylim(-500, 1500)
ggplot(df_mid, aes(x = x, y = y)) + geom_point(size = 0.5, alpha = 0.2) + theme_cus() + xlim(-500, 1500) + ylim(-500, 1500)
ggplot(df_high, aes(x = x, y = y)) + geom_point(size = 0.5, alpha = 0.2) + theme_cus() + xlim(-500, 1500) + ylim(-500, 1500)

# for 2D spin model only
ggplot(df_low, aes(x = angle0, y = angle1)) + geom_point(size = 0.5, alpha = 0.2) + theme_cus()
ggplot(df_mid, aes(x = angle0, y = angle1)) + geom_point(size = 0.5, alpha = 0.2) + theme_cus()
ggplot(df_high, aes(x = angle0, y = angle1)) + geom_point(size = 0.5, alpha = 0.2) + theme_cus()

## Additional figures specific to 2D spin model
## Load data
df <- read.csv("decision_geometry_spin2d.csv")
head(df)

df_low <- df[which(df$temp == 0.0),]
df_mid <- df[which(df$temp == 0.05),]
df_high <- df[which(df$temp == 0.1),]

ggplot(df_low, aes(x = cue_reached)) + geom_histogram() + theme_cus()
ggplot(df_mid, aes(x = cue_reached)) + geom_histogram() + theme_cus()
ggplot(df_high, aes(x = cue_reached)) + geom_histogram() + theme_cus()
