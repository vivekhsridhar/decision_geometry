rm(list = ls())

## Load packages
library(ggplot2)
library(viridis)
library(reshape)
library(fields)
library(arm)

## Set working directory
setwd("/Users/viveksridhar/Documents/Code/multi-choice_decision_geometry/multi-choice_decision_geometry/output/Data/Agent paths/Symmetric")
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
df <- read.csv("decision_geometry_5.csv")
head(df)

## Clean data
df_clean <- df$replicate[which(df$good_run == 0)]
df <- df[!df$replicate %in% df_clean,]
df$pos_x <- df$pos_x - 500
df$pos_y <- df$pos_y - 500
df$cue_reached <- df$cue_reached + 1
df$cue_reached <- df$cue_reached[!!df$cue_reached][cumsum(!!df$cue_reached)+1]
df$cue_reached <- df$cue_reached - 2
df$rotated_x <- 0
df$rotated_y <- 0
df <- df[-nrow(df),]
summary(df)

ggplot(df, aes(x = pos_x, y = pos_y)) + geom_point(alpha = 0.05) + 
  xlim(-500,500) + ylim(-500,500) + theme_cus()

n_cues = length(unique(df$cue_reached))
for(i in 1:nrow(df)) {
  rotation <- matrix(c(cos(2*pi/n_cues * df$cue_reached[i]), -sin(2*pi/n_cues * df$cue_reached[i]), sin(2*pi/n_cues * df$cue_reached[i]), cos(2*pi/n_cues * df$cue_reached[i])), ncol = 2, byrow = TRUE)
  tmp <- c(df$pos_x[i], df$pos_y[i])
  rotated <- tmp%*%rotation
  df$rotated_x[i] <- rotated[1]
  df$rotated_y[i] <- rotated[2]
}

ggplot(df, aes(x = rotated_x, y = rotated_y)) + geom_point(alpha = 0.05) + 
  xlim(-10,500) + ylim(-200,200) + theme_cus()

df$rotated_x <- floor(df$rotated_x)
df$rotated_y <- abs(df$rotated_y)
t <- aggregate(df$rotated_y, list(df$rotated_x), median)

ggplot(t, aes(x = Group.1, y = x)) + geom_line() + 
  xlim(10,500) + ylim(-200,200) + theme_cus()

## Load data
df <- read.csv("group_properties.csv")
head(df)

df <- melt(tapply(df$accuracy, list(df$group_size, df$n_informed), mean))
df <- df[complete.cases(df), ]
names(df) <- c("group_size", "n_informed", "Accuracy")
df$prop <- df$n_informed/df$group_size
df$`Group Size` <- as.factor(df$group_size)
  
ggplot(df, aes(x = prop, y = Accuracy, colour = `Group Size`)) +
  geom_point() + geom_line() + theme_cus() + xlim(0,1) + ylim(0,1) + 
  xlab("Proportion of informed individuals")

## Load data
df <- read.csv("dGeom_fig3.csv")
head(df)

## Clean data
df_clean <- df$replicate[which(df$good_run == 0)]
df <- df[!df$replicate %in% df_clean,]
df$angle0 <- floor(df$angle0)
df$angle1 <- floor(df$angle1)
summary(df)

x <- rep(1:180,each=180)
y <- rep(1:180,180)
mat <- as.data.frame(cbind(x, y))
mat$count <- 0

for (i in 1:length(df$angle0)) {
  idx = df$angle0[i] * 180 + df$angle1[i]
  mat$count[idx] = mat$count[idx] + 1
}

maxs <- tapply(mat$count, mat$x, max)
for (i in 1:length(maxs)) {
  if (maxs[i] != 0) mat$count[which(mat$x == i)] = mat$count[which(mat$x == i)] / maxs[i]
}
  
ggplot(mat, aes(x = x, y = y, fill = count)) + geom_tile() + 
  scale_fill_viridis() + xlim(61,180) + theme_cus()

#Set working directory
setwd("/Users/viveksridhar/Documents/Code/multi-choice_decision_geometry/multi-choice_decision_geometry/output/")
dir <- getwd()

## Load data
df <- read.csv("cue_reached.csv")
head(df)

# Reformat data
n_cues = length(unique(df$cue_reached))
df <- aggregate(df$cue_reached, list(df$group_size, df$n_informed), length)
head(df)
names(df) <- c("Group Size", "n_informed", "number")
df$`P(Group splitting)` <- 1 - df$number / max(df$number)
df$n_uninformed <- df$`Group Size` -  n_cues * df$n_informed
df$`Prop. uninformed` <- as.factor(df$n_uninformed / df$`Group Size`)
head(df)

ggplot(df, aes(x = `Group Size`, y = `P(Group splitting)`, colour = `Prop. uninformed`)) + 
  geom_point() + geom_line() + theme_cus() + ylim(0, 1)

## Load data
df <- read.csv("cue_reached.csv")
head(df)

# Reformat and clean data
df$numbers <- 1
df <- aggregate(df$numbers, list(df$cue_reached, df$group_size, df$n_informed), length)
names(df) <- c("Cue ID", "Group Size", "n_informed", "number")
head(df)

n_cues = length(unique(df$`Cue ID`))
df$n_uninformed <- df$`Group Size` -  n_cues * df$n_informed
df$`Prop. uninformed` <- as.factor(df$n_uninformed / df$`Group Size`)

df <- df[order(df$`Group Size`, df$n_informed, df$`Cue ID`),]
sums <- aggregate(df$number, list(df$`Group Size`, df$n_informed), sum)
sums <- sums[order(sums$Group.1, sums$Group.2),]

df$total <- rep(sums$x, each=3)
df$`Group Size` <- as.factor(df$`Group Size`)
summary(df)

ggplot(df, aes(x = `Cue ID`, y = number/total, shape = `Prop. uninformed`, colour = `Group Size`)) + 
  geom_point() + geom_line() + theme_cus() + 
  ylim(0, 1) + ylab("Prop. of times cue was chosen")

df1 <- aggregate(df$number, list(df$`Prop. uninformed`, df$`Cue ID`), sum)
names(df1) <- c("Prop. uninformed", "Cue ID", "number")
df1$total <- rep(aggregate(df1$number, list(df1$`Prop. uninformed`), sum)$x, 3)
ggplot(df1, aes(x = `Cue ID`, y = number/total, colour = `Prop. uninformed`)) + 
  geom_point() + geom_line() + theme_cus() + 
  ylim(0, 1) + ylab("Prop. of times cue was chosen")

df2 <- aggregate(df$number, list(df$`Group Size`, df$`Cue ID`), sum)
names(df2) <- c("Group Size", "Cue ID", "number")
df2$`Group Size` <- as.factor(df2$`Group Size`)
df2$total <- rep(aggregate(df2$number, list(df2$`Group Size`), sum)$x, 3)
ggplot(df2, aes(x = `Cue ID`, y = number/total, colour = `Group Size`)) + 
  geom_point() + geom_line() + theme_cus() + 
  ylim(0, 1) + ylab("Prop. of times cue was chosen")

