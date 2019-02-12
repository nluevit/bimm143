#' ---
#'title: "Class 05 R Graphics Intro"
#'author: Nasha Luevit
#'date: Jan 22nd, 2019
#' ---

# Class 05 R Graphics Intro
#' this is some text and I can have **bold** and **italic** and `code`

# my first boxplot 
x <- rnorm(1000,0)
boxplot(x)

summary(x)
hist(x)

boxplot(x, horizontal = TRUE)

# Hands On Session

#import weight chart into R
#read.table("weight_chart.txt", header = TRUE)
#read.table("weight_chart.txt", header = FALSE)
weight <- read.table("bimm143_05_rstats/weight_chart.txt", header = TRUE)
#graph weight chart with changes
#plot(weight)
plot(weight, typ = "o", pch = 15, cex = 1.5, lwd = 2, ylim = c(2,10), xlab = "Age (months)", ylab = "Weight (kg)", main = "Baby Weight with Age")

#import feature counts
#read.table("feature_counts.txt", header = TRUE, sep = "\t")
feature_counts <- read.table("bimm143_05_rstats/feature_counts.txt", header = TRUE, sep = "\t")
#plot feature counts
barplot(height = feature_counts[,2])
barplot(feature_counts$Count)
barplot(feature_counts$Count, horiz = TRUE, names.arg = feature_counts$Feature, main = "Number of Features in the Mouse GRCm38 Genome", las = 1, xlim = c(0,80000))
#change margins so we can see the labels 
par()$mar
par(mar = c(3.1, 11.1, 4.1, 2))
barplot(feature_counts$Count, horiz = TRUE, names.arg = feature_counts$Feature, main = "Number of Features in the Mouse GRCm38 Genome", las = 1, xlim = c(0,80000))
#import male female counts
#read.table("male_female_counts.txt", header = TRUE, sep = "\t")
male_female_counts <- read.table("bimm143_05_rstats/male_female_counts.txt", header = TRUE, sep = "\t")
#plot male female counts
barplot(male_female_counts$Count, names.arg = male_female_counts$Sample, ylab = "Counts", col = rainbow(10))
barplot(male_female_counts$Count, names.arg = male_female_counts$Sample, ylab = "Counts", col = rainbow(nrow(male_female_counts)))
#plot with color by gender
barplot(male_female_counts$Count, names.arg = male_female_counts$Sample, ylab = "Counts", col = c("blue2", "pink2"))

#import up down expression
#read.table("up_down_expression.txt", header = TRUE, sep = "\t")
up_down_expression <- read.table("bimm143_05_rstats/up_down_expression.txt", header = TRUE, sep = "\t")
#plot up down expression
plot.default(up_down_expression)
## determine how many genes are up, down, or unchanging
table(up_down_expression$State)
## plot w color by up, down, unchanging 
plot.default(up_down_expression$Condition1, up_down_expression$Condition2, ylab = "Expression Condition 2", xlab = "Expression Condition 1", col = up_down_expression$State)
palette()
levels(up_down_expression$State)
palette(c("blue", "gray", "red"))
plot.default(up_down_expression$Condition1, up_down_expression$Condition2, ylab = "Expression Condition 2", xlab = "Expression Condition 1", col = up_down_expression$State)

