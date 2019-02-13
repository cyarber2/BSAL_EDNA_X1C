rm(list=ls())
setwd("~/Desktop/MSResearch/BSAL_EDNA_X1C")

library("readxl")

Plate_11_raw <- read_excel(path="Bsal_Plate_11_11132018_resave.xls",
                           sheet=1,
                           skip=19)
Plate_12_raw <- read_excel(path="Bsal_Plate_12_11152018_resave.xls",
                           sheet=1,
                           skip=19)
Plate_13_raw <- read_excel(path="Bsal_Plate_13_12122018_resave.xls",
                           sheet=1,
                           skip=19)
Plate_14_raw <- read_excel(path="Bsal_Plate_14_12122018_resave.xls",
                           sheet=1,
                           skip=19)

Plate_17_raw <- read_excel(path="Bsal_Plate_17_01092019_resave.xls",
                           sheet=1,
                           skip=19)
Plate_17 <- Plate_17_raw[!grepl(pattern="X2.", x=Plate_17_raw$Sample), ]
Plate_17$Sample <- gsub(x=Plate_17$Sample, pattern="X1C.", replacement="")

Bsal <- rbind(Plate_11_raw, Plate_12_raw, Plate_13_raw, Plate_14_raw, Plate_17)

Bsal <- subset(Bsal, select=c(1,4,5,6,7), subset=Target=="Bsal")
colnames(Bsal) <- c("Well", "Target", "Sample", "Starting_Quantity", "Cq")
Bsal <- Bsal[!grepl(pattern="gBlock", x=Bsal$Sample), ]
Bsal <- Bsal[!grepl(pattern="NEG", x=Bsal$Sample, ignore.case=T), ]
Bsal.agg <- aggregate(x=Bsal$Starting_Quantity,
                      by=list(Sample=Bsal$Sample),
                      FUN=mean,
                      na.rm=T)
colnames(Bsal.agg)[2] <- "Estimated_Copies"
# Manually remove 1st run when samples are rerun for concensus
Bsal.agg <- Bsal.agg[-c(7, 18, 21, 24, 28, 33, 35, 37, 65, 83), ]
o <- unlist(regmatches(Bsal.agg$Sample, gregexpr("[[:digit:]]+", Bsal.agg$Sample)))
my.order <- order(as.numeric(o))
Bsal.agg <- Bsal.agg[my.order, ]
Bsal.agg$Sample <- gsub(x=Bsal.agg$Sample, pattern=".C", replace="", fixed=T)

key <- read_excel(path="X1C_Key.xls",
                  sheet=1)

final <- merge(Bsal.agg, key, by="Sample")
flag <- which(!final$Status_ByHand)
final[flag, "Estimated_Copies"] <- 0

z.power <- c(1, 1.5, 2, 2.5, 3)
prob <- c(0.6, 0.85, 0.95, 1, 0.95)

plot(x=z.power, y=prob,
     xlim=c(-0.25, 7.25),
     ylim=c(0, 1),
     xlab="Number of Zoospores Processed (log10 scale)",
     ylab="Detection Probability")

load("X1B_FinalFilterData.Rdata")

composite.positives <- c(final$Status_ByHand, f.final$Positive)
composite.powers <- c(final$Zoospore_Power, f.final$power.z)

LOD <- data.frame(Pos=composite.positives, Z=as.numeric(composite.powers))
LOD <- LOD[which(LOD$Z!="control"), ]

# Changed the graph so it shows proportion of samples positive
library(tidyverse)

# First group and calculate proportion and number of samples
LOD %>%
  group_by(Z) %>%
  summarise(Proportion = mean(Pos),
            N = n()) %>%
  ggplot(data=., aes(x=Z, y=Proportion)) + # feed this grouped data into ggplot
  geom_smooth(data=LOD, aes(x=Z, y=Pos), # unfortunately, need to orginal ungrouped
              inherit.aes = FALSE,       # data to get regression to work
              method=glm, method.args=list(family=binomial(link="probit"))) +
  geom_point(aes(size = N)) + # add in the points
  geom_hline(yintercept=0.95) + # add a horizontal line
  geom_segment(aes(x=log10(199), xend=log10(199), y=0.95, yend=0)) + # extend vertical segment at LOD
  scale_size_area(breaks = c(5,20,25)) +
  xlab("Number of Zoospores Filtered (log10)") +
  ylab("Detection Probability")

# Previous version
ggplot(data=LOD, aes(x=Z, y=Pos)) +
  stat_sum() +
  xlab("Number of Zoospores Filtered (log10)") +
  ylab("Detection Probability") +
  geom_smooth(method=glm, method.args=list(family=binomial(link="probit"))) +
  scale_size_area()

LODProbit <- glm(Pos ~ Z,
                data=LOD,
                family=binomial(link="probit"))
x.index <- seq(0, 3, length.out=5000)
predLOD <- predict(LODProbit,
                   data.frame(Z=x.index),
                   type="response")
LOD_95 <- x.index[which.min(abs(predLOD - 0.95)) ]

# 95% LOD for this Bsal assay is 10^LOD_95 ~ 199 zoospores



