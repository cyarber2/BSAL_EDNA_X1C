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





