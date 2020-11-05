### ED Figure 4a
load("RData/Dataset1_tsneData.RData", verbose = T)
tmp <- tsneData[!duplicated(tsneData$TumorID), ]
tumorCol <- tmp$TumorCol
names(tumorCol) <- tmp$TumorID

library(ggplot2)
tumorID <- factor(tsneData$TumorID)
d <- ggplot(tsneData, aes(x = tSNE_1, y = tSNE_2, color = tumorID)) + 
  geom_point(size = 1.5) +
  scale_color_manual(values = tumorCol[levels(tumorID)]) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  labs(x = "tSNE 1", y = "tSNE 2") +
  theme_classic(base_size = 15) + theme(legend.position = "none")
d


### ED Figure 4b
load("RData/Dataset2_tsneData.RData", verbose = T)
tmp <- tsneData[!duplicated(tsneData$TumorID), ]
tumorCol <- tmp$TumorCol
names(tumorCol) <- tmp$TumorID

library(ggplot2)
tumorID <- factor(tsneData$TumorID)
d <- ggplot(tsneData, aes(x = tSNE_1, y = tSNE_2, color = tumorID)) + 
  geom_point(size = 1.5) +
  scale_color_manual(values = tumorCol[levels(tumorID)]) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  labs(x = "tSNE 1", y = "tSNE 2") +
  theme_classic(base_size = 15) + theme(legend.position = "none")
d


### ED Figure 4c
load("RData/Dataset3_tsneData.RData", verbose = T)
tmp <- tsneData[!duplicated(tsneData$TumorID), ]
tumorCol <- tmp$TumorCol
names(tumorCol) <- tmp$TumorID

library(ggplot2)
tumorID <- factor(tsneData$TumorID)
d <- ggplot(tsneData, aes(x = tSNE_1, y = tSNE_2, color = tumorID)) + 
  geom_point(size = 1.5) +
  scale_color_manual(values = tumorCol[levels(tumorID)]) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  labs(x = "tSNE 1", y = "tSNE 2") +
  theme_classic(base_size = 15) + theme(legend.position = "none")
d


### ED Figure 4d
load("RData/Dataset1_tsneData.RData", verbose = T)
tmp <- tsneData[!duplicated(tsneData$concordanceClassCol), ]
classCol <- tmp$concordanceClassCol
names(classCol) <- classCol

library(ggplot2)
class <- factor(tsneData$concordanceClassCol)
d <- ggplot(tsneData, aes(x = tSNE_1, y = tSNE_2, color = class)) + 
  geom_point(size = 1.5) +
  scale_color_manual(values = classCol[levels(class)]) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  labs(x = "tSNE 1", y = "tSNE 2") +
  theme_classic(base_size = 15) + theme(legend.position = "none")
d


### ED Figure 4e
load("RData/Dataset2_tsneData.RData", verbose = T)
tmp <- tsneData[!duplicated(tsneData$concordanceClassCol), ]
classCol <- tmp$concordanceClassCol
names(classCol) <- classCol

library(ggplot2)
class <- factor(tsneData$concordanceClassCol)
d <- ggplot(tsneData, aes(x = tSNE_1, y = tSNE_2, color = class)) + 
  geom_point(size = 1.5) +
  scale_color_manual(values = classCol[levels(class)]) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  labs(x = "tSNE 1", y = "tSNE 2") +
  theme_classic(base_size = 15) + theme(legend.position = "none")
d


### ED Figure 4f
load("RData/Dataset3_tsneData.RData", verbose = T)
tmp <- tsneData[!duplicated(tsneData$concordanceClassCol), ]
classCol <- tmp$concordanceClassCol
names(classCol) <- classCol

library(ggplot2)
class <- factor(tsneData$concordanceClassCol)
d <- ggplot(tsneData, aes(x = tSNE_1, y = tSNE_2, color = class)) + 
  geom_point(size = 1.5) +
  scale_color_manual(values = classCol[levels(class)]) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  labs(x = "tSNE 1", y = "tSNE 2") +
  theme_classic(base_size = 15) + theme(legend.position = "none")
d