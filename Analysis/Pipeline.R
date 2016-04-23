#####################################################################
###                                                               ###
###             Agreement ambiguity analysis pipeline             ###
###             Christopher Hammerly                              ###
###             UMass Amherst - 04.22.16                          ###
###                                                               ###
#####################################################################

### Condition key for convenient reference: order is N1 N2 V; S = singular, P = plural
# agr-amb-a = SPP
# agr-amb-b = SPS
# agr-amb-c = PSP
# agr-amb-d = PSS

#### Load libraries

library(MASS)
library(dplyr)
library(ggplot2)
library(tikzDevice)
library(downloader)

#### Set working directory
#### Eventually, I will add a github link so people will be able to download the data directly and run analyses seamlessly.
#### Here is a template for doing that:
#url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/msleep_ggplot2.csv"
#filename <- "msleep_ggplot2.csv"
#if (!file.exists(filename)) download(url,filename)

setwd("/Users/chrishammerly")

#### Load data

mycols <- c("Subject","MD5","TrialType","Number","Element","Experiment","Item","Sentence","Response","X","RT")
results <- read.csv('/Users/chrishammerly/Desktop/agr_amb_dist/results/fake_results.csv',
                    header = 0, 
                    sep = ",", 
                    comment.char = "#",
                    col.names=mycols)

#### segregate ambiguity data from fillers and other experiment

data <- subset(results , Experiment %in% c('agr-amb-a','agr-amb-b','agr-amb-c','agr-amb-d'))

### Adding columns for the factors: attachment (low, high) by attractor number (singular, plural)

data$attachment <- ifelse(data$Experiment=='agr-amb-b' | data$Experiment=='agr-amb-c','high','low')
data$attractor <- ifelse(data$Experiment=='agr-amb-a' | data$Experiment=='agr-amb-c','singular','plural')


#### segregate acceptability judgments and reading time data

data.acceptability <- subset(data, TrialType == 'Question')
data.RT <- subset(data, TrialType == 'DashedSentence')

#### change the values for the judgments and RT from characters to numbers 

data.acceptability$Response <- as.numeric(as.character(data.acceptability$Response))
data.RT$X <- as.numeric(as.character(data.RT$X))

#### Basic descriptive stats

mean(data.acceptability$Response)
tapply(data.acceptability$Response, data.acceptability$attachment, mean)
tapply(data.acceptability$Response, data.acceptability$attachment, sd)
tapply(data.acceptability$Response, data.acceptability$attractor, mean)
tapply(data.acceptability$Response, data.acceptability$attractor, sd)

#### Hypothesis testing

aov(Response ~ attachment*attractor, data=data.acceptability)

t.test(data.acceptability$Response~data.acceptability$attachment)
t.test(data.acceptability$Response~data.acceptability$attractor)



