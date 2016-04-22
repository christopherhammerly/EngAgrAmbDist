#####################################################
###     Agreement ambiguity analysis pipeline     ###
###     Christopher Hammerly                      ###
###     UMass Amherst - 04.21.16                  ###
#####################################################

#### Load libraries

library(MASS)
library(tikzDevice)


#### Set working directory

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

#### segregate acceptability judgments and reading time data

data.acceptability <- subset(agr.amb, TrialType == 'Question')
data.RT <- subset(agr.amb, TrialType == 'DashedSentence')




