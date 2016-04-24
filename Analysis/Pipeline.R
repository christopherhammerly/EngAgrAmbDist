#####################################################################
###                                                               ###
###             Agreement ambiguity analysis pipeline             ###
###             Christopher Hammerly                              ###
###             UMass Amherst - 04.22.16                          ###
###                                                               ###
#####################################################################

### Condition key for convenient reference: order is N1 N2 V; S = singular, P = plural
# agr-amb-a = SPP = low, singular
# agr-amb-b = SPS = high, plural
# agr-amb-c = PSP = high, singular
# agr-amb-d = PSS = low, plural

#### Load libraries

library(MASS)
library(dplyr)
library(ggplot2)
library(tikzDevice)
library(downloader)
library(ez)

#### Set working directory
#### Eventually, I will add a github link so people will be able to download the data directly and run analyses seamlessly.
#### Here is a template for doing that:
#url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/msleep_ggplot2.csv"
#filename <- "msleep_ggplot2.csv"
#if (!file.exists(filename)) download(url,filename)

setwd("/Users/chrishammerly")

#### Load data

mycols <- c("Subject","MD5","TrialType","Number","Element","Experiment","Item","Sentence","Response","X","RT")
results <- read.csv('/Users/chrishammerly/EngAgrAmbDist/Data/results.csv',
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

subj.by.cond <- data.acceptability %>% 
    group_by(Subject, Experiment) %>% 
    summarise(mean = mean(Response))

cond.summ <- subj.by.cond %>%
  group_by(Experiment) %>%
  summarise(mean_cond = mean(mean),
            SEM = sd(mean)/sqrt(n_distinct(Subject)))

subj.by.factor <- data.acceptability %>% 
  group_by(Subject, attractor, attachment) %>% 
  summarise(mean = mean(Response))

factor.summ <- subj.by.factor %>%
  group_by(attachment, attractor) %>%
  summarise(mean_cond = mean(mean),
            SEM = sd(mean)/sqrt(n_distinct(Subject)))

attractor.summ <- subj.by.factor %>%
  group_by(attractor) %>%
  summarise(mean_cond = mean(mean),
            SEM = sd(mean)/sqrt(n_distinct(Subject)))

attachment.summ <- subj.by.factor %>%
  group_by(attachment) %>%
  summarise(mean_cond = mean(mean),
            SEM = sd(mean)/sqrt(n_distinct(Subject)))

#### Hypothesis testing

# something seems off about this one
ezANOVA(data = data.frame(subj.by.cond), dv = mean, wid = Subject, within = Experiment)

#maybe this is right?
ezANOVA(data = data.frame(subj.by.factor), dv = mean, wid = Subject, within = c(attractor, attachment))

# here are some attempts with the old method, but it comes up with different results than exANOVA
anova <- aov(Response ~ attachment*attractor, data=data.acceptability)
summary(anova)

anova2 <- aov(mean ~ attachment*attractor, data=subj.by.factor)
summary(anova2)

