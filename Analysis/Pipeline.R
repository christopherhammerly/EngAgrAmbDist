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

library(dplyr)
library(ggplot2)
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

mycols <- c("Subject","MD5","TrialType","Number","Element","Experiment","Item", "question", "response","null","RT","null")
results <- read.csv('/Users/chrishammerly/EngAgrAmbDist/Data/results.csv',
                    header = 0, 
                    sep = ",",
                    quote = "",
                    comment.char = "#",
                    col.names=mycols,
                    fill = TRUE)

#### segregate ambiguity data from fillers and other experiment

data <- subset(results , Experiment %in% c('agr-amb-a','agr-amb-b','agr-amb-c','agr-amb-d'))

### Adding columns for the factors: attachment (low, high) by attractor number (singular, plural)

data$attachment <- ifelse(data$Experiment=='agr-amb-b' | data$Experiment=='agr-amb-c','high','low')
data$attractor <- ifelse(data$Experiment=='agr-amb-a' | data$Experiment=='agr-amb-c','singular','plural')


#### segregate acceptability judgments and reading time data

data.acceptability <- subset(data, TrialType == 'Question')
data.RT <- subset(data, TrialType == 'DashedSentence')

print(data.acceptability %>% summarise(number = n_distinct(Subject)))

#### change the values for the judgments and RT from characters to numbers 

data.acceptability$response <- as.numeric(as.character(data.acceptability$response))
data.RT$X <- as.numeric(as.character(data.RT$X))

#### Basic descriptive stats

subj.by.cond <- data.acceptability %>% 
    group_by(Subject, Experiment) %>% 
    summarise(mean.resp = mean(response),
              SEM = sd(response)/sqrt(n_distinct(Experiment)))

cond.summ <- subj.by.cond %>%
  group_by(Experiment) %>%
  summarise(mean_cond = mean(mean),
            SEM = sd(mean)/sqrt(n_distinct(Subject)))

subj.by.factor <- data.acceptability %>% 
  group_by(Subject, attractor, attachment) %>% 
  summarise(mean = mean(response))

factor.summ <- subj.by.factor %>%
  group_by(attachment, attractor) %>%
  summarise(mean_cond = mean(mean),
            SEM = sd(mean)/sqrt(n_distinct(Subject)))

#### Hypothesis testing

ezANOVA(data = data.acceptability, dv = response, wid = Subject, within = c(attractor, attachment))

    #$ANOVA
    #                Effect DFn DFd          F   p p<.05                ges
    #2            attractor   1  47 2.62525336 0.1118668       6.829641e-03
    #3           attachment   1  47 0.01401451 0.9062688       3.681569e-05
    #4 attractor:attachment   1  47 0.05863566 0.8097190       1.022591e-04

