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
library(lme4)

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

data.acceptability <- droplevels(subset(data, TrialType == 'Question'))

data.RT <- subset(data, TrialType == 'DashedSentence')
names(data.RT) <- c("Subject","MD5","TrialType","Number","Element","Experiment","Item", "region", "fragment","RT","null","sentence", "attachment", "attractor")

#######################################
###                                 ###
###     Acceptability analysis      ###
###                                 ###
#######################################

#### a sanity check to ensure the correct number of subjects have been processed through
print(data.acceptability %>% summarise(number = n_distinct(Subject)))

#### change the values for the judgments from characters to numbers 

data.acceptability$response <- as.numeric(as.character(data.acceptability$response))
data.acceptability$RT <- as.numeric(as.character(data.acceptability$RT))


data.acceptability$z <- ave(data.acceptability$response, data.acceptability$Subject, FUN = scale)

#### Basic descriptive stats

subj.by.cond <- data.acceptability %>% 
    group_by(Subject, Experiment) %>% 
    summarise(mean.resp = mean(response), 
              mean.z = mean(z))

cond.summ <- subj.by.cond %>%
  group_by(Experiment) %>%
  summarise(mean_cond = mean(mean.resp),
            SEM = sd(mean.resp)/sqrt(n_distinct(Subject)),
            average.z = mean(mean.z))

subj.by.factor <- data.acceptability %>% 
  group_by(Subject, attractor, attachment) %>% 
  summarise(mean = mean(response))

factor.summ <- subj.by.factor %>%
  group_by(attachment, attractor) %>%
  summarise(mean_cond = mean(mean),
            SEM = sd(mean)/sqrt(n_distinct(Subject)))

subj.by.factor.high <- subset(subj.by.factor, attachment == 'high') %>%
  group_by(attachment) %>%
  summarise(mean_cond = mean(mean),
            SEM = sd(mean)/sqrt(n_distinct(Subject)))

subj.by.factor.low <- subset(subj.by.factor, attachment == 'low') %>%
  group_by(attachment) %>%
  summarise(mean_cond = mean(mean),
            SEM = sd(mean)/sqrt(n_distinct(Subject)))


#### Hypothesis testing

ezANOVA(data = data.acceptability, dv = response, wid = Subject, within = c(attractor, attachment))

    #$ANOVA
    #                Effect DFn DFd          F   p p<.05                ges
    #2            attractor   1  47 2.62525336 0.1118668       6.829641e-03
    #3           attachment   1  47 0.01401451 0.9062688       3.681569e-05
    #4 attractor:attachment   1  47 0.05863566 0.8097190       1.022591e-04

#tapply(region2$RT, list(region1$Item, region1$Experiment), mean)
#tapply(region1$RT, list(region1$Subject, region1$Experiment), mean)

###

truehist(data.acceptability$RT)

subj.by.cond.judge.RT <- data.acceptability %>% 
  group_by(Subject, Experiment) %>% 
  summarise(mean.RT = mean(RT))

cond.summ.judge.RT <- subj.by.cond.judge.RT %>%
  group_by(Experiment) %>%
  summarise(mean_cond = mean(mean.RT),
            SEM = sd(mean.RT)/sqrt(n_distinct(Subject)))

tapply(data.acceptability$RT, list(data.acceptability$Experiment), mean)

ezANOVA(data = data.acceptability, dv = RT, wid = Subject, within = c(attractor, attachment))



#########################################
###                                   ###
###     Raw Reading time analysis     ###
###                                   ###
#########################################

#### Sanity check for subject number

print(data.RT %>% summarise(number = n_distinct(Subject)))

#### Ensure RT data is numeric

data.RT$RT <- as.numeric(as.character(data.RT$RT))

#### Descriptive stats

RT.subj.by.cond <- data.RT %>%
  group_by(Subject, Experiment, region) %>%
  summarise(average = mean(RT))

RT.cond.summ <- RT.subj.by.cond %>%
  group_by(Experiment, region) %>%
  summarise(mean = mean(average),
            SEM = sd(average)/sqrt(n_distinct(Subject)))

#### ANOVA by region
region1 <- droplevels(subset(data.RT, region == "1"))
region2 <- droplevels(subset(data.RT, region == "2"))
region3 <- droplevels(subset(data.RT, region == "3"))
region4 <- droplevels(subset(data.RT, region == "4"))

ezANOVA(data = region1, dv = RT, wid = Subject, within = c(attractor, attachment))

ezANOVA(data = region2, dv = RT, wid = Subject, within = c(attractor, attachment))

ezANOVA(data = region3, dv = RT, wid = Subject, within = c(attractor, attachment))

ezANOVA(data = region4, dv = RT, wid = Subject, within = c(attractor, attachment))

#### Plotting RT data

pdf('RTplotRaw.pdf')
ggplot(subset(RT.cond.summ,Experiment %in% c("agr-amb-a","agr-amb-b","agr-amb-c","agr-amb-d")),aes(x=region,y=mean,color=Experiment,base=6,group=Experiment))+ 
  labs(y="Reading time",x="Region",group=1) +geom_point(stat = "identity",size=1)+
  geom_errorbar(aes(ymax = mean+SEM,ymin=mean-SEM,width=0.05))+ 
  theme(text = element_text(size=10))+stat_identity(geom="line")+
  scale_colour_manual(values = c("steelblue2","firebrick4","navyblue","firebrick2"),
                      name="Conditions",
                      labels= c("high, singular","low, singular", "high, plural", "low, plural"),
                      breaks= c("agr-amb-c",'agr-amb-a','agr-amb-b','agr-amb-d'))
dev.off()

#####################################
###                               ###
###      Residual RT analysis     ###
###                               ###
#####################################

data.model <- subset(results, TrialType == 'DashedSentence')
names(data.model) <- c("Subject","MD5","TrialType","Number","Element","Experiment","Item", "region", "fragment","RT","null","sentence")

data.model$fragment <- as.character(data.model$fragment)

data.model$char.length <- nchar(data.model$fragment)

data.model$char.length <- as.numeric(as.character(data.model$char.length))
data.model$RT <- as.numeric(as.character(data.model$RT))

resid.model <- lmer(RT ~ char.length + (1|Subject), data.model)

data.model$resid <- resid(resid.model)

data.model.cond <- data.model %>%
  group_by(Subject, Experiment, region) %>%
  summarise(average = mean(resid))

data.model.cond.summ <- data.model.cond %>%
  group_by(Experiment, region) %>%
  summarise(mean = mean(average),
            SEM = sd(average)/sqrt(n_distinct(Subject)))

#### ANOVA PREP...THIS NEEDS TO BE MAJORLY CLEANED UP SO ALL OF THIS ADDING OF COLUMNS IS FRONT-ENDED

data.model.RT <- droplevels(subset(data.model , Experiment %in% c('agr-amb-a','agr-amb-b','agr-amb-c','agr-amb-d')))
data.model.RT$attachment <- ifelse(data.model.RT$Experiment=='agr-amb-b' | data.model.RT$Experiment=='agr-amb-c','high','low')
data.model.RT$attractor <- ifelse(data.model.RT$Experiment=='agr-amb-a' | data.model.RT$Experiment=='agr-amb-c','singular','plural')

region1.resid <- droplevels(subset(data.model.RT, region == "1"))
region2.resid <- droplevels(subset(data.model.RT, region == "2"))
region3.resid <- droplevels(subset(data.model.RT, region == "3"))
region4.resid <- droplevels(subset(data.model.RT, region == "4"))

#### ACTUAL ANOVAs

ezANOVA(data = region1.resid, dv = resid, wid = Subject, within = c(attractor, attachment))

ezANOVA(data = region2.resid, dv = resid, wid = Subject, within = c(attractor, attachment))

ezANOVA(data = region3.resid, dv = resid, wid = Subject, within = c(attractor, attachment))

ezANOVA(data = region4.resid, dv = resid, wid = Subject, within = c(attractor, attachment))

pdf('RTplotResidual.pdf')
ggplot(subset(data.model.cond.summ,Experiment %in% c("agr-amb-a","agr-amb-b","agr-amb-c","agr-amb-d")),aes(x=region,y=mean,color=Experiment,base=6,group=Experiment))+ 
  labs(y="Reading time",x="Region",group=1) +geom_point(stat = "identity",size=1)+
  geom_errorbar(aes(ymax = mean+SEM,ymin=mean-SEM,width=0.05))+ 
  theme(text = element_text(size=10))+stat_identity(geom="line")+
  scale_colour_manual(values = c("steelblue2","firebrick4","navyblue","firebrick2"),
                      name="Conditions",
                      labels= c("high, singular","low, singular", "high, plural", "low, plural"),
                      breaks= c("agr-amb-c",'agr-amb-a','agr-amb-b','agr-amb-d'))
dev.off()

#                        labels= c("high, singular","low, singular", "high, plural", "low, plural"),
#                        breaks= c("agr-amb-c",'agr-amb-a','agr-amb-b','agr-amb-d'))



