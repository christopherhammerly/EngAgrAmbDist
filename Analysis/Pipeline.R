#####################################################################
###                                                               ###
###             Agreement ambiguity analysis pipeline             ###
###             Christopher Hammerly                              ###
###             UMass Amherst - 05.12.16                          ###
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
library(gridExtra)

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

###################
###  JUDGMENTS  ###
###################

#### change the values for the judgments from characters to numbers 

data.acceptability$response <- as.numeric(as.character(data.acceptability$response))
data.acceptability$RT <- as.numeric(as.character(data.acceptability$RT))

pdf('agr-amb-alljudgments.pdf')
ggplot(data.acceptability,aes(x=response))+geom_histogram(binwidth=1)+xlim(1,7)+ggtitle("JUDGMENT DISTRIBUTION")
dev.off()

pdf('agr-amb-alljudgementRT.pdf')
ggplot(data.acceptability,aes(x=RT))+geom_histogram(binwidth=50)+xlim(0,7000)+ggtitle("JUDGEMENT RT DISTRIBUTION")
dev.off()


#### Filter our the abnormally long responses

data.acceptability <- filter(data.acceptability, RT < 4000)

#### Basic descriptive stats

subj.by.cond <- data.acceptability %>% 
    group_by(Subject, attachment, attractor) %>% 
    summarise(mean.resp = mean(response))

cond.summ <- subj.by.cond %>%
  group_by(attachment, attractor) %>%
  summarise(mean_cond = mean(mean.resp),
            SEM = sd(mean.resp)/sqrt(n_distinct(Subject)))

#### calculating difference scores and confidence intervals

differences.by.subj <- data.acceptability %>% 
  group_by(Subject,attractor) %>%
  filter(attachment=='high') %>%
  summarise(high.mean = mean(response))

differences.by.subj <- data.acceptability %>% 
  group_by(Subject,attractor) %>%
  filter(attachment=='low') %>%
  summarise(low.mean = mean(response)) %>%
  right_join(differences.by.subj)

differences.summary <- differences.by.subj                                  %>%
  group_by(attractor)                                                       %>%
  summarise(
    N = n_distinct(Subject),
    mean.diff = mean(low.mean-high.mean),
    mean.sem = sd(low.mean-high.mean)/sqrt(n_distinct(Subject))) %>%
  mutate(ci.lower.mean = mean.diff - qt(.975,df=N-1)*mean.sem,
         ci.upper.mean = mean.diff + qt(.975,df=N-1)*mean.sem
  )

differences.by.subj.num <- data.acceptability %>% 
  group_by(Subject,attachment) %>%
  filter(attractor=='singular') %>%
  summarise(singular.mean = mean(response))

differences.by.subj.num <- data.acceptability %>% 
  group_by(Subject,attachment) %>%
  filter(attractor=='plural') %>%
  summarise(plural.mean = mean(response)) %>%
  right_join(differences.by.subj.num)

differences.summary.num <- differences.by.subj.num                          %>%
  group_by(attachment)                                                     %>%
  summarise(
    N = n_distinct(Subject),
    mean.diff = mean(singular.mean-plural.mean),
    mean.sem = sd(singular.mean-plural.mean)/sqrt(n_distinct(Subject))) %>%
  mutate(ci.lower.mean = mean.diff - qt(.975,df=N-1)*mean.sem,
         ci.upper.mean = mean.diff + qt(.975,df=N-1)*mean.sem
  )

#### Combined plots of difference scores and raw ratings

diff.plot <- ggplot(differences.summary,aes(x=attractor,y=mean.diff,fill=attractor))+geom_bar(position='dodge',stat = "identity")+geom_linerange(stat='identity',mapping=aes(ymax = ci.lower.mean,ymin=ci.upper.mean,width=.25))+ scale_fill_manual(values=c("navy","maroon","grey"))+xlab("ATTRACTOR TYPE")+ylab("HIGH ATTACHMENT PENALTY")
diff.plot.num <- ggplot(differences.summary.num,aes(x=attachment,y=mean.diff,fill=attachment))+geom_bar(position='dodge',stat = "identity")+geom_linerange(stat='identity',mapping=aes(ymax = ci.lower.mean,ymin=ci.upper.mean,width=.25))+ scale_fill_manual(values=c("navy","maroon","grey"))+xlab("ATTACHMENT TYPE")+ylab("PLURAL PENALTY")
raw.plot <- ggplot(cond.summ,aes(x=attractor,y=mean_cond,fill=attachment))+geom_bar(position='dodge',stat = "identity")+geom_linerange(stat='identity',position = position_dodge(width = 0.9),mapping=aes(ymax = mean_cond+SEM,ymin=mean_cond-SEM)) + scale_fill_manual(values=c("navy","maroon"))+xlab("ATTRACTOR TYPE")+ylab("MEAN RATING")+ylim(0,7)

pdf('agr-amb-judgmentplot.pdf')
grid.arrange(raw.plot, diff.plot, diff.plot.num,ncol=3)
dev.off()

#### Some extra descriptive stuff

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

###      $ANOVA
###      Effect DFn DFd           F         p p<.05          ges
###      2            attractor   1  47 2.354599742 0.1316187       5.907001e-03
###      3           attachment   1  47 0.016589218 0.8980662       4.595903e-05
###      4 attractor:attachment   1  47 0.006684451 0.9351860       1.177561e-05


###################
###     RT      ###
###################

#### Descriptive 

subj.by.cond.judge.RT <- data.acceptability %>% 
  group_by(Subject, attractor, attachment) %>% 
  summarise(mean.RT = mean(RT))

cond.summ.judge.RT <- subj.by.cond.judge.RT %>%
  group_by(attractor, attachment) %>%
  summarise(mean_cond = mean(mean.RT),
            SEM = sd(mean.RT)/sqrt(n_distinct(Subject)))

#### Plots 

pdf('agr-amb-judgmentRTplot.pdf')
ggplot(cond.summ.judge.RT,aes(x=attractor,y=mean_cond,fill=attachment))+geom_bar(position='dodge',stat = "identity")+geom_linerange(stat='identity',position = position_dodge(width = 0.9),mapping=aes(ymax = mean_cond+SEM,ymin=mean_cond-SEM)) + scale_fill_manual(values=c("navy","maroon"))+xlab("ATTRACTOR TYPE")+ylab("MEAN RATING")+ylim(0,1500)+ggtitle("JUDGEMENT RT")
dev.off()

#### Hypothesis testing

ezANOVA(data = data.acceptability, dv = RT, wid = Subject, within = c(attractor, attachment))

###       $ANOVA
###       Effect DFn DFd        F          p p<.05         ges
###       2            attractor   1  47 2.383681 0.12931551       0.003595689
###       3           attachment   1  47 1.237401 0.27163004       0.001683808
###       4 attractor:attachment   1  47 3.657623 0.06191823       0.006346966


########### BY ITEM STUFF ##############

### by item judgments

item.by.cond.judge <- data.acceptability %>% 
  group_by(Item, attractor, attachment) %>% 
  summarise(mean.judge = mean(response),
            SEM = sd(response)/sqrt(n_distinct(Subject)))

pdf('agr-amb-item-judge-attach.pdf')
ggplot(item.by.cond.judge,aes(x=Item,y=mean.judge,fill=attachment))+geom_bar(position='dodge',stat = "identity")+ scale_fill_manual(values=c("navy","maroon"))+xlab("ITEM")+ylab("MEAN RATING")+ylim(0,7)+ggtitle("JUDGEMENT BY ITEM")
dev.off()

pdf('agr-amb-item-judge-attract.pdf')
ggplot(item.by.cond.judge,aes(x=Item,y=mean.judge,fill=attractor))+geom_bar(position='dodge',stat = "identity")+ scale_fill_manual(values=c("navy","maroon"))+xlab("ITEM")+ylab("MEAN RATING")+ylim(0,7)+ggtitle("JUDGEMENT BY ITEM")
dev.off()

### by item judgment reaction times

item.by.cond.judge.RT <- data.acceptability %>% 
  group_by(Item, attractor, attachment) %>% 
  summarise(mean.RT = mean(RT))

pdf('agr-amb-item-judgeRT.pdf')
ggplot(item.by.cond.judge.RT,aes(x=Item,y=mean.RT,fill=attachment))+geom_bar(position='dodge',stat = "identity")+ scale_fill_manual(values=c("navy","maroon"))+xlab("ITEM")+ylab("MEAN RATING")+ylim(0,1500)+ggtitle("JUDGEMENT RT BY ITEM")
dev.off()


#########################################
###                                   ###
###          RAW SPR ANALYSIS         ###
###                                   ###
#########################################

#### Ensure RT data is numeric

data.RT$RT <- as.numeric(as.character(data.RT$RT))

pdf('agr-amb-SPR-RT.pdf')
ggplot(data.RT,aes(x=RT))+geom_histogram(binwidth=50)+xlim(0,7000)+ggtitle("SPR RT DISTRIBUTION")
dev.off()

#### Filter our the abnormally long responses

data.RT <- filter(data.RT, RT < 5000)

#### Descriptive stats

RT.subj.by.cond <- data.RT %>%
  group_by(Subject, Experiment, attractor, attachment, region) %>%
  summarise(average = mean(RT))

RT.cond.summ <- RT.subj.by.cond %>%
  group_by(attractor, Experiment, attachment, region) %>%
  summarise(mean = mean(average),
            SEM = sd(average)/sqrt(n_distinct(Subject)))

#### Plotting RT data

pdf('agr-amb-SPR-RTplotRaw.pdf')
ggplot(subset(RT.cond.summ,Experiment %in% c("agr-amb-a","agr-amb-b","agr-amb-c","agr-amb-d")),aes(x=region,y=mean,color=Experiment,base=6,group=Experiment))+ 
  labs(y="Reading time",x="Region",group=1) +geom_point(stat = "identity",size=1)+
  geom_errorbar(aes(ymax = mean+SEM,ymin=mean-SEM,width=0.05))+ 
  theme(text = element_text(size=10))+stat_identity(geom="line")+
  scale_colour_manual(values = c("steelblue2","firebrick4","navyblue","firebrick2"),
                      name="Conditions",
                      labels= c("high, singular","low, singular", "high, plural", "low, plural"),
                      breaks= c("agr-amb-c",'agr-amb-a','agr-amb-b','agr-amb-d'))+
  ggtitle("SPR RAW RT")
dev.off()

#### ANOVA by region
region1 <- droplevels(subset(data.RT, region == "1"))
region2 <- droplevels(subset(data.RT, region == "2"))
region3 <- droplevels(subset(data.RT, region == "3"))
region4 <- droplevels(subset(data.RT, region == "4"))

ezANOVA(data = region1, dv = RT, wid = Subject, within = c(attractor, attachment))

###       $ANOVA
###       Effect DFn DFd         F         p p<.05          ges
###       2            attractor   1  47 0.4355169 0.5125154       0.0006952119
###       3           attachment   1  47 0.6308108 0.4310502       0.0010781077
###       4 attractor:attachment   1  47 3.4796906 0.0683754       0.0042759098

ezANOVA(data = region2, dv = RT, wid = Subject, within = c(attractor, attachment))

###       $ANOVA
###       Effect DFn DFd         F         p p<.05          ges
###       2            attractor   1  47 0.1691656 0.6827240       0.0002059359
###       3           attachment   1  47 0.7979166 0.3762682       0.0005181058
###       4 attractor:attachment   1  47 0.5255524 0.4720767       0.0006412581

ezANOVA(data = region3, dv = RT, wid = Subject, within = c(attractor, attachment))

###       $ANOVA
###       Effect DFn DFd            F         p p<.05          ges
###       2            attractor   1  47 2.5939346706 0.1139701       2.046357e-03
###       3           attachment   1  47 0.0001588828 0.9899964       1.252671e-07
###       4 attractor:attachment   1  47 0.0630814220 0.8027852       3.287930e-05

ezANOVA(data = region4, dv = RT, wid = Subject, within = c(attractor, attachment))

###       $ANOVA
###       Effect DFn DFd         F         p p<.05          ges
###       2            attractor   1  47 0.5927148 0.4452237       0.0007829055
###       3           attachment   1  47 1.6133420 0.2102749       0.0022916065
###       4 attractor:attachment   1  47 0.2903031 0.5925694       0.0003194378

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

### Create the model

resid.model <- lmer(RT ~ char.length + (1|Subject), data.model)

data.model$resid <- resid(resid.model)

### Get rid of outliers

data.model <- filter(data.model, RT < 5000)

### data frames to interpret the model results

data.model.cond <- data.model %>%
  group_by(Subject, Experiment, region) %>%
  summarise(average = mean(resid))

data.model.cond.summ <- data.model.cond %>%
  group_by(Experiment, region) %>%
  summarise(mean = mean(average),
            SEM = sd(average)/sqrt(n_distinct(Subject)))

#### ANOVA PREP

data.model.RT <- droplevels(subset(data.model , Experiment %in% c('agr-amb-a','agr-amb-b','agr-amb-c','agr-amb-d')))
data.model.RT$attachment <- ifelse(data.model.RT$Experiment=='agr-amb-b' | data.model.RT$Experiment=='agr-amb-c','high','low')
data.model.RT$attractor <- ifelse(data.model.RT$Experiment=='agr-amb-a' | data.model.RT$Experiment=='agr-amb-c','singular','plural')

region1.resid <- droplevels(subset(data.model.RT, region == "1"))
region2.resid <- droplevels(subset(data.model.RT, region == "2"))
region3.resid <- droplevels(subset(data.model.RT, region == "3"))
region4.resid <- droplevels(subset(data.model.RT, region == "4"))

### Plotting residualized RT data

pdf('agr-amb-SPR-RTplotResidual.pdf')
ggplot(subset(data.model.cond.summ,Experiment %in% c("agr-amb-a","agr-amb-b","agr-amb-c","agr-amb-d")),aes(x=region,y=mean,color=Experiment,base=6,group=Experiment))+ 
  labs(y="Reading time",x="Region",group=1) +geom_point(stat = "identity",size=1)+
  geom_errorbar(aes(ymax = mean+SEM,ymin=mean-SEM,width=0.05))+ 
  theme(text = element_text(size=10))+stat_identity(geom="line")+
  scale_colour_manual(values = c("steelblue2","firebrick4","navyblue","firebrick2"),
                      name="Conditions",
                      labels= c("high, singular","low, singular", "high, plural", "low, plural"),
                      breaks= c("agr-amb-c",'agr-amb-a','agr-amb-b','agr-amb-d'))+
  ggtitle("SPR RESIDUALIZED RT")
dev.off()

#### Inferential stats

ezANOVA(data = region1.resid, dv = resid, wid = Subject, within = c(attractor, attachment))

###       $ANOVA
###       Effect DFn DFd          F         p p<.05          ges
###       2            attractor   1  47 0.42064215 0.5197733       6.778155e-04
###       3           attachment   1  47 0.58152619 0.4495252       1.087430e-03
###       4 attractor:attachment   1  47 0.06196518 0.8045010       7.597401e-05

ezANOVA(data = region2.resid, dv = resid, wid = Subject, within = c(attractor, attachment))

###       $ANOVA
###       Effect DFn DFd         F         p p<.05          ges
###       2            attractor   1  47 0.2329751 0.6315672       0.0002434144
###       3           attachment   1  47 0.9358488 0.3382996       0.0006242376
###       4 attractor:attachment   1  47 0.1959143 0.6600695       0.0002743053

ezANOVA(data = region3.resid, dv = resid, wid = Subject, within = c(attractor, attachment))

###       $ANOVA
###       Effect DFn DFd            F         p p<.05          ges
###       2            attractor   1  47 0.9557999182 0.3332529       8.128818e-04
###       3           attachment   1  47 0.0001789849 0.9893824       1.493456e-07
###       4 attractor:attachment   1  47 0.1217080151 0.7287474       6.570922e-05

ezANOVA(data = region4.resid, dv = resid, wid = Subject, within = c(attractor, attachment))

###       $ANOVA
###       Effect DFn DFd         F         p p<.05          ges
###       2            attractor   1  47 0.6153849 0.4367036       0.0008373330
###       3           attachment   1  47 1.7279633 0.1950500       0.0023893650
###       4 attractor:attachment   1  47 0.2642216 0.6096430       0.0003367803
