# load packages
library(tidyverse)
library(haven)
library(kableExtra)
library(stargazer)
library(TAM)
library(spatstat)
library(SmartEDA)
library(stats)

#----- READING DATA --------
# read in 5% level adjustment factors from Lee et al. (2021)
dat05 <- readRDS("./adj_05.rda")
# read in 1% level adjustment factors from Lee et al. (2021)
dat01 <- readRDS("./adj_01.rda")

# read in my assembled dataset
papers <- read.csv("./papers_dataset.csv") %>% mutate(
  first_stage_f = as.numeric(gsub(",", "", first_stage_f)),
  iv_sd = as.numeric(gsub(",", "", iv_sd)),
  iv_beta = as.numeric(gsub(",", "", iv_beta)),
  journal = ifelse(grepl("macro", journal, ignore.case = T), "Macroeconomics",
                   ifelse(grepl("micro", journal, ignore.case = T), "Microeconomics",
                          ifelse(grepl("applied", journal, ignore.case = T), "Applied Economics", "Economic Policy")))
) %>% filter(journal != "")

# filter specifications with my required fields
papers_dat <- filter(papers, !is.na(first_stage_f)  &  !is.na(iv_sd)) %>%
  select(journal, year, title, strategy, first_stage_f, iv_beta, iv_sd, n_authors, n_female, top_uni, topic) %>%
  mutate(female_author_share = n_female / n_authors,
         has_female_author = n_female >= 1,
         topic_aggregate = NA) %>%
  mutate(n_female_authors = n_female, 
         n_authors_top_uni = top_uni,
         .keep = "unused") 

aggregates <- c("Health", "Labor", "Politics", "Tax", "Micro", "Education", "Law", "Macro")
for (row in 1:nrow(papers_dat)) {
  topics <- str_split(papers_dat[row, 'topic'], "; ", simplify = T)
  for (topic in topics) {
    if (topic %in% aggregates) {
      papers_dat[row,]['topic_aggregate'] <- topic
      break
    }
  }
}
papers_dat <- subset(papers_dat, select = -topic)
#----- ASSIGNING ADJUSTMENT --------
# Adjustment tfse
adj_factor_index <- function(fstat, adj_dat) {
  return(findInterval(fstat, adj_dat$f))
}

# return adjusted SE
tf_se <- function(fstat, se, adj_dat) {
  if (fstat <= min(adj_dat$f)) {
    adj_factor <- max(adj_dat$adj)
  } else if (fstat >= max(adj_dat$f)) {
    adj_factor <- min(adj_dat$adj)
  } else {
    index <- adj_factor_index(fstat, adj_dat)
    if (adj_dat$f[index] == fstat) {
      adj_factor <- adj_dat$adj[index]
    } else {
      adj_factor <- adj_dat$adj[index+1] + 
        (adj_dat$f[index+1] - fstat)/(adj_dat$f[index+1] - adj_dat$f[index]) *
        (adj_dat$adj[index] - adj_dat$adj[index+1])
    }
  }
  return(se * adj_factor)
}

stopifnot(abs(tf_se(10, 1, dat05) - 1.751) < 10e-3)

#--------- ADJUSTED SE & CONFIDENCE INTERVALS AT 5% and 1%-------
# for every spec, get adjusted SE at 0.05 level
papers_dat$tfse_05 <- NA
papers_dat$tfse_01 <- NA
for (row in 1:nrow(papers_dat)) {
  tf_05 <- tf_se(papers_dat[row,]['first_stage_f'], papers_dat[row,]['iv_sd'], dat05)
  tf_01 <- tf_se(papers_dat[row,]['first_stage_f'], papers_dat[row,]['iv_sd'], dat01)
  papers_dat[row,]['tfse_05'] <- tf_05
  papers_dat[row,]['tfse_01'] <- tf_01
}

z05 <- 1.96
z01 <- 2.576

papers_dat <- papers_dat %>% mutate(
  CI05_includes_0 = iv_beta - z05*iv_sd <= 0 & iv_beta + z05*iv_sd >= 0,
  CI01_includes_0 = iv_beta - z01*iv_sd <= 0 & iv_beta + z01*iv_sd >= 0,
  tfCI05_includes_0 = iv_beta - z05*tfse_05 <= 0 & iv_beta + z05*tfse_05 >= 0,
  tfCI01_includes_0 = iv_beta - z01*tfse_01 <= 0 & iv_beta + z01*tfse_01 >= 0
)

######-------- BUILDING CONFIDENCE INTERVALS--------
papers_dat <- papers_dat %>% 
  mutate(ci05_low = iv_beta - z05*iv_sd,
         ci05_up = iv_beta + z05*iv_sd,
         ci01_low = iv_beta - z01*iv_sd,
         ci01_up = iv_beta + z01*iv_sd,
         ci05_low_adj = iv_beta - z05*tfse_05,
         ci05_up_adj = iv_beta + z05*tfse_05,
         ci01_low_adj = iv_beta - z01*tfse_01,
         ci01_up_adj = iv_beta + z01*tfse_01,
         f_over_10 = ifelse(first_stage_f >= 10, 1, 0))


######-------- COUNTING HOW MANY PAPERS CHANGE ---------
sum(papers_dat$tfCI05_includes_0 == TRUE & papers_dat$CI05_includes_0 == FALSE)
sum(papers_dat$CI05_includes_0 == FALSE)

# get new significant stats at 0.01 level
sum(papers_dat$tfCI01_includes_0 == TRUE & papers_dat$CI01_includes_0 == FALSE)
sum(papers_dat$CI01_includes_0 == FALSE)

# CHECK: Out of 437 specifications, 262 are significant (at the 0.05 level),
# 65 of which ($65/262 \approx 25\%$) are no longer significantly different from zero
#at 0.05 level after SE adjustment. Similarly, 187 are significant (at the 0.01 level),
#56 of which ($56/187 \approx 30\%$) are no longer significantly different from zero
# at the 0.01 level after adjustment.

######--------category for significance



######---------- GRAPHS----------
# weights
papers_titles <- papers_dat %>% 
  group_by(title) %>% 
  summarise(specifications = n(),
            share = n()/nrow(papers_dat)) %>% 
  mutate(weight = 1/specifications)

papers_dat <- left_join(papers_dat, papers_titles, by = "title")
papers_dat$proportion <- papers_dat$weight/sum(papers_dat$weight)

# adjustment factor
papers_dat <- papers_dat %>% 
  mutate(adjustment_inverse = 1/(tfse_05/iv_sd),
         adjustment = (tfse_05/iv_sd),
         adjustment01 = (tfse_01/iv_sd),
         adjustment_inverse01 = 1/(tfse_01/iv_sd))
summary(papers_dat$adjustment)

# GRAPH 1: Distribution of f-stats
summary(papers_dat$first_stage_f, weight = weight) # 25%: 14.88; 50%: 31.19; 75%: 241.77

require(scales)
a <- papers_dat %>% 
  ggplot(mapping = aes(x = first_stage_f, y = stat(count / sum(count)), weight = weight)) +
  geom_histogram(bins = 40, color = "lightgrey", fill = "darkgrey") +
  scale_x_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)),
                     ) +
  xlab("First-stage F-statistic") + ylab("Proportion") +
  geom_vline(aes(xintercept = 14.88), colour="black", linetype = "dashed") +
  geom_vline(aes(xintercept = 31.19), colour="black", linetype = "dashed") +
  geom_vline(aes(xintercept = 241.77), colour="black", linetype = "dashed") + 
  geom_vline(aes(xintercept = 142.6), colour="black") +
  annotate(geom = "text",
           label = c("25th", 142.6, "50th", "75th"),
           x = c(10, 90, 45, 360),
           y = c(0.12,0.14,0.12,0.12),
           angle = 0,
           vjust = 1,
           size=2.5) +
  theme_classic()
a

# Graph 2: distribution of inverse adjustment factors
summary(papers_dat$adjustment_inverse, weight = weight) # 25%: 0.6829; 50%: 0.84; 75%: 1
require(scales)
b <- papers_dat %>% 
  ggplot(mapping = aes(x = adjustment_inverse, y = stat(count / sum(count)), weight = weight)) +
  geom_histogram(bins = 40, color = "lightgrey", fill = "darkgrey") +
  xlab("1.96/srqt[c_.05(F)]") + ylab("Proportion") +
  geom_vline(aes(xintercept = 0.6986143), colour="black", linetype = "dashed") +
  geom_vline(aes(xintercept = 0.9489362), colour="black", linetype = "dashed") +
  geom_vline(aes(xintercept = 1), colour="black", linetype = "dashed") +
  scale_y_continuous(limits = c(0,0.5)) +
  annotate(geom = "text",
           label = c("25th", "50th", "75th"),
           x = c(0.66, 0.92, 0.98),
           y = c(0.2,0.2,0.25),
           angle = 0,
           vjust = 1,
           size=2.5) +
  theme_classic()
b

# Graph 2.2: distribution of inverse adjustment factors 1%
summary(papers_dat$adjustment_inverse, weight = weight) # 25%: 0.4596993; 50%: 0.7887594  ; 75%: 0.9418363 
require(scales)
b.2 <- papers_dat %>% 
  ggplot(mapping = aes(x = adjustment_inverse01, y = stat(count / sum(count)), weight = weight)) +
  geom_histogram(bins = 40, color = "lightgrey", fill = "darkgrey") +
  xlab("2.58/srqt[c_.01(F)]") + ylab("Proportion") +
  geom_vline(aes(xintercept = 0.4596993), colour="black", linetype = "dashed") +
  geom_vline(aes(xintercept = 0.7887594), colour="black", linetype = "dashed") +
  geom_vline(aes(xintercept = 0.9418363), colour="black", linetype = "dashed") +
  scale_y_continuous(limits = c(0,0.5)) +
  annotate(geom = "text",
           label = c("25th", "50th", "75th"),
           x = c(0.43, 0.80, 0.92),
           y = c(0.2,0.2,0.2),
           angle = 0,
           vjust = 1,
           size=2.5) +
  theme_classic()
b.2

# Graph 3: figure 6 in Lee et al (2021):
papers_figure6 <- papers_dat %>% 
  select(c(5,6,7,14,15,16,17,18,19,28, 31, 34, 35)) %>% 
  mutate(f_normalized = (first_stage_f/10)/(1+first_stage_f/10),
         t_normalized = ((iv_beta^2/iv_sd^2))/(1+(iv_beta^2/iv_sd^2)),
         )


### smoothed tf-function
dat05 <-dat05 %>%
  mutate(minimum_t = 1.96*dat05$adj,
         f_normalized = (f/10)/(1+f/10),
         t_normalized = (minimum_t^2/1.96^2)/(1+(minimum_t^2/1.96^2)))

dat01 <-dat01 %>%
  mutate(minimum_t = 2.73*dat05$adj,
         f_normalized = (f/10)/(1+f/10),
         t_normalized = (minimum_t^2/1.96^2)/(1+(minimum_t^2/1.96^2)))

## SETTING COLORS W SIGNIFICANCE
papers_figure6$significance = 0
papers_figure6$significance[papers_figure6$f_over_10==1 &
                              abs(papers_figure6$iv_beta/papers_figure6$iv_sd) > 1.96 &
                              papers_figure6$CI05_includes_0 == 0 &
                              papers_figure6$tfCI05_includes_0 == 1] <- 1
papers_figure6$significance[papers_figure6$f_over_10==1 &
                              abs(papers_figure6$iv_beta/papers_figure6$iv_sd) > 1.96 &
                              papers_figure6$CI05_includes_0 == 0 &
                              papers_figure6$tfCI05_includes_0 == 0 &
                              papers_figure6$tfCI01_includes_0 == 1] <- 2
papers_figure6$significance[papers_figure6$f_over_10==1 &
                              abs(papers_figure6$iv_beta/papers_figure6$iv_sd) > 1.96 &
                              papers_figure6$CI05_includes_0 == 0 &
                              papers_figure6$tfCI05_includes_0 == 0 &
                              papers_figure6$tfCI01_includes_0 == 0] <- 3

summary(papers_figure6$significance[papers_figure6$f_over_10==1 | abs(papers_figure6$iv_beta/papers_figure6$iv_sd) > 1.96])
summary(papers_figure6$significance[papers_figure6$f_normalized==1/2 & abs(papers_figure6$iv_beta/papers_figure6$iv_sd) > 1.96])


papers_figure6_filtered <- filter(papers_figure6, !is.na(tfCI05_includes_0))
papers_figure6_filtered2 <- filter(papers_figure6, !is.na(significance))

c <- ggplot() +
  geom_point(data = papers_figure6, aes(x = f_normalized, y = t_normalized, color = as.factor(significance),
                 size = weight), shape = 1, stroke = 0.8) + 
  scale_color_manual(values = c("black", "blue", "#8B008B", "darkred")) +
  geom_smooth(data = dat05, aes(x = f_normalized, y = t_normalized), color="black") +
  geom_segment(aes(x=0.91281,xend=1,y=0.505, yend=0.505), color= "black", size =0.85) +
  geom_smooth(data = dat01, aes(x = f_normalized, y = t_normalized), color="darkgrey") +
  geom_segment(aes(x=0.96188,xend=1,y=0.65986, yend=0.664), color= "darkgrey", size =0.8) +
  geom_vline(aes(xintercept = 0.27754), colour="grey", linetype = "dashed") + 
  geom_vline(aes(xintercept = 0.39963), colour="grey", linetype = "dashed") + 
  geom_vline(aes(xintercept = 0.5), colour="grey", linetype = "dashed") + 
  geom_vline(aes(xintercept = 0.91281), colour="grey", linetype = "dashed") +
  geom_hline(aes(yintercept = 0.5), colour="grey", linetype = "dashed") + 
  geom_hline(aes(yintercept = 0.65986), colour="grey", linetype = "dashed") + 
  ylab(bquote(t^2~"statistic")) +
  xlab("F statistic") + 
  scale_x_continuous(
    limits = c(-0.05, 1.05),
    expand = c(0, 0), # The horizontal axis does not extend to either side
    breaks = c(0, 0.27754, 0.39963, 1/2, 0.91281,1),  # Set custom break locations
    labels = c("0", "1.96^2", "2.58^2", "10","104.6", expression(infinity)) # And custom labels on those breaks!
  )  +
  scale_y_continuous(
    limits = c(-0.05, 1.05),
    breaks = c(0, 1/2, 0.65986, 1), 
    labels = c("0", "1.96^2", "2.73^2", expression(infinity)),
    expand = c(0, 0)
  )  + theme_classic() + theme(legend.position ="none")
c

### smoothed tf-function
dat05 <-dat05 %>%
  mutate(minimum_t = 1.96*dat05$adj,
         f_normalized = (f/10)/(1+f/10),
         t_normalized = (minimum_t^2/1.96^2)/(1+(minimum_t^2/1.96^2)))
d <- dat05 %>% 
  ggplot() + 
  geom_smooth(aes(x = f_normalized, y = t_normalized), color="black") + 
  scale_x_continuous(limits = c(0,1)) + 
  scale_y_continuous(limits = c(0,1)) + 
  geom_segment(aes(x=0.91281,xend=1,y=0.505, yend=0.505), color= "black", size =0.7)
d


## d-stat table weighted
####------------ LENGTH OF CONFIDENCE INTERVALS
papers_dat <- papers_dat %>% 
  mutate(CI05_length = ci05_up-ci05_low,
         CI01_length = ci01_up-ci01_low,
         tfCI05_length = ci05_up_adj-ci05_low_adj,
         tfCI01_length = ci01_up_adj-ci01_low_adj
         )

table_manual <- 
  round(data.frame(min = c(min(papers_dat$tfse_05), min(papers_dat$tfse_01), min(papers_dat$adjustment), min(papers_dat$CI05_includes_0),
                     min(papers_dat$CI01_includes_0),min(papers_dat$tfCI05_includes_0),min(papers_dat$tfCI01_includes_0),
                     min(papers_dat$CI05_includes_0),min(papers_dat$CI05_length),min(papers_dat$CI01_length),
                                                              min(papers_dat$tfCI05_length),min(papers_dat$tfCI01_length)),
             perc25 = c(weighted_quantile(papers_dat$tfse_05, papers_dat$weight, probs = 0.25), weighted_quantile(papers_dat$tfse_01, papers_dat$weight, probs = 0.25),
                        weighted_quantile(papers_dat$adjustment, papers_dat$weight, probs = 0.25), weighted_quantile(papers_dat$CI05_includes_0, papers_dat$weight, probs = 0.25),
                        weighted_quantile(papers_dat$CI01_includes_0, papers_dat$weight, probs = 0.25),weighted_quantile(papers_dat$tfCI05_includes_0, papers_dat$weight, probs = 0.25),
                        weighted_quantile(papers_dat$tfCI01_includes_0, papers_dat$weight, probs = 0.25),
                        weighted_quantile(papers_dat$CI05_includes_0, papers_dat$weight, probs = 0.25),weighted_quantile(papers_dat$CI05_length, papers_dat$weight, probs = 0.25),
                        weighted_quantile(papers_dat$CI01_length, papers_dat$weight, probs = 0.25),
                        weighted_quantile(papers_dat$tfCI05_length, papers_dat$weight, probs = 0.25),weighted_quantile(papers_dat$tfCI01_length, papers_dat$weight, probs = 0.25)),
             perc50 = c(weighted_quantile(papers_dat$tfse_05, papers_dat$weight, probs = 0.5), weighted_quantile(papers_dat$tfse_01, papers_dat$weight, probs = 0.5),
                        weighted_quantile(papers_dat$adjustment, papers_dat$weight, probs = 0.5), weighted_quantile(papers_dat$CI05_includes_0, papers_dat$weight, probs = 0.5),
                        weighted_quantile(papers_dat$CI01_includes_0, papers_dat$weight, probs = 0.5),weighted_quantile(papers_dat$tfCI05_includes_0, papers_dat$weight, probs = 0.5),
                        weighted_quantile(papers_dat$tfCI01_includes_0, papers_dat$weight, probs = 0.5),
                        weighted_quantile(papers_dat$CI05_includes_0, papers_dat$weight, probs = 0.5),weighted_quantile(papers_dat$CI05_length, papers_dat$weight, probs = 0.5),
                        weighted_quantile(papers_dat$CI01_length, papers_dat$weight, probs = 0.5),
                        weighted_quantile(papers_dat$tfCI05_length, papers_dat$weight, probs = 0.5),weighted_quantile(papers_dat$tfCI01_length, papers_dat$weight, probs = 0.5)),
             mean = c(weighted_mean(papers_dat$tfse_05, papers_dat$weight), weighted_mean(papers_dat$tfse_01, papers_dat$weight),
                        weighted_mean(papers_dat$adjustment, papers_dat$weight), weighted_mean(papers_dat$CI05_includes_0, papers_dat$weight),
                        weighted_mean(papers_dat$CI01_includes_0, papers_dat$weight),weighted_mean(papers_dat$tfCI05_includes_0, papers_dat$weight),
                        weighted_mean(papers_dat$tfCI01_includes_0, papers_dat$weight),
                        weighted_mean(papers_dat$CI05_includes_0, papers_dat$weight),weighted_mean(papers_dat$CI05_length, papers_dat$weight),
                        weighted_mean(papers_dat$CI01_length, papers_dat$weight),
                        weighted_mean(papers_dat$tfCI05_length, papers_dat$weight),weighted_mean(papers_dat$tfCI01_length, papers_dat$weight)),
             perc75 = c(weighted_quantile(papers_dat$tfse_05, papers_dat$weight, probs = 0.75), weighted_quantile(papers_dat$tfse_01, papers_dat$weight, probs = 0.75),
                        weighted_quantile(papers_dat$adjustment, papers_dat$weight, probs = 0.75), weighted_quantile(papers_dat$CI05_includes_0, papers_dat$weight, probs = 0.75),
                        weighted_quantile(papers_dat$CI01_includes_0, papers_dat$weight, probs = 0.75),weighted_quantile(papers_dat$tfCI05_includes_0, papers_dat$weight, probs = 0.75),
                        weighted_quantile(papers_dat$tfCI01_includes_0, papers_dat$weight, probs = 0.75),
                        weighted_quantile(papers_dat$CI05_includes_0, papers_dat$weight, probs = 0.75),weighted_quantile(papers_dat$CI05_length, papers_dat$weight, probs = 0.75),
                        weighted_quantile(papers_dat$CI01_length, papers_dat$weight, probs = 0.75),
                        weighted_quantile(papers_dat$tfCI05_length, papers_dat$weight, probs = 0.75),weighted_quantile(papers_dat$tfCI01_length, papers_dat$weight, probs = 0.75)),
             max = c(max(papers_dat$tfse_05), max(papers_dat$tfse_01), max(papers_dat$adjustment), max(papers_dat$CI05_includes_0),
               max(papers_dat$CI01_includes_0),max(papers_dat$tfCI05_includes_0),max(papers_dat$tfCI01_includes_0),
               max(papers_dat$CI05_includes_0),max(papers_dat$CI05_length),max(papers_dat$CI01_length),
               max(papers_dat$tfCI05_length),max(papers_dat$tfCI01_length)), 
        sd = c(weighted_sd(papers_dat$tfse_05, papers_dat$weight), weighted_sd(papers_dat$tfse_01, papers_dat$weight),
               weighted_sd(papers_dat$adjustment, papers_dat$weight), weighted_sd(papers_dat$CI05_includes_0, papers_dat$weight),
               weighted_sd(papers_dat$CI01_includes_0, papers_dat$weight),weighted_sd(papers_dat$tfCI05_includes_0, papers_dat$weight),
               weighted_sd(papers_dat$tfCI01_includes_0, papers_dat$weight),
               weighted_sd(papers_dat$CI05_includes_0, papers_dat$weight),weighted_sd(papers_dat$CI05_length, papers_dat$weight),
               weighted_sd(papers_dat$CI01_length, papers_dat$weight),
               weighted_sd(papers_dat$tfCI05_length, papers_dat$weight),weighted_sd(papers_dat$tfCI01_length, papers_dat$weight))),
        digits = 2)
           
sd <- t(as.data.frame(table_manual$sd))

# d-stats table CATEGORICAL
options(width = 150)
dstats_category <- as.data.frame(ExpCustomStat(papers_dat,Cvar = c("journal","year","strategy","topic_aggregate"), stat = c('Count','Prop'),gpby=FALSE))
cat1 <- as.data.frame(dstats_category$Prop)

# Table 2 in lee et al (2021):
papers_dat <- papers_dat %>% 
  mutate(tab1 = ifelse(first_stage_f <=10 & (iv_beta/iv_sd)^2 >= 1.96^2, 1, 0),
         tab2 = ifelse(first_stage_f > 10 & (iv_beta/iv_sd)^2 >= 1.96^2, 1, 0),
         tab3 = ifelse(first_stage_f <=10 & (iv_beta/iv_sd)^2 < 1.96^2, 1, 0),
         tab4 = ifelse(first_stage_f > 10 & (iv_beta/iv_sd)^2 < 1.96^2, 1, 0))
sum(papers_dat$tab1==1)
sum(papers_dat$tab2==1)
sum(papers_dat$tab3==1)
sum(papers_dat$tab4==1)

weighted.mean(papers_dat$tab1==1, papers_dat$weight)
weighted.mean(papers_dat$tab2==1, papers_dat$weight)
weighted.mean(papers_dat$tab3==1, papers_dat$weight)
weighted.mean(papers_dat$tab4==1, papers_dat$weight)


### NUMERICAL VARIABLES D-STATS
table <- round(data.frame(min = c(min(papers_dat$n_authors),
                                  min(papers_dat$n_female),
                                  min(papers_dat$top_uni),
                                  min(papers_dat$first_stage_f),
                                  min(papers_dat$iv_beta),
                                  min(papers_dat$iv_sd),
                                  min(papers_dat$specifications)),
                          quantile1 = c(weighted_quantile(papers_dat$n_authors, w=papers_dat$weight, probs=.25),
                                        weighted_quantile(papers_dat$n_female, w=papers_dat$weight, probs=.25),
                                        weighted_quantile(papers_dat$top_uni, w=papers_dat$weight, probs=.25),
                                        weighted_quantile(papers_dat$first_stage_f, w=papers_dat$weight, probs=.25),
                                        weighted_quantile(papers_dat$iv_beta, w=papers_dat$weight, probs=.25),
                                        weighted_quantile(papers_dat$iv_sd, papers_dat$weight, probs=.25),
                                        weighted_quantile(papers_dat$specifications, w=papers_dat$weight, probs=.25)),
                          
                          mean = c(weighted.mean(papers_dat$n_authors, papers_dat$weight),
                                   weighted.mean(papers_dat$n_female, papers_dat$weight),
                                   weighted.mean(papers_dat$top_uni, papers_dat$weight),
                                   weighted.mean(papers_dat$first_stage_f, papers_dat$weight),
                                   weighted.mean(papers_dat$iv_beta, papers_dat$weight),
                                   weighted.mean(papers_dat$iv_sd, papers_dat$weight),
                                   weighted.mean(papers_dat$specifications, papers_dat$weight)),
                          median = c(weighted.quantile(papers_dat$n_authors, w = papers_figure6$weight, prob = 0.5),
                                     weighted.quantile(papers_dat$n_female, w = papers_figure6$weight, prob = 0.5),
                                     weighted.quantile(papers_dat$top_uni, w = papers_figure6$weight, prob = 0.5),
                                     weighted.quantile(papers_dat$first_stage_f, w = papers_figure6$weight, prob = 0.5),
                                     weighted.quantile(papers_dat$iv_beta, w = papers_figure6$weight, prob = 0.5),
                                     weighted.quantile(papers_dat$iv_sd, w = papers_figure6$weight, prob = 0.5),
                                     weighted.quantile(papers_dat$specifications, w = papers_figure6$weight, prob = 0.5)),
                          
                          
                          quantile3 = c(weighted.quantile(papers_dat$n_authors, w = papers_figure6$weight, prob = 0.75),
                                        weighted.quantile(papers_dat$n_female, w = papers_figure6$weight, prob = 0.75),
                                        weighted.quantile(papers_dat$top_uni, w = papers_figure6$weight, prob = 0.75),
                                        weighted.quantile(papers_dat$first_stage_f, w = papers_figure6$weight, prob = 0.75),
                                        weighted.quantile(papers_dat$iv_beta, w = papers_figure6$weight, prob = 0.75),
                                        weighted.quantile(papers_dat$iv_sd, w = papers_figure6$weight, prob = 0.75),
                                        weighted.quantile(papers_dat$specifications, w = papers_figure6$weight, prob = 0.75)),
                          
                          max = c(max(papers_dat$n_authors),
                                  max(papers_dat$n_female),
                                  max(papers_dat$top_uni),
                                  max(papers_dat$first_stage_f),
                                  max(papers_dat$iv_beta),
                                  max(papers_dat$iv_sd),
                                  max(papers_dat$specifications)),
                          
                          sd = c(weighted_sd(papers_dat$n_authors, w=papers_dat$weight),
                                 weighted_sd(papers_dat$n_female, w=papers_dat$weight),
                                 weighted_sd(papers_dat$top_uni, w=papers_dat$weight),
                                 weighted_sd(papers_dat$first_stage_f, w=papers_dat$weight),
                                 weighted_sd(papers_dat$iv_beta, w=papers_dat$weight),
                                 weighted_sd(papers_dat$iv_sd, w=papers_dat$weight),
                                 weighted_sd(papers_dat$specifications, w=papers_dat$weight))), digits = 2)


papers_dat$n_authors_top_uni <- ifelse(papers_dat$n_authors_top_uni == 0, 0, 1)

## categorical table weighted
table_cat1 <- papers_dat %>% 
  group_by(title, strategy) %>% 
  summarise(n = n())

table_cat1.1 <- table_cat1 %>% 
  group_by(strategy) %>% 
  summarise(n = n(), share = n()/37)
sum(table_cat1.1$share)

#### REGRESSIONS:
papers_dat <- papers_dat %>% 
  mutate(trend = year-2010,
         t_squared = (iv_beta/iv_sd)^2,
         adj = tfse_05/iv_sd,
         adj1 = tfse_01/iv_sd)
# without year dummies
logit1 <- glm(CI05_includes_0~factor(strategy)+first_stage_f+n_authors+has_female_author+n_authors_top_uni+
                factor(topic_aggregate)+specifications, papers_dat,
              weights = weight, family = binomial(link = "logit"))

logit2 <- glm(CI05_includes_0~factor(strategy)+first_stage_f+n_authors+has_female_author+n_authors_top_uni+
                factor(topic_aggregate)+specifications + factor(year), papers_dat,
              weights = weight, family = binomial(link = "logit"))

logit3 <- glm(CI05_includes_0~factor(strategy)+first_stage_f+n_authors+has_female_author+n_authors_top_uni+
                factor(topic_aggregate)+specifications + trend, papers_dat,
              weights = weight, family = binomial(link = "logit"))

logit4 <- glm(CI05_includes_0~factor(strategy)+first_stage_f+n_authors+has_female_author+n_authors_top_uni+
                factor(topic_aggregate)+specifications + factor(year) + trend, papers_dat,
              weights = weight, family = binomial(link = "logit"))

stargazer(logit1, logit2, logit3, logit4)


## MAIN REGRESSIONS
reg1 <- lm(adj~t_squared+first_stage_f+factor(strategy)+n_authors+has_female_author+n_authors_top_uni+
             factor(topic_aggregate)+specifications, papers_dat,
           weights = weight)
reg2 <- lm(adj~t_squared+first_stage_f+factor(strategy)+n_authors+has_female_author+n_authors_top_uni+
             factor(topic_aggregate)+specifications+trend, papers_dat,
           weights = weight)
reg3 <- lm(adj~t_squared+first_stage_f+factor(strategy)+n_authors+has_female_author+n_authors_top_uni+
             factor(topic_aggregate)+specifications+factor(year), papers_dat,
           weights = weight)
reg4 <- lm(adj~t_squared+first_stage_f+factor(strategy)+n_authors+has_female_author+n_authors_top_uni+
             factor(topic_aggregate)+specifications+trend+factor(year), papers_dat,
           weights = weight)

stargazer(reg1, reg2, reg3, reg4,
          type = 'latex',
          
          header=FALSE, # to get rid of r package output text
          
          single.row = TRUE, # to put coefficients and standard errors on same line
          
          no.space = TRUE, # to remove the spaces after each line of coefficients
          
          column.sep.width = "1pt", # to reduce column width
          
          font.size = "small", digits = 2, df = FALSE,
          dep.var.labels = "\textit{tF}-Adjustment Factor at the 5% Significance Level")

predict1 <- predict(reg1)
summary(as.data.frame(predict1))

reg5 <- lm(adj1~t_squared+first_stage_f+factor(strategy)+n_authors+has_female_author+n_authors_top_uni+
             factor(topic_aggregate)+specifications, papers_dat,
           weights = weight)
reg6 <- lm(adj1~t_squared+first_stage_f+factor(strategy)+n_authors+has_female_author+n_authors_top_uni+
             factor(topic_aggregate)+specifications+trend, papers_dat,
           weights = weight)
reg7 <- lm(adj1~t_squared+first_stage_f+factor(strategy)+n_authors+has_female_author+n_authors_top_uni+
             factor(topic_aggregate)+specifications+factor(year), papers_dat,
           weights = weight)
reg8 <- lm(adj1~t_squared+first_stage_f+factor(strategy)+n_authors+has_female_author+n_authors_top_uni+
             factor(topic_aggregate)+specifications+trend+factor(year), papers_dat,
           weights = weight)

stargazer(reg5, reg6, reg7, reg8,
          type = 'latex',
          
          header=FALSE, # to get rid of r package output text
          
          single.row = TRUE, # to put coefficients and standard errors on same line
          
          no.space = TRUE, # to remove the spaces after each line of coefficients
          
          column.sep.width = "1pt", # to reduce column width
          
          font.size = "small", digits = 2, df = FALSE,
          dep.var.labels = "\textit{tF}-Adjustment Factor at the 1% Significance Level")
           