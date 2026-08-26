# Script for Testing Spatial Frequencies Task - AL
# v0.1 : AL 12.06.2024


##########################
#### Load data frames ####
##########################

# clear work space
rm(list = ls())

# set working directory
setwd('path/')

# loading / installing Packages
requiredPackages <- c('plyr','dplyr','tidyr','sjPlot' ,'ggplot2','r2glmm'
                      ,'lme4', 'lmerTest' ,'multcomp', 'broom', 'RLRsim', 'esquisse',
                      'grid', 'gridExtra', 'scales','car','data.table','emmeans', 'tidyverse','ggrain',
                      'scales','car','sm','forestplot','cowplot') #psycho
idx <- which(!requiredPackages %in% installed.packages())
if (length(idx)>0) {
  install.packages(requiredPackages[idx])
}
for (pkg in requiredPackages) {
  require(pkg, character.only = TRUE)
}

library('ggsignif')

# load files
files <- list.files('online_frequencies/') # for exp 1
#files <- list.files('online_frequencies_2/') # for exp 2

csv_files <- grep ('.csv',files)
files <- files[csv_files]

# initialize an empty data.table to store all df
df.all <- data.table()

# loop through each file name in my file list

  for (filename in paste0('online_frequencies/',files)) { # exp 1
  #for (filename in paste0('online_frequencies_2/',files)) { # exp 2
  
  df <- fread(filename, na.strings=c('','NA'))
  
  columns_to_extract <- c('participant','condition','target','responseMain.keys','responseMain.rt',
                          'responseSimple.keys', 'responseSimple.rt',
                          'slider_panas.response','slider_panas.rt','words','valence',
                          'slider_mathys1.response','slider_mathys1.rt','item_type','direction',
                          'slider_mathys2.response','slider_mathys2.rt','words_mathys',
                          'mathys3.text',"left_side","right_side")#,'image_1.tStart','image_2.tStart','mask.tStart')
  
  # convert keys if vb_positive is 1
  if (df$vb_positive[1] == 1) {
    df$responseMain.keys <- ifelse(df$responseMain.keys == 'v', 0, 1)
    df$responseSimple.keys <- ifelse(df$responseSimple.keys == 'v', 0, 1)
    
  } else {
    df$responseMain.keys <- ifelse(df$responseMain.keys == 'b', 0, 1)
    df$responseSimple.keys <- ifelse(df$responseSimple.keys == 'b', 0, 1)
  }
  
  # convert to numeric
  df$responseMain.keys <- as.numeric(df$responseMain.keys)
  df$responseSimple.keys <- as.numeric(df$responseSimple.keys)
  
  # select desired columns
  selected_data <- df[, ..columns_to_extract]
  
  # combine selected data with df
  df.all <- rbind(df.all, selected_data)
}

# read demographics
demographics <- read.csv('demographics_freq/demo_freq_v1.csv') # exp 1
#demographics <- read.csv('demographics_freq/demo_freq_v2.csv') # exp 2

# # check timings
# df.all$image1_pres <- df.all$image_2.tStart - df.all$image_1.tStart
# df.all$image2_pres <- df.all$mask.tStart - df.all$image_2.tStart


########################
### Data definition ####
########################

# rename condition in the response main to be able to combine everything
# get target data
# block 1
df.simple <- df.all[!is.na(df.all$responseSimple.keys),]
df.simple$Block <- 1
# block 2
df.double <- df.all[!is.na(df.all$responseMain.keys),]
df.double$Block <- 2

#specific to the double
df.double$condition <- plyr::revalue(df.double$condition, c("LSF"="CtF","HSF"="FtC"))

#same column for both responses
df.simple$Response <- df.simple$responseSimple.keys
df.simple$RT <- df.simple$responseSimple.rt
df.double$Response <- df.double$responseMain.keys
df.double$RT <- df.double$responseMain.rt

# combine them
df.target <- merge(df.simple,df.double, all=T)

# rename columns
colnames(df.target)[colnames(df.target) == 'target'] <- 'Target'
colnames(df.target)[colnames(df.target) == 'condition'] <- 'Condition'
colnames(df.target)[colnames(df.target) == 'participant'] <- 'ID'

# rename values in columns
df.target$Target <- plyr::revalue(df.target$Target, c("sad"="Sad","happy"="Happy","morph"="Morph"))


# define factors
df.target$Target <- factor(df.target$Target)
df.target$Condition <- factor(df.target$Condition)
df.target$ID <- factor(df.target$ID)

# convert RT in ms
df.target$RT <- df.target$RT*1000

#add a correct column for happy and sad targets
df.target$correct.happy <- ifelse(df.target$Target=='Happy' & df.target$Response == 1, 1, 0)
df.target$correct.sad <- ifelse(df.target$Target=='Sad' & df.target$Response == 0, 1, 0)



# check aberrant RT per targets per participant
df.hs <- df.target[!df.target$Target=='Morph',]
df.hs.out <- aggregate(RT ~ID, data=df.hs, mean)
df.hs.out$Z.RT <- scale(df.hs.out$RT)
hs.out <- df.hs.out[df.hs.out$Z.RT > 2 | df.hs.out$Z.RT < (-2),]

out.ppt <- hs.out

# delete outliers from the df
df.target <- df.target[!df.target$ID %in% out.ppt$ID,]
#df.target %>% count(ID) # check how many ppt I have -> 28

### outliers rt per trial per target
####################################

# get unique id per line to get rid of trials later
df.target <- df.target %>%
  mutate(row_extract = row_number())

# exclude trials with RT > 3s (meaning technical problem)
df.target <- df.target[!df.target$RT > 3000,]

# plot RT distribution
ggplot(df.target) + aes(x=RT) + geom_histogram() + facet_grid(~Target)+
  labs(x='Reaction time',
       y='Count')



# trim trials above 3 SD from the mean in block 1 and 2 
detect_rt_outliers <- function(df) {
  df %>%
    group_by(Target, Block, Condition) %>%
    mutate(Z.RT = scale(RT)) %>%
    filter(abs(Z.RT) > 3) %>%
    ungroup()
}

trials.out <- detect_rt_outliers(df.target)

# exclude aberrant trials from the analysis
df.target <- df.target[!df.target$row_extract %in% trials.out$row_extract,]

# replot histogram
ggplot(df.target) + aes(x=RT) + geom_histogram() + facet_wrap(~Block+Target)+
  labs(x='Reaction time',
       y='Count')



#########################################
#### Sort questionnaire data for PCA ####
#########################################

# ## make the df for the PCA
# uncomment to create the df
# # #panas
# df.panas <- df.all[!is.na(df.all$slider_panas.response),]
# extract_panas <- c('slider_panas.response','participant','words')
# df.panas <- df.panas[,..extract_panas]
# panas_pca <- pivot_wider(df.panas, names_from = words, values_from = slider_panas.response)
# panas_pca <- panas_pca[!(panas_pca$participant %in% out.ppt$ID),]
# 
# # #mathys
# df.mathys1 <- df.all[!is.na(df.all$slider_mathys1.response),]
# extract_mathys1 <- c('participant','slider_mathys1.response','item_type','direction','left_side','right_side')
# df.mathys1 <- df.mathys1[,..extract_mathys1]
# df.mathys1$slider_mathys1.response[df.mathys1$direction == "inverse"] <- 10 - df.mathys1$slider_mathys1.response[df.mathys1$direction == "inverse"]
# df.mathys1 <- df.mathys1 %>% mutate(item = ifelse(direction == "normal", right_side, left_side))
# extract_mathys1 <- c('participant','slider_mathys1.response','item')
# df.mathys1 <- df.mathys1[,..extract_mathys1]
# mathys1_pca <- pivot_wider(df.mathys1, names_from = item, values_from = slider_mathys1.response)
# mathys1_pca <- mathys1_pca[!(mathys1_pca$participant %in% out.ppt$ID),]
# 
# df.mathys2 <- df.all[!is.na(df.all$slider_mathys2.response),]
# extract_mathys2 <- c('participant','slider_mathys2.response','words_mathys')
# df.mathys2 <- df.mathys2[,..extract_mathys2]
# mathys2_pca <- pivot_wider(df.mathys2, names_from = words_mathys, values_from = slider_mathys2.response)
# mathys2_pca <- mathys2_pca[!(mathys2_pca$participant %in% out.ppt$ID),]
# 
# # # merge them all
# pca_all <- merge(panas_pca, mathys1_pca, by='participant')
# pca_all <- merge(pca_all, mathys2_pca, by='participant')
# 
# colnames(pca_all)[colnames(pca_all) == 'participant'] <- 'participant_id'
# 
# write.table(pca_all, 'pca/quest_freq_v2.csv')

pca <- read.table('pca/results/pca_score_freq_v1.csv') # exp 1
#pca <- read.table('pca/results/pca_score_freq_v2.csv') # exp 2
colnames(pca)[colnames(pca) == 'participant_id'] <- 'ID'

# demographics
demographics$Age <- as.numeric(demographics$Age)
demographics$Ethnicity.simplified <- factor(demographics$Ethnicity.simplified)
demographics$Time.taken <- demographics$Time.taken / 60
demographics$ID <- demographics$Participant.id
demographics$ID <- as.factor(demographics$ID)
demographics <- demographics[demographics$Status == 'APPROVED',]

summary_stats_demographics <- demographics %>%
  summarise(
    n = n(),
    n_female = sum(Sex == "Female", na.rm = TRUE),
    pct_female = n_female / n * 100,
    mean.Age = mean(Age, na.rm = T),
    sd.Age = sd(Age),
    mean.Time = mean(Time.taken, na.rm = T),
    sd.Time = sd(Time.taken))

demographics <- demographics[!demographics$ID %in% out.ppt$ID,]
df.target <- merge(df.target, demographics, by='ID')

################################################################################
################################################################################
##################################
#### Graphical representation ####
##################################
################################################################################
################################################################################

fig_outdir <- 'path/Figures/'


df.target$Condition <- factor(df.target$Condition, levels = c("HSF","LSF","FtC","CtF"))
df.target$Block <- factor(df.target$Block,
                          levels = c(1, 2),
                          labels = c("Block 1", "Block 2"))


# aggregate all scores per ID, Target, Block and Condition
agg.response <- df.target %>%
  group_by(ID,Target,Block,Condition) %>%
  summarise(
    mean = mean(Response, na.rm = T),
    sd =  sd(Response, na.rm = T),
    se = sd(Response) / sqrt(n())
  )

mean.resp <- df.target %>%
  group_by(Target,Block,Condition) %>%
  summarise(
    mean = mean(Response, na.rm = T),
    sd =  sd(Response, na.rm = T),
    se = sd(Response) / sqrt(n())
  )


agg.rt <- df.target %>%
  group_by(ID,Target,Block,Condition) %>%
  summarise(
    mean = mean(RT, na.rm = T),
    se = sd(RT) / sqrt(n())
  )


# scale age
df.target$Age <- scale(df.target$Age, scale = FALSE)

# only keep correct happy and sad RT
df.happy <- df.target[df.target$Target=='Happy',]
df.sad <- df.target[df.target$Target=='Sad',]
df.happy.corr <- df.happy[df.happy$correct.happy==1,]
df.sad.corr <- df.sad[df.sad$correct.sad==1,]
df.hs <- rbind(df.happy.corr, df.sad.corr, fill = TRUE)

# morphs based on positive or negative responses
df.morph <- df.target[df.target$Target=='Morph',]
df.morph$Target <- ifelse(df.morph$Response==0, "Morph -", "Morph +")

# merge
df.rt <- merge(df.hs, df.morph, all=T)

df.rt$Block <- as.factor(df.rt$Block)
# not normally distributed, try to log transform
df.rt$logRT <- log(df.rt$RT)


agg.rt <- df.rt %>%
  group_by(ID,Target,Block,Condition) %>%
  summarise(
    mean = mean(RT, na.rm = T),
    se = sd(RT) / sqrt(n())
  )



# faire des plots par block

# block 1
b1Palette <- c('#00416A','#86bbd8')

# cat
agg.response.1 <- agg.response[agg.response$Block=="Block 1",]

cat.b1 <- 
  ggplot(agg.response.1) + 
  aes(y = mean, x = Target, fill = Condition, color = Condition) +
  geom_boxplot(alpha = 0.7, width = 0.7, outlier.shape = NA, 
               position = position_dodge(width = 0.85)) +
  geom_jitter(aes(color = Condition),
              alpha = 0.6,
              size = 1.4,
              position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.6))+
  labs(x='Target Emotion', y="Emotional Rating", fill="Filtering", color='Filtering')+
  stat_summary(fun = mean, shape = 23, size = 0.8, stroke = 0.8, 
               color = "black", alpha = 0.8, position = position_dodge(width = 0.85)) +  
  theme_blank() +  
  scale_color_manual(values = b1Palette) +
  scale_fill_manual(values = b1Palette) +
  ggtitle("Emotional Rating to Targets by Filtering")+
  theme(plot.title = element_text(size = 20),
        axis.title.x = element_text(size = 20),   # x‐axis label
        axis.title.y = element_text(size = 20),   # y‐axis label
        axis.text.x  = element_text(size = 18, color="grey30"),   # x‐axis tick labels
        axis.text.y  = element_text(size = 26, color="grey30"),
        strip.text.x = element_text(size = 18),
        legend.title = element_text(size = 20),   # "Condition"
        legend.text  = element_text(size = 18)) +
  
  #scale_y_continuous(breaks = c(0, 1), labels = c("Negative", "Positive")) +
  # scale_y_continuous(breaks = c(0, 1), labels = c("-", "+")) +
  # scale_y_continuous(breaks = c(0.25, 0.75), labels = c("↓", "↑")) +
  
  # scale_y_continuous(
  #   breaks = c(0, 0.25, 0.5, 0.75, 1),
  #   labels = c("- ↓", "", "", "", "+ ↑")
  # ) +
  
  scale_y_continuous(
    breaks = c(0, 0.15, 0.35, 0.5, 0.65, 0.85,1),
    labels = c("0", "-", "↓", "", "↑", "+","1")
  ) +
  
  labs(tag = "A") +
  theme(
    plot.tag.position = c(0.02, 0.98),  # x, y in [0,1], (0,1) = top-left
    plot.tag = element_text(size = 21, face = "bold"))+
  geom_hline(yintercept = 0.5,
             linetype   = "dashed",
             color      = "grey40",
             size       = 0.45)+
  # # exp 1
  # # Happy
  # annotate("segment", x = 0.8, xend = 1.2, y = 1.1, yend = 1.1, color = "black", size = 0.6) +
  # annotate("text", x = 1, y = 1.12, label = "***", size = 7) +
  # # Sad
  # annotate("segment", x = 2.8, xend = 3.2, y = 1.05, yend = 1.05, color = "black", size = 0.6) +
  # annotate("text", x = 3, y = 1.07, label = "***", size = 7)
  #   
  #exp 2
  # Happy
  annotate("segment", x = 0.8, xend = 1.2, y = 1.1, yend = 1.1, color = "black", size = 0.6) +
  annotate("text", x = 1, y = 1.12, label = "**", size = 7) +
  # Morphs
  annotate("segment", x = 1.8, xend = 2.2, y = 1.09, yend = 1.09, color = "black", size = 0.6) +
  annotate("text", x = 2, y = 1.11, label = "***", size = 7)




ggsave(file.path(fig_outdir,"Resp_B1_1st_Exp.png"), 
       plot = cat.b1, width = 8, height = 6, units = "in")

ggsave(file.path(fig_outdir,"Resp_B1_2nd_Exp.png"), 
       plot = cat.b1, width = 8, height = 6, units = "in")




model_block1 <- glmer(Response ~ Target * Condition + Age + Sex + (1 | ID),
                      data = df.target %>% filter(Block == "Block 1"),
                      family = binomial(link='logit'))
summary(model_block1)
Anova(model_block1, type=3)

# main effects
print(emmeans(model_block1, pairwise ~ Target))
print(emmeans(model_block1, pairwise ~ Condition))

# planned comparisons
print(emmeans(model_block1, pairwise ~ Condition | Target))

# check if morphs are rated as more positive or negative overall
df_morph_block_1 <- df.target %>% filter(Block == "Block 1", Target == "Morph")

# count % of positive resp
k <- sum(df_morph_block_1$Response == 1)   # si codage binaire 1=positif, 0=négatif
n <- nrow(df_morph_block_1)

# binomial test against 0.5
binom.test(k, n, p = 0.5)


# RT
agg.rt$Target <- factor(agg.rt$Target, levels = c("Happy", "Morph +","Sad","Morph -"))

agg.rt.1 <- agg.rt[agg.rt$Block=="Block 1",]



rt.b1 <- 
  ggplot(agg.rt.1) + aes(y=mean, x= Target, fill=Condition, color=Condition)+
  geom_violin(alpha=0.65, position = position_dodge(width = 1.1))+
  labs(x = "Target Emotion", y = "log(RT)", 
       title = 'RT to Targets by Filtering', fill='Filtering',color='Filtering')+
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, 
               color = "black", fatten = 2, position = position_dodge(width = 1.1)) +  
  geom_jitter(aes(color = Condition),
              alpha = 1,
              size = 1.6,
              shape = 21,
              position = position_jitterdodge(jitter.width = 0.7, dodge.width = 1.1))+
  theme_blank()+
  theme(
    plot.title = element_text(size = 21, margin = margin(b=35)),
    axis.title.x = element_text(size = 20),   # x‐axis label
    axis.title.y = element_text(size = 20),   # y‐axis label
    axis.text.x  = element_text(size = 18,color='grey30'),   # x‐axis tick labels
    axis.text.y  = element_text(size = 18,color='grey30'),
    strip.text.x = element_text(size = 18),
    legend.title = element_text(size = 20),   # "Condition"
    legend.text  = element_text(size = 18)
  ) +
  scale_color_manual(values=b1Palette)+
  scale_fill_manual(values = b1Palette) +
  
  labs(tag = "B") +
  theme(
    plot.tag.position = c(0.02, 0.98),  # x, y in [0,1], (0,1) = top-left
    plot.tag = element_text(size = 21, face = "bold"))+ 
  
  # exp 1   
  # Happy
  annotate("segment", x = 0.8, xend = 1.2, y = 1240, yend = 1240, color = "black", size = 0.6) +
  annotate("text", x = 1, y = 1260, label = "***", size = 7) +
  
  # Sad
  annotate("segment", x = 2.8, xend = 3.2, y = 1330, yend = 1330, color = "black", size = 0.6) +
  annotate("text", x = 3, y = 1350, label = "***", size = 7)


# exp 2
# # Happy 
# annotate("segment", x = 0.8, xend = 1.2, y = 920, yend = 920, color = "black", size = 0.6) +
# annotate("text", x = 1, y = 940, label = "***", size = 7) +
# 
# # Morph Negative
# annotate("segment", x = 3.8, xend = 4.2, y = 1110, yend = 1110, color = "black", size = 0.6) +
# annotate("text", x = 4, y = 1130, label = "*", size = 7)



ggsave(file.path(fig_outdir,"RT_B1_1st_Exp.png"), 
       plot = rt.b1, width = 8, height = 6, units = "in")

ggsave(file.path(fig_outdir,"RT_B1_2nd_Exp.png"), 
       plot = rt.b1, width = 8, height = 6, units = "in")



# stats

model_rt_block1 <- lmer(logRT ~ Target * Condition + Age + Sex + (1 | ID), data = df.rt %>% 
                          filter(Block == "Block 1"))#, !(ID %in% agg.rt.O$ID)))
summary(model_rt_block1)
Anova(model_rt_block1, type = 3)

# main effects
print(emmeans(model_rt_block1, pairwise ~ Target))
print(emmeans(model_rt_block1, pairwise ~ Condition))

# planned comparisons
print(emmeans(model_rt_block1, pairwise ~ Condition | Target))

aggregate(RT ~ Target, data = df.rt %>% filter(Block == "Block 1"), mean)
aggregate(RT ~ Condition, data = df.rt %>% filter(Block == "Block 1"), mean)


# descriptive statistics
desc_stat_RTs <- df.rt %>%
  group_by(Target, Condition) %>%
  summarise(
    meanRT = mean(RT, na.rm = TRUE),
    sdRT   = sd(RT, na.rm = TRUE))



# block 2 
# faire des plots par block
b2Palette <- c('#7A5B8F', "#80B192")

# cat
agg.response.2 <- agg.response[agg.response$Block=="Block 2",]

cat.b2 <- 
  ggplot(agg.response.2) + 
  aes(y = mean, x = Target, fill = Condition, color = Condition) +
  geom_boxplot(alpha = 0.7, width = 0.7, outlier.shape = NA, 
               position = position_dodge(width = 0.85)) +
  geom_jitter(aes(color = Condition),
              alpha = 0.6,
              size = 1.4,
              position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.6))+
  labs(x='Target Emotion', y="Response",fill="Sequence", color='Sequence')+
  stat_summary(fun = mean, shape = 23, size = 0.8, stroke = 0.8, 
               color = "black", alpha = 0.8, position = position_dodge(width = 0.85)) +  
  theme_blank() +  
  scale_color_manual(values = b2Palette) +
  scale_fill_manual(values = b2Palette) +
  ggtitle("Response to Targets by Sequence")+
  theme(plot.title = element_text(size = 21),
        axis.title.x = element_text(size = 20),   # x‐axis label
        axis.title.y = element_text(size = 20),   # y‐axis label
        axis.text.x  = element_text(size = 18, color='grey30'),   # x‐axis tick labels
        axis.text.y  = element_text(size = 26, color="grey30"),
        strip.text.x = element_text(size = 18),
        legend.title = element_text(size = 20),   # "Condition"
        legend.text  = element_text(size = 18)) +
  
  labs(tag = "A") +
  theme(
    plot.tag.position = c(0.02, 0.98),  # x, y in [0,1], (0,1) = top-left
    plot.tag = element_text(size = 21, face = "bold"))+
  
  scale_y_continuous(
    breaks = c(0, 0.15, 0.35, 0.5, 0.65, 0.85,1),
    labels = c("0", "-", "↓", "", "↑", "+","1")
  ) +
  
  geom_hline(yintercept = 0.5,
             linetype   = "dashed",
             color      = "grey40",
             size       = 0.45)+
  #exp 2:
  # Happy
  annotate("segment", x = 0.8, xend = 1.2, y = 1.1, yend = 1.1, color = "black", size = 0.6) +
  annotate("text", x = 1, y = 1.12, label = "*", size = 7) +
  # Morph +
  annotate("segment", x = 1.8, xend = 2.2, y = 1.05, yend = 1.05, color = "black", size = 0.6) +
  annotate("text", x = 2, y = 1.07, label = "*", size = 7)





ggsave(file.path(fig_outdir,"Resp_B2_1st_Exp.png"), 
       plot = cat.b2, width = 8, height = 6, units = "in")

ggsave(file.path(fig_outdir,"Resp_B2_2nd_Exp.png"), 
       plot = cat.b2, width = 8, height = 6, units = "in")


# stats
model_block2 <- glmer(Response ~ Target * Condition  + Age + Sex +(1 | ID),
                      data = df.target %>% filter(Block == "Block 2"),
                      family = binomial(link='logit'))

Anova(model_block2, type=3)

# main effects
print(emmeans(model_block2, pairwise ~ Target))
print(emmeans(model_block2, pairwise ~ Condition))

# planned comparisons
print(emmeans(model_block2, pairwise ~ Condition | Target))


# rt

agg.rt.2 <- agg.rt[agg.rt$Block=="Block 2",]

rt.b2 <- 
  ggplot(agg.rt.2) + aes(y=mean, x= Target, fill=Condition, color=Condition)+
  geom_violin(alpha=0.65, position = position_dodge(width = 1.1))+
  labs(x = "Target Emotion", y = "RT (ms)", title = 'RT to Targets by Sequence',
       color='Sequence',fill='Sequence')+
  #facet_wrap( ~Block, strip.position = 'bottom')+
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, 
               color = "black", fatten = 2, position = position_dodge(width = 1.1)) +  
  geom_jitter(aes(color = Condition),
              alpha = 1,
              size = 1.6,
              shape = 21,
              #stroke = 0.1,
              #color='black',
              position = position_jitterdodge(jitter.width = 0.7, dodge.width = 1.1))+
  theme_blank()+
  theme(
    plot.title = element_text(size = 21, margin = margin(b=35)),
    axis.title.x = element_text(size = 20),   # x‐axis label
    axis.title.y = element_text(size = 20),   # y‐axis label
    axis.text.x  = element_text(size = 18,color='grey30'),   # x‐axis tick labels
    axis.text.y  = element_text(size = 18,color='grey30'),
    strip.text.x = element_text(size = 18),
    legend.title = element_text(size = 20),   # "Condition"
    legend.text  = element_text(size = 18)
  ) +
  scale_color_manual(values=b2Palette)+
  scale_fill_manual(values = b2Palette) +
  
  labs(tag = "B") +
  theme(
    plot.tag.position = c(0.02, 0.98),  # x, y in [0,1], (0,1) = top-left
    plot.tag = element_text(size = 21, face = "bold")) +
  
  
  # exp 1
  # Happy
  annotate("segment", x = 0.8, xend = 1.2, y = 1000, yend = 1000, color = "black", size = 0.6) +
  annotate("text", x = 1, y = 1020, label = "***", size = 7) +
  
  # Sad
  annotate("segment", x = 2.8, xend = 3.2, y = 940, yend = 940, color = "black", size = 0.6) +
  annotate("text", x = 3, y = 960, label = "***", size = 7) +
  
  # Morph Positive
  annotate("segment", x = 1.8, xend = 2.2, y = 1240, yend = 1240, color = "black", size = 0.6) +
  annotate("text", x = 2, y = 1260, label = "***", size = 7)



# # exp 2
# # Happy 
# annotate("segment", x = 0.8, xend = 1.2, y = 940, yend = 940, color = "black", size = 0.6) +
# annotate("text", x = 1, y = 960, label = "*", size = 7) +
# 
# # Morph Positive 
# annotate("segment", x = 1.8, xend = 2.2, y = 1240, yend = 1240, color = "black", size = 0.6) +
# annotate("text", x = 2, y = 1260, label = "***", size = 7)
# 



ggsave(file.path(fig_outdir,"RT_B2_1st_Exp.png"), 
       plot = rt.b2, width = 8, height = 6, units = "in")

ggsave(file.path(fig_outdir,"RT_B2_2nd_Exp.png"), 
       plot = rt.b2, width = 8, height = 6, units = "in")


# block 2
# model_rt_block2 <- lmer(RT ~ Target * Condition + Age + Sex + (1 | ID),
#                         data = df.rt %>% filter(Block == "Block 2"))
# summary(model_rt_block2)
# Anova(model_rt_block2, type = 3)
# print(emmeans(model_rt_block2, pairwise ~ Condition | Target))
# 
# res <- residuals(model_rt_block2)
# qqnorm(res)
# qqline(res, col = "red")
# shapiro.test(res)

model_rt_block2 <- lmer(logRT ~ Target * Condition + Age + Sex + (1 | ID), data = df.rt %>% 
                          filter(Block == "Block 2"))
summary(model_rt_block2)
Anova(model_rt_block2, type = 3)

# main effects
print(emmeans(model_rt_block2, pairwise ~ Target))
print(emmeans(model_rt_block2, pairwise ~ Condition))


# planned comparisons
print(emmeans(model_rt_block2, pairwise ~ Condition | Target))

# in exp 1, effect of age -> exploratory analyses???

# exp 2 also effect of age! explore?


#####___________________________________________________________________________

# just mood and morphs
agg.morph <- aggregate(df.morph, Response ~ ID + Block + Condition + Age + Sex, mean)
agg.morph <- merge(agg.morph,pca,by='ID')


ggplot(agg.morph) +
  aes(x = Response, y = PC1) +
  geom_point(colour = "#112446") +
  theme_minimal() +
  facet_wrap(~Block*Condition)+
  geom_smooth(method='lm')



mood.df <- left_join(df.morph, pca, by = "ID")
mood.df$ID <- as.factor(mood.df$ID)
mood.df$Block <- as.factor(mood.df$Block)

ggplot(mood.df, aes(x = PC1, y = Response)) + 
  geom_point() + stat_smooth(method = "glm",method.args = list(family = binomial(link='logit')), se = TRUE, color='#2f4858') + 
  xlab("Mood (PC1)") + ylab("Probability of Positive Answer") + 
  facet_wrap(~Block)+
  ggtitle("Probability of positive answer to morphs by Mood") +
  theme_blank()

model_block1 <- glmer(Response ~ PC1 + Age + Sex + (1 | ID), 
                      data = mood.df %>% filter(Block == "Block 1"),
                      family = binomial(link='logit'))
summary(model_block1)

# sum <- mood.df %>%
#   group_by(Block,Response) %>%
#   summarise(RT.mean = mean(RT, na.rm = T),
#             RT.sd = sd(RT, na.rm = T))

model_block2 <- glmer(Response ~ PC1 + Age + Sex + (1 | ID), 
                      data = mood.df %>% filter(Block == "Block 2"),
                      family = binomial(link='logit'))
summary(model_block2)


# facilitation of RT to + or - depending on mood
mood.df$RTlog <- log(mood.df$RT)

model_rt_block1 <- lmer(RTlog ~ PC1 * Response + Age + Sex + (1 | ID),
                        data = mood.df %>% filter(Block == "Block 1"))
summary(model_rt_block1)


# sex and responses
sum <- mood.df %>%
  group_by(Block,Response) %>% # or sex
  summarise(RT.mean = mean(RT, na.rm = T),
            RT.sd = sd(RT, na.rm = T))

# age
ggplot(data = mood.df %>% filter(Block == "Block 1"), aes(x = Age, y = RT)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Effet de l'âge sur les temps de réaction",
       x = "Âge",
       y = "Temps de réaction (logRT)") +
  theme_minimal()



model_rt_block2 <- lmer(RTlog ~ PC1 * Response + Age + Sex + (1 | ID),
                        data = mood.df %>% filter(Block == "Block 2"))
summary(model_rt_block2)
#Anova(model_rt_block2, type = 3)
#print(emmeans(model_rt_block1, pairwise ~ Condition | Target))


ggplot(mood.df %>% filter(Block == "Block 1"),
       aes(x = PC1, y = RTlog, color = as.factor(Response))) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "Mood score (PC1)",
    y = "Log RT",
    color = "Response",
    title = "Mood (PC1) × Response interaction on log RT"
  ) +
  theme_minimal()


ggplot(mood.df %>% filter(Block == "Block 2"),
       aes(x = Age, y = RT, color = as.factor(Response))) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "Mood score (PC1)",
    y = "Log RT",
    color = "Response",
    title = "Mood (PC1) × Response interaction on log RT"
  ) +
  theme_minimal()




# plot pca results

# define palette
pal <- c("#f26419",'#f6ae2d','#592720')

# PC1

# sort from higher to lower
pca <- pca %>%
  mutate(Feature = reorder(feature, score)) %>%
  arrange(desc(score))

# threshold - only show items with feature contribution below -0.15 and above 0.15
pca <- pca[pca$score > 0.15 | pca$score < -0.15,]

# plot
pc2_exp2 <- ggplot(pca) +
  aes(x = score, y = feature, fill = Scale) +
  geom_bar(stat = "summary", alpha=0.7, width = 0.7) +
  labs(x='PC2 scores', y="Items")+
  scale_fill_manual(values=pal)+
  theme_blank()+
  ggtitle("Items loadings on PC2") +
  theme(plot.title = element_text(size = 41, face='bold'),
        axis.title.x = element_text(size = 39,color='black'),   # x‐axis label
        axis.title.y = element_text(size = 39,color='black'),   # y‐axis label
        axis.text.x  = element_text(size = 36,color='black'),   # x‐axis tick labels
        axis.text.y  = element_text(size = 28, color='#393E46'),
        legend.text = element_text(size = 38,color='black'),
        legend.title = element_text(size = 39,color='black')) 


ggsave(file.path(fig_outdir,"PC1_exp1.png"), 
       plot = pc1_exp1, width = 20, height = 15, units = "in")

ggsave(file.path(fig_outdir,"PC1_exp2.png"), 
       plot = pc1_exp2, width = 20, height = 15, units = "in")


ggsave(file.path(fig_outdir,"sup_mat/PC2_exp1.png"), 
       plot = pc2_exp1, width = 20, height = 15, units = "in")

ggsave(file.path(fig_outdir,"sup_mat/PC2_exp2.png"), 
       plot = pc2_exp2, width = 20, height = 15, units = "in")
