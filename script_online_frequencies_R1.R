# =============================================================================
# Emotional Face Perception — Spatial Frequencies Task (online data)
# Behavioural analysis script — Experiments 1 & 2
#
# Reproduces the emotional-rating and reaction-time analyses and the figures
# reported in the manuscript.
#
# CONTENTS
#   0. Setup .................  packages + user-configurable paths
#   1. Load raw data .........  read PsychoPy csv exports, assemble one table
#   2. Preprocessing .........  block/response coding, factors, RT outlier removal
#   3. Questionnaires / PCA ..  run-once block to build the PCA input (kept commented)
#   4. Demographics ..........  summary + merge into the trial table
#   5. Per-experiment export .  write exp1.csv / exp2.csv, then pool both
#   6. Aggregations ..........  subject-level means used by the figures
#   7. Block 1 — Response ....  figure + models (single experiment / pooled)
#   8. Block 1 — RT ..........  figure + models (single experiment / pooled)
#   9. Block 2 — Response ....  figure + models (single experiment / pooled)
#  10. Block 2 — RT ..........  figure + models (single experiment / pooled)
#  11. Mood (PCA) analyses ...  morph responses & RT as a function of PC1
#  12. PCA loading figures ...  item contributions on PC1 / PC2
#
# WORKFLOW
#   Sections 1–5 are run ONCE PER EXPERIMENT. Set `experiment` (section 0) to
#   1 or 2, run down to the export step to create exp1.csv / exp2.csv, then the
#   pooled models read both files back in for the between-experiment comparison.
#
#   Each figure carries manually placed significance brackets that differ
#   between experiments; the Exp 1 and Exp 2 variants are kept as clearly
#   labelled blocks — uncomment the one matching the experiment being plotted.
#
# All local paths live in section 0 — edit them to match your setup.
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------

# Clear the workspace
rm(list = ls())

# ---- Paths (edit to match your local setup) --------------------------------
# The working directory is assumed to be the project root; everything below is
# relative to it. Uncomment and set if you prefer an explicit path.
# setwd("path/to/project")

# Which experiment to preprocess in sections 1–5: 1 or 2
experiment <- 2

# Experiment-specific inputs (folder / file names differ between Exp 1 and Exp 2)
raw_dir       <- if (experiment == 1) "online_frequencies/"              else "online_frequencies_2/"
demo_file     <- if (experiment == 1) "demographics_freq/demo_freq_v1.csv" else "demographics_freq/demo_freq_v2.csv"
pca_file      <- if (experiment == 1) "pca/results/pca_score_freq_v1.csv"  else "pca/results/pca_score_freq_v2.csv"
pca_quest_out <- if (experiment == 1) "pca/quest_freq_v1.csv"             else "pca/quest_freq_v2.csv"
pc2_file      <- if (experiment == 1) "pca/results/freq_v1PC2.csv"        else "pca/results/freq_v2PC2.csv"

# Outputs
csv_outdir <- "output/"    # intermediate per-experiment tables (exp1.csv / exp2.csv)
fig_outdir <- "figures/"   # figures

# ---- Packages --------------------------------------------------------------
requiredPackages <- c('plyr', 'dplyr', 'tidyr', 'sjPlot', 'ggplot2', 'r2glmm',
                      'lme4', 'lmerTest', 'multcomp', 'broom', 'RLRsim', 'esquisse',
                      'grid', 'gridExtra', 'scales', 'car', 'data.table', 'emmeans',
                      'tidyverse', 'ggrain', 'sm', 'forestplot', 'cowplot')

idx <- which(!requiredPackages %in% installed.packages())
if (length(idx) > 0) {
  install.packages(requiredPackages[idx])
}
for (pkg in requiredPackages) {
  require(pkg, character.only = TRUE)
}
library('ggsignif')

# If dplyr verbs get masked by plyr/data.table, force the dplyr versions:
# conflicts_prefer(dplyr::mutate)
# conflicts_prefer(dplyr::summarise)


# -----------------------------------------------------------------------------
# 1. Load raw data
# -----------------------------------------------------------------------------

# List the csv exports for the current experiment
files     <- list.files(raw_dir)
csv_files <- grep('.csv', files)
files     <- files[csv_files]

# Columns kept from each PsychoPy export
columns_to_extract <- c('participant', 'condition', 'target',
                        'responseMain.keys', 'responseMain.rt',
                        'responseSimple.keys', 'responseSimple.rt',
                        'slider_panas.response', 'slider_panas.rt', 'words', 'valence',
                        'slider_mathys1.response', 'slider_mathys1.rt', 'item_type', 'direction',
                        'slider_mathys2.response', 'slider_mathys2.rt', 'words_mathys',
                        'mathys3.text', 'left_side', 'right_side')

# Empty table to accumulate all participants
df.all <- data.table()

# Read every file, recode responses, keep the columns of interest, stack them
for (filename in file.path(raw_dir, files)) {

  df <- fread(filename, na.strings = c('', 'NA'))

  # Response keys are counterbalanced across participants (vb_positive flag).
  # Recode to a common scheme: 0 = negative response, 1 = positive response.
  if (df$vb_positive[1] == 1) {
    df$responseMain.keys   <- ifelse(df$responseMain.keys   == 'v', 0, 1)
    df$responseSimple.keys <- ifelse(df$responseSimple.keys == 'v', 0, 1)
  } else {
    df$responseMain.keys   <- ifelse(df$responseMain.keys   == 'b', 0, 1)
    df$responseSimple.keys <- ifelse(df$responseSimple.keys == 'b', 0, 1)
  }

  df$responseMain.keys   <- as.numeric(df$responseMain.keys)
  df$responseSimple.keys <- as.numeric(df$responseSimple.keys)

  selected_data <- df[, ..columns_to_extract]
  df.all <- rbind(df.all, selected_data)
}

# Demographics (Prolific export)
demographics <- read.csv(demo_file)


# -----------------------------------------------------------------------------
# 2. Preprocessing / data definition
# -----------------------------------------------------------------------------

# The two blocks are stored in different response columns. Split them, tag the
# block, then stack back together into a single long table.

# Block 1: single-image trials (responseSimple)
df.simple       <- df.all[!is.na(df.all$responseSimple.keys), ]
df.simple$Block <- 1

# Block 2: sequence trials (responseMain)
df.double       <- df.all[!is.na(df.all$responseMain.keys), ]
df.double$Block <- 2

# In Block 2, the filtering condition encodes a sequence direction:
#   LSF -> Coarse-to-Fine (CtF), HSF -> Fine-to-Coarse (FtC)
df.double$condition <- plyr::revalue(df.double$condition, c("LSF" = "CtF", "HSF" = "FtC"))

# Harmonise response / RT column names across blocks
df.simple$Response <- df.simple$responseSimple.keys
df.simple$RT       <- df.simple$responseSimple.rt
df.double$Response <- df.double$responseMain.keys
df.double$RT       <- df.double$responseMain.rt

# Recombine both blocks
df.target <- merge(df.simple, df.double, all = TRUE)

# Rename key columns
colnames(df.target)[colnames(df.target) == 'target']      <- 'Target'
colnames(df.target)[colnames(df.target) == 'condition']   <- 'Condition'
colnames(df.target)[colnames(df.target) == 'participant'] <- 'ID'

# Tidy target labels
df.target$Target <- plyr::revalue(df.target$Target, c("sad" = "Sad", "happy" = "Happy", "morph" = "Morph"))

# Factors
df.target$Target    <- factor(df.target$Target)
df.target$Condition <- factor(df.target$Condition)
df.target$ID        <- factor(df.target$ID)

# RT to milliseconds
df.target$RT <- df.target$RT * 1000

# Accuracy flags for the unambiguous (Happy / Sad) targets
df.target$correct.happy <- ifelse(df.target$Target == 'Happy' & df.target$Response == 1, 1, 0)
df.target$correct.sad   <- ifelse(df.target$Target == 'Sad'   & df.target$Response == 0, 1, 0)


# ---- Participant-level RT outliers -----------------------------------------
# Flag participants whose mean RT (Happy/Sad trials) is > 2 SD from the group mean.
df.hs        <- df.target[!df.target$Target == 'Morph', ]
df.hs.out    <- aggregate(RT ~ ID, data = df.hs, mean)
df.hs.out$Z.RT <- scale(df.hs.out$RT)
hs.out       <- df.hs.out[df.hs.out$Z.RT > 2 | df.hs.out$Z.RT < (-2), ]
out.ppt      <- hs.out

# Drop flagged participants
df.target <- df.target[!df.target$ID %in% out.ppt$ID, ]
# df.target %>% count(ID)   # sanity check on remaining N


# ---- Trial-level RT outliers -----------------------------------------------
# Unique row id so specific trials can be removed later
df.target <- df.target %>% mutate(row_extract = row_number())

# Trials slower than 3 s are treated as technical problems and removed
df.target <- df.target[!df.target$RT > 3000, ]

# RT distribution before trimming
ggplot(df.target) + aes(x = RT) + geom_histogram() + facet_grid(~ Target) +
  labs(x = 'Reaction time', y = 'Count')

# Remove trials > 3 SD from the cell mean (Target x Block x Condition)
detect_rt_outliers <- function(df) {
  df %>%
    group_by(Target, Block, Condition) %>%
    mutate(Z.RT = scale(RT)) %>%
    filter(abs(Z.RT) > 3) %>%
    ungroup()
}
trials.out <- detect_rt_outliers(df.target)
df.target  <- df.target[!df.target$row_extract %in% trials.out$row_extract, ]

# RT distribution after trimming
ggplot(df.target) + aes(x = RT) + geom_histogram() + facet_wrap(~ Block + Target) +
  labs(x = 'Reaction time', y = 'Count')


# -----------------------------------------------------------------------------
# 3. Questionnaire data -> PCA input   (run once to build the file)
# -----------------------------------------------------------------------------
# The block below reshapes the PANAS and Mathys questionnaire responses into a
# wide, participant-by-item matrix and writes it out. It only needs to be run
# once per experiment to produce the PCA input; the PCA itself is computed
# elsewhere and read back below. Kept commented on purpose.
#
# # PANAS
# df.panas <- df.all[!is.na(df.all$slider_panas.response), ]
# extract_panas <- c('slider_panas.response', 'participant', 'words')
# df.panas <- df.panas[, ..extract_panas]
# panas_pca <- pivot_wider(df.panas, names_from = words, values_from = slider_panas.response)
# panas_pca <- panas_pca[!(panas_pca$participant %in% out.ppt$ID), ]
#
# # Mathys (part 1) — reverse-scored items flipped, item label resolved
# df.mathys1 <- df.all[!is.na(df.all$slider_mathys1.response), ]
# extract_mathys1 <- c('participant', 'slider_mathys1.response', 'item_type', 'direction', 'left_side', 'right_side')
# df.mathys1 <- df.mathys1[, ..extract_mathys1]
# df.mathys1$slider_mathys1.response[df.mathys1$direction == "inverse"] <-
#   10 - df.mathys1$slider_mathys1.response[df.mathys1$direction == "inverse"]
# df.mathys1 <- df.mathys1 %>% mutate(item = ifelse(direction == "normal", right_side, left_side))
# extract_mathys1 <- c('participant', 'slider_mathys1.response', 'item')
# df.mathys1 <- df.mathys1[, ..extract_mathys1]
# mathys1_pca <- pivot_wider(df.mathys1, names_from = item, values_from = slider_mathys1.response)
# mathys1_pca <- mathys1_pca[!(mathys1_pca$participant %in% out.ppt$ID), ]
#
# # Mathys (part 2)
# df.mathys2 <- df.all[!is.na(df.all$slider_mathys2.response), ]
# extract_mathys2 <- c('participant', 'slider_mathys2.response', 'words_mathys')
# df.mathys2 <- df.mathys2[, ..extract_mathys2]
# mathys2_pca <- pivot_wider(df.mathys2, names_from = words_mathys, values_from = slider_mathys2.response)
# mathys2_pca <- mathys2_pca[!(mathys2_pca$participant %in% out.ppt$ID), ]
#
# # Merge all questionnaires and export
# pca_all <- merge(panas_pca, mathys1_pca, by = 'participant')
# pca_all <- merge(pca_all, mathys2_pca, by = 'participant')
# colnames(pca_all)[colnames(pca_all) == 'participant'] <- 'participant_id'
# write.table(pca_all, pca_quest_out)

# PCA scores computed externally, read back in
pca <- read.table(pca_file)
colnames(pca)[colnames(pca) == 'participant_id'] <- 'ID'


# -----------------------------------------------------------------------------
# 4. Demographics
# -----------------------------------------------------------------------------

demographics$Age                 <- as.numeric(demographics$Age)
demographics$Ethnicity.simplified <- factor(demographics$Ethnicity.simplified)
demographics$Time.taken          <- demographics$Time.taken / 60   # seconds -> minutes
demographics$ID                  <- as.factor(demographics$Participant.id)
demographics                     <- demographics[demographics$Status == 'APPROVED', ]

# Sample description
summary_stats_demographics <- demographics %>%
  summarise(
    n          = n(),
    n_female   = sum(Sex == "Female", na.rm = TRUE),
    pct_female = n_female / n * 100,
    mean.Age   = mean(Age, na.rm = TRUE),
    sd.Age     = sd(Age),
    mean.Time  = mean(Time.taken, na.rm = TRUE),
    sd.Time    = sd(Time.taken)
  )

# # Ethnicity breakdown by group
# eth_counts <- demographics %>%
#   count(Group.x, Ethnicity.simplified, name = "n") %>%
#   group_by(Group.x) %>%
#   mutate(pct = n / sum(n) * 100) %>%
#   pivot_wider(names_from = Ethnicity.simplified,
#               values_from = c(n, pct), values_fill = 0)

# Drop flagged participants from demographics, then join onto the trial table
demographics <- demographics[!demographics$ID %in% out.ppt$ID, ]
df.target    <- merge(df.target, demographics, by = 'ID')


# -----------------------------------------------------------------------------
# 5. Per-experiment export, then pool both experiments
# -----------------------------------------------------------------------------

# Fix condition / block ordering for plotting and modelling
df.target$Condition <- factor(df.target$Condition, levels = c("HSF", "LSF", "FtC", "CtF"))
df.target$Block     <- factor(df.target$Block, levels = c(1, 2), labels = c("Block 1", "Block 2"))

# Tag the experiment and write this experiment's table.
# Run the script once with experiment = 1 and once with experiment = 2 so that
# BOTH exp1.csv and exp2.csv exist before the pooled step below.
df.target$Exp <- experiment
write.csv(df.target, file.path(csv_outdir, paste0("exp", experiment, ".csv")))

# Pool the two experiments for all analyses / figures that follow
exp1 <- read.csv(file.path(csv_outdir, "exp1.csv"))
exp2 <- read.csv(file.path(csv_outdir, "exp2.csv"))
df.target <- merge(exp1, exp2, all = TRUE)


# -----------------------------------------------------------------------------
# 6. Aggregations used by the figures
# -----------------------------------------------------------------------------

# Subject-level mean rating (per Target x Block x Condition x Exp)
agg.response <- df.target %>%
  group_by(ID, Target, Block, Condition, Exp) %>%
  summarise(mean = mean(Response, na.rm = TRUE),
            sd   = sd(Response, na.rm = TRUE),
            se   = sd(Response) / sqrt(n()))

# Cell-level mean rating
mean.resp <- df.target %>%
  group_by(Target, Block, Condition, Exp) %>%
  summarise(mean = mean(Response, na.rm = TRUE),
            sd   = sd(Response, na.rm = TRUE),
            se   = sd(Response) / sqrt(n()))

# Centre Age (used as covariate in the models below)
df.target$Age <- scale(df.target$Age, scale = FALSE)

# Keep only correctly answered Happy / Sad trials for the RT analyses
df.happy      <- df.target[df.target$Target == 'Happy', ]
df.sad        <- df.target[df.target$Target == 'Sad', ]
df.happy.corr <- df.happy[df.happy$correct.happy == 1, ]
df.sad.corr   <- df.sad[df.sad$correct.sad == 1, ]
df.hs         <- rbind(df.happy.corr, df.sad.corr, fill = TRUE)

# Split morphs by the response given (perceived negative / positive)
df.morph        <- df.target[df.target$Target == 'Morph', ]
df.morph$Target <- ifelse(df.morph$Response == 0, "Morph -", "Morph +")

# RT analysis table (correct Happy/Sad + morphs), plus log-RT (RTs are skewed)
df.rt       <- merge(df.hs, df.morph, all = TRUE)
df.rt$Block <- as.factor(df.rt$Block)
df.rt$logRT <- log(df.rt$RT)

# Subject-level mean RT (per Target x Block x Condition x Exp)
agg.rt <- df.rt %>%
  group_by(ID, Target, Block, Condition, Exp) %>%
  summarise(mean = mean(RT, na.rm = TRUE),
            se   = sd(RT) / sqrt(n()))

# Colour palettes: Block 1 (single image) vs Block 2 (sequence)
b1Palette <- c('#00416A', '#86bbd8')
b2Palette <- c('#7A5B8F', "#80B192")


# =============================================================================
# 7. BLOCK 1 — Emotional rating (Response)
# =============================================================================

agg.response.1 <- agg.response[agg.response$Block == "Block 1", ]

# ---- Figure: rating by target and filtering (single experiment) ------------
# Set filter(Exp == 1) or filter(Exp == 2) and the matching brackets below.
cat.b1 <-
  ggplot(agg.response.1 %>% filter(Exp == 1)) +
  aes(y = mean, x = Target, fill = Condition, color = Condition) +
  geom_boxplot(alpha = 0.7, width = 0.7, outlier.shape = NA,
               position = position_dodge(width = 0.85)) +
  geom_jitter(aes(color = Condition), alpha = 0.6, size = 1.4,
              position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.6)) +
  labs(x = 'Target Emotion', y = "Emotional Rating", fill = "Filtering", color = 'Filtering') +
  stat_summary(fun = mean, shape = 23, size = 0.8, stroke = 0.8,
               color = "black", alpha = 0.8, position = position_dodge(width = 0.85)) +
  theme_blank() +
  scale_color_manual(values = b1Palette) +
  scale_fill_manual(values = b1Palette) +
  ggtitle("Emotional Rating to Targets by Filtering") +
  theme(plot.title   = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 18, color = "grey30"),
        axis.text.y  = element_text(size = 26, color = "grey30"),
        strip.text.x = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.text  = element_text(size = 18)) +
  # y axis annotated with negative/positive poles
  scale_y_continuous(breaks = c(0, 0.15, 0.35, 0.5, 0.65, 0.85, 1),
                     labels = c("0", "-", "\u2193", "", "\u2191", "+", "1")) +
  labs(tag = "A") +
  theme(plot.tag.position = c(0.02, 0.98),
        plot.tag = element_text(size = 21, face = "bold")) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40", size = 0.45) +

  ## --- Significance brackets: Experiment 1 ---
  # annotate("segment", x = 0.8, xend = 1.2, y = 1.1,  yend = 1.1,  color = "black", size = 0.6) +  # Happy
  # annotate("text",    x = 1,   y = 1.12, label = "***", size = 7) +
  # annotate("segment", x = 2.8, xend = 3.2, y = 1.05, yend = 1.05, color = "black", size = 0.6) +  # Sad
  # annotate("text",    x = 3,   y = 1.07, label = "***", size = 7)

  ## --- Significance brackets: Experiment 2 ---
  annotate("segment", x = 0.8, xend = 1.2, y = 1.1,  yend = 1.1,  color = "black", size = 0.6) +   # Happy
  annotate("text",    x = 1,   y = 1.12, label = "**",  size = 7) +
  annotate("segment", x = 1.8, xend = 2.2, y = 1.09, yend = 1.09, color = "black", size = 0.6) +   # Morphs
  annotate("text",    x = 2,   y = 1.11, label = "***", size = 7)

# NOTE: saves the plot currently in `cat.b1`. Re-run this section with the other
# experiment (and its brackets) to populate the second file.
ggsave(file.path(fig_outdir, "Resp_B1_1st_Exp.png"), plot = cat.b1, width = 8, height = 6, units = "in")
ggsave(file.path(fig_outdir, "Resp_B1_2nd_Exp.png"), plot = cat.b1, width = 8, height = 6, units = "in")

# ---- Model: rating, single experiment --------------------------------------
model_block1 <- glmer(Response ~ Target * Condition + Age + Sex + (1 | ID),
                      data = df.target %>% filter(Block == "Block 1"),
                      family = binomial(link = 'logit'))
summary(model_block1)
Anova(model_block1, type = 3)

# Main effects
print(emmeans(model_block1, pairwise ~ Target))
print(emmeans(model_block1, pairwise ~ Condition))

# Planned comparisons: filtering within each target
print(emmeans(model_block1, pairwise ~ Condition | Target))

# Are morphs rated positive vs negative overall? (binomial test against chance)
df_morph_block_1 <- df.target %>% filter(Block == "Block 1", Target == "Morph")
k <- sum(df_morph_block_1$Response == 1)   # 1 = positive response
n <- nrow(df_morph_block_1)
binom.test(k, n, p = 0.5)


# ---- Model: rating, pooled (Experiment 1 vs 2) -----------------------------
# Exploratory both-experiments plot (faceted by Exp)
ggplot(agg.response.1) +
  aes(y = mean, x = Target, fill = Condition, color = Condition) +
  geom_boxplot(alpha = 0.7, width = 0.7, outlier.shape = NA,
               position = position_dodge(width = 0.85)) +
  geom_jitter(aes(color = Condition), alpha = 0.6, size = 1.4,
              position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.6)) +
  facet_wrap(~ Exp) +
  labs(x = 'Target Emotion', y = "Emotional Rating", fill = "Filtering", color = 'Filtering') +
  stat_summary(fun = mean, shape = 23, size = 0.8, stroke = 0.8,
               color = "black", alpha = 0.8, position = position_dodge(width = 0.85)) +
  theme_blank() +
  scale_color_manual(values = b1Palette) +
  scale_fill_manual(values = b1Palette)

model_block1 <- glmer(Response ~ Target * Condition * Exp + Age + Sex + (1 | Exp:ID),
                      data = df.target %>% filter(Block == "Block 1"),
                      family = binomial(link = 'logit'),
                      control = glmerControl(optimizer = "bobyqa",
                                             optCtrl = list(maxfun = 2e5)))
summary(model_block1)
Anova(model_block1, type = 3)

print(emmeans(model_block1, pairwise ~ Target))
print(emmeans(model_block1, pairwise ~ Condition))

# HSF vs LSF between experiments, per target
print(emmeans(model_block1, pairwise ~ Exp | Target * Condition))

# Difference-of-differences: (filtering effect) x (experiment) per target
emmeans_obj <- emmeans(model_block1, ~ Condition * Exp | Target)
print(contrast(emmeans_obj, interaction = c(Condition = "pairwise", Exp = "pairwise")))


# =============================================================================
# 8. BLOCK 1 — Reaction time
# =============================================================================

agg.rt$Target <- factor(agg.rt$Target, levels = c("Happy", "Morph +", "Sad", "Morph -"))
agg.rt.1      <- agg.rt[agg.rt$Block == "Block 1", ]

# ---- Figure: RT by target and filtering (single experiment) ----------------
rt.b1 <-
  ggplot(agg.rt.1) + aes(y = mean, x = Target, fill = Condition, color = Condition) +
  geom_violin(alpha = 0.65, position = position_dodge(width = 1.1)) +
  labs(x = "Target Emotion", y = "log(RT)", title = 'RT to Targets by Filtering',
       fill = 'Filtering', color = 'Filtering') +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black",
               fatten = 2, position = position_dodge(width = 1.1)) +
  geom_jitter(aes(color = Condition), alpha = 1, size = 1.6, shape = 21,
              position = position_jitterdodge(jitter.width = 0.7, dodge.width = 1.1)) +
  theme_blank() +
  theme(plot.title   = element_text(size = 21, margin = margin(b = 35)),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 18, color = 'grey30'),
        axis.text.y  = element_text(size = 18, color = 'grey30'),
        strip.text.x = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.text  = element_text(size = 18)) +
  scale_color_manual(values = b1Palette) +
  scale_fill_manual(values = b1Palette) +
  labs(tag = "B") +
  theme(plot.tag.position = c(0.02, 0.98),
        plot.tag = element_text(size = 21, face = "bold")) +

  ## --- Significance brackets: Experiment 1 ---
  annotate("segment", x = 0.8, xend = 1.2, y = 1240, yend = 1240, color = "black", size = 0.6) +   # Happy
  annotate("text",    x = 1,   y = 1260, label = "***", size = 7) +
  annotate("segment", x = 2.8, xend = 3.2, y = 1330, yend = 1330, color = "black", size = 0.6) +   # Sad
  annotate("text",    x = 3,   y = 1350, label = "***", size = 7)

  ## --- Significance brackets: Experiment 2 ---
  # annotate("segment", x = 0.8, xend = 1.2, y = 920,  yend = 920,  color = "black", size = 0.6) +  # Happy
  # annotate("text",    x = 1,   y = 940,  label = "***", size = 7) +
  # annotate("segment", x = 3.8, xend = 4.2, y = 1110, yend = 1110, color = "black", size = 0.6) +  # Morph -
  # annotate("text",    x = 4,   y = 1130, label = "*",   size = 7)

ggsave(file.path(fig_outdir, "RT_B1_1st_Exp.png"), plot = rt.b1, width = 8, height = 6, units = "in")
ggsave(file.path(fig_outdir, "RT_B1_2nd_Exp.png"), plot = rt.b1, width = 8, height = 6, units = "in")

# ---- Model: RT, single experiment ------------------------------------------
model_rt_block1 <- lmer(logRT ~ Target * Condition + Age + Sex + (1 | ID),
                        data = df.rt %>% filter(Block == "Block 1"))
summary(model_rt_block1)
Anova(model_rt_block1, type = 3)

print(emmeans(model_rt_block1, pairwise ~ Target))
print(emmeans(model_rt_block1, pairwise ~ Condition))
print(emmeans(model_rt_block1, pairwise ~ Condition | Target))   # planned comparisons

# Descriptives
aggregate(RT ~ Target,    data = df.rt %>% filter(Block == "Block 1"), mean)
aggregate(RT ~ Condition, data = df.rt %>% filter(Block == "Block 1"), mean)
desc_stat_RTs <- df.rt %>%
  group_by(Target, Condition) %>%
  summarise(meanRT = mean(RT, na.rm = TRUE),
            sdRT   = sd(RT, na.rm = TRUE))

# ---- Model: RT, pooled (Experiment 1 vs 2) ---------------------------------
# Exploratory both-experiments plot (faceted by Exp)
ggplot(agg.rt.1) + aes(y = mean, x = Target, fill = Condition, color = Condition) +
  geom_violin(alpha = 0.65, position = position_dodge(width = 1.1)) +
  labs(x = "Target Emotion", y = "log(RT)", title = 'RT to Targets by Filtering',
       fill = 'Filtering', color = 'Filtering') +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black",
               fatten = 2, position = position_dodge(width = 1.1)) +
  facet_wrap(~ Exp) +
  geom_jitter(aes(color = Condition), alpha = 1, size = 1.6, shape = 21,
              position = position_jitterdodge(jitter.width = 0.7, dodge.width = 1.1)) +
  theme_blank() +
  scale_color_manual(values = b1Palette) +
  scale_fill_manual(values = b1Palette)

model_rt_block1 <- lmer(logRT ~ Target * Condition * Exp + Age + Sex + (1 | ID),
                        data = df.rt %>% filter(Block == "Block 1"))
summary(model_rt_block1)
Anova(model_rt_block1, type = 3)

print(emmeans(model_rt_block1, pairwise ~ Target))
print(emmeans(model_rt_block1, pairwise ~ Condition))
print(emmeans(model_rt_block1, pairwise ~ Exp | Target * Condition))

emmeans_obj <- emmeans(model_rt_block1, ~ Condition * Exp | Target)
print(contrast(emmeans_obj, interaction = c(Condition = "pairwise", Exp = "pairwise")))


# =============================================================================
# 9. BLOCK 2 — Emotional rating (Response)
# =============================================================================

agg.response.2 <- agg.response[agg.response$Block == "Block 2", ]

# ---- Figure: rating by target and sequence (single experiment) -------------
cat.b2 <-
  ggplot(agg.response.2) +
  aes(y = mean, x = Target, fill = Condition, color = Condition) +
  geom_boxplot(alpha = 0.7, width = 0.7, outlier.shape = NA,
               position = position_dodge(width = 0.85)) +
  geom_jitter(aes(color = Condition), alpha = 0.6, size = 1.4,
              position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.6)) +
  labs(x = 'Target Emotion', y = "Response", fill = "Sequence", color = 'Sequence') +
  stat_summary(fun = mean, shape = 23, size = 0.8, stroke = 0.8,
               color = "black", alpha = 0.8, position = position_dodge(width = 0.85)) +
  theme_blank() +
  scale_color_manual(values = b2Palette) +
  scale_fill_manual(values = b2Palette) +
  ggtitle("Response to Targets by Sequence") +
  theme(plot.title   = element_text(size = 21),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 18, color = 'grey30'),
        axis.text.y  = element_text(size = 26, color = "grey30"),
        strip.text.x = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.text  = element_text(size = 18)) +
  labs(tag = "A") +
  theme(plot.tag.position = c(0.02, 0.98),
        plot.tag = element_text(size = 21, face = "bold")) +
  scale_y_continuous(breaks = c(0, 0.15, 0.35, 0.5, 0.65, 0.85, 1),
                     labels = c("0", "-", "\u2193", "", "\u2191", "+", "1")) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40", size = 0.45) +

  ## --- Significance brackets: Experiment 2 ---
  annotate("segment", x = 0.8, xend = 1.2, y = 1.1,  yend = 1.1,  color = "black", size = 0.6) +   # Happy
  annotate("text",    x = 1,   y = 1.12, label = "*", size = 7) +
  annotate("segment", x = 1.8, xend = 2.2, y = 1.05, yend = 1.05, color = "black", size = 0.6) +   # Morph +
  annotate("text",    x = 2,   y = 1.07, label = "*", size = 7)

ggsave(file.path(fig_outdir, "Resp_B2_1st_Exp.png"), plot = cat.b2, width = 8, height = 6, units = "in")
ggsave(file.path(fig_outdir, "Resp_B2_2nd_Exp.png"), plot = cat.b2, width = 8, height = 6, units = "in")

# ---- Model: rating, single experiment --------------------------------------
model_block2 <- glmer(Response ~ Target * Condition + Age + Sex + (1 | ID),
                      data = df.target %>% filter(Block == "Block 2"),
                      family = binomial(link = 'logit'))
Anova(model_block2, type = 3)

print(emmeans(model_block2, pairwise ~ Target))
print(emmeans(model_block2, pairwise ~ Condition))
print(emmeans(model_block2, pairwise ~ Condition | Target))   # planned comparisons

# ---- Model: rating, pooled (Experiment 1 vs 2) -----------------------------
# Exploratory both-experiments plot (faceted by Exp)
ggplot(agg.response.2) +
  aes(y = mean, x = Target, fill = Condition, color = Condition) +
  geom_boxplot(alpha = 0.7, width = 0.7, outlier.shape = NA,
               position = position_dodge(width = 0.85)) +
  geom_jitter(aes(color = Condition), alpha = 0.6, size = 1.4,
              position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.6)) +
  facet_wrap(~ Exp) +
  labs(x = 'Target Emotion', y = "Emotional Rating", fill = "Filtering", color = 'Filtering') +
  stat_summary(fun = mean, shape = 23, size = 0.8, stroke = 0.8,
               color = "black", alpha = 0.8, position = position_dodge(width = 0.85)) +
  theme_blank() +
  scale_color_manual(values = b2Palette) +
  scale_fill_manual(values = b2Palette)

model_block2 <- glmer(Response ~ Target * Condition * Exp + Age + Sex + (1 | Exp:ID),
                      data = df.target %>% filter(Block == "Block 2"),
                      family = binomial(link = 'logit'),
                      control = glmerControl(optimizer = "bobyqa",
                                             optCtrl = list(maxfun = 2e5)))
summary(model_block2)
Anova(model_block2, type = 3)

print(emmeans(model_block2, pairwise ~ Target))
print(emmeans(model_block2, pairwise ~ Condition))
print(emmeans(model_block2, pairwise ~ Exp | Target * Condition))

emmeans_obj <- emmeans(model_block2, ~ Condition * Exp | Target)
print(contrast(emmeans_obj, interaction = c(Condition = "pairwise", Exp = "pairwise")))


# =============================================================================
# 10. BLOCK 2 — Reaction time
# =============================================================================

agg.rt.2 <- agg.rt[agg.rt$Block == "Block 2", ]

# ---- Figure: RT by target and sequence (single experiment) -----------------
rt.b2 <-
  ggplot(agg.rt.2) + aes(y = mean, x = Target, fill = Condition, color = Condition) +
  geom_violin(alpha = 0.65, position = position_dodge(width = 1.1)) +
  labs(x = "Target Emotion", y = "RT (ms)", title = 'RT to Targets by Sequence',
       color = 'Sequence', fill = 'Sequence') +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black",
               fatten = 2, position = position_dodge(width = 1.1)) +
  geom_jitter(aes(color = Condition), alpha = 1, size = 1.6, shape = 21,
              position = position_jitterdodge(jitter.width = 0.7, dodge.width = 1.1)) +
  theme_blank() +
  theme(plot.title   = element_text(size = 21, margin = margin(b = 35)),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x  = element_text(size = 18, color = 'grey30'),
        axis.text.y  = element_text(size = 18, color = 'grey30'),
        strip.text.x = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.text  = element_text(size = 18)) +
  scale_color_manual(values = b2Palette) +
  scale_fill_manual(values = b2Palette) +
  labs(tag = "B") +
  theme(plot.tag.position = c(0.02, 0.98),
        plot.tag = element_text(size = 21, face = "bold")) +

  ## --- Significance brackets: Experiment 1 ---
  annotate("segment", x = 0.8, xend = 1.2, y = 1000, yend = 1000, color = "black", size = 0.6) +   # Happy
  annotate("text",    x = 1,   y = 1020, label = "***", size = 7) +
  annotate("segment", x = 2.8, xend = 3.2, y = 940,  yend = 940,  color = "black", size = 0.6) +   # Sad
  annotate("text",    x = 3,   y = 960,  label = "***", size = 7) +
  annotate("segment", x = 1.8, xend = 2.2, y = 1240, yend = 1240, color = "black", size = 0.6) +   # Morph +
  annotate("text",    x = 2,   y = 1260, label = "***", size = 7)

  ## --- Significance brackets: Experiment 2 ---
  # annotate("segment", x = 0.8, xend = 1.2, y = 940,  yend = 940,  color = "black", size = 0.6) +  # Happy
  # annotate("text",    x = 1,   y = 960,  label = "*",   size = 7) +
  # annotate("segment", x = 1.8, xend = 2.2, y = 1240, yend = 1240, color = "black", size = 0.6) +  # Morph +
  # annotate("text",    x = 2,   y = 1260, label = "***", size = 7)

ggsave(file.path(fig_outdir, "RT_B2_1st_Exp.png"), plot = rt.b2, width = 8, height = 6, units = "in")
ggsave(file.path(fig_outdir, "RT_B2_2nd_Exp.png"), plot = rt.b2, width = 8, height = 6, units = "in")

# ---- Model: RT, single experiment ------------------------------------------
model_rt_block2 <- lmer(logRT ~ Target * Condition + Age + Sex + (1 | ID),
                        data = df.rt %>% filter(Block == "Block 2"))
summary(model_rt_block2)
Anova(model_rt_block2, type = 3)

print(emmeans(model_rt_block2, pairwise ~ Target))
print(emmeans(model_rt_block2, pairwise ~ Condition))
print(emmeans(model_rt_block2, pairwise ~ Condition | Target))   # planned comparisons

# ---- Model: RT, pooled (Experiment 1 vs 2) ---------------------------------
# Exploratory both-experiments plot (faceted by Exp)
ggplot(agg.rt.2) + aes(y = mean, x = Target, fill = Condition, color = Condition) +
  geom_violin(alpha = 0.65, position = position_dodge(width = 1.1)) +
  labs(x = "Target Emotion", y = "log(RT)", title = 'RT to Targets by Filtering',
       fill = 'Filtering', color = 'Filtering') +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black",
               fatten = 2, position = position_dodge(width = 1.1)) +
  facet_wrap(~ Exp) +
  geom_jitter(aes(color = Condition), alpha = 1, size = 1.6, shape = 21,
              position = position_jitterdodge(jitter.width = 0.7, dodge.width = 1.1)) +
  theme_blank() +
  scale_color_manual(values = b2Palette) +
  scale_fill_manual(values = b2Palette)

model_rt_block2 <- lmer(logRT ~ Target * Condition * Exp + Age + Sex + (1 | ID),
                        data = df.rt %>% filter(Block == "Block 2"))
summary(model_rt_block2)
Anova(model_rt_block2, type = 3)

print(emmeans(model_rt_block2, pairwise ~ Target))
print(emmeans(model_rt_block2, pairwise ~ Condition))
print(emmeans(model_rt_block2, pairwise ~ Exp | Target * Condition))

emmeans_obj <- emmeans(model_rt_block2, ~ Condition * Exp | Target)
print(contrast(emmeans_obj, interaction = c(Condition = "pairwise", Exp = "pairwise")))


# =============================================================================
# 11. Mood (PCA) analyses — morph responses & RT as a function of PC1
# =============================================================================

pca$ID <- pca$participant_id

# Morph responses aggregated with mood scores
agg.morph <- aggregate(df.morph, Response ~ ID + Block + Condition + Age + Sex, mean)
agg.morph <- merge(agg.morph, pca, by = 'ID')

ggplot(agg.morph) +
  aes(x = Response, y = PC1) +
  geom_point(colour = "#112446") +
  theme_minimal() +
  facet_wrap(~ Block * Condition) +
  geom_smooth(method = 'lm')

# Trial-level morph data joined with mood scores
mood.df       <- left_join(df.morph, pca, by = "ID")
mood.df$ID    <- as.factor(mood.df$ID)
mood.df$Block <- as.factor(mood.df$Block)

# Probability of a positive answer to morphs as a function of mood (PC1)
ggplot(mood.df, aes(x = PC1, y = Response)) +
  geom_point() +
  stat_smooth(method = "glm", method.args = list(family = binomial(link = 'logit')),
              se = TRUE, color = '#2f4858') +
  xlab("Mood (PC1)") + ylab("Probability of Positive Answer") +
  facet_wrap(~ Block) +
  ggtitle("Probability of positive answer to morphs by Mood") +
  theme_blank()

# Mood -> morph response, per block
model_block1 <- glmer(Response ~ PC1 + Age + Sex + (1 | ID),
                      data = mood.df %>% filter(Block == "Block 1"),
                      family = binomial(link = 'logit'))
summary(model_block1)

model_block2 <- glmer(Response ~ PC1 + Age + Sex + (1 | ID),
                      data = mood.df %>% filter(Block == "Block 2"),
                      family = binomial(link = 'logit'))
summary(model_block2)

# Does mood facilitate RT to positive / negative morphs?
mood.df$RTlog <- log(mood.df$RT)

model_rt_block1 <- lmer(RTlog ~ PC1 * Response + Age + Sex + (1 | ID),
                        data = mood.df %>% filter(Block == "Block 1"))
summary(model_rt_block1)

# RT summary by response (or by Sex — change the grouping variable)
sum <- mood.df %>%
  group_by(Block, Response) %>%
  summarise(RT.mean = mean(RT, na.rm = TRUE),
            RT.sd   = sd(RT, na.rm = TRUE))

# Age vs RT (Block 1)
ggplot(data = mood.df %>% filter(Block == "Block 1"), aes(x = Age, y = RT)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Effect of age on reaction time", x = "Age", y = "Reaction time (log RT)") +
  theme_minimal()

model_rt_block2 <- lmer(RTlog ~ PC1 * Response + Age + Sex + (1 | ID),
                        data = mood.df %>% filter(Block == "Block 2"))
summary(model_rt_block2)

# Mood (PC1) x Response interaction on log RT (Block 2)
ggplot(mood.df %>% filter(Block == "Block 2"),
       aes(x = PC1, y = RTlog, color = as.factor(Response))) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Mood score (PC1)", y = "Log RT", color = "Response",
       title = "Mood (PC1) x Response interaction on log RT") +
  theme_minimal()

# Age x Response on RT (Block 2)
ggplot(mood.df %>% filter(Block == "Block 2"),
       aes(x = Age, y = RT, color = as.factor(Response))) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Age", y = "RT", color = "Response",
       title = "Age x Response interaction on RT") +
  theme_minimal()


# =============================================================================
# 12. PCA loading figures — item contributions on PC1 / PC2
# =============================================================================
# Loadings are read from the externally computed PCA output.
# (`pc1_exp1`, `pc1_exp2`, `pc2_exp1` are the corresponding PC1/other-experiment
#  plot objects built with the same recipe as `pc2_exp2` below.)

pca <- read.csv2(pc2_file)

# Palette for the three questionnaire scales
pal <- c("#f26419", '#f6ae2d', '#592720')

# Optionally put PANAS first in the legend
# pca$Scale <- relevel(as.factor(pca$Scale), ref = "PANAS")

# Order items by loading and keep only |loading| > 0.15
pca <- pca %>%
  mutate(Feature = reorder(feature, score)) %>%
  arrange(desc(score))
pca <- pca[pca$score > 0.15 | pca$score < -0.15, ]

pc2_exp2 <- ggplot(pca) +
  aes(x = score, y = feature, fill = Scale) +
  geom_bar(stat = "summary", alpha = 0.7, width = 0.7) +
  labs(x = 'PC2 scores', y = "Items") +
  scale_fill_manual(values = pal) +
  theme_blank() +
  ggtitle("Items loadings on PC2") +
  theme(plot.title   = element_text(size = 41, face = 'bold'),
        axis.title.x = element_text(size = 39, color = 'black'),
        axis.title.y = element_text(size = 39, color = 'black'),
        axis.text.x  = element_text(size = 36, color = 'black'),
        axis.text.y  = element_text(size = 28, color = '#393E46'),
        legend.text  = element_text(size = 38, color = 'black'),
        legend.title = element_text(size = 39, color = 'black'))

ggsave(file.path(fig_outdir, "PC1_exp1.png"), plot = pc1_exp1, width = 20, height = 15, units = "in")
ggsave(file.path(fig_outdir, "PC1_exp2.png"), plot = pc1_exp2, width = 20, height = 15, units = "in")
ggsave(file.path(fig_outdir, "sup_mat/PC2_exp1.png"), plot = pc2_exp1, width = 20, height = 15, units = "in")
ggsave(file.path(fig_outdir, "sup_mat/PC2_exp2.png"), plot = pc2_exp2, width = 20, height = 15, units = "in")
