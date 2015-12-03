library('seqinr')
library('stringr')
library('dplyr')
library('reshape2')
library('ggplot2')

patient_data <- read.fasta("INUTIL-Trial-sample.fst") %>%
  names() %>%
  str_split('_((CD4|VL)/|day)')

# Name columns extracted from raw FASTA

ptdf <- data_frame(
  patient = lapply(patient_data, '[[', 1) %>%
    unlist %>%
    str_extract('[^Pt].+'),
  day = lapply(patient_data, '[[', 2) %>% unlist,
  CD4 = lapply(patient_data, '[[', 3) %>% unlist,
  VL = lapply(patient_data, '[[', 4) %>% unlist,
  VL_log = lapply(patient_data, '[[', 4) %>% unlist %>% as.numeric %>% log10
)

# Function to later make a new column with: VL and CD4 change over time course

TimeSeriesDiffs <- function (pt.id, stat, df = ptdf.long) {
  # get change in stat ('CD4' or 'VL') for a given patient ID ('01'...'17')
  df %>%
    filter(patient == pt.id) %>%
    group_by(variable) %>%
    summarise(delta = diff(value)) %>%
    filter(variable == stat) %>%
    select(delta) %>%
    as.numeric %>%
    return
}

deltas <- lapply(patient.ids,
       function(pt.id) {
         data_frame(patient = pt.id,
                    cd4.del = TimeSeriesDiffs(pt.id, 'CD4'),
                    vl.del = TimeSeriesDiffs(pt.id, 'VL_log'))
         }) %>%
  rbind_all


DeltaVal <- function(pt.id, stat, df=deltas) {
#   df %>%
#     filter(patient == pt.id) %>%
#     # basic subset operation - seems no better way than anon func!
#     (function(x) { x[stat %>%
#                        tolower %>%
#                        paste0('.del')]
#     }
#     ) %>%
#     as.numeric
  col <- stat %>% tolower %>% paste0('.del')
  df[which(df$patient == pt.id),][col] %>%
    na.omit %>%
    as.numeric %>%
    return
}

# make the day explicitly ordered categoric [factor], "Start" and "End" levels
ptdf.long <- melt(
  ptdf, id.vars = c('patient','day')) %>%
  mutate(cat_day = ifelse(day == 0, 'Start', 'End') %>%
           factor(levels = c("Start","End")),
         value = value %>% as.numeric)

# ugly ye olde base R table join operation to delta values, oh well

cd4.delta.vals <- lapply(ptdf.long$patient %>% unique, function(x){
  DeltaVal(x, 'CD4')
}) %>% unlist
df.sub.cd4 <- ptdf.long[which(ptdf.long$variable == 'CD4'),]
df.sub.cd4.deltavals <- lapply(cd4.delta.vals, function(x){
  rep(x,2)
}) %>% unlist
df.sub.cd4$delta <- df.sub.cd4.deltavals

vl.delta.vals <- lapply(ptdf.long$patient %>% unique, function(x){
  DeltaVal(x, 'VL')
}) %>% unlist
df.sub.vl <- ptdf.long[which(ptdf.long$variable == 'VL_log'),]
df.sub.vl.deltavals <- lapply(vl.delta.vals, function(x){
  rep(x,2)
}) %>% unlist
df.sub.vl$delta <-df.sub.vl.deltavals

# Add a "codelta" column for corresponding value in VL/CD4 and vice versa
df.sub.cd4$codelta <- df.sub.vl.deltavals
df.sub.vl$codelta <- df.sub.cd4.deltavals

ptdf.full <- rbind(df.sub.cd4, df.sub.vl)

# Plot CD4 change
cd4.ptdf.full <- ptdf.full %>% filter(variable == 'CD4')
cd4.plot <- ggplot(cd4.ptdf.full,
  aes(cat_day,value,group=patient,color=codelta)) +
  geom_line() +
  scale_x_discrete('Time point', expand=c(0,0)) +
  scale_y_continuous('CD4 count') +
  ggtitle('HIV CD4 count in cohort over approx. 4 weeks') +
  theme_light() +
  theme(panel.grid=element_blank()) +
  scale_colour_gradient(limits=range(cd4.ptdf.full$codelta),
                        low="#ff3232", high="#47c147",
                        name="log"~Delta~"viral load")
# low is green high is red here

# Plot viral load change
vl.ptdf.full <- ptdf.full %>% filter(variable == 'VL_log')
vl.plot <- ggplot(vl.ptdf.full,
                  aes(cat_day,value,group=patient,color=codelta)) +
  geom_line() +
  scale_x_discrete('Time point', expand=c(0,0)) +
  scale_y_continuous() +
  ylab(bquote('Viral load ('~log[10]~')')) +
  ggtitle('HIV viral load in cohort over approx. 4 weeks') +
  theme_light() +
  theme(panel.grid=element_blank()) +
  scale_colour_gradient(limits=range(vl.ptdf.full$codelta),
                        low="#47c147", high="#ff3232",
                        name=~Delta~"CD4 count")
# low is red high is green here


### cutoff of 4 from graph

ptdf.long %>%
  filter(cat_day == "End",
         variable == "VL_log") %>%
  filter(value > 4)

### matches patients with G48V mutation