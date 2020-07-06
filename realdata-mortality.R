options(stringsAsFactors = F)
library(ggplot2)
library(dplyr)
library(tidyr)

setwd("~/Documents/kiwis/changepointdet/covid/")

dt = read.csv("eurostat/demo_r_mweek3_1_Data.csv")
dt = mutate(dt, ValueNum = sub(":", "", Value)) %>%
  mutate(ValueNum = as.numeric(sub(",", "", ValueNum)))
dt

dt = filter(dt, !GEO  %in% c("Germany (until 1990 former territory of the FRG)", "United Kingdom"))
dt = filter(dt, GEO %in% c("France", "Italy", "Spain", "Sweden", "Switzerland"))

pl = filter(dt, AGE!="TOTAL")
ggplot(pl, aes(x=seq_along(TIME), y=ValueNum, col=AGE)) + 
  geom_point() + facet_wrap(~GEO, scales="free_y") +
  theme_minimal()

filter(dt, AGE=="TOTAL") %>%
  ggplot(aes(x=seq_along(TIME), y=ValueNum, col=AGE)) + 
  geom_point() + facet_wrap(~GEO, scales="free_y") +
  theme_minimal()

filter(dt, GEO=="Spain") %>%
  ggplot(aes(x=seq_along(TIME), y=ValueNum)) + 
  geom_point() + facet_wrap(~AGE, scales="free_y") +
  theme_minimal()

tsdf = filter(dt, GEO=="Spain", AGE==AGE=="55-59" | AGE=="Y60-64") %>%
  group_by(TIME) %>%
  summarize(ValueNum = mean(ValueNum))
  
tsdf = separate(tsdf, "TIME", c("YR", "WK"), sep="W", remove=F) %>%
  mutate(YR=as.numeric(YR), WK=as.numeric(WK))

# two missing values at the end:
tail(tsdf)
which(is.na(tsdf$ValueNum))
# also filter to an exact 4 yr period (2017W22 to 2020W21)
tsdf = filter(tsdf, !is.na(ValueNum), YR>2017 | WK>21)
table(tsdf$YR)

ts = tsdf$ValueNum
plot(ts)

setwd("~/Documents/gitrep/changepoint-detection/")
source("CLEAN-algorithm.R")

## Apply Algorithm 1:
pen = autoset_penalty(ts)
maxpeaklen = 10
SD = sd(ts[1:52])
alg1 = shortsgd(ts, MAXLOOKBACK=maxpeaklen, PEN=pen, SD=SD)
print(alg1$segs)

finaltheta = alg1$wt[length(alg1$wt)]
finaltheta

## Apply Algorithm 2:
pen = autoset_penalty(ts)
maxpeaklen = 10
alg2 = fulldetector_prune(ts, theta0=finaltheta, MAXLOOKBACK=maxpeaklen, PEN=pen, PEN2=pen, BURNIN=2, SD=SD, prune=2)
print(alg2$segs)

alg2res = data.frame(alg2$segs)
colnames(alg2res) = c("start", "end", "theta", "segtype")
rownames(alg2res) = NULL
alg2res$segtype = ifelse(alg2res$segtype==1, "S", "N")

# fix time for plotting
library(lubridate)
tsdf$date = parse_date_time(paste(tsdf$YR, tsdf$WK, "1"), "%Y %W %w")
alg2res$stdate = min(tsdf$date) + dweeks(alg2res$start-1)
alg2res$edate = min(tsdf$date) + dweeks(alg2res$end-1)

p2 = ggplot(tsdf) +
  geom_point(aes(x=as.Date(date), y=ValueNum), size=0.5) +
  theme_bw()+ theme(text=element_text(size=14), axis.text.x=element_text(size=12, angle=30, vjust=0.8, hjust=1)) +
  scale_color_manual(values=c("S"="blue", "N"="gold"), name="detections", labels=c("nuisance", "signal")) +
  ylab("weekly deaths") + xlab(NULL) + # ggtitle("Weekly death count in Spain over 2017-2020, ages 55-64") +
  geom_segment(aes(as.Date(stdate), 220, xend=as.Date(edate), yend=220, color=segtype, size=segtype),
              data=alg2res) +
  scale_x_date(date_breaks = "3 months", date_labels="%Y-%m", minor_breaks = NULL) +
  scale_size_manual(values=c("false positive"=2, correct=0.7, "N"=7, "S"=5.5)) +
  guides(size="none", color=guide_legend(order=2, override.aes=list(size=7)))
print(p2)

ggsave(paste0("~/Documents/gitrep/drafts/changepoint-method/results-appl/", "det-mortality.png"), width=20, height=8.5, units="cm")
