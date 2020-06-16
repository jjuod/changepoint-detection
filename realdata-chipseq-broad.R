# Real data analysis: Chipseq Broad

options(stringsAsFactors = F)
library(ggplot2)
library(dplyr)

setwd("~/Documents/gitrep/changepoint-detection/")
source("CLEAN-algorithm.R")

indir = "~/Documents/kiwis/changepointdet/chipseq/BROAD/"
bedC = read.table(paste0(indir, "filtered_control.bedGraph"))
bedH = read.table(paste0(indir, "filtered_h3k27ac.bedGraph"))

plot(bedC$V2, bedC$V4)
plot(bedH$V2, bedH$V4)

bedCf = filter(bedC, V2>120100e3, V3<=120700e3)
bedHf = filter(bedH, V2>120100e3, V3<=120700e3)

chr = "chr1"

# align different start positions
actualstart = max(min(bedCf$V2), min(bedHf$V2))
actualend = min(max(bedCf$V2), max(bedHf$V2))
bedCf = filter(bedCf, V2>=actualstart, V3<=actualend)
bedHf = filter(bedHf, V2>=actualstart, V3<=actualend)

# sections with no reads appear as gaps - need to fill those in with 0 counts
toappend = data.frame()
for (i in 2:nrow(bedCf)){
  if(bedCf$V2[i]-bedCf$V3[i-1]>0){
    neededpos = seq(bedCf$V3[i-1], bedCf$V2[i]-1, 25)
    toappend = bind_rows(toappend, data.frame(V1=chr, V2=neededpos, V3=neededpos+25, V4=0))
  }
}
bedCf = bind_rows(bedCf, toappend) %>% arrange(V2)
# same for treated
toappend = data.frame()
for (i in 2:nrow(bedHf)){
  if(bedHf$V2[i]-bedHf$V3[i-1]>0){
    neededpos = seq(bedHf$V3[i-1], bedHf$V2[i]-1, 25)
    toappend = bind_rows(toappend, data.frame(V1=chr, V2=neededpos, V3=neededpos+25, V4=0))
  }
}
bedHf = bind_rows(bedHf, toappend) %>% arrange(V2)

range(lead(bedCf$V2)-bedCf$V3, na.rm = T)
range(lead(bedHf$V2)-bedHf$V3, na.rm = T)
head(bedCf)
head(bedHf)

# convert to long vectors (for each bp)
bedCv = rep(bedCf$V4, bedCf$V3-bedCf$V2)
bedHv = rep(bedHf$V4, bedHf$V3-bedHf$V2)
length(bedCv); length(bedHv)

# downsample to a manageable number
DSfactor = 500
bedCds = data.frame(pos = round((seq_along(bedCv)+actualstart)/DSfactor)*DSfactor, m=bedCv) %>%
  group_by(pos) %>%
  summarize(m = mean(m))
bedHds = data.frame(pos = round((seq_along(bedHv)+actualstart)/DSfactor)*DSfactor, m=bedHv) %>%
  group_by(pos) %>%
  summarize(m = mean(m))

ggplot(bind_rows(bedCds,bedHds, .id="trt")) + geom_point(aes(pos, m, col=trt), alpha=0.5)

# Apply Algorithm 1
ts = bedHds$m
pen = autoset_penalty(ts)
maxpeaklen = round(50e3/DSfactor)
alg1 = shortsgd(ts, MAXLOOKBACK=maxpeaklen, PEN=pen, SD=sd(ts))
alg1$segs

# final theta estimate:
finaltheta = alg1$wt[length(alg1$wt)]
print(finaltheta)

# we are interested only in segments of increased mean:
alg1incr = alg1$segs[alg1$segs[,3]>finaltheta,,drop=F]

# empty dataframe with a "segment" at 0-1 to allow plotting
alg1res = data.frame(chrom=chr, chromStart=actualstart, chromEnd=actualstart+1)
if(nrow(alg1incr)>0){
  alg1res = data.frame(chrom=chr,
                       chromStart=alg1incr[,1]*DSfactor + actualstart,
                       chromEnd=alg1incr[,2]*DSfactor + actualstart+1, 
                       row.names = NULL)  
}
print(alg1res)

p1 = ggplot(bind_rows(bedCds,bedHds, .id="trt")) +
  geom_point(aes(pos/1000, m, col=trt), alpha=0.5) +
  geom_segment(aes(chromStart/1000, -3, xend=chromEnd/1000, yend=-3),
               data=alg1res, color="deepskyblue", size=3)+
  xlab("position on chromosome, kbp") + ylab("coverage") +
  theme_bw()
print(p1)


## Apply Algorithm 2:
ts = bedHds$m
pen = autoset_penalty(ts)
maxpeaklen = round(50e3/DSfactor)
alg2 = fulldetector_prune(ts, theta0=finaltheta, MAXLOOKBACK=maxpeaklen, PEN=pen, PEN2=pen, BURNIN=2, SD=sd(ts), prune=2)
print(alg2$segs)

# we are interested only in segments of increased mean:
alg2incr = alg2$segs[alg2$segs[,3]>finaltheta,, drop=F]

alg2res = data.frame(chrom=chr,
                     chromStart=alg2incr[,1]*DSfactor + actualstart,
                     chromEnd=alg2incr[,2]*DSfactor + actualstart+DSfactor, 
                     segtype=ifelse(alg2incr[,4]==1, "S", "N"),
                     row.names = NULL)
print(alg2res)

p2 = ggplot(bind_rows(bedCds,bedHds, .id="trt")) +
  geom_point(aes(pos/1000, m, col=trt), alpha=0.5) +
  geom_segment(aes(chromStart/1000, -3, xend=chromEnd/1000, yend=-3),
               data=alg2res[alg2res$segtype=="S",], color="deepskyblue", size=3)+
  xlab("position on chromosome, kbp") + ylab("coverage") +
  theme_bw()
if("N" %in% alg2res$segtype){
  p2 = p2 + geom_segment(aes(chromStart/1000, -4, xend=chromEnd/1000, yend=-4),
                         data=alg2res[alg2res$segtype=="N",], color="purple", size=3)
}
print(p2)

# merged plot
toplot1 = mutate(alg1res, segtype="S", treatment=2)
toplot2 = mutate(alg2res, treatment=2)
toplotL = data.frame(label=c("Standard", "Proposed"), x=120725, y=c(-3, -7), treatment=2)
p3 = bind_rows(bedCds,bedHds, .id="treatment") %>%
  ggplot() +
  geom_point(aes(pos/1000, m), alpha=0.5) +
  geom_segment(aes(chromStart/1000, -3, xend=chromEnd/1000, yend=-3, color=segtype),
               data=toplot1, size=7)+
  geom_segment(aes(chromStart/1000, -7, xend=chromEnd/1000, yend=-7, color=segtype),
               data=toplot2, size=7) +
  geom_text(aes(x=x, y=y, label=label), data=toplotL) + 
  facet_grid(treatment~., labeller = labeller(treatment=c("1"="control", "2"="H3K27ac"))) + 
  xlab("position on chromosome, kbp") + ylab("coverage") +
  xlim(c(120100, 120750)) + ylim(c(-8, 27)) +
  scale_color_manual(values=c("S"="blue", "N"="gold")) +
  theme_bw()
print(p3)  

# save output
outprefix = "../drafts/changepoint-method/results-appl/broad-chr1-120mbp-ds500-"
write.table(alg1res, paste0(outprefix, "alg1.tsv"), quote=F, row.names=F, sep="\t")
save(alg1, file=paste0(outprefix, "alg1.RData"))
ggsave(paste0(outprefix, "alg1.png"), plot=p1, width=14, height=8, units="cm")
write.table(alg2res, paste0(outprefix, "alg2.tsv"), quote=F, row.names=F, sep="\t")
save(alg2, file=paste0(outprefix, "alg2.RData"))
ggsave(paste0(outprefix, "alg2.png"), plot=p2, width=14, height=8, units="cm")
ggsave(paste0(outprefix, "merged.png"), plot=p3, width=18, height=15, units="cm")
