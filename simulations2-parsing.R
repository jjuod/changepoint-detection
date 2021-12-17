#### ----- CLEAN SIMULATIONS 2 PARSING -------
options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(ggplot2)
setwd("~/Documents/gitrep/changepoint-detection/")
source("MAIN-algorithm.R")

## RUN THE SIMULATIONS:
# Recommend running that outside Rstudio, as older versions of it sometimes
# don't manage to do garbage collection efficiently
#system("Rscript simulations2-separate.R")

# Set the simulation settings
ntotest2 = c(30, 60, 100, 150, 220)
niter2 = 1000

# Create a dataframe of true chp positions
trueposdf = data.frame()
for(n in ntotest2){
    truepos = floor(c(0.3*n+1, 0.5*n))
    df1 = data.frame(scen=1, n=rep(n, each=length(truepos)), truepos=truepos)
    # same positions but var changes:
    df4 = data.frame(scen=4, n=rep(n, each=length(truepos)), truepos=truepos)
    
    truepos = floor(c(0.5*n+1, 0.6*n, 0.7*n+1, 0.8*n))
    df2 = data.frame(scen=2, n=rep(n, each=length(truepos)), truepos=truepos)
    
    l1 = floor(n*0.1)
    l2 = floor(n*0.05)
    l3 = ceiling(n*0.05)
    # Note that the starts and ends are interleaved to maintain pairs
    truepos = c(l1+1 + 0:8*(l2+l3),  # starts
                l1 + l2 + 0:8*(l2+l3))  # ends
    truepos = sort(truepos)  # interleave starts-ends
    df3 = data.frame(scen=3, n=rep(n, each=length(truepos)), truepos=truepos)
    
    trueposdf = bind_rows(trueposdf, rbind(df1, df2, df3, df4))
}

# Reports the estimation error for each reported chp, checking for start/end match
get_dist_errors_est = function(pos1, pos2, scen, n, tpdf, segtype){
    d = rep(Inf, length(pos1)*2)
    if (segtype[1]=="nuis"){ return(d) }
    truepos = tpdf$truepos[tpdf$scen==scen[1] & tpdf$n==n[1]]
    trueposS = truepos[seq(1, length(truepos), 2)]
    trueposE = truepos[seq(2, length(truepos), 2)]
    
    # for each true segment, calculate distance to chp, in points:
    for(pairix in seq(length(pos1))){
        chp1 = pos1[pairix]
        chp2 = pos2[pairix]
        # pick the closest true chp for each estimated one
        d1 = min(abs(trueposS-chp1))
        d2 = min(abs(trueposE-chp2))
        d[2*pairix-1] = d1
        d[2*pairix] = d2
    }
    # d<thr <=> there was a TP for this est chp within thr
    return(d)
}
# we'll also need a per-true point detection function:
get_dist_errors2_ordered = function(pos1, pos2, scen, n, tpdf){
    truepos = tpdf$truepos[tpdf$scen==scen[1] & tpdf$n==n[1]]
    
    # for each true segment, calculate distance to chp, in points:
    d = rep(Inf, length(truepos))
    for(pairix in seq(1, length(truepos), 2)){
        chps = truepos[pairix]
        chpe = truepos[pairix+1]
        # pick the closest estimated chp for start and end columns
        d1 = min(abs(pos1-chps))
        d2 = min(abs(pos2-chpe))
        d[pairix] = d1
        d[pairix+1] = d2
    }
    # d<thr <=> there was a TP for this true chp within thr
    return(d)
}


##  READ OUTPUT
allsegs2 = read.table("../drafts/changepoint-method/results-sim-new/sim2-detections.tsv", h=T)
allsegs2_old = read.table("../drafts/changepoint-method/results-sim/sim2-detections.tsv", h=T)
allsegs2_old$i = allsegs2_old$i + 500
# Old NOT positions were stored incorrectly:
allsegs2_old$V1[allsegs2_old$alg=="NOT"] = allsegs2_old$V1[allsegs2_old$alg=="NOT"] +1
allsegs2_old$V2[allsegs2_old$alg=="NOT"] = allsegs2_old$V2[allsegs2_old$alg=="NOT"] +1
allsegs2 = bind_rows(allsegs2, allsegs2_old)

# NOT needs post-processing as it does not distinguish background and segments
# so we arbitrarily define any segs w/ mean within 1 SD of theta as background:
allsegs2_allscens = filter(allsegs2, alg!="NOT" | abs(V3)>1)

allsegs2_scen4 = filter(allsegs2_allscens, scen==4)
allsegs2 = filter(allsegs2_allscens, scen %in% 1:3)


## ----- SUPPL TABLE: pruning makes no difference ----- 
dfF = filter(allsegs2, alg=="full")
dfP = filter(allsegs2, alg=="pruned")
diffis = bind_rows(anti_join(dfF, dfP, by=c("NPOINTS", "i", "segtype", "scen", "V1", "V2")),
                   anti_join(dfP, dfF, by=c("NPOINTS", "i", "segtype", "scen", "V1", "V2")))

# all iterations that differed
# (skip scenario 3 as that one has more differences)
semi_join(allsegs2, diffis, by=c("NPOINTS", "scen", "i")) %>%
    filter(alg %in% c("full", "pruned"), scen %in% 1:2) %>%
    mutate(segtype=ifelse(segtype=="nuis", "N", "S")) %>%
    .[,c("alg", "scen", "NPOINTS", "i", "segtype", "V1", "V2")]

# how frequently something differed?
diffis_u = unique(diffis[,c("NPOINTS", "scen", "i")])
table(diffis_u$scen)
nrow(diffis_u) / niter2 / length(ntotest2) / 3  # 3 scenarios


## -----------------------------------------------------
# TPRs (by true changepoints)
# *** (Supplement Table) ***
disterrors_bytrue = allsegs2 %>%
    filter(segtype=="seg") %>%
    group_by(alg, scen, NPOINTS, i) %>%
    summarize(disterr=get_dist_errors2_ordered(V1, V2, scen, NPOINTS, trueposdf)) %>%
    mutate(disterr=disterr/NPOINTS) %>%
    summarize(tps = max(disterr)<0.05) %>%
    summarize(tps = sum(tps)/niter2)  # per penalty level
st3 = disterrors_bytrue %>%
    spread(key="alg", value="tps") %>%
    arrange(scen)
st3[,c("scen", "NPOINTS", "pruned", "anomaly", "aPELT", "sparse", "NOT")]


# localization errors for each detected segment:
disterrors_byest = allsegs2 %>%
    group_by(alg, segtype, scen, NPOINTS, i) %>%
    summarize(disterr=get_dist_errors_est(V1, V2, scen, NPOINTS, trueposdf, segtype),
              nsegs=n())
# convert distance errors to fractions of full data
disterrors_byest = mutate(disterrors_byest, disterr=disterr/NPOINTS)
# actually summarize per iteration:
worst_errs_byest = group_by(disterrors_byest, NPOINTS, alg, segtype, scen, i) %>%
    summarize(hits = sum(disterr<0.05), nsegs = max(nsegs))

# MEAN ABS ERROR of detected SIGNAL chps
# *** (Supplement Table) ***
mae_table = filter(disterrors_byest, disterr<0.05, segtype=="seg") %>%
    group_by(NPOINTS, alg, scen) %>%
    summarize(mae=mean(disterr*NPOINTS)) %>%
    spread(key="alg", value="mae") %>%
    arrange(scen)
mae_table[,c("scen", "NPOINTS", "pruned", "anomaly", "aPELT", "sparse", "NOT")]


# precision / PPV / 1-FDR (signal segments only)
# and also nsegs
# *** (Main Table) *** 
ppv_table = worst_errs_byest %>%
    group_by(NPOINTS, alg, segtype, scen) %>%
    summarize(totalhits=sum(hits), totaln=sum(nsegs)*2, ppv=totalhits/totaln)
t2 = filter(ppv_table, segtype=="seg", alg!="full")[,c(1,2,4,7)] %>%
    spread(key="alg", value="ppv") %>%
    arrange(scen)
t2[,c("scen", "NPOINTS", "pruned", "anomaly", "aPELT", "sparse", "NOT")]


# Segment number bias
# *** (Main Figure) ***
# attach the expected number of true segments
# (NOTE: it is the number of segments, while other columns
# are numbers of chps, so multiply by 2 as needed)
distrnsegs = ppv_table %>%
    filter(alg!="full", segtype=="seg" | alg=="pruned") %>%  # no nuisance segs in standard chp algs
    mutate(ntrue = ifelse(scen==3, ifelse(segtype=="seg", 9, 0),
                   ifelse(scen==2, ifelse(segtype=="seg", 2, 1),
                          1))) %>%
    mutate(meannseg = totaln/2/niter2, bias = meannseg-ntrue) %>%
    ungroup

distrnsegs %>%
    mutate(algseg=factor(ifelse(segtype=="seg", alg, "nuisance"),
                levels=c("anomaly", "aPELT", "sparse", "NOT", "pruned", "nuisance"),
                labels=c("anomaly", "aPELT", "sparse", "not", "proposed", "proposed (nuisance)"))) %>%
    ggplot(aes(x=NPOINTS)) +
    geom_hline(aes(yintercept=0), col="grey30") +
    geom_line(aes(y=bias, col=algseg, lty=algseg), lwd=1) +
    geom_point(aes(y=bias, shape=algseg, col=algseg, alpha=algseg), size=3) +
    scale_alpha_manual(name="detections:", values=c("anomaly"=0, "aPELT"=1, "sparse"=1,
                                        "not"=1, "proposed"=0, "proposed (nuisance)"=0)) +
    scale_shape_manual(name="detections:", values=c("anomaly"=17, "aPELT"=4, "sparse"=15,
                                        "not"=20, "proposed"=1, "proposed (nuisance)"=1)) +
    scale_linetype_manual(values=c("anomaly"=1, "aPELT"=1, "sparse"=1, "not"=1, "proposed"=1,
                                   "proposed (nuisance)"=3), name="detections:") +
    scale_color_manual(name="detections:", values =c("dodgerblue2", "blue4", "mediumturquoise",
                                                     "skyblue1", "coral1", "coral1")) +
    facet_wrap(~scen, labeller=labeller(scen=c("1"="scenario 1", "2"="scenario 2", "3"="scenario 3")),
               scales="free_y") +
    theme_bw() + xlab("n") + ylab(expression(E(~hat(k)-k))) +
    theme(legend.position = "bottom", text=element_text(size=14),
          legend.text=element_text(size=13), plot.margin = unit(c(0.1,0.3,0,0.1), "cm"))
ggsave("../drafts/changepoint-method/v21dec/fig2.eps", width=17, height=10, units="cm")

#  ---------------
# Compare theta_i for the first true segment:
# (using the signal segment w/ most overlap)
# *** (Main Text) ***
best1seg = filter(allsegs2, scen==1) %>%
    group_by(alg, segtype, NPOINTS, i) %>%
    mutate(firstend = pmin(V2, floor(0.5*NPOINTS)),
           laststart = pmax(V1, floor(0.3*NPOINTS+1))) %>%
    mutate(overl = pmax(firstend-laststart, 0)) %>%
    top_n(1, rank(overl, ties.method = "r"))
extract_theta = function(V1, V2, V3, segtype){
    if("nuis" %in% segtype & "seg" %in% segtype){
        ni = segtype=="nuis"
        si = segtype=="seg"
        # if segment overlaps a nuisance:
        if(V1[ni]<=V1[si] & V2[ni]>=V2[si]){
            thetaN = V3[segtype=="nuis"]
            thetaS = V3[segtype=="seg"] - thetaN  
        } else{
            thetaS = V3[si] 
        }
    } else if ("seg" %in% segtype){
        thetaS = V3[segtype=="seg"]
    } else {
        thetaS = NA
    }
    return(thetaS)
}
thetadf = group_by(best1seg, alg, NPOINTS, i) %>%
    summarize(thetaS = extract_theta(V1, V2, V3, segtype)) %>%
    mutate(thetaS = ifelse(alg=="anomaly", sqrt(thetaS), thetaS)) %>%  # anomaly reports squared estimate
    summarize(m=mean(thetaS, na.rm=T), s=sd(thetaS, na.rm=T), n=sum(!is.na(thetaS)))
thetadf %>%
    print.data.frame
filter(thetadf, NPOINTS==220)


# TODO finish

# SCENARIO 4 sensitivity analysis results
# *** (Supplement Table) ***
niter2 = 500
disterrors_byest = allsegs2_scen4 %>%
    group_by(alg, segtype, scen, NPOINTS, i) %>%
    summarize(disterr=get_dist_errors_est(V1, V2, scen, NPOINTS, trueposdf, segtype),
              nsegs=n())
# convert distance errors to fractions of full data
disterrors_byest = mutate(disterrors_byest, disterr=disterr/NPOINTS)
# actually summarize per iteration:
worst_errs_byest = group_by(disterrors_byest, NPOINTS, alg, segtype, i) %>%
    summarize(hits = sum(disterr<0.05), nsegs = max(nsegs))

mae_table = filter(disterrors_byest, disterr<0.05, segtype=="seg") %>%
    group_by(NPOINTS, alg) %>%
    summarize(mae=mean(disterr*NPOINTS)) %>%
    spread(key="alg", value="mae")
mae_table[,c("NPOINTS", "pruned", "anomaly", "aPELT", "sparse", "NOT")]

ppv_table = worst_errs_byest %>%
    group_by(NPOINTS, alg, segtype) %>%
    summarize(totalhits=sum(hits), totaln=sum(nsegs)*2)
# attach the expected number of true segments and calculate bias, tpr
distrnsegs = ppv_table %>%
    filter(alg!="full", segtype=="seg" | alg=="pruned") %>%
    mutate(ntrue = 1) %>%
    mutate(meannseg=totaln/2/niter2, bias=meannseg-ntrue) %>%
    ungroup
distrnsegs %>%
    select(one_of(c("NPOINTS", "alg", "segtype", "meannseg"))) %>%
    spread(key="alg", value="meannseg")

disterrors_bytrue = allsegs2_scen4 %>%
    filter(segtype=="seg") %>%
    group_by(alg, NPOINTS, i) %>%
    summarize(disterr=get_dist_errors2_ordered(V1, V2, 1, NPOINTS, trueposdf)) %>%
    mutate(disterr=disterr/NPOINTS) %>%
    summarize(tps = max(disterr)<0.05) %>%
    summarize(tps = sum(tps)/niter2)  # per penalty level
disterrors_bytrue %>%
    spread(key="alg", value="tps")

