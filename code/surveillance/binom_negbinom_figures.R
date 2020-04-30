# All for M15+
# For high or low (best/worst case scenario)
# primaquine vs. no primaquine
# 2020 and 2025
# Evaluate probability of detecting at least one case across values of n 

library(ggplot2)

dat <- read.csv("./aggregated_results.csv")

ns <- seq(100,10000, by=100)

long.dat.b <- dat[rep(1:nrow(dat), each=length(ns)),]

long.dat.b$ns <- rep(ns, times=nrow(long.dat.b)/length(ns))

long.dat.b$probs <- 1-dbinom(0,size=long.dat.b$ns, prob = long.dat.b$Prevalence.in.M15.)

long.dat.b$Baseline.incidence <- factor(long.dat.b$Baseline.incidence, levels = c("low","high"), ordered = T)



# All variables
ggplot(long.dat.b, aes(x=ns, y=probs, colour=Baseline.incidence)) + 
  scale_y_continuous("Probability of detection") +
  scale_x_continuous("Number of samples") +
  geom_point() + geom_line() + scale_colour_discrete() + facet_grid(Province~Year+Intervention.scenario)

ggsave(paste0("./prob_detect1_allprovinces.pdf"), width = 11.69, height = 8.27, units = "in")

single.prov <- "Pailin"
ggplot(long.dat.b[long.dat.b$Province==single.prov,], aes(x=ns, y=probs, colour=Baseline.incidence)) + 
  scale_y_continuous("Probability of detection") +
  scale_x_continuous("Number of samples") +
  geom_point(size=1.5) + geom_line(size=1) + scale_colour_discrete() + facet_grid(Intervention.scenario~Year)

ggsave(paste0("./prob_detect1_",single.prov,".pdf"), width = 11.69, height = 8.27, units = "in")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Part 2: What N is needed to have >0.95 or >0.99 probability of detecting at least 1 case

# dat <- read.csv("~/Dropbox/Consulting/R. Hickson/aggregated_results.csv")


ps <- c(0.95, 0.99)

long.dat.nb <- dat[rep(1:nrow(dat), each=length(ps)),]

long.dat.nb$ps <- rep(ps, times=nrow(long.dat.nb)/length(ps))


long.dat.nb$min.Ns <- qgeom(long.dat.nb$ps, prob = long.dat.nb$Prevalence.in.M15.)

long.dat.nb$xlabs <- paste(long.dat.nb$Baseline.incidence, long.dat.nb$Intervention.scenario, sep="\n")

long.dat.nb$Baseline.incidence <- factor(long.dat.nb$Baseline.incidence, levels = c("low","high"), ordered = T)
long.dat.nb$xlabs <- factor(long.dat.nb$xlabs, 
                            levels = c("low\nno_prim","low\nprim","high\nno_prim","high\nprim"),
                            labels = c("Low Incidence\nNo Primaquine","Low Incidence\nPrimaquine",
                                       "High Incidence\nNo Primaquine","High Incidence\nPrimaquine"),
                            ordered = T)

long.dat.nb$Province <- factor(long.dat.nb$Province, levels = levels(long.dat.nb$Province), 
                               labels = gsub(pattern = "_", replacement = "\n", x = levels(long.dat.nb$Province)))

# ggplot(long.dat.nb[long.dat.nb$ps==0.95,], aes(x=xlabs, y=min.Ns, fill=Baseline.incidence)) + 
#   geom_bar(stat="identity") + 
#   scale_y_continuous("Target Sample Size", trans="log10")+
#   scale_x_discrete("Scenario") +
#   facet_grid(Province ~ Year) +
#   geom_text(aes(label=min.Ns, x=xlabs, y=5)) +
#   scale_fill_discrete() +
#   theme(legend.position = "none", axis.text.x = element_text(angle=90), text=element_text(size=12))
# 
# ggsave("./SampleSize_95pc_detectcase.pdf", width = 11.69, height = 8.27, units = "in")



ggplot(long.dat.nb[long.dat.nb$ps==0.99,], aes(x=xlabs, y=min.Ns, fill=Baseline.incidence)) + 
  geom_bar(stat="identity") + 
  scale_y_continuous("Target Sample Size", trans="log10")+
  scale_x_discrete("Scenario") +
  facet_grid(Province ~ Year) +
  geom_text(aes(label=min.Ns, x=xlabs, y=5)) +
  scale_fill_discrete() +
  theme(legend.position = "none", axis.text.x = element_text(angle=90), text=element_text(size=12))

ggsave("./SampleSize_99pc_detectcase.pdf", width = 11.69, height = 8.27, units = "in")



