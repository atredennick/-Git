#Code to replicate results from Tredennick et al. 2013 PLoS One

#This code, and any updates on it, can also be found on GitHub: 
#http://github.com/atredennick/-Git/tree/master/Savanna_Allometry

### NOTICE!! ### NOTICE!! ### NOTICE!! ### NOTICE!! ### NOTICE!! ### NOTICE!! ### NOTICE!! ### NOTICE!! ####
#                                                                                                          #
# To run this code you must first install JAGS on your computer.                                           #
# Binary packages for Windows and OS X can be found here: http://sourceforge.net/projects/mcmc-jags/files/ # 
# Information on JAGS can be found here: http://mcmc-jags.sourceforge.net/                                 #
# There is a line of code that will automatically send you to the link for JAGS source files                #
# Comment out (#) the "browseURL" line after JAGS has been installed                                       #
#                                                                                                          #
############################################################################################################

#Code written by A. Tredennick (atredenn@gmail.com)

rm(list=ls())

#First, install JAGS from this link, comment out (#) after install
# browseURL("http://sourceforge.net/projects/mcmc-jags/files/", browser = getOption("browser"), encodeIfNeeded = FALSE)

#Find the JAGS HB model file path on your computer
HB.model <- file.choose()

#Install necessary packages // comment in/out if you need these
# install.packages('rjags')
# install.packates('coda')
# install.packages('rdryad')
# install.packages('ggplot2')
# install.packages('scales')
# install.packages('gridExtra')

#Load packages
library(rjags)
library(coda)
library(rdryad)
library(ggplot2)
library(scales)
library(gridExtra)

##Bring in data
##Data is available on Dryad by searching by DOI
data.file = file.choose()
data=read.csv(data.file)

names(data) = tolower(names(data))

##Set up matrix to store LogDiameter, LogLength, LogMass, LogWoodMass and LogLeafMass. 
Y = matrix(nrow=length(data[,1]), ncol=4)

#Log data and put in dependent variable data matrix
Y[,1] = log10(data$length)
Y[,2] = log10(data$total_mass)
Y[,3] = log10(data$aggregated_wood_wt)
Y[,4] = log10(data$aggregated_leaf_wt)

LogDiameter = log10(data$diameter)

##Set up positive definite matrix for derivition of Omega precision matrix (R)
R = structure(.Data=c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1), .Dim=c(4,4))

## Assign "N's" for indexing	
N=length(data$species)
Nsp = 3 #3 species
Ntree=c(10,5,10) #10 trees of demi, 5 of coge, and 10 of cogl

#Assign a prediction vector for posterior predictions
x.pred<-log10(seq(2,20,by=0.1)) #prediction vector
y.pred <- matrix(nrow=length(x.pred), ncol=4) #holder for posterior predictions

#Create matrices to hold derived quantities of interest
beta.spp.diff <- matrix(nrow=3,ncol=4) #stores species differences from HB model
alpha.spp.diff <- matrix(nrow=3,ncol=4) #stores species differences from HB model


##Set-up for call to JAGS
data.J = list(R=R, N=N, Nsp=Nsp, Y=Y, SP=data$species, tree=data$tree, 
              Ntree=Ntree, LogDiameter=LogDiameter, x.pred=x.pred, 
              beta.spp.diff=beta.spp.diff, alpha.spp.diff=alpha.spp.diff,
              y.pred=y.pred)

inits=list(
	list(Omega = R),
	list(Omega = R),
	list(Omega = R)
)

#The 'vars' variable holds a list of possible nodes that can be tracked by JAGS during the MCMC
#Pick and choose which ones you want output for.
#The first assignment for vars (which is commented out) contains the entire list
#The second vars assignment (which can be added to) just as the posterior Y predictions tracked to create Figure 1

vars <- c("mu.alpha", "mu.beta", "alpha.species", "beta.species", 
          "Sigma", "p.mu", "p.fit",  "spp.var", "ind.var", 
          "betalength.coge", "betamass.coge", "betawoodm.coge", "betaleafm.coge", 
          "betalength.demi", "betamass.demi", "betawoodm.demi", "betaleafm.demi", 
          "betalength.cogl", "betamass.cogl", "betawoodm.cogl", "betaleafm.cogl", 
          "beta.spp.diff", "alpha.spp.diff", "y.pred")

# vars <- c("y.pred")

#For production runs: 50,000 adapt; 200,000 burn-in; 1,000,000 samples; thin chains by 10
#Number of chains set by length of 'inits' list
#Change these to lower values (specifically sample.n) for test runs
adapt.n <- 50000
burn.n <- 200000
samples.n <- 1000000
thin.n <- 10

#Call JAGS HB model
jm <- jags.model(HB.model, data=data.J, inits, n.chains=length(inits), n.adapt=adapt.n)

#Burn in
update(jm, n.iter=burn.n)

#Run MCMC and track nodes as listed in 'vars'
zm <- coda.samples(jm, variable.names=vars, n.iter=samples.n, thin=thin.n)

#summarize MCMC list as statistics and quantiles
zm.stats <- summary(zm)$stat
zm.quants <- summary(zm)$quantile

#Create quantile and statistical matrix from posterior predictions for Figure 1
y.quants <- zm.quants[grep("y.pred", rownames(zm.quants)),]
y.stats <- zm.stats[grep("y.pred", rownames(zm.stats)),]

#Check convergence with Heidel diagnostic
# h <- heidel.diag(zm)
# print, h


##########################################################
######################## FIGURE 1 ########################
##########################################################
#Create Figure 1: fitted 'global', interspecific allometries

#set-up dataframes with predicted data and CIs
#back-transform data from Log10 for dataframes, can be plotted in log or not in ggplot
x.pred<-seq(2,20,by=0.1)

df.length <- data.frame(x=x.pred,
                        y=as.numeric(10^y.stats[1:181,1]),
                        y.low=as.numeric(10^y.quants[1:181,1]),
                        y.hi=as.numeric(10^y.quants[1:181,5]))
df.agmass <- data.frame(x=x.pred,
                        y=as.numeric(10^y.stats[182:362,1]),
                        y.low=as.numeric(10^y.quants[182:362,1]),
                        y.hi=as.numeric(10^y.quants[182:362,5]))
df.smass <- data.frame(x=x.pred,
                       y=as.numeric(10^y.stats[363:543,1]),
                       y.low=as.numeric(10^y.quants[363:543,1]),
                       y.hi=as.numeric(10^y.quants[363:543,5]))
df.lmass <- data.frame(x=x.pred,
                       y=as.numeric(10^y.stats[544:724,1]),
                       y.low=as.numeric(10^y.quants[544:724,1]),
                       y.hi=as.numeric(10^y.quants[544:724,5]))

#set transparency value
alpha = 0.5


#legend plot
legend.plot <- ggplot() + 
  geom_point(data=data, aes(x=data$diameter, y=data$length, shape=spp, color=spp), size=3) +
  labs(shape="") +
  scale_colour_manual(name = "",
                      labels = c("coge", "cogl", "demi"),
                      values = c("red", "orange", "blue")) +   
  scale_shape_manual(name = "",
                     labels = c("coge", "cogl", "demi"),
                     values = c(15, 17, 19)) +
  theme_bw() +
  opts(axis.title.x = theme_text(size=20),
       axis.title.y = theme_text(size=20, angle=90), 
       axis.text.x = theme_text(size=14), 
       axis.text.y = theme_text(size=14), 
       panel.grid.major = theme_blank(),
       panel.grid.minor = theme_blank(),
       legend.position = "right",
       legend.key = theme_blank()   
  )

#start making plots
length = ggplot() + 
  geom_ribbon(data=df.length, aes(x=x, ymin=y.low, ymax=y.hi), alpha=alpha) +
  geom_point(data=data, aes(x=data$diameter, y=data$length, shape=spp, color=spp), alpha=alpha, size=3) +
  geom_line(data=df.length, aes(x=x, y=y), size=1.2) +
  xlab("Diameter (cm)") + ylab("Length (cm)") + 
  theme_bw() +
  scale_fill_hue(l=40) +
  #   scale_y_continuous(limits = c(0,5000)) +
  #   scale_x_continuous(limits = c(0,20)) +
  labs(shape="") +
  geom_text(aes(2,20000,label = "A")) +
  scale_colour_manual(name = "",
                      labels = c("coge", "cogl", "demi"),
                      values = c("red", "orange", "blue")) +   
  scale_shape_manual(name = "",
                     labels = c("coge", "cogl", "demi"),
                     values = c(15, 17, 19)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  opts(axis.title.x = theme_text(size=16),
       axis.title.y = theme_text(size=16, angle=90), 
       axis.text.x = theme_text(size=12), 
       axis.text.y = theme_text(size=12), 
       panel.grid.major = theme_blank(),
       panel.grid.minor = theme_blank(),
       legend.position = "right",
       legend.key = theme_blank()   
  )

agmass = ggplot() + 
  geom_ribbon(data=df.agmass, aes(x=x, ymin=y.low/1000, ymax=y.hi/1000), alpha=alpha) +
  geom_point(data=data, aes(x=data$diameter, y=data$total_mass/1000, shape=spp, color=spp), alpha=alpha, size=3) +
  geom_line(data=df.agmass, aes(x=x, y=y/1000), size=1.2) +
  xlab("Diameter (cm)") + ylab("AG Mass (kg)") + 
  theme_bw() +
  scale_fill_hue(l=40) +
  labs(shape="") +
  geom_text(aes(2,8000,label = "B")) +
  scale_colour_manual(name = "",
                      labels = c("coge", "cogl", "demi"),
                      values = c("red", "orange", "blue")) +   
  scale_shape_manual(name = "",
                     labels = c("coge", "cogl", "demi"),
                     values = c(15, 17, 19)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  opts(axis.title.x = theme_text(size=16),
       axis.title.y = theme_text(size=16, angle=90), 
       axis.text.x = theme_text(size=12), 
       axis.text.y = theme_text(size=12),  
       panel.grid.major = theme_blank(),
       panel.grid.minor = theme_blank(),
       legend.position = "right",
       legend.key = theme_blank()   
  )


smass = ggplot() + 
  geom_ribbon(data=df.smass, aes(x=x, ymin=y.low/1000, ymax=y.hi/1000), alpha=alpha) +
  geom_point(data=data, aes(x=data$diameter, y=data$aggregated_wood_w/1000, shape=spp, color=spp), alpha=alpha, size=3) +
  geom_line(data=df.smass, aes(x=x, y=y/1000), size=1.2) +
  xlab("Diameter (cm)") + ylab("Stem Mass (kg)") + 
  theme_bw() +
  scale_fill_hue(l=40) +
  #   scale_y_continuous(limits = c(0,5000)) +
  #   scale_x_continuous(limits = c(0,20)) +
  labs(shape="") +
  geom_text(aes(2,8000,label = "C")) +
  scale_colour_manual(name = "",
                      labels = c("coge", "cogl", "demi"),
                      values = c("red", "orange", "blue")) +   
  scale_shape_manual(name = "",
                     labels = c("coge", "cogl", "demi"),
                     values = c(15, 17, 19)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  opts(axis.title.x = theme_text(size=16),
       axis.title.y = theme_text(size=16, angle=90), 
       axis.text.x = theme_text(size=12), 
       axis.text.y = theme_text(size=12),  
       panel.grid.major = theme_blank(),
       panel.grid.minor = theme_blank(),
       legend.position = "right",
       legend.key = theme_blank()   
  )



lmass = ggplot() + 
  geom_ribbon(data=df.lmass, aes(x=x, ymin=y.low/1000, ymax=y.hi/1000), alpha=alpha) +
  geom_point(data=data, aes(x=data$diameter, y=data$aggregated_leaf_wt/1000, shape=spp, color=spp), alpha=alpha, size=3) +
  geom_line(data=df.lmass, aes(x=x, y=y/1000), size=1.2) +
  xlab("Diameter (cm)") + ylab("Leaf Mass (kg)") + 
  theme_bw() +
  scale_fill_hue(l=40) +
  #   scale_y_continuous(limits = c(0,5000)) +
  #   scale_x_continuous(limits = c(0,20)) +
  labs(shape="") +
  geom_text(aes(2,8000,label = "D")) +
  scale_colour_manual(name = "",
                      labels = c("coge", "cogl", "demi"),
                      values = c("red", "orange", "blue")) +   
  scale_shape_manual(name = "",
                     labels = c("coge", "cogl", "demi"),
                     values = c(15, 17, 19)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  opts(axis.title.x = theme_text(size=16),
       axis.title.y = theme_text(size=16, angle=90), 
       axis.text.x = theme_text(size=12), 
       axis.text.y = theme_text(size=12), 
       panel.grid.major = theme_blank(),
       panel.grid.minor = theme_blank(),
       legend.position = "right",
       legend.key = theme_blank()   
  )


g1 <- ggplotGrob(legend.plot)
legend <- editGrob(getGrob(g1, gPath("guide-box"), grep=TRUE), vp=viewport())

# Arrange the four charts
grid.arrange(arrangeGrob(length + opts(legend.position = "none"), 
                         agmass + opts(legend.position = "none"),
                         smass + opts(legend.position = "none"),
                         lmass + opts(legend.position = "none"),
                         nrow=2),
             legend, 
             widths=unit.c(unit(1, "npc") - legend$width, legend$width), 
             nrow=1)



##########################################################
######################## FIGURE 2 ########################
##########################################################
#Create Figure 2: scaling exponents at all HB levels for each trait

#Get line indexes for sub-setting posterior samples
l.coge.in = grep("betalength.coge", rownames(zm.quants))
m.coge.in = grep("betamass.coge", rownames(zm.quants))
w.coge.in = grep("betawoodm.coge", rownames(zm.quants))
s.coge.in = grep("betaleafm.coge", rownames(zm.quants))

l.demi.in = grep("betalength.demi", rownames(zm.quants))
m.demi.in = grep("betamass.demi", rownames(zm.quants))
w.demi.in = grep("betawoodm.demi", rownames(zm.quants))
s.demi.in = grep("betaleafm.demi", rownames(zm.quants))

l.cogl.in = grep("betalength.cogl", rownames(zm.quants))
m.cogl.in = grep("betamass.cogl", rownames(zm.quants))
w.cogl.in = grep("betawoodm.cogl", rownames(zm.quants))
s.cogl.in = grep("betaleafm.cogl", rownames(zm.quants))

l.sp = grep("beta.species", rownames(zm.quants))

lbeta=grep("mu.beta", rownames(zm.quants))

#Create x-values for spacing in plots
x=seq(1,38,1)
x2=seq(1,10,1)
x1=seq(12,16,1)
x3=seq(18,27,1)
x4=seq(30,32,1)
x5=37
lcol="grey55"
col=c("blue", "red", "orange")


#Set 2x2 plot window
par(mfrow=c(2,2), mar=c(3.1,4.1,4.1,2.1), oma=c(0,0,0,0), 
    mgp = c(1.8, 0.2, 0), tck=0.02, mar=c(1,3,1,1))

#Length
plot(x, seq(0.2,1.6,(1.4/(max(x)-1))), cex=0.000000001, pch=21, 
     bg="white", col="white", axes=F, xlab="", ylab="Length scaling exponent", 
     cex.lab=1.4)
abline(h=0.67,lwd=1.5) #MST prediction
abline(h=1.00, lty="dashed", lwd=1.5) #GEOM prediction
abline(h=0.5, lty="dotted", lwd=1.5) #STRESS prediction
for(i in 1:5) lines(rep(x1[i],2), zm.quants[l.coge.in[i],c(1,5)], col=lcol)
points(x1, zm.stats[l.coge.in,1], pch=22, bg="red", cex=1)
for(i in 1:10) lines(rep(x2[i],2), zm.quants[l.demi.in[i],c(1,5)], col=lcol)
points(x2, zm.stats[l.demi.in,1], pch=21, bg="blue", cex=1)
for(i in 1:10) lines(rep(x3[i],2), zm.quants[l.cogl.in[i],c(1,5)], col=lcol)
points(x3, zm.stats[l.cogl.in,1], pch=24, bg="orange", cex=1)
for (i in 1:3) lines(rep(x4[i],2), zm.quants[l.sp[i], c(1,5)], col=lcol)
points(x4, zm.stats[l.sp[1:3],1], pch=c(21,22,24), bg=col, cex=1.3)
lines(rep(x5,2), zm.quants[lbeta[1],c(1,5)], col=lcol)
points(x5, zm.stat[lbeta[1],1], pch=23, cex=2, col="black", bg="grey")
box()
abline(v=c(28.5, 34), lty=4, col="grey")
axis(2,at=c(seq(0,4,0.2)), las=1)
legend(3,1.6,legend=c("demi", "coge", "cogl", "'global'"), pch=c(21,22,24,23), pt.bg=c("blue", "red", "orange", "grey"), cex=c(1,1,1,1), bty="n")
legend(12,1.6,c("MST","GEOM","STRESS"), lty=c("solid","dashed","dotted"), cex=1,lwd=1.3, bty="n")
mtext("A", side=2, las=1, line=-1.4, at=1.55, cex=1)

#Aboveground mass
plot(x, seq(2.2,3.2,(1/(max(x)-1))), cex=0.000000001, pch=21, bg="white", col="white", axes=F, xlab="", ylab="AG mass scaling exponent", cex.lab=1.4)
abline(h=2.67,lwd=1.5) #MST prediction
abline(h=3.00, lty="dashed", lwd=1.5) #GEOM prediction
abline(h=2.5, lty="dotted", lwd=1.5) #STRESS prediciton
for(i in 1:5) lines(rep(x1[i],2), zm.quants[m.coge.in[i],c(1,5)], col=lcol)
points(x1, zm.stats[m.coge.in,1], pch=22, bg="red", cex=1)
for(i in 1:10) lines(rep(x2[i],2), zm.quants[m.demi.in[i],c(1,5)], col=lcol)
points(x2, zm.stats[m.demi.in,1], pch=21, bg="blue", cex=1)
for(i in 1:10) lines(rep(x3[i],2), zm.quants[m.cogl.in[i],c(1,5)], col=lcol)
points(x3, zm.stats[m.cogl.in,1], pch=24, bg="orange", cex=1)
for (i in 1:3) lines(rep(x4[i],2), zm.quants[l.sp[(i+3)], c(1,5)], col=lcol)
points(x4, zm.stats[l.sp[4:6],1], pch=c(21,22,24), bg=c("blue", "red", "orange"), cex=1.3)
lines(rep(x5,2), zm.quants[lbeta[2],c(1,5)], col=lcol)
points(x5, zm.stats[lbeta[2],1], pch=23, cex=2, col="black", bg="grey")
box()
abline(v=c(28.5, 34), lty=4, col="grey")
axis(2,at=c(seq(0,5,0.2)), las=1)
mtext("B", side=2, las=1, line=-1.4, at=3.15, cex=1)

#Stem mass
plot(x, seq(2.2,3.4,(1.2/(max(x)-1))), cex=0.000000001, pch=21, bg="white", col="white", axes=F, xlab="", ylab="Stem mass scaling exponent", cex.lab=1.4)
abline(h=2.67,lwd=1.5) #MST prediction
abline(h=3.00, lty="dashed", lwd=1.5) #GEOM prediction
for(i in 1:5) lines(rep(x1[i],2), zm.quants[w.coge.in[i],c(1,5)], col=lcol)
points(x1, zm.stats[w.coge.in,1], pch=22, bg="red", cex=1)
for(i in 1:10) lines(rep(x2[i],2), zm.quants[w.demi.in[i],c(1,5)], col=lcol)
points(x2, zm.stats[w.demi.in,1], pch=21, bg="blue", cex=1)
for(i in 1:10) lines(rep(x3[i],2), zm.quants[w.cogl.in[i],c(1,5)], col=lcol)
points(x3, zm.stats[w.cogl.in,1], pch=24, bg="orange", cex=1)
for (i in 1:3) lines(rep(x4[i],2), zm.quants[l.sp[(i+6)], c(1,5)], col=lcol)
points(x4, zm.stats[l.sp[7:9],1], pch=c(21,22,24), bg=c("blue", "red", "orange"), cex=1.3)
lines(rep(x5,2), zm.quants[lbeta[3],c(1,5)], col=lcol)
points(x5, zm.stats[lbeta[3],1], pch=23, cex=2, col="black", bg="grey")
box()
abline(v=c(28.5, 34), lty=4, col="grey")
axis(2,at=c(seq(0,5,0.2)), las=1)
mtext("C", side=2, las=1, line=-1.4, at=3.35, cex=1)

#Leaf mass
plot(x, seq(1,2.8,(1.8/(max(x)-1))), cex=0.000000001, pch=21, bg="white", col="white", axes=F, xlab="", ylab="Leaf mass scaling exponent", cex.lab=1.4)
abline(h=2,lwd=1.5)
for(i in 1:5) lines(rep(x1[i],2), zm.quants[s.coge.in[i],c(1,5)], col=lcol)
points(x1, zm.stats[s.coge.in,1], pch=22, bg="red", cex=1)
for(i in 1:10) lines(rep(x2[i],2), zm.quants[s.demi.in[i],c(1,5)], col=lcol)
points(x2, zm.stats[s.demi.in,1], pch=21, bg="blue", cex=1)
for(i in 1:10) lines(rep(x3[i],2), zm.quants[s.cogl.in[i],c(1,5)], col=lcol)
points(x3, zm.stats[s.cogl.in,1], pch=24, bg="orange", cex=1)
for (i in 1:3) lines(rep(x4[i],2), zm.quants[l.sp[(i+9)], c(1,5)], col=lcol)
points(x4, zm.stats[l.sp[10:12],1], pch=c(21,22,24), bg=c("blue", "red", "orange"), cex=1.3)
lines(rep(x5,2), zm.quants[lbeta[4],c(1,5)], col=lcol)
points(x5, zm.stats[lbeta[4],1], pch=23, cex=2, col="black", bg="grey")
box()
abline(v=c(28.5, 34), lty=4, col="grey")
axis(2,at=c(seq(0,5,0.2)), las=1)
mtext("D", side=2, las=1, line=-1.4, at=2.75, cex=)


##########################################################
######################## FIGURE 3 ########################
##########################################################
#Create Figure 3: posterior species differences of length scaling exponent (density plot)
#Combine chains as a dataframe
df <- as.data.frame(rbind(zm[[1]], zm[[2]], zm[[3]]))
df.diff <- grep("beta.spp.diff", colnames(df))

#Calculate probabilities of difference not being 0
diff=numeric(3)
diff[1] = 1-ecdf(df[,df.diff[1]])(0)
diff[2] = 1-ecdf(df[,df.diff[2]])(0)
diff[3] = 1-ecdf(df[,df.diff[3]])(0)

par(mar=c(3.1,4.1,2,2.1), oma=c(0,0,0,0),mgp = c(1.8, 0.2, 0),
    tck=0.02, mfrow=c(1,1))
plot(density(df[,df.diff[1]], adjust=2), lwd=1.5, xlim=c(-0.2,0.5),
     ylim=c(0,5.2),xlab="Exponent Difference", ylab="Poseterior Density", 
     las=1, main="", cex.lab=1.5)
lines(density(df[,df.diff[2]], adjust=2), col="blue", lwd=1.5)
lines(density(df[,df.diff[3]], adjust=2), col="red", lwd=1.5)
abline(v=0, lty="dashed")
legend(0.08,5.5,legend=c(paste("demi-coge (", round(diff[1]*100), "%)"), 
                         paste("demi-cogl (", round(diff[2]*100), "%)"), 
                         paste("coge-cogl (", round(diff[3]*100), "%)")), 
       col=c("black", "blue", "red"), bty="n", lty=1, lwd=1.5)



##########################################################
######################## FIGURE 4 ########################
##########################################################
#Create Figure 4: posterior means and CIs of normalizing constants at species-level

#Get index for species-level normalizing constants
l.sp = grep("alpha.species", rownames(zm.quants))

#Set up x-values for plot spacing
x=seq(1,18,1)
x1=seq(1,3,1)
x2=seq(6,8,1)
x3=seq(11,13,1)
x4=seq(16,18,1)


par(mar=c(3.1,4.1,4.1,2.1), oma=c(0,0,0,0),mgp = c(1.8, 0.2, 0),
    tck=0.02, mfrow=c(1,1))
plot(x, seq(1,2,(1/(max(x)-1))), cex=0.000000001, pch=21, bg="white", 
     col="white", axes=F, ylab="Normalizing Constant", las=1, xlab="Trait", cex.lab=1.3)
for (i in 1:3) lines(rep(x1[i],2), zm.quants[l.sp[i], c(1,5)], col="black")
points(x1, zm.stats[l.sp[1:3],1], pch=c(21,22,24), col="black", 
       bg=c("blue", "red", "orange"),  ylim=c(0,2), cex=1)
for (i in 1:3) lines(rep(x2[i],2), zm.quants[l.sp[i+3], c(1,5)], col="black")
points(x2, zm.stats[l.sp[4:6],1], pch=c(21,22,24), col="black", 
       bg=c("blue", "red", "orange"), ylim=c(0,2), cex=1)
for (i in 1:3) lines(rep(x3[i],2), zm.quants[l.sp[i+6], c(1,5)], col="black")
points(x3, zm.stats[l.sp[7:9],1], pch=c(21,22,24), col="black", 
       bg=c("blue", "red", "orange"), ylim=c(0,2), cex=1)
for (i in 1:3) lines(rep(x4[i],2), zm.quants[l.sp[i+9], c(1,5)], col="black")
points(x4, zm.stats[l.sp[10:12],1], pch=c(21,22,24),col="black", 
       bg=c("blue", "red", "orange"), ylim=c(0,2), cex=1)
axis(2,at=seq(0,2,0.2), las=1)
axis(1, at=c(2,7,12,17), labels=c("length", "ag mass", "stem mass", "leaf mass"), cex.lab=0.3)
abline(v=4.5, lty="dotted")
abline(v=9.5, lty="dotted")
abline(v=14.5, lty="dotted")
box()
legend(4.4,2,legend=c("demi", "coge", "cogl"), pch=c(21,22,24), pt.bg=c("blue", "red", "orange"), cex=1, bty="o", bg="white", border="white", box.col="white", ncol=3)







