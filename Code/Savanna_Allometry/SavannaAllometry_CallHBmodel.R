
### NOTICE!! ### NOTICE!! ### NOTICE!! ### NOTICE!! ### NOTICE!! ### NOTICE!! ### NOTICE!! ### NOTICE!! ####
#                                                                                                          #
# To run this code you must first install JAGS on your computer.                                           #
# Binary packages for Windows and OS X can be found here: http://sourceforge.net/projects/mcmc-jags/files/ # 
# Information on JAGS can be found here: http://mcmc-jags.sourceforge.net/                                 #
# There is a line of code that will automatically send you to the link for JAGS source files                #
# Comment out (#) the "browseURL" line after JAGS has been installed                                       #
#                                                                                                          #
############################################################################################################

rm(list=ls())

#First, install JAGS from this link, comment out (#) after install
browseURL("http://sourceforge.net/projects/mcmc-jags/files/", browser = getOption("browser"), encodeIfNeeded = FALSE)

#Find the JAGS HB model file path on your computer
##ANDREW -- can I do this from a link????
HB.model <- file.choose()

#Install necessary packages
install.packages('rjags')
install.packates('coda')
install.packages('rdryad')
install.packages('ggplot2')

#Load packages
library(rjags)
library(coda)
library(rdryad)
library(ggplot2)

##Bring in data from Dryad
data.link <- "hyperlink/here/"
data <- dryad_getfile(data.link)

names(data) = tolower(names(data))
data = data[data$aggregated_leaf_wt>-9000,] #take out bad data

##Set up matrix to store LogDiameter, LogLength, LogMass, LogWoodMass and LogLeafMass. 
Y = matrix(nrow=length(data[,1]), ncol=4)

#Log data and put in dependent variable data matrix
Y[,1] = log10(data$length)
Y[,2] = log10(data$total_mass)
Y[,3] = log10(data$aggregated_wood_wt)
Y[,4] = log10(data$aggregated_leaf_wt)

LogDiameter = log10(data$diameter)

##Set up positive definite matrix for derivition of Omega precision matrix
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

# vars <- c("mu.alpha", "mu.beta", "alpha.species", "beta.species", 
#             "Sigma", "p.mu", "p.fit",  "spp.var", "ind.var", 
#               "betalength.coge", "betamass.coge", "betawoodm.coge", "betaleafm.coge", 
#               "betalength.demi", "betamass.demi", "betawoodm.demi", "betaleafm.demi", 
#               "betalength.cogl", "betamass.cogl", "betawoodm.cogl", "betaleafm.cogl", 
#               "beta.spp.diff", "alpha.spp.diff", "y.pred")

vars <- c("y.pred")

#for production runs: 50,000 adapt; 200,000 burn-in; 1,000,000 samples; thin chains by 10
#Change these to lower values (specifically sample.n) for test runs
adapt.n <- 50000
burn.n <- 200000
sample.n <- 1000000
thin.n <- 10

#Call JAGS HB model
jm <- jags.model(HB.model, data=data.J, inits, n.chains=length(inits), n.adapt=adapt.n)

#Burn in
update(jm, n.iter=burn.n)

#Run MCMC and track nodes as listed in 'vars'
zm <- coda.samples(jm, variable.names=vars, n.iter=samples.n, thin=thin.n)

#Create quantile matrix from posterior samples
y.quants <- summary(zm)$quantile

#Create statistical matrix (mean, sd, etc.) from posterior samples
y.stats <- summary(zm)$stat

#Check convergence with Heidel diagnostic
h <- heidel.diag(zm)
print, h


##########################################################
######################## FIGURE 1 ########################
##########################################################
#Create Figure 1: fitted 'global', interspecific allometries

#set-up dataframes with predicted data and CIs
#back-transform data from Log10 for dataframes, can be plotted in log or not in ggplot
df.length <- data.frame(x=x.pred,
                        y=as.numeric(10^y.stats[1:181,2]),
                        y.low=as.numeric(10^y.quants[1:181,2]),
                        y.hi=as.numeric(10^y.quants[1:181,6]))
df.agmass <- data.frame(x=x.pred,
                        y=as.numeric(10^y.stats[182:362,2]),
                        y.low=as.numeric(10^y.quants[182:362,2]),
                        y.hi=as.numeric(10^y.quants[182:362,6]))
df.smass <- data.frame(x=x.pred,
                       y=as.numeric(10^y.stats[363:543,2]),
                       y.low=as.numeric(10^y.quants[363:543,2]),
                       y.hi=as.numeric(10^y.quants[363:543,6]))
df.lmass <- data.frame(x=x.pred,
                       y=as.numeric(10^y.stats[544:724,2]),
                       y.low=as.numeric(10^y.quants[544:724,2]),
                       y.hi=as.numeric(10^y.quants[544:724,6]))

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
  #   scale_y_continuous(limits = c(0,5000)) +
  #   scale_x_continuous(limits = c(0,20)) +
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
library(gridExtra)


grid.arrange(arrangeGrob(length + opts(legend.position = "none"), 
                         agmass + opts(legend.position = "none"),
                         smass + opts(legend.position = "none"),
                         lmass + opts(legend.position = "none"),
                         nrow=2),
             legend, 
             widths=unit.c(unit(1, "npc") - legend$width, legend$width), 
             nrow=1)









