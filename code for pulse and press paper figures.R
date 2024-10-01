rm(list=ls())

####======Figures 2 and 4 for "From disturbances to nonlinear fitness"=====####
#== Histogram from sensitivity, J0 and second derivatives
# Phaseolus lunatus, Lima Bean

# Function to get second derivatives of population growth rate lambda
# The output is a matrix of size n^2 cross n^2, where n is dimension of original MPM
second_derivative_on_lamda0 <- function(input_mat){
  lamda0 <- Re(eigen(input_mat)$values[1]) # the largest eigenvalue
  #lamda1 <- Re(eigen(input_mat)$values[2]) # the second largest eigenvalue
  v <- Re(eigen(t(input_mat))$vectors[,1]) #first left eigenvector
  u <- Re(eigen(input_mat)$vectors[,1]) #first right eigenvector
  u <- u/sum(u)
  v <- v/c(t(v)%*%u)
  Q0 <- u%*%t(v)
  Q1 <- (1/lamda0)*(input_mat - lamda0*Q0)
  I <- diag(nrow(input_mat))
  # J0 <- solve(I-Q1) %*% (I-Q0)/lamda0 # follow equation 12 in "From disturbances to nonlinear fitness and back"
  J0 <- get_J0(input_mat)
  n <- nrow(input_mat)
  d2_vector <- rep(NA, n^2)
  e2_vector <- rep(NA, n^2)
  i <- 1
  # enumerate by column: from column first and then row, elements here are input_mat[p,q] and input_mat[k,l]
  for(q in 1:n){
    for (p in 1:n) {
      for (l in 1:n) {
        for (k in 1:n) {
          d2_vector[i] <- (v[p]*J0[q,k]*u[l] + v[k]*J0[l,p]*u[q])
          #print(paste0("(",p,", ",q,")(",k,", ",l,")"))
          i <- i+1
        }
      }
    }
  }
  #d2mat <- d2vector_to_d2matrix(d2_vector)
  d2mat <- matrix(d2_vector, n^2, n^2, byrow = T)
  return(d2mat)
}

# Function to get Transient response matrix: J0
get_J0 <- function(input_mat){
  lamda0 <- Re(eigen(input_mat)$values[1]) # the largest eigenvalue
  #lamda1 <- Re(eigen(input_mat)$values[2]) # the second largest eigenvalue
  v <- Re(eigen(t(input_mat))$vectors[,1]) #first left eigenvector
  u <- Re(eigen(input_mat)$vectors[,1]) #first right eigenvector
  u <- u/sum(u)
  v <- v/c(t(v)%*%u)
  Q0 <- u%*%t(v)
  Q1 <- (1/lamda0)*(input_mat - lamda0*Q0)
  I <- diag(nrow(input_mat))
  J0 <- solve(I-Q1) %*% (I-Q0)/lamda0 # follow equation 48 on overleaf
  return(J0)
}

# Function to get first order sensitivity
get_sensitivity <- function(input_mat){
  lamda0 <- Re(eigen(input_mat)$values[1]) # the largest eigenvalue
  #lamda1 <- Re(eigen(input_mat)$values[2]) # the second largest eigenvalue
  v <- Re(eigen(t(input_mat))$vectors[,1]) #first left eigenvector
  u <- Re(eigen(input_mat)$vectors[,1]) #first right eigenvector
  u <- u/sum(u)
  v <- v/c(t(v)%*%u)
  sens <- v%*%t(u)
  return(sens)
}

# We use a matrix from Compadre database. The plant is Phaselous Lunatus (Lima bean)
# Save all the files in the same working directory including the mpm1 matrix!!!
mpm1 <- read.csv("mpm1.csv", header = T, colClasses=c("NULL",NA,NA,NA,NA,NA,NA))
sde <- second_derivative_on_lamda0(mpm1)
J0_mat <- get_J0(mpm1)
sen <- get_sensitivity(mpm1)

# Figure 4 in the manuscript 
pdf("Hists_sens-J0-2ndder1.pdf", width = 6, height = 8)
nf <- layout(matrix(c(1,2,3), ncol=1))
par(mar = c(5.1, 5.1, 4.1, 2.1))
hist(sde, main = NULL,
     ylab = "# of values",
     xlab = "Second Derivative",
     cex.lab = 2, cex.axis=1.5)
abline(v=0, lwd = 3, lty = 1)
hist(J0_mat, main = NULL,
     ylab = "# of values",
     xlab = "TRM",
     cex.lab = 2, cex.axis=1.5)
abline(v=0, lwd = 3, lty = 1)
#segments(0,-0.2,3.5,-0.2, lwd = 5)
hist(sen, main = NULL,
     ylab = "# of values",
     xlab = "Sensitivity",
     cex.lab = 2, cex.axis=1.5)
abline(v=0, lwd = 3)
dev.off()

#== SSD dynamic on a give stage ==
mpm1 <- read.csv("mpm1.csv", header = T, colClasses=c("NULL",NA,NA,NA,NA,NA,NA))
perMat <- matrix(0, nrow = 6, ncol = 6)

g <- perMat[1,5] <- 0.02*mpm1[1,5]
f <- perMat[6,5] <- 0.02*mpm1[6,5]

ts <- 20
#pulse
bmat <- as.matrix(mpm1)
bee <- eigen(bmat)
lam0 <- Re(bee$values[1])
u0 <- Re(bee$vectors[,1])
u0 <- u0/sum(u0)
pulse_ssd <- u0
bmat <- as.matrix(mpm1+perMat)
u0 <- bmat%*%u0/sum(bmat%*%u0)
pulse_ssd <- cbind(pulse_ssd,u0)
bmat <- as.matrix(mpm1)

for (i in 3:ts) {
  u0 <- bmat%*%u0/sum(bmat%*%u0)
  pulse_ssd <- cbind(pulse_ssd,u0)
}

#press
bmat <- as.matrix(mpm1)
bee <- eigen(bmat)
lam0 <- Re(bee$values[1])
u0 <- Re(bee$vectors[,1])
u0 <- u0/sum(u0)
press_ssd <- u0
bmat <- as.matrix(mpm1+perMat)
for (i in 2:ts) {
  u0 <- bmat%*%u0/sum(bmat%*%u0)
  press_ssd <- cbind(press_ssd,u0)
}

# Figure 2 in the manuscript
pdf("Stage_dynamics_together0.02_largefont.pdf", width = 10, height = 5)
nf <- layout(matrix(c(1,2), ncol=2, byrow = T))
par(mar = c(5.1, 5.1, 3.1, 2.1))
plot(0:(ts-1), pulse_ssd[1,]*10, type = "l",
     xlab = "Time Steps",
     ylab = expression("The proportion of stage 1, (x"*10^{-1}*")"),
     #main = "Pulse disturbance",
     cex.lab = 1.5, cex.axis=1.5)
abline(h=pulse_ssd[1,1]*10, lty = 2)
lines(0:(ts-1), press_ssd[1,]*10, col = "red")
plot(pulse_ssd[5,-1]*1000, pulse_ssd[1,-1]*10, type = "p", pch = 16,
     xlim = range(c(press_ssd[5,]*1000,pulse_ssd[5,]*1000)),
     ylim = range(c(press_ssd[1,]*10,pulse_ssd[1,]*10)),
     xlab = expression("The proportion of stage 5, (x"*10^{-3}*")"),
     ylab = expression("The proportion of stage 1, (x"*10^{-1}*")"),
     cex.lab = 1.5, cex.axis=1.5)

points(press_ssd[5,-1]*1000, press_ssd[1,-1]*10, type = "p", pch = 16, col = "red")
s <- 2:(ts-1)
for (i in s) {
  arrows(press_ssd[5,i]*1000, press_ssd[1,i]*10, press_ssd[5,i+1]*1000, press_ssd[1,i+1]*10, length=min(0.16,0.16/(i*0.3)), col = "red")
  arrows(pulse_ssd[5,i]*1000, pulse_ssd[1,i]*10, pulse_ssd[5,i+1]*1000, pulse_ssd[1,i+1]*10, length=min(0.2,0.2/(i*0.3)))
}
text(pulse_ssd[5,2]*1000+0.015, pulse_ssd[1,2]*10, expression(t==1), cex = 2)
points(pulse_ssd[5,1], pulse_ssd[1,1]*10, pch = 17, cex = 3)
segments(pulse_ssd[5,1]*1000, pulse_ssd[1,1]*10, pulse_ssd[5,2]*1000, pulse_ssd[1,2]*10,lty = 1)
segments(press_ssd[5,1]*1000, press_ssd[1,1]*10, press_ssd[5,2]*1000, press_ssd[1,2]*10,lty = 2, col = "red")
dev.off()

###############################################
# Code for Figure 5 in the manuscript
###############################################
install.packages(c("ape",
                   "phytools",
                   "reshape",
                   "caper",
                   "tidyverse",
                   "picante",
                   "ggtext",
                   "ggpubr"), dependencies = T)
library(ape)
library(phytools)
library(reshape)
library(caper)
library(picante)
library(ggtext)

# Check your working directory
getwd()

# Save your working directory to where the data files are saved/downloaded
setwd("/Users/harmanjaggi/Documents/Research/Damping/Damping codes")

full_uni <- read.csv2("./match_list_uncorrected_newversion v2.csv",
                      sep = ",")
all_data <- read.csv2("./full animal and plant data v2.csv")
length(unique(all_data$SpeciesAccepted))
final_tree_read<-read.tree(file="./final_tree_corrected_newversion v2.tre")
final_tree_read$tip.label <-  gsub("_", " ",final_tree_read$tip.label)
species_tree = as.vector(unique(final_tree_read$tip.label))

species_name_correct <- as.vector(all_data$SpeciesAccepted)


library(tidyverse)
## check the variance of each variable within species
check <- all_data %>%
  dplyr::group_by(SpeciesAccepted, db_sep)%>%
  summarise(number = n(),
            Tc.var = var(Tc),
            sigma.var = var(sigma),
            eig.j0.var = var(eig.j0),
            damping.time.var = var(damping.time),
            Tc.mean = mean(Tc),
            sigma.mean = mean(sigma),
            damping.time.mean = mean(damping.time),
            Tc.cv = sqrt(Tc.var)/Tc.mean,
            sigma.cv = sqrt(sigma.var)/sigma.mean,
            damping.time.cv = sqrt(damping.time.var)/damping.time.mean)
# summary(check$Tc.mean)
# summary(check$sigma.mean)
# summary(check$Tc.var)
# summary(check$sigma.var)
# summary(check$Tc.cv)
# summary(check$sigma.cv)
# summary(check$damping.time.cv)

# hist(log(check$number))
# summary(check$number)

# unique(all_data$Class)
# length(unique(filter(all_data, Class %in% c("Actinopterygii", "Aves", "Mammalia", "Reptilia"))$SpeciesAccepted))

matrix_select<-all_data %>%
  group_by(SpeciesAccepted, db_sep)%>%
  mutate(number = n())

## use the median value of damping time to choose one matrix for each species
median_index = function(x) {
  lx = length(x)
  if (lx %% 2 == 1) {
    return(median(x))
  }
  return(median(c(x,0))) # for even number, add 0 in the end
}

median_data_sep <- all_data %>%
  dplyr::select(SpeciesAccepted, Tc, sigma, eig.j0, damping.time, db_sep, db_taxa,db_source,
                Class, omega, alpha, damping.cal, damping.approx, Order, ProjectionInterval)%>%
  group_by(SpeciesAccepted)%>%
  mutate(sigma.median = median_index(sigma))%>%
  distinct(SpeciesAccepted,.keep_all=TRUE)%>%
  mutate(logTc = log(Tc),
         logsigma = log(sigma),
         logeigj0 = log(eig.j0),
         logdamping.time = log(damping.time),
         logreptime = log(omega - alpha))%>%
  filter(SpeciesAccepted %in% species_tree)%>%
  as.data.frame()

length(unique(median_data_sep$SpeciesAccepted))
summary(median_data_sep$ProjectionInterval) # NAs are for the GMO dataset,for which the number should be 1
table(median_data_sep$ProjectionInterval)
table(median_data_sep$db_source,
      median_data_sep$ProjectionInterval)

median_data_sep %>%
  group_by(db_source)%>%
  summarise(count = n_distinct(ProjectionInterval))

median_data_sep %>%
  group_by(db_source)%>%
  summarise(count = n_distinct(SpeciesAccepted))

median_data_sep %>%
  group_by(db_source)%>%
  summarise(count = n_distinct(SpeciesAccepted))

# Check if the tree is rooted
tree <- final_tree_read

is.rooted(tree)

# If it is binary
is.binary(tree)

# Make the node labels unique
tree<- makeNodeLabel(tree)

# PGLS analyses -------------------------------------------------------------------------------------------------
#### PGLS function ####
# need to have a grp variable in input data
PGLS_fun = function(data, variable, data_group, outputlist){
  output_combine=data.frame()
  db_combine=data.frame()

  for (d in 1:(length(data_group))) {
    data_filter <- filter(data, grp %in% data_group[d])
    comp_data <- comparative.data(phy = final_tree_read,
                                  data = data_filter,
                                  names.col = SpeciesAccepted,
                                  vcv = TRUE,
                                  vcv.dim = 3)
    rowname_data <- rownames(comp_data$data)
    count=0
    for (expl in 1:(length(variable)-1)) {
      for (resp in (expl+1):(length(variable))) {
        count=count+1
        # take log
        explanatory=as.numeric(as.matrix(comp_data$data[variable[expl]]))
        response=as.numeric(as.matrix(comp_data$data[variable[resp]]))
        mod=pgls(explanatory ~ response, comp_data, lambda='ML', param.CI = 0.95)

        output=c(data_group[d],
                 paste(variable[expl],"~",variable[resp]),
                 summary(mod)$coeff[1,1],
                 summary(mod)$coeff[2,],
                 mod$param[1],
                 mod$param[2],
                 mod$param.CI$lambda$ci.val[1],
                 mod$param.CI$lambda$ci.val[2],
                 mod$mlVals[2],
                 summary(mod)$r.squared,
                 unique(data_filter$db_taxa))

        db = comp_data[["data"]] %>%
          mutate(residual = mod$residual,
                 variable = paste(variable[expl],"~",variable[resp]))

        print(paste("PGLS model: ",variable[expl]," ~ ",variable[resp],sep=""))

        db_combine = rbind(db_combine,db)
        output_combine = rbind(output_combine,output)
      }
    }
  }
  colnames(output_combine) <- c(outputlist[-length(outputlist)])
  output_combine[,"Padjusted"]=p.adjust(output_combine[,"P"],"BH")
  output_combine <- as.data.frame(output_combine)%>%
    mutate_at(c(3:13,15), as.numeric)%>%
    arrange(grp)%>%
    as.data.frame()

  return(list(output_combine, db_combine))
}


variable <- c("logdamping.time", "logeigj0", "logTc")

median_data_sep = median_data_sep %>%
  mutate(grp = db_sep)
data_group <- unique(median_data_sep$grp)

# We create an empty lists where we will store the results
outputlist=c("grp","variable","intercept",
             "slope","SE","t","P",
             "Kappa",
             "Pagel lambda","lambda CI lower","lambda CI upper",
             "LambdaOptimization",
             "R2","db_taxa", "Padjusted")

result = PGLS_fun(median_data_sep, variable, data_group, outputlist)

saveTable = result[[1]]%>%
  mutate(slope_lower = slope - SE,
         slope_upper = slope + SE)

names(saveTable)
saveTable %>% group_by(db_taxa) %>% summarise(`Pagel lambda`)
residual = result[[2]]
# Save the results
# write.csv(saveTable, file = "./phylogenetic tree/PGLsList_uncorrected by seperate groups (median).csv") # PGLS raw correlations

pos = median_data_sep %>% group_by(db_sep) %>%
  summarise(Tc.pos = 0.8 * min(log(Tc)),
            sigma.pos = 0.9 * max(log(sigma)),
            eigj0.pos = 0.9 * max(log(eig.j0)),
            damp.pos = 0.9 * max(log(damping.time)),
            sigma.pos2 = sigma.pos - (max(log(sigma) - min(log(sigma))))*0.15,
            damp.pos2 = damp.pos - (max(log(damping.time) - min(log(damping.time))))*0.15)

#### Tc and Eig of J0 #####
summ <- saveTable %>%
  group_by(grp) %>%
  filter(variable %in% "logeigj0 ~ logTc")%>%
  mutate_if(is.numeric, round, digits=2)%>%
  mutate(intercept2 = case_when(
    intercept >=0 ~ paste("+", intercept),
    TRUE ~ paste("-",as.character(-intercept))
  ))%>%
  mutate(p_report = case_when(
    Padjusted < 0.001 ~ "***",
    0.001< Padjusted &  Padjusted<= 0.01 ~ "**",
    0.001< Padjusted  &  Padjusted<= 0.05 ~ "*"))

f_labels <- data.frame(
  db_sep = summ$grp,
  label1 = c(paste0("y = ",summ$slope, "x ",summ$intercept2)),
  label2 = c(paste0("R<sup>2</sup> = ",summ$R2))
  # label2 = c(paste0("R<sup>2</sup> = ",summ$R2,", ", summ$p_report))
)

unique(median_data_sep$db_sep)

minj <- min(log(median_data_sep$eig.j0))
maxj <- max(log(median_data_sep$eig.j0))
mint <- min(log(median_data_sep$Tc))
maxt<-max(log(median_data_sep$Tc))

names(median_data_sep)

# Figure 5 in the paper is a composite of plot1, plot2, and plot3
# The Figure correlates generation time Tc with dominant eigenvalue for J0 matrix

plot1 <-  ggplot(filter(median_data_sep, db_sep=="Animal by age"),
                     aes(log(Tc), log(eig.j0)))+
  geom_point(size=4, alpha = 0.8, color="orange2")+
  # stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
  #              parse = TRUE, color = "blue", 
  #              rr.digits = 2, coef.digits = 2, size = 5) +
  
  # scale_shape_manual(values=myshape)+
  theme_bw()+
  ylim(c(minj,maxj))+
  xlim(c(mint,maxt))+
  labs(title = "COMADRE Age: Animals")+
  theme(axis.text.x = element_text(color = "grey20", size = 17),
        axis.text.y = element_text(color = "grey20", size = 17), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        # axis.title.x = element_text(color = "grey20", size = 22),
        # axis.title.y = element_text(color = "grey20", size = 22),
        plot.title = element_text(size=18),
        legend.text=element_text(size=18),
        strip.text.x = element_text(size = 15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
    geom_smooth(method='lm', formula= y~x, se=F, color=c("darkblue"), size=0.2, linetype="dashed")
  # facet_wrap(. ~db_sep, nrow = 3, scales = "free_y")+
  # geom_richtext(data = f_labels, aes(label = label1), size = 5,
  #               x = -2, y = pos$sigma.pos-2,
  #               hjust = 0,
  #               color = "black",
  #               fill = NA,
  #               label.colour = NA)+
  # geom_richtext(data = f_labels, aes(label = label2), size = 5,
  #               x = -2, y = pos$sigma.pos2-3,
  #               hjust = 0,
  #               color = "black",
  #               fill = NA,
  #               label.colour = NA)

plot2 <-  ggplot(filter(median_data_sep, db_sep=="Animal by stage"),
                 aes(log(Tc), log(eig.j0)))+
  geom_point(size=4, alpha = 0.8, color="cyan3")+
  # stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
  #              parse = TRUE, color = "blue", 
  #              rr.digits = 2, coef.digits = 2, size = 5) +
  
  # scale_shape_manual(values=myshape)+
  theme_bw()+
  ylim(c(minj,maxj))+
  xlim(c(mint,maxt))+
  labs(title = "COMADRE Stage: Animals")+
  theme(axis.text.x = element_text(color = "grey20", size = 17),
        axis.text.y = element_text(color = "grey20", size = 17), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        # axis.title.x = element_text(color = "grey20", size = 22),
        # axis.title.y = element_text(color = "grey20", size = 22),
        plot.title = element_text(size=18),
        legend.text=element_text(size=18),
        strip.text.x = element_text(size = 15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  geom_smooth(method='lm', formula= y~x, se=F, color=c("darkblue"), size=0.2, linetype="dashed")

plot3 <-  ggplot(filter(median_data_sep, db_sep=="Plant by stage"),
                 aes(log(Tc), log(eig.j0)))+
  geom_point(size=4, alpha = 0.8, color="darkolivegreen3")+
  # stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
  #              parse = TRUE, color = "blue", 
  #              rr.digits = 2, coef.digits = 2, size = 5) +
  
  # scale_shape_manual(values=myshape)+
  theme_bw()+
  ylim(c(minj,maxj))+
  xlim(c(mint,maxt))+
  labs(title = "COMPADRE Stage: Plants")+
  theme(axis.text.x = element_text(color = "grey20", size = 17),
        axis.text.y = element_text(color = "grey20", size = 17), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        # axis.title.x = element_text(color = "grey20", size = 22),
        # axis.title.y = element_text(color = "grey20", size = 22),
        plot.title = element_text(size=18),
        legend.text=element_text(size=18),
        strip.text.x = element_text(size = 15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  geom_smooth(method='lm', formula= y~x, se=F, color=c("darkblue"), size=0.2, linetype="dashed")

j0_phylo <- ggpubr::ggarrange(plot1, plot2, plot3, ncol=3)

ggpubr::annotate_figure(j0_phylo, left = grid::textGrob(expression("log of Dominant Eigenvalue of J"["0"]),
                                          rot = 90, vjust = 0.5, gp = grid::gpar(cex = 1.7)),
                bottom = grid::textGrob("log Generation time ", gp = grid::gpar(cex = 1.7)))

# For displaying the slopes and R^2 values, please uncomment the stat_poly_eq() lines and rerun by loading library(ggpmisc)!
