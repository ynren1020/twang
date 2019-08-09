####################2019-08-09########################
##simulate data and plot log x and logy scatterplot###
##reference to the figure 32.2 p.483 of book #########
##Computational Exome and Genome Analysis  ###########
######################################################

## Show the different plot types
x <- 0:12
y <- sin(pi/5 * x)
op <- par(mfrow = c(3,3), mar = .1+ c(2,2,3,1))
for (tp in c("p","l","b",  "c","o","h",  "s","S","n")) {
  plot(y ~ x, type = tp, main = paste0("plot(*, type = \"", tp, "\")"))
  if(tp == "S") {
    lines(x, y, type = "s", col = "red", lty = 2)
    mtext("lines(*, type = \"s\", ...)", col = "red", cex = 0.8)
  }
}
par(op)





##--- Log-Log Plot  with  custom axes
lx <- seq(1, 5, length = 41)
yl <- expression(e^{-frac(1,2) * {log[10](x)}^2})
y <- exp(-.5*lx^2)
op <- par(mfrow = c(2,1), mar = par("mar")-c(1,0,2,0), mgp = c(2, .7, 0))
plot(10^lx, y, log = "xy", type = "l", col = "purple",
     main = "Log-Log plot", ylab = yl, xlab = "x")
plot(10^lx, y, log = "xy", type = "o", pch = ".", col = "forestgreen",
     main = "Log-Log plot with custom axes", ylab = yl, xlab = "x",
     axes = FALSE, frame.plot = TRUE)
my.at <- 10^(1:5)
axis(1, at = my.at, labels = formatC(my.at, format = "fg"))
e.y <- -5:-1 ; at.y <- 10^e.y
axis(2, at = at.y, col.axis = "red", las = 1,
     labels = as.expression(lapply(e.y, function(E) bquote(10^.(E)))))
par(op)


##ref to figure##
x1<-seq(exp(-4.5),exp(-2.5),length=500)
x2<-seq(exp(-2.5),exp(-2),length=20)
x3<-seq(exp(-6),exp(-4.5),length=80)
x4<-seq(exp(-9),exp(-6),length=5)

x<-c(x4,x3,x1,x2)


y1<-seq(exp(-5),exp(0),length=600)
y2<-seq(exp(-20),exp(-5),length=5)
y<-c(y1,y2)
ly<--log(y)
lx<-log(x)

sly<-sample(ly,605)
slx<-sample(lx,605)

plot(slx,sly)


######3
y1<-seq(0,0.3,length=50)
x1<-c(-9,-7,-6.7,-6.67,-6.65,-6.63,-6.61,-6.6,-6.55,-6.5,-6.4,-6.45,-6.3,-6.2,-6.1,-6,-5.9,-5.8,-5.7,-5.6,seq(-5.6,-4.2,length=30))

y2<-c(5,6,9,13,19)
x2<-c(-2.5,-3.7,-4,-4.1,-2.9)

x3<-sample(seq(-5.8,-5.6,length=30),50,replace = TRUE)
y3<-sample(seq(0,0.3,length=30),50,replace = TRUE)

x4<-sample(seq(-5.6,-5.3,length=50),230,replace = TRUE)
y4<-c(sample(seq(0,0.45,length=50),130,replace = TRUE),sample(seq(0.45,0.7,length=20),70,replace = TRUE),sample(seq(0.7,1.5,length=10),30,replace = TRUE))

x5<-sample(seq(-5.3,-5.0,length=50),450,replace = TRUE)
y5<-c(sample(seq(0,1.0,length=50),280,replace = TRUE),sample(seq(1,2,length=30),130,replace = TRUE),sample(seq(2,2.5,length=30),40,replace = TRUE))

x6<-sample(seq(-5.0,-4.7,length=50),900,replace = TRUE)
y6<-c(sample(seq(0,1.2,length=50),500,replace = TRUE),sample(seq(1.2,3,length=30),350,replace = TRUE),sample(seq(3,4.5,length=30),50,replace = TRUE))

x7<-sample(seq(-4.7,-4.4,length=50),550,replace = TRUE)
y7<-c(sample(seq(0.2,1.7,length=50),270,replace = TRUE),sample(seq(1.7,3,length=20),200,replace = TRUE),sample(seq(3,5,length=10),80,replace = TRUE))

x8<-sample(seq(-4.4,-4.0,length=50),300,replace = TRUE)
y8<-c(sample(seq(0.4,1.9,length=50),170,replace = TRUE),sample(seq(1.9,4.5,length=20),130,replace = TRUE))

x9<-sample(seq(-4.0,-3.5,length=50),80,replace = TRUE)
y9<-c(sample(seq(0.45,2,length=50),60,replace = TRUE),sample(seq(2,4.9,length=20),20,replace = TRUE))

x10<-sample(seq(-3.5,-2.5,length=30),30,replace = TRUE)
y10<-sample(seq(0.5,2.5,length=30),30,replace = TRUE)




#x3<-sample(seq(-5.8,-3,length=300),1550,replace = TRUE)
#y3<-c(sample(seq(0,3,length=250),1500,replace=TRUE),sample(seq(0,5,length=50),50))
#x4<-sample(seq(-3.2,-2,length=50),50,replace=TRUE)
#y4<-sample(seq(0,5,length=50),50,replace=TRUE)

plot(c(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10),c(y1,y2,y3,y4,y5,y6,y7,y8,y9,y10))



##beta distribution##use this one##
set.seed(1)
p_sim <- rbeta(1000, 5, 1)
y_sim<-dbeta(p_sim,5, 1)
plot(p_sim,y_sim)
dat<-data.frame(p_sim=p_sim,y_sim=y_sim)
library(ggplot2)
jitter <- position_jitter(width = 0.01*p_sim, height = 0.1*y_sim)
ggplot(dat, aes(p_sim, y_sim)) +
  #geom_jitter()+
  geom_point(position = jitter, color = "red", aes(p_sim + 0.2, y_sim + 0.2))+
  theme_classic()
ggsave("simulate.beta.jitter.pdf",width = 6,height = 6,dpi = 300)

logit=function(x){log(x/(1-x))}
de=1/(1+exp(-rnorm(10^3,mean=logit(p_sim),sd=2)))
#hist(p_sim, freq = FALSE)
#curve(dbeta(x, 5, 1),
#      add = TRUE, col = "red", 
#      lwd = 2)

#plot(p_sim)
##normal##
c <- 1
d <- 2

a <- -2
b <- 3.5

ll <- pnorm(a, c, d)
ul <- pnorm(b, c, d)

x <- qnorm( runif(3000, ll, ul), c, d )
hist(x)
range(x)
mean(x)
sd(x)
plot(x, type='l')

##exp##
d<-rexp(50)
