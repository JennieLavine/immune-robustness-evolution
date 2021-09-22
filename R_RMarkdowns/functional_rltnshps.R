
######
# Relationships between viral genome size and immune power
#####3
par(lwd=3)
curve(1/(9*x+1), ylim=c(0,1), xlim=c(0,1), 
      xlab = "proportion of max disruption", 
      ylab='proportion of max immune power',
      col='blue')
curve(1/(9*(x-(11.25/9.225))+1) + 1.125, add=T, lty='dashed', col='cyan')
curve(-0.9*x+1, add=T, col='blue')
curve(1 - (x^4 / (x^4 +0.1)), add=T, lty='dashed', col='cyan')
legend('bottomleft', lty=1:2, col=c('blue', 'cyan'), legend=c('non-robust', 'robust'))


#######
#Relationship between viral genome size and viral replication rate
########3

curve(1.5- (x/(x+0.4)), xlim=c(0,1),
      xlab= "proportion of max # additional genes",
      ylab="viral replication rate")

