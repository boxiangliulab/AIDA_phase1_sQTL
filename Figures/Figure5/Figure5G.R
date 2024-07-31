library(ggplot2)
library(ggpubr)

spliz1.subset <- read.delim('data/HNRNPLL.CD4TcmTemTnaive.counts.txt')
table(spliz1.subset$HNRNPLL_expression)
#not expressed     expressed 
#214504         53064 

spliz1.subset.0.scZ <- subset(spliz1.subset, expr_0_1 == 'not expressed')$scZ
spliz1.subset.1.scZ <- subset(spliz1.subset, expr_0_1 == 'expressed')$scZ

median(spliz1.subset.0.scZ)
median(spliz1.subset.1.scZ)
mean(spliz1.subset.0.scZ)
mean(spliz1.subset.1.scZ)

CI0 <- confint(lm(spliz1.subset.0.scZ ~ 1, as.data.frame(spliz1.subset.0.scZ)), level=0.95)[2] - mean(spliz1.subset.0.scZ)
CI1 <- confint(lm(spliz1.subset.1.scZ ~ 1, as.data.frame(spliz1.subset.1.scZ)), level=0.95)[2] - mean(spliz1.subset.1.scZ)
data <- data.frame(
  name=factor(c('not expressed', 'expressed'), levels = c('not expressed', 'expressed')),
  value=c(mean(spliz1.subset.0.scZ), mean(spliz1.subset.1.scZ)),
  sd=c(sd(spliz1.subset.0.scZ)/sqrt(length(spliz1.subset.0.scZ)), sd(spliz1.subset.1.scZ)/sqrt(length(spliz1.subset.1.scZ))),
  CI=c(CI0, CI1)
)

ggplot(data) +
  geom_errorbar( aes(x=name, ymin=value-CI, ymax=value+CI), color="#807dba", width=0) + theme_classic() +
  geom_point( aes(x=name, y=value), stat="identity", color="#807dba") +
  ylab('spliZ') + xlab('hnRNPLL expression')
ggsave('Fig.5g.pdf', width=2.5, height=2.2)


t1 <- t.test(spliz1.subset.0.scZ, spliz1.subset.1.scZ)
t1$p.value
#[1] 3.068184e-197
t1$statistic
#        t 
#-30.04473 
t1$parameter
#df 
#81578.31 
