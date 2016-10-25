n1 = 10000 # number of samples at each p
n2 = 20 # number of steps along p

p1 = seq(0,1, length.out=n2)
p2 = p1

ps <- numeric(n1*n2)
s1 <- p1
s2 <- s1

for(j in 1:n2){
	print(j)
	for(i in 1:n1){
		ps[(j-1)*n1 + i] <- p1[j]
		s1[(j-1)*n1 + i] <- mean(rbinom(58,1,p1[j]))
		s2[(j-1)*n1 + i] <- mean(rbinom(48,1,p2[j]))
	}
}

plot(s1, s2, xlim=c(0,1), ylim=c(0,1))

plot(s1, abs(s1-s2), xlim=c(0,1), ylim=c(0,1))
plot(s1, abs(s1-s2), xlim=c(0,1), ylim=c(0,1), pch=16, cex=0.5, col=rgb(0,0,0,0.2))

plot(jitter(s1), abs(s1-s2), xlim=c(0,1), ylim=c(0,1))


i <- s1 > 0 & s1 < 0.1
d1 <- density((abs(s1[i] - s2[i])), adjust=2)
i <- s1 > 0.2 & s1 < 0.3
d2 <- density((abs(s1[i] - s2[i])), adjust=2)
i <- s1 > 0.4 & s1 < 0.5
d3 <- density((abs(s1[i] - s2[i])), adjust=2)
i <- s1 > 0.6 & s1 < 0.7
d4 <- density((abs(s1[i] - s2[i])), adjust=2)
i <- s1 > 0.8 & s1 < 0.9
d5 <- density((abs(s1[i] - s2[i])), adjust=2)
i <- s1 > 0.9 & s1 < 1
d6 <- density((abs(s1[i] - s2[i])), adjust=2)



plot(d1)
lines(d2, col='red')
lines(d3, col='blue')
lines(d4, col='green')
lines(d5, col='purple')
lines(d6, col='orange')