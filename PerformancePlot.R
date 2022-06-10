time=c(5.98,6.15,7.93,28.1,153.3)
N=c(10,100,1000,10000,58673)

plot(lm(time~N))

ggplot()+
  geom_smooth(aes(N,time),method = "loess",)+
  geom_point(aes(N,time))

'''
Coefficients:
  (Intercept)            N  
5.137157     0.002519  
'''

