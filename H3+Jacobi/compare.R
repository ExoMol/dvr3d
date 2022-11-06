data1=read.table("dipole.data.com")
data2=read.table("tmp2")

odata1=data1[order(abs(data1$V5)),]
odata2=data2[order(abs(data2$V4)),]

freq1=odata1$V5
freq2=odata2$V4

coeff1=odata1$V10
coeff2=odata2$V3


#(label1)
#plot(abs(freq1),coeff1,type = "l")
#points(freq2,coeff2,col=2)
