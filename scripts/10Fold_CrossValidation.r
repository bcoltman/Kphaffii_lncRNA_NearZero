library(ROCR)

args <- commandArgs(TRUE)
CPAT_feature <- args[1]
output_dir <- args[2]

system(paste("mkdir ", output_dir, sep = ""))

data=read.table(file=CPAT_feature,header=T,sep="\t",quote="", fill=FALSE)
attach(data)
names(data)

len = nrow(data)
indices = seq(1,len)
splits = floor(len/10)

all = sample(indices,len, replace=FALSE) 


#total 20000
d1 = all[seq(1,splits)]
d2 = all[seq(1+splits,splits*2)]
d3 = all[seq(1+splits*2,splits*3)]
d4 = all[seq(1+splits*3,splits*4)]
d5 = all[seq(1+splits*4,splits*5)]
d6 = all[seq(1+splits*5,splits*6)]
d7 = all[seq(1+splits*6,splits*7)]
d8 = all[seq(1+splits*7,splits*8)]
d9 = all[seq(1+splits*8,splits*9)]
d10 = all[seq(1+splits*9,len-1)]

#d1
vlabel = Label[-d1]
vmrna = mRNA[-d1]
vorf = ORF[-d1]
vfickett = Fickett[-d1]
vhexamer = Hexamer[-d1]

mylogit <- glm(vlabel ~ vmrna + vorf + vfickett + vhexamer, family=binomial(link="logit"), na.action=na.pass)
test <- data.frame(vmrna = mRNA[d1], vorf = ORF[d1],vfickett = Fickett[d1], vhexamer = Hexamer[d1], vlabel=Label[d1])
test$prob <- predict(mylogit,newdata=test,type="response")
output = cbind("mRNA"=test$vmrna, "ORF"=test$vorf, "Fickett"=test$vfickett, "Hexamer" = test$vhexamer, "Label"=test$vlabel,"Prob"=test$prob)
write.table(output,file= paste(output_dir, "test1.xls", sep = "/"),quote=F,sep="\t",row.names=ID[d1])



#d2
vlabel = Label[-d2]
vmrna = mRNA[-d2]
vorf = ORF[-d2]
vfickett = Fickett[-d2]
vhexamer = Hexamer[-d2]

mylogit <- glm(vlabel ~ vmrna + vorf + vfickett + vhexamer, family=binomial(link="logit"), na.action=na.pass)
test <- data.frame(vmrna = mRNA[d2], vorf = ORF[d2],vfickett = Fickett[d2], vhexamer = Hexamer[d2], vlabel=Label[d2])
test$prob <- predict(mylogit,newdata=test,type="response")
output = cbind("mRNA"=test$vmrna, "ORF"=test$vorf, "Fickett"=test$vfickett, "Hexamer" = test$vhexamer, "Label"=test$vlabel,"Prob"=test$prob)
write.table(output,file=paste(output_dir, "test2.xls", sep = "/"),quote=F,sep="\t",row.names=ID[d2])



#d3
vlabel = Label[-d3]
vmrna = mRNA[-d3]
vorf = ORF[-d3]
vfickett = Fickett[-d3]
vhexamer = Hexamer[-d3]

mylogit <- glm(vlabel ~ vmrna + vorf + vfickett + vhexamer, family=binomial(link="logit"), na.action=na.pass)
test <- data.frame(vmrna = mRNA[d3], vorf = ORF[d3],vfickett = Fickett[d3], vhexamer = Hexamer[d3], vlabel=Label[d3])
test$prob <- predict(mylogit,newdata=test,type="response")
output = cbind("mRNA"=test$vmrna, "ORF"=test$vorf, "Fickett"=test$vfickett, "Hexamer" = test$vhexamer, "Label"=test$vlabel,"Prob"=test$prob)
write.table(output,file=paste(output_dir, "test3.xls", sep = "/"),quote=F,sep="\t",row.names=ID[d3])



#d4
vlabel = Label[-d4]
vmrna = mRNA[-d4]
vorf = ORF[-d4]
vfickett = Fickett[-d4]
vhexamer = Hexamer[-d4]

mylogit <- glm(vlabel ~ vmrna + vorf + vfickett + vhexamer, family=binomial(link="logit"), na.action=na.pass)
test <- data.frame(vmrna = mRNA[d4], vorf = ORF[d4],vfickett = Fickett[d4], vhexamer = Hexamer[d4], vlabel=Label[d4])
test$prob <- predict(mylogit,newdata=test,type="response")
output = cbind("mRNA"=test$vmrna, "ORF"=test$vorf, "Fickett"=test$vfickett, "Hexamer" = test$vhexamer, "Label"=test$vlabel,"Prob"=test$prob)
write.table(output,file=paste(output_dir, "test4.xls", sep = "/"),quote=F,sep="\t",row.names=ID[d4])



#d5
vlabel = Label[-d5]
vmrna = mRNA[-d5]
vorf = ORF[-d5]
vfickett = Fickett[-d5]
vhexamer = Hexamer[-d5]

mylogit <- glm(vlabel ~ vmrna + vorf + vfickett + vhexamer, family=binomial(link="logit"), na.action=na.pass)
test <- data.frame(vmrna = mRNA[d5], vorf = ORF[d5],vfickett = Fickett[d5], vhexamer = Hexamer[d5], vlabel=Label[d5])
test$prob <- predict(mylogit,newdata=test,type="response")
output = cbind("mRNA"=test$vmrna, "ORF"=test$vorf, "Fickett"=test$vfickett, "Hexamer" = test$vhexamer, "Label"=test$vlabel,"Prob"=test$prob)
write.table(output,file=paste(output_dir, "test5.xls", sep = "/"),quote=F,sep="\t",row.names=ID[d5])



#d6
vlabel = Label[-d6]
vmrna = mRNA[-d6]
vorf = ORF[-d6]
vfickett = Fickett[-d6]
vhexamer = Hexamer[-d6]

mylogit <- glm(vlabel ~ vmrna + vorf + vfickett + vhexamer, family=binomial(link="logit"), na.action=na.pass)
test <- data.frame(vmrna = mRNA[d6], vorf = ORF[d6],vfickett = Fickett[d6], vhexamer = Hexamer[d6], vlabel=Label[d6])
test$prob <- predict(mylogit,newdata=test,type="response")
output = cbind("mRNA"=test$vmrna, "ORF"=test$vorf, "Fickett"=test$vfickett, "Hexamer" = test$vhexamer, "Label"=test$vlabel,"Prob"=test$prob)
write.table(output,file=paste(output_dir, "test6.xls", sep = "/"),quote=F,sep="\t",row.names=ID[d6])



#d7
vlabel = Label[-d7]
vmrna = mRNA[-d7]
vorf = ORF[-d7]
vfickett = Fickett[-d7]
vhexamer = Hexamer[-d7]

mylogit <- glm(vlabel ~ vmrna + vorf + vfickett + vhexamer, family=binomial(link="logit"), na.action=na.pass)
test <- data.frame(vmrna = mRNA[d7], vorf = ORF[d7],vfickett = Fickett[d7], vhexamer = Hexamer[d7], vlabel=Label[d7])
test$prob <- predict(mylogit,newdata=test,type="response")
output = cbind("mRNA"=test$vmrna, "ORF"=test$vorf, "Fickett"=test$vfickett, "Hexamer" = test$vhexamer, "Label"=test$vlabel,"Prob"=test$prob)
write.table(output,file=paste(output_dir, "test7.xls", sep = "/"),quote=F,sep="\t",row.names=ID[d7])



#d8
vlabel = Label[-d8]
vmrna = mRNA[-d8]
vorf = ORF[-d8]
vfickett = Fickett[-d8]
vhexamer = Hexamer[-d8]

mylogit <- glm(vlabel ~ vmrna + vorf + vfickett + vhexamer, family=binomial(link="logit"), na.action=na.pass)
test <- data.frame(vmrna = mRNA[d8], vorf = ORF[d8],vfickett = Fickett[d8], vhexamer = Hexamer[d8], vlabel=Label[d8])
test$prob <- predict(mylogit,newdata=test,type="response")
output = cbind("mRNA"=test$vmrna, "ORF"=test$vorf, "Fickett"=test$vfickett, "Hexamer" = test$vhexamer, "Label"=test$vlabel,"Prob"=test$prob)
write.table(output,file=paste(output_dir, "test8.xls", sep = "/"),quote=F,sep="\t",row.names=ID[d8])



#d9
vlabel = Label[-d9]
vmrna = mRNA[-d9]
vorf = ORF[-d9]
vfickett = Fickett[-d9]
vhexamer = Hexamer[-d9]

mylogit <- glm(vlabel ~ vmrna + vorf + vfickett + vhexamer, family=binomial(link="logit"), na.action=na.pass)
test <- data.frame(vmrna = mRNA[d9], vorf = ORF[d9],vfickett = Fickett[d9], vhexamer = Hexamer[d9], vlabel=Label[d9])
test$prob <- predict(mylogit,newdata=test,type="response")
output = cbind("mRNA"=test$vmrna, "ORF"=test$vorf, "Fickett"=test$vfickett, "Hexamer" = test$vhexamer, "Label"=test$vlabel,"Prob"=test$prob)
write.table(output,file=paste(output_dir, "test9.xls", sep = "/"),quote=F,sep="\t",row.names=ID[d9])



#d10
vlabel = Label[-d10]
vmrna = mRNA[-d10]
vorf = ORF[-d10]
vfickett = Fickett[-d10]
vhexamer = Hexamer[-d10]

mylogit <- glm(vlabel ~ vmrna + vorf + vfickett + vhexamer, family=binomial(link="logit"), na.action=na.pass)
test <- data.frame(vmrna = mRNA[d10], vorf = ORF[d10],vfickett = Fickett[d10], vhexamer = Hexamer[d10], vlabel=Label[d10])
test$prob <- predict(mylogit,newdata=test,type="response")
output = cbind("mRNA"=test$vmrna, "ORF"=test$vorf, "Fickett"=test$vfickett, "Hexamer" = test$vhexamer, "Label"=test$vlabel,"Prob"=test$prob)
write.table(output,file=paste(output_dir, "test10.xls", sep = "/"),quote=F,sep="\t",row.names=ID[d10])



#ROC
test1=read.table(file=paste(output_dir, "test1.xls", sep = "/"),header=T,sep="\t",quote="", fill=FALSE, row.names=NULL)
test2=read.table(file=paste(output_dir, "test2.xls", sep = "/"),header=T,sep="\t",quote="", fill=FALSE, row.names=NULL)
test3=read.table(file=paste(output_dir, "test3.xls", sep = "/"),header=T,sep="\t",quote="", fill=FALSE, row.names=NULL)
test4=read.table(file=paste(output_dir, "test4.xls", sep = "/"),header=T,sep="\t",quote="", fill=FALSE, row.names=NULL)
test5=read.table(file=paste(output_dir, "test5.xls", sep = "/"),header=T,sep="\t",quote="", fill=FALSE, row.names=NULL)
test6=read.table(file=paste(output_dir, "test6.xls", sep = "/"),header=T,sep="\t",quote="", fill=FALSE, row.names=NULL)
test7=read.table(file=paste(output_dir, "test7.xls", sep = "/"),header=T,sep="\t",quote="", fill=FALSE, row.names=NULL)
test8=read.table(file=paste(output_dir, "test8.xls", sep = "/"),header=T,sep="\t",quote="", fill=FALSE, row.names=NULL)
test9=read.table(file=paste(output_dir, "test9.xls", sep = "/"),header=T,sep="\t",quote="", fill=FALSE, row.names=NULL)
test10=read.table(file=paste(output_dir, "test10.xls", sep = "/"),header=T,sep="\t",quote="", fill=FALSE, row.names=NULL)


Response = list(test1$Prob,test2$Prob,test3$Prob,test4$Prob,test5$Prob,test6$Prob,test7$Prob,test8$Prob,test9$Prob,test10$Prob)
Labls = list(test1$Label,test2$Label,test3$Label,test4$Label,test5$Label,test6$Label,test7$Label,test8$Label,test9$Label,test10$Label)
ROCR_data = list(predictions=Response,Labels=Labls)
pred <- prediction(ROCR_data$predictions, ROCR_data$Labels)
#perf <- performance(pred,"auc")
#avergae AUC = 0.9927


pdf(paste(output_dir, "Figure.pdf", sep = "/"))
par(mfrow=c(2,2),mar=c(5,4,2,2),cex.axis=1.2, cex.lab=1.2)
#ROC curve
#pdf("Human_10fold.ROC.pdf")
perf <- performance(pred,"tpr","fpr")
plot(perf,col="blue",lty=3,xlab="1-Specificity",ylab="Sensitivity",ylim=c(0.7,1),xlim=c(0,0.3),main="",cex.axis=1.5,cex.label=1.5)	#AUC = 0.9927 
plot(perf,lwd=2,avg="vertical",add=TRUE,col="red",xlab="1-specificity",ylab="sensitivity",main="",cex.axis=1.2,cex.label=1.2) 
abline(v=0,lty="dashed",lwd=0.5)
abline(h=1.0,lty="dashed",lwd=0.5)
abline(v=0.05,lty="dashed",lwd=0.5)
abline(h=0.95,lty="dashed",lwd=0.5)
#dev.off()

#precision
#pdf("Human_10fold.precision_vs_recall.pdf")
d=performance(pred,measure="prec", x.measure="rec")
plot(d,col="blue",lty=3,xlab="Recall (TPR)",ylab="Precision (PPV)",xlim=c(0.7,1),ylim=c(0.7,1),cex.axis=1.2,cex.label=1.2)
plot(d,lwd=2,avg="vertical",col="red",xlab="Recall (TPR)",ylab="Precision (PPV)",add=T,cex.axis=1.2,cex.label=1.2)
abline(v=1.0,lty="dashed",lwd=0.5)
abline(h=1.0,lty="dashed",lwd=0.5)
abline(v=0.95,lty="dashed",lwd=0.5)
abline(h=0.95,lty="dashed",lwd=0.5)
#dev.off()


#Accuracy
#pdf("Human_10fold.Accuracy.pdf")
perf <- performance(pred,"acc")
plot(perf,col="blue",lty=3,xlab="Coding probability cutoff",ylab="Accuracy",ylim=c(0.7,1),cex.axis=1.2,cex.label=1.2) 
plot(perf,lwd=2,avg="vertical",add=TRUE,col="red",cex.axis=1.2,cex.label=1.2) 
abline(h=1,lty="dashed",lwd=0.5)
abline(h=0.95,lty="dashed",lwd=0.5)
#dev.off()


#sensitivity vs specificity
pred <- prediction(ROCR_data$predictions, ROCR_data$Labels)
S <- performance(pred,measure="sens")
P <- performance(pred,measure="spec")


#both.eq <- which.min(abs(S@y.values[[1]]-P@y.values[[1]]))
#best.sum <- which.max(S@y.values[[1]]+P@y.values[[1]])
#cpc_cutoff_min <- S@x.values[[1]][both.eq]
#cpc_cutoff_max <- S@x.values[[1]][best.sum]
#tmp_results <- (S@y.values-P@y.values)
#print(matrix(unlist(S@y.values),ncol=10))

roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
opt.cut = function(perf, pred){
    cut.ind = mapply(FUN=function(x, y, p){
        d = (x - 0)^2 + (y-1)^2
        ind = which(d == min(d))
        c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
            cutoff = p[[ind]])
    }, perf@x.values, perf@y.values, pred@cutoffs)
}
av_cutoff <- mean(opt.cut(roc.perf, pred)["cutoff",])
av_acc <- mean(opt.cut(roc.perf, pred)[c("sensitivity","specificity"),])


#pdf("Human_10fold_sens_vs_spec.pdf")
plot(S,col="blue",lty=3,ylab="Performance",xlab="Coding Probability Cutoff",ylim=c(0.8,1),cex.axis=1.2,cex.label=1.2) 
plot(S,lwd=2,avg="vertical",add=TRUE,col="blue") 
plot(P,col="red",lty=3, add=TRUE,) 
plot(P,lwd=2,avg="vertical",add=TRUE,col="red") 

print(paste("Cutoff:", av_cutoff))
abline(h=av_acc,lty="dashed",lwd=0.5)
abline(v=av_cutoff,lty="dashed",lwd=0.5)
#abline(h=best.sum,lty="dashed",lwd=0.5)
#abline(v=cpc_cutoff_max,lty="dashed",lwd=0.5)

legend(0.4,0.85,col=c("blue","red"),lwd=2,legend=c("Sensitivity","Specificity"))

