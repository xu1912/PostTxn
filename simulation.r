library("lme4")
library("lmerTest")
source("gen_exp.R")

n_1=25
n_2=25
prem_g1=1
prem_g2=2
truem_g1=1
truem_g2=0
k=5000

	n_prob_intron=3
	n_prob_exon=5
tr_m=matrix(0,nrow=k, ncol = n_prob_intron + n_prob_exon)
tr_m2=matrix(0,nrow=k, ncol = 4)

for (i in 1:k){

	dataT=gen_exp(n_1,prem_g1,truem_g1, n_2, prem_g2, truem_g2)
	dataT$group=as.factor(dataT$group)
	dataT$exon_t=as.factor(dataT$exon_t)

	d1=dataT[which(dataT$group==1),]
      d2=dataT[which(dataT$group==2),]

      tr_x=0
      for (j in levels(dataT$prb_x)){
		tr_x=tr_x+1
		tr_m[i,tr_x]=sum(t.test(d1[which(d1$prb_x==j),1], d2[which(d2$prb_x==j),1])$p.value<0.05)

		#dataTM=dataT[which(dataT$prb_x==j),]
		#dresult=lm(gexp~group,data=dataTM)
      }

	ddresult=lmer(gexp~group*exon_t+(1|id_x),data=dataT)
	summary(ddresult)
	tr_m2[i, 1]=sum(summary(ddresult)$coefficients[1,5]<0.05)
	tr_m2[i, 2]=sum(summary(ddresult)$coefficients[2,5]<0.05)
	tr_m2[i, 3]=sum(summary(ddresult)$coefficients[3,5]<0.05)
	tr_m2[i, 4]=sum(summary(ddresult)$coefficients[4,5]<0.05)

}

colMeans(tr_m)
colMeans(tr_m2)
