##Generate gene expression data for number of "n_id" samples with mean level of precursor to be "prem", mature RNA to be "truem"

gen_exp=function(n1=NULL, prem_g1=NULL, truem_g1=0, n2=NULL, prem_g2=NULL, truem_g2=0){

	if (is.null(prem_g1) || is.null(prem_g2) || is.null(n1) || is.null(n2)){
		print("ERR!")
		return()
	}

	n=n1+n2
	v_id=5
	v_err=2

	n_prob_intron=5
	n_prob_exon=3

	re_id=rnorm(n,mean=0,sd=sqrt(v_id))
	re_err=rnorm(n*(n_prob_intron+n_prob_exon), mean=0, sd=sqrt(v_err))
	
	##Group 1
	re_err_intron=re_err[1:(n1*n_prob_intron)]
	re_id_g1=re_id[1:n1]

	idx=0
	prb_intron=matrix(0,nrow=n1,ncol=n_prob_intron)
	for (i in 1:n1){
		for (j in 1:n_prob_intron){
			idx=idx+1
			prb_intron[i,j]=prem_g1+re_id_g1[i]+re_err_intron[idx]
		}
	}

	dimnames(prb_intron) <- list(rownames(prb_intron, do.NULL = FALSE, prefix = "id_"),
                          colnames(prb_intron, do.NULL = FALSE, prefix = "prob_intron"))

	d_intron=data.frame(prb_intron)

	re_err_exon=re_err[(1+n1*n_prob_intron):(n1*(n_prob_intron+n_prob_exon))]
	idx=0
	prb_exon=matrix(0,nrow=n1,ncol=n_prob_exon)
	for (i in 1:n1){
		for (j in 1:n_prob_exon){
			idx=idx+1
			prb_exon[i,j]=truem_g1+prem_g1+re_id_g1[i]+re_err_exon[idx]
		}
	}

	dimnames(prb_exon) <- list(rownames(prb_exon, do.NULL = FALSE, prefix = "id_"),
                          colnames(prb_exon, do.NULL = FALSE, prefix = "prob_exon"))

	gexp=rep(0, n1*n_prob_intron)
	id_x=rep("", n1*n_prob_intron)
	prb_x=rep("", n1*n_prob_intron)
	exon_t=rep(0, n1*n_prob_intron)
	group=rep(1, n1*n_prob_intron)

	for(t in 1:n1){

		for (tt in 1:n_prob_intron){
			gexp[(t-1)*n_prob_intron+tt]=prb_intron[t,tt]
			id_x[(t-1)*n_prob_intron+tt]=rownames(prb_intron)[t]
			prb_x[(t-1)*n_prob_intron+tt]=colnames(prb_intron)[tt]
		}

	}

	dd=data.frame(gexp,id_x,prb_x,exon_t,group)

	gexp=rep(0, n1*n_prob_exon)
	id_x=rep("", n1*n_prob_exon)
	prb_x=rep("", n1*n_prob_exon)
	exon_t=rep(1, n1*n_prob_exon)
	group=rep(1, n1*n_prob_exon)

	for(t in 1:n1){

		for (tt in 1:n_prob_exon){
			gexp[(t-1)*n_prob_exon+tt]=prb_exon[t,tt]
			id_x[(t-1)*n_prob_exon+tt]=rownames(prb_exon)[t]
			prb_x[(t-1)*n_prob_exon+tt]=colnames(prb_exon)[tt]
		}

	}

	dd_exon=data.frame(gexp,id_x,prb_x,exon_t,group)

	dd_g1=rbind(dd,dd_exon)

	##Group 2
	re_err_intron=re_err[(1+n1*(n_prob_intron+n_prob_exon)):(n1*(n_prob_intron+n_prob_exon)+n2*n_prob_intron)]
	re_id_g2=re_id[(1+n1):n]

	idx=0
	prb_intron=matrix(0,nrow=n2,ncol=n_prob_intron)
	for (i in 1:n2){
		for (j in 1:n_prob_intron){
			idx=idx+1
			prb_intron[i,j]=prem_g2+re_id_g2[i]+re_err_intron[idx]
		}
	}

	dimnames(prb_intron) <- list(rownames(prb_intron, do.NULL = FALSE, prefix = "id_"),
                          colnames(prb_intron, do.NULL = FALSE, prefix = "prob_intron"))

	d_intron=data.frame(prb_intron)

	re_err_exon=re_err[(1+n1*(n_prob_intron+n_prob_exon)+n2*n_prob_intron):(n*(n_prob_intron+n_prob_exon))]

	idx=0
	prb_exon=matrix(0,nrow=n2,ncol=n_prob_exon)
	for (i in 1:n2){
		for (j in 1:n_prob_exon){
			idx=idx+1
			prb_exon[i,j]=truem_g2+prem_g2+re_id_g2[i]+re_err_exon[idx]
		}
	}

	dimnames(prb_exon) <- list(rownames(prb_exon, do.NULL = FALSE, prefix = "id_"),
                          colnames(prb_exon, do.NULL = FALSE, prefix = "prob_exon"))

	gexp=rep(0, n2*n_prob_intron)
	id_x=rep("", n2*n_prob_intron)
	prb_x=rep("", n2*n_prob_intron)
	exon_t=rep(0, n2*n_prob_intron)
	group=rep(2, n2*n_prob_intron)

	for(t in 1:n2){

		for (tt in 1:n_prob_intron){
			gexp[(t-1)*n_prob_intron+tt]=prb_intron[t,tt]
			id_x[(t-1)*n_prob_intron+tt]=paste("id_", t+n1, sep="")
			prb_x[(t-1)*n_prob_intron+tt]=colnames(prb_intron)[tt]
		}

	}

	dd=data.frame(gexp,id_x,prb_x,exon_t,group)

	gexp=rep(0, n2*n_prob_exon)
	id_x=rep("", n2*n_prob_exon)
	prb_x=rep("", n2*n_prob_exon)
	exon_t=rep(1, n2*n_prob_exon)
	group=rep(2, n2*n_prob_exon)

	for(t in 1:n2){

		for (tt in 1:n_prob_exon){
			gexp[(t-1)*n_prob_exon+tt]=prb_exon[t,tt]
			id_x[(t-1)*n_prob_exon+tt]=paste("id_", t+n1, sep="")
			prb_x[(t-1)*n_prob_exon+tt]=colnames(prb_exon)[tt]
		}

	}

	dd_exon=data.frame(gexp,id_x,prb_x,exon_t,group)

	dd=rbind(dd_g1,dd,dd_exon)

	return(dd)

}
