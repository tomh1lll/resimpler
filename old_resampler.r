
reps = 1
count = 1
stat_resample = data.frame(matrix(ncol = 9, nrow = 0))
for (x in c(1:1000)){
for (i in levels(innSubset$cat2)){
	ss = droplevels(innSubset[innSubset$cat2 == i,])
	fst_v = c()
	tp_v = c()
	td_v = c()
	alpha_v = c()
	dos_v = c()
	snipre_v = c()
	snipreb_v = c()
	for (j in ss$gene){
		chr = as.character(droplevels(ss$chr[ss$gene == j]))
		start = ss$start[ss$gene == j]
		end = ss$end[ss$gene == j]
		size = ss$end[ss$gene == j]-ss$start[ss$gene == j]
		ssO = droplevels(innOther[innOther$chr == chr & innOther$start > start -100000 & innOther$end < end + 100000,])
		ssO = na.omit(ssO)
		if (length(ssO$chr) > 0){
			fst_v = c(fst_v,ss$Fst1[ss$gene == j] - as.numeric(sample(ssO$Fst1,size = 1)))
			tp_v = c(tp_v,ss$tP[ss$gene == j] - as.numeric(sample(ssO$tP,size = 1)))
			td_v = c(td_v,ss$tD[ss$gene == j] - as.numeric(sample(ssO$tD,size = 1)))
			alpha_v = c(alpha_v,ss$alpha[ss$gene == j] - as.numeric(sample(ssO$alpha,size = 1)))
			dos_v = c(dos_v,ss$DoS[ss$gene == j] - as.numeric(sample(ssO$DoS,size = 1)))
			snipre_v = c(snipre_v,ss$snipre[ss$gene == j] - as.numeric(sample(ssO$snipre,size = 1)))
			snipreb_v = c(snipreb_v,ss$snipre.b[ss$gene == j] - as.numeric(sample(ssO$snipre.b,size = 1)))
		}
	}
	cat_rep = c(as.character(reps), as.character(i), as.character(mean(fst_v)) ,as.character(mean(tp_v)) ,as.character(mean(td_v)) ,as.character(mean(alpha_v)) ,as.character(mean(dos_v)) ,as.character(mean(snipre_v)) ,as.character(mean(snipreb_v)))
	stat_resample[count,] = cat_rep
	count = count + 1
}
reps = reps + 1
if (reps%%10 == 0){cat(reps,"... ")}
}

innS = read.delim("file:///E:/DinnPopGen/Dinn_total_stats.txt")
innSubset = droplevels(innS[innS$cat2 != "Other",])
innOther = droplevels(innS[innS$cat2 == "Other",])

innSubset2 = innSubset
for (i in levels(innSubset$cat2)){
	ss = droplevels(innSubset[innSubset$cat2 == i,])
	for (j in ss$gene){
		chr = as.character(droplevels(ss$chr[ss$gene == j]))
		start = ss$start[ss$gene == j]
		end = ss$end[ss$gene == j]
		size = ss$end[ss$gene == j]-ss$start[ss$gene == j]
		ssO = droplevels(innOther[innOther$chr == chr & innOther$start > start -100000 & innOther$end < end + 100000,])
		innSubset2 = rbind(innSubset2,ssO)
	}
}

count = 1
stat_resample_bypop = data.frame(matrix(ncol = 7, nrow = 0))
for (p in c("CH")){
	reps = 1
	Dinn_stat_bygene4_ip = droplevels(Dinn_stat_bygene4_i[Dinn_stat_bygene4_i$pop == p,])
	Dinn_stat_bygene4_op = droplevels(Dinn_stat_bygene4_o[Dinn_stat_bygene4_o$pop == p,])
	for (x in c(1:1000)){
	for (i in levels(Dinn_stat_bygene4_ip$cat2)){
		ss = droplevels(Dinn_stat_bygene4_ip[Dinn_stat_bygene4_ip$cat2 == i,])
		fst_1 = c()
		fst_2 = c()
		tp_v = c()
		td_v = c()
		for (j in ss$gene){
			chr = as.character(droplevels(ss$chr[ss$gene == j]))
			start = ss$start[ss$gene == j]
			end = ss$end[ss$gene == j]
			ssO = droplevels(Dinn_stat_bygene4_op[Dinn_stat_bygene4_op$chr == chr & Dinn_stat_bygene4_op$start > start -100000 & Dinn_stat_bygene4_op$end < end + 100000,])
			ssO = na.omit(ssO)
			if (length(ssO$chr) > 0){
				fst_1 = c(fst_1,ss$Fst1[ss$gene == j] - as.numeric(sample(ssO$Fst1,size = 1)))
				fst_2 = c(fst_2,ss$Fst2[ss$gene == j] - as.numeric(sample(ssO$Fst2,size = 1)))
				tp_v = c(tp_v,ss$tP[ss$gene == j] - as.numeric(sample(ssO$tP,size = 1)))
				td_v = c(td_v,ss$tD[ss$gene == j] - as.numeric(sample(ssO$tD,size = 1)))
			}
		}
	cat_rep = c(as.character(reps), as.character(i), as.character(p),as.character(mean(fst_1)),as.character(mean(fst_2)) ,as.character(mean(tp_v)) ,as.character(mean(td_v)))
	stat_resample_bypop[count,] = cat_rep
	count = count + 1
}
reps = reps + 1
if (reps%%10 == 0){cat(reps,"... ")}
}
}




count = 1
ptm <- proc.time()
stat_resample_MKv = data.frame(matrix(ncol = 7, nrow = 0))
for (v in c("i","m")){
	reps = 1
	imm_subs = droplevels(totalMK[totalMK$species == v & totalMK$Category != "Other",])
	nor_subs = droplevels(totalMK[totalMK$species == v & totalMK$Category == "Other",])
	for (x in c(1:1000)){
	for (i in levels(imm_subs$Category)){
		ss = droplevels(imm_subs[imm_subs$Category == i,])
		alpha_v = c()
		dos_v = c()
		snipre_v = c()
		snipreb_v = c()
		for (j in ss$GeneID){
			chr = as.character(droplevels(ss$chr[ss$GeneID == j]))
			start = ss$start[ss$GeneID == j]
			end = ss$end[ss$GeneID == j]
			size = ss$end[ss$GeneID == j]-ss$start[ss$GeneID == j]
			ssO = droplevels(nor_subs[nor_subs$chr == chr & nor_subs$start > start -100000 & nor_subs$end < end + 100000,])
			ssO = na.omit(ssO)
			if (length(ssO$chr) > 0){
				alpha_v = c(alpha_v,ss$basic_alpha[ss$GeneID == j] - as.numeric(sample(ssO$basic_alpha,size = 1)))
				dos_v = c(dos_v,ss$basic_dos[ss$GeneID == j] - as.numeric(sample(ssO$basic_dos,size = 1)))
				snipre_v = c(snipre_v,ss$SnIPRE.est[ss$GeneID == j] - as.numeric(sample(ssO$SnIPRE.est,size = 1)))
				snipreb_v = c(snipreb_v,ss$SnIPRE.Rest[ss$GeneID == j] - as.numeric(sample(ssO$SnIPRE.Rest,size = 1)))
			}
		}
		cat_rep = c(as.character(v),as.character(reps), as.character(i) ,as.character(mean(alpha_v)) ,as.character(mean(dos_v)) ,as.character(mean(snipre_v)) ,as.character(mean(snipreb_v)))
		stat_resample_MKv[count,] = cat_rep
		count = count + 1
	}
	reps = reps + 1
	if (reps%%100 == 0){cat(reps,"... ")}
	}
}
proc.time() - ptm
View(stat_resample_MKv)



count = 1
ptm <- proc.time()
stat_resample_MKv = data.frame(matrix(ncol = 7, nrow = 0))
for (v in c("m")){
	reps = 1
	imm_subs = droplevels(totalMK[totalMK$species == v & totalMK$Category != "Other",])
	nor_subs = droplevels(totalMK[totalMK$species == v & totalMK$Category == "Other",])
	for (x in c(1:1000)){
	for (i in c("IMD")){
		ss = droplevels(imm_subs[imm_subs$Category == i,])
		alpha_v = c()
		dos_v = c()
		snipre_v = c()
		snipreb_v = c()
		for (j in ss$GeneID){
			chr = as.character(droplevels(ss$chr[ss$GeneID == j]))
			start = ss$start[ss$GeneID == j]
			end = ss$end[ss$GeneID == j]
			size = ss$end[ss$GeneID == j]-ss$start[ss$GeneID == j]
			ssO = droplevels(nor_subs[nor_subs$chr == chr & nor_subs$start > start -100000 & nor_subs$end < end + 100000,])
			ssO = na.omit(ssO)
			if (length(ssO$chr) > 0){
				alpha_v = c(alpha_v,ss$basic_alpha[ss$GeneID == j] - as.numeric(sample(ssO$basic_alpha,size = 1)))
				dos_v = c(dos_v,ss$basic_dos[ss$GeneID == j] - as.numeric(sample(ssO$basic_dos,size = 1)))
				snipre_v = c(snipre_v,ss$SnIPRE.est[ss$GeneID == j] - as.numeric(sample(ssO$SnIPRE.est,size = 1)))
				snipreb_v = c(snipreb_v,ss$SnIPRE.Rest[ss$GeneID == j] - as.numeric(sample(ssO$SnIPRE.Rest,size = 1)))
			}
		}
		cat_rep = c(as.character(v),as.character(reps), as.character(i) ,as.character(mean(alpha_v)) ,as.character(mean(dos_v)) ,as.character(mean(snipre_v)) ,as.character(mean(snipreb_v)))
		stat_resample_MKv[count,] = cat_rep
		count = count + 1
	}
	reps = reps + 1
	if (reps%%100 == 0){cat(reps,"... ")}
	}
}
proc.time() - ptm


args=commandArgs(trailingOnly = TRUE)
#args="/GoogleDrive/UncklessLab/AMPBalancingSelection_NEW/data/mauritianaPool/Dmau_merged.fixed.txt"
stats<-read.table(args[1],header=FALSE)

stats$pi<-as.numeric(stats$pi)
stats$theta<-as.numeric(stats$theta)
stats$D<-as.numeric(stats$D)

names(stats)<-c("chr","start","end","gene","length","gene2","pi","theta","D")

amps=read.table("/GoogleDrive/UncklessLab/AMPBalancingSelection_NEW/data/mauritianaPool/Dmau.AMP.list",header=FALSE)$V1

#amps2<-amps[-which(amps %in% c("CG18107-PA_cds", "IM1-PA_cds",     "CG15068-PA_cds", "IM2-PA_cds",
                   "CG15067-PA_cds","IM23-PA_cds", "IM3-PA_cds","CG16836-PA_cds","CG15065-PA_cds"))]

#amps<-amps2
AMP.data<-stats[which(stats$gene %in% amps),]

cat(length(AMP.data$gene)," AMPs before filtering\n",sep="")

AMP.data=AMP.data[!is.na(AMP.data$D),]

cat(length(AMP.data$gene)," AMPs after filtering\n",sep="")

###Window
win=100000
###size dif
size=10
###runs
runs = 1000

Dlist=c()
Plist=c()
Tlist=c()

for(r in 1:runs) {
if(r %in% seq(0,runs,runs/10)) {
  cat(r,"...",sep="")
}
D.diff=c()
Pi.diff=c()
Theta.diff=c()

for(i in 1:length(AMP.data$gene)) {
  test=stats[intersect(intersect(which(stats$chr==AMP.data$chr[i]),which(stats$start>AMP.data$start[i]-win)),which(stats$start<AMP.data$start[i]+win)),]
  flag=0
  while(flag==0) {
    pick=test[sample(1:dim(test)[[1]],1),]
    if(pick$length[1]<size*AMP.data$length[i]) {
        flag=1
    }
    if(pick$gene[1] == AMP.data$gene[i] || is.na(pick$D[1])) {
        flag=0
    }
  }
  D.diff=append(D.diff,AMP.data$D[i]-pick$D[1])
  Pi.diff=append(Pi.diff,AMP.data$pi[i]-pick$pi[1])
  Theta.diff=append(Theta.diff,AMP.data$theta[i]-pick$theta[1])
}

Dlist=append(Dlist,mean(D.diff))
Plist=append(Plist,mean(Pi.diff))
Tlist=append(Tlist,mean(Theta.diff))

}

cat("\n")

df=data.frame(rep=1:r,D=Dlist,Pi=Plist,Theta=Tlist)
write.table(df,file=paste(args[1],".diffs",sep=""),quote=FALSE,row.names=FALSE)

hist(df$D)
hist(df$Pi)
hist(df$Theta)
