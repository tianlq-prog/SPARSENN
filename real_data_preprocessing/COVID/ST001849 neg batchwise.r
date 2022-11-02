library(apLCMS)

folder<-"D:\\new_met_data\\ST001849_polar_neg"
setwd(folder)

data(metabolite.table)
data(adduct.table)
known.table.common.neg<-make.known.table(metabolite.table, adduct.table[c(34:39,44),])
known.table.common.neg<-known.table.common.neg[known.table.common.neg[,6]>10,]

file.pattern=".mzML"
n.nodes=60
baseline.correct=0
baseline.correct.noise.percentile=0.05
min.bw=NA
max.bw=NA
shape.model="bi-Gaussian"
mz.tol=1e-5
align.mz.tol=1e-5
batch.align.mz.tol=1e-5
align.chr.tol=100
batch.align.chr.tol=150
pre.process=FALSE
recover.mz.range=NA
recover.chr.range=NA
use.observed.range=TRUE
sd.cut=c(0.25,60)
sigma.ratio.lim=c(0.2,5)
peak.estim.method="moment"
max.align.mz.diff=0.015
match.tol.ppm=NA
new.feature.min.count=4
recover.min.count=4
subs=NA
component.eliminate=0.01
moment.power=1


files<-dir(pattern=".mzML")
load("polar_neg_batch_info.Rdata")
info<-a[,2:3]
info<-info[order(info[,1]),]
sum(files == info[,1])

###########################



combos<-matrix(c(8,12,0.4,0.5),ncol=2)

for(moment.power in c(2))
{
	for(combos.i in c(2,1))
	{
		min.run<-combos[combos.i,1]
		min.pres<-combos[combos.i,2]

		for(min.within.batch.prop.report in c(0.3, 0.5))
		{

		  for(min.within.batch.prop.detect in c(min.within.batch.prop.report/1.5))#, min.within.batch.prop.detect+(1-min.within.batch.prop.detect)/2))
		  {
		  
			for(min.batch.prop in c(0.3, 0.5))
			{
				this.name<-paste("negative batchwise HMDB", min.run, min.pres, moment.power, min.within.batch.prop.detect, min.within.batch.prop.report, min.batch.prop, "new.bin")
				is.there<-dir(pattern=this.name)
				if(length(is.there)>0)
				{
					message(this.name, " done.")
				}else{
					b<-two.step.hybrid(folder=folder, info=info,  min.within.batch.prop.detect = min.within.batch.prop.detect, min.within.batch.prop.report =min.within.batch.prop.report, min.batch.prop = min.batch.prop, batch.align.mz.tol=batch.align.mz.tol, batch.align.chr.tol=batch.align.chr.tol,n.nodes=n.nodes, file.pattern=file.pattern,known.table=known.table.common.neg, sd.cut=sd.cut,sigma.ratio.lim=sigma.ratio.lim, component.eliminate=component.eliminate, moment.power=moment.power, min.pres=min.pres, min.run=min.run, mz.tol=mz.tol, baseline.correct.noise.percentile=baseline.correct.noise.percentile, align.mz.tol=align.mz.tol, align.chr.tol=align.chr.tol, max.align.mz.diff=max.align.mz.diff, recover.mz.range=recover.mz.range, recover.chr.range=recover.chr.range,use.observed.range=use.observed.range, shape.model=shape.model,new.feature.min.count=new.feature.min.count, recover.min.count=recover.min.count)

						g<-b$final.ftrs
						save(g, file=this.name)
				}
			}
		  }
		}
	}
}
