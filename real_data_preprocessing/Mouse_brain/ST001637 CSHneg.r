folder<-"D:\\ST001637_Brain atlas files_CSHneg"
setwd(folder)
library(apLCMS)
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

data(metabolite.table)
data(adduct.table)
known.table.common.neg<-make.known.table(metabolite.table, adduct.table[c(34:39,44),])
known.table.common.neg<-known.table.common.neg[known.table.common.neg[,6]>10,]

file.pattern=".mzML"
n.nodes=32
mz.tol=1e-5
baseline.correct=0
min.bw=NA
max.bw=NA
shape.model="bi-Gaussian"
align.mz.tol=1e-5
align.chr.tol=NA
pre.process=FALSE
recover.mz.range=NA
recover.chr.range=NA
use.observed.range=TRUE
sd.cut=c(0.125,60)
sigma.ratio.lim=c(0.2,5)
peak.estim.method="moment"
max.align.mz.diff=0.015
match.tol.ppm=NA
new.feature.min.count=300
recover.min.count=4
baseline.correct.noise.percentile=0
component.eliminate = 0.01
moment.power = 1

###########################


combos<-matrix(c(5,10,0.1,0.15),ncol=2)

for(mz.tol in c(1e-4, 5e-5, 5e-4))
{
    align.mz.tol=mz.tol
for(moment.power in c(1))
{
for(combos.i in c(1,2))
{
min.run<-combos[combos.i,1]
min.pres<-combos[combos.i,2]

for(min.exp in round(length(dir(pattern=file.pattern))*c(0.2,0.4, 0.6,0.8)))
    {
        this.name<-paste("negative HMDB", mz.tol, min.run, min.pres, moment.power, min.exp, ".bin")
        is.there<-dir(pattern=this.name)
		
        if(length(is.there)>0)
        {
            message(this.name, " done.")
        }else{
                   b<-semi.sup(folder,n.nodes=n.nodes, file.pattern=file.pattern,known.table=known.table.common.neg,  sd.cut=sd.cut,sigma.ratio.lim=sigma.ratio.lim, min.pres=min.pres, min.run=min.run, min.exp=min.exp, mz.tol=mz.tol, baseline.correct=baseline.correct, baseline.correct.noise.percentile=baseline.correct.noise.percentile, align.mz.tol=align.mz.tol, align.chr.tol=align.chr.tol, max.align.mz.diff=max.align.mz.diff, recover.mz.range=recover.mz.range, recover.chr.range=recover.chr.range,use.observed.range=use.observed.range,moment.power=moment.power, shape.model=shape.model,new.feature.min.count=new.feature.min.count, recover.min.count=recover.min.count)
        
                g<-b$final.ftrs
                save(g, file=this.name)
        }
    }
  }
}
}

