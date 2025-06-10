# read Experiment 1 data #

# load tab delimited data file
x=read.delim("Experiment_1_data.txt")

# recode a few things
x$RT=x$tlag/65  # time coded in 'clicks' of the computer
names(x)[1]="subj"
x$subj=factor(x$subj)

# recode base price and overage as numeric (without $ sign)
x$b1=as.numeric(sapply(strsplit(as.character(x$b1),"\\$"),function(y) y[2]))
x$b2=as.numeric(sapply(strsplit(as.character(x$b2),"\\$"),function(y) y[2]))
x$b3=as.numeric(sapply(strsplit(as.character(x$b3),"\\$"),function(y) y[2]))
x$c1=as.numeric(sapply(strsplit(as.character(x$c1),"\\$"),function(y) y[2]))
x$c2=as.numeric(sapply(strsplit(as.character(x$c2),"\\$"),function(y) y[2]))
x$c3=as.numeric(sapply(strsplit(as.character(x$c3),"\\$"),function(y) y[2]))
# recode usage to reflect minutes OVER the cutoff price (500)
x$use=x$use-500
x$use[x$use<0]=0

# generate objective value for each option:
#    base price + charge for cents/minute over 500 minutes * minutes over 500
x$S1=x$b1+(x$c1*x$use)
x$S2=x$b2+(x$c2*x$use)
x$S3=x$b3+(x$c3*x$use)

# fix a few features
# remove very slow trials
x$badtrial=x$badsubj=FALSE
x$badtrial[x$RT>30]=TRUE

# subjects 868, 996, 3188 are missing data in one round - remove
x$badsubj[x$subj=="868"]=TRUE
x$badsubj[x$subj=="996"]=TRUE
x$badsubj[x$subj=="3188"]=TRUE

# exclude subjects with accuracy below chance given number of trials they completed
#   justified because each participant completed different n so fixed cutoff not appropriate
for(s in levels(x$subj)) {
  tmp1=subset(x,subj==s & rnd<=2) ; tmp2=subset(x,subj==s & rnd>2)
  # exclude subject if accuracy too low
  # check both binomial density. some subjects were considerably below chance so had p<.05
  #   density for under chance - insert manual cut as well
  # separately for each phase
  if((dbinom(sum(tmp1$cor),nrow(tmp1),.33)>.05) | (dbinom(sum(tmp2$cor),nrow(tmp2),.33)>.05) |
    (mean(tmp1$cor)<.33) | (mean(tmp2$cor)<.33)) x$badsubj[x$subj==s]=TRUE
}

originaldata=x

# removing slow and fast responses
#   exclude trial if the difference in RT between the slowest and second slowest RT is x times larger
#   (same for fastest and second fastest).  than the difference between the second and
#   third slowest RTs. continue until exclusion condition==FALSE. do this recursively for the last
#   (first) sorted RTs over window.for.tail (window.for.leading.edge). x is defined as
#   abs.val.upper (abs.val.lower)

abs.val.lower=.8  # for leading edge of the distribution
abs.val.upper=5  # for tail of the distribution
window.for.leading.edge=5  # number of trials minus 1 indicating the window over which to test for fast outliers
window.for.tail=2 # number of trials minus 1 indicating the window over which to test for slow outliers

for(s in levels(x$subj)) {
  tmp=subset(x,subj==s)
  # initialise cutoff values at 0 and Inf (keep all trials from subject)
  cutoff.tail.rt=Inf
  cutoff.leading.edge.rt=0
  sorted.rt=sort(tmp$RT)
  diff.rt=diff(sorted.rt)
  sorted.n=length(sorted.rt)
  diff.n=length(diff.rt)
  tail.rt=diff.rt[(diff.n-window.for.tail):diff.n]>abs.val.upper
  if(any(tail.rt)) {
    cutoff.tail.rt=sorted.rt[sorted.n-((window.for.tail+1)-min(which(tail.rt)))]
  }
  leading.edge.rt=diff.rt[1:window.for.leading.edge]>abs.val.lower
  if(any(leading.edge.rt)) {
    cutoff.leading.edge.rt=sorted.rt[max(which(leading.edge.rt))]
  }
  tmp$badtrial[tmp$RT<=cutoff.leading.edge.rt | tmp$RT>=cutoff.tail.rt]=TRUE
  x[x$subj==s,]=tmp
}

data=subset(x,!(badtrial |badsubj))
# re-factor subject column
data$subj=factor(as.character(data$subj))
