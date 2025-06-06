# read Experiment 2 data - dots #

# load tab delimited data file
x=read.delim("Experiment_2_data_dots.txt")

# recode a few things
x$RT=x$tlag/65  # time coded in 'clicks' of the computer
names(x)[1]="subj"
x$subj=factor(x$subj)

# estimate cumulative normal link function for the drift rates for option 1 (<50 dots) vs 2 (>50 dots)

# re-code some features for easier modelling
x$R=3-x$choice

# fix a few features
# remove very fast and very slow trials
x$badtrial=x$badsubj=FALSE
x$badtrial[x$RT<.5 | x$RT>20]=TRUE

# exclude subjects with very low accuracy
for(s in levels(x$subj)) {
  tmp=subset(x,subj==s)
  # exclude subject if accuracy too low
  if(mean(tmp$cor)<.55) x$badsubj[x$subj==s]=TRUE
}

originaldata=x

dots.data=subset(x,!(badtrial |badsubj))
# re-factor subject column
dots.data$subj=factor(as.character(dots.data$subj))
