# prepare data to fit LBA to Exp1 data #
source("read-data-Exp1.r")

# work with smaller data file throughout fitting
data$phase <- factor((data$rnd>2)+1)
data$cond <- factor(data$cond)
  levels(data$cond) <- c("FF","FS","SS","SF")
originaldata <- data
# first pass: use objective value of each option (base price + charge for cents/minute over 500 minutes * minutes over 500)
data <- data[,c("subj","cond","phase","rnd","S1","S2","S3","choice","cor","RT")]
  names(data)[8:9] <- c("R","correct")

# final exclusion: remove subjects who completed fewer than 5 trials in a phase (rounds 1/2 or 3/4)
keep.subj <- with(data, tapply(phase, subj, function(x) { all(table(x) > 5) }))
data <- data[(data$subj %in% names(keep.subj)[keep.subj]),]
data$subj <- droplevels(data$subj)

# split data by SAT condition
all.data <- list()
all.data$FF <- subset(data, cond == "FF")
all.data$FS <- subset(data, cond == "FS")
all.data$SF <- subset(data, cond == "SF")
all.data$SS <- subset(data, cond == "SS")

for(i in 1:length(all.data)) {
  # remove unused factor levels and then reorder
  all.data[[i]]$subj <- droplevels(all.data[[i]]$subj)
  all.data[[i]]$subject <- NA
  s.count <- 1
  all.data[[i]]$subject[1] <- s.count
  for(j in 2:nrow(all.data[[i]])) {
    if(all.data[[i]]$subj[j] != all.data[[i]]$subj[j-1]) {
      s.count <- s.count +1
    }
    all.data[[i]]$subject[j] <- s.count
  }
  all.data[[i]]$subject <- factor(all.data[[i]]$subject)
}

# assume contaminant mixture process
# ps[1]=% assumed contamination, ps[2]=RT range in data
p.contaminant <- .05
