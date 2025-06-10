# prepare data to fit LBA to Exp2 data #
source("read-data-Exp2-pricing.r")
source("read-data-Exp2-dots.r")

# final prep for store.data - smaller file for fitting
store.data$phase <- factor((store.data$rnd>2)+1)
store.data$cond <- factor(store.data$cond)
  levels(store.data$cond)=c("FS","SS")
original.store.data <- store.data
store.data$task <- "store"
store.data=store.data[,c("subj","task","cond","phase","rnd","sumstore1","sumstore2","sumstore3","sumstore4","choice","cor","RT")]
  names(store.data)[6:11]=c("S1","S2","S3","S4","R","correct")

# final exclusion: remove subjects who completed fewer than 5 trials in a phase (rounds 1/2 or 3/4)
keep.subj <- with(store.data, tapply(phase, subj, function(x) { all(table(x) > 5) }))
store.data <- store.data[(store.data$subj %in% names(keep.subj)[keep.subj]),]
store.data$subj <- droplevels(store.data$subj)

# and some final prep for dots.data
dots.data$cond <- factor(dots.data$cond)
  levels(dots.data$cond)=c("FS","SS")
dots.data$task <- "dots"
dots.data=dots.data[,c("subj","task","cond","numdot","R","cor","RT")]
  names(dots.data)[6]="correct"

data <- bind_rows(store.data, dots.data)

# some subjects only have data in one task, so remove from joint model of both tasks
tabs <- table(data$subj, data$task)
tmp <- apply(tabs, 1, function(x) { all(x>0) })
keep.subj <- tabs[tmp,] %>% rownames
data <- data %>%
    filter(subj %in% keep.subj) %>%
    mutate(subj=as.numeric(subj)) %>%
    arrange(subj) %>%
    mutate(subj=factor(subj), task=factor(task))

# split data by SAT condition
all.data <- list()
all.data$FS <- subset(data, cond == "FS")
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
