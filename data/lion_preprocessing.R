library(dplyr)

## data source:
# https://www.movebank.org/cms/webapp?gwt_fragment=page=studies,path=study3809257699

lions <- read.csv("./data/African lions in Central Kalahari Botswana.csv")
nrow(lions)

# drop outliers
lions <- lions[lions$manually.marked.outlier != "true", ]
lions <- lions[lions$algorithm.marked.outlier != "true", ]

unique(lions$tag.local.identifier)

goodtags <- c("GSM08448", "GSM08449", "AL152", "AL153", "AL158", "AL154", "AL156", "AL155")
idx <- which(lions$tag.local.identifier %in% goodtags)

lions2 <- lions[idx, ]
nrow(lions2)
plot(lions2$location.long, lions2$location.lat, asp = 1)

# lions2$time <- substr(lions2$timestamp, 1, 13)
lions2$time <- lions2$timestamp
lions2$time <- as.POSIXct(lions2$time, format="%Y-%m-%d %H:%M:%S", tz = "UTC")

lions3 <- lions2[c("time", "location.long", "location.lat")]
lions3$ID <- lions2$tag.local.identifier
lions3$aniID <- lions2$individual.local.identifier

fill_hourly <- function(d, window = 15*60) {
  # d: columns time, location.long, location.lat

  # compute nearest hour
  round_hour <- as.POSIXct(round(as.numeric(d$time)/3600)*3600, origin="1970-01-01", tz="UTC")
  diff_sec <- abs(as.numeric(d$time - round_hour))

  d$round_hour <- round_hour
  d$diff_sec <- diff_sec

  # keep only points within window
  d <- d[diff_sec <= window, ]

  # if multiple points per hour, keep closest
  d <- d[order(d$round_hour, d$diff_sec), ]
  d <- d[!duplicated(d$round_hour), ]

  # create complete hourly sequence
  full_hours <- seq(from = min(d$round_hour),
                    to   = max(d$round_hour),
                    by   = "hour")

  # merge with full grid
  out <- merge(data.frame(time = full_hours),
               d[, c("round_hour", "location.long", "location.lat")],
               by.x = "time", by.y = "round_hour", all.x = TRUE)

  out
}

data.l <- split(lions3, lions3$ID)

for(id in names(data.l)) {
  d <- data.l[[id]]
  aid <- d$aniID[1]

  data.l[[id]] <- fill_hourly(d)
  data.l[[id]]$ID <- id
  data.l[[id]]$aniID <- aid
}

lions_hourly <- do.call(rbind, data.l)

data <- momentuHMM::prepData(lions_hourly, type = "LL",
                             coordNames = c("location.long", "location.lat"))

# interpolate missing x and y per individual

data <- data %>% arrange(ID, time)

data <- data %>%
  group_by(ID) %>%
  mutate(
    x_int = approx(seq_along(x), x, xout = seq_along(x))$y,
    y_int = approx(seq_along(y), y, xout = seq_along(y))$y
  ) %>%
  ungroup()

data <- as.data.frame(data)

sum(is.na(data$x_int))
sum(is.na(data$y_int))

# add hour
data$hour <- as.POSIXlt(data$time)$hour + 1

# save preprocessed data
saveRDS(data, "./data/lions.rds")
