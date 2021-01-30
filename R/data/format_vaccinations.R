## For now, simply generating some dummy data based on the early provisional NIMS data shared by Andre.
require(tidyverse)
require(cubelyr)
require(lubridate)

r.even <- function(vaccs, len) rmultinom(1, vaccs, rep(1, len))
regions <- c("East_of_England", "London", "Midlands", "North_East", "North_West", "South_East", "South_West")


## ====== DATA ON FIRST VACCINATION

vaccs <- rbind(r.even(269, 20), r.even(11760, 20), r.even(9327, 10), r.even(5746, 5), r.even(5681, 5), r.even(4153, 5), r.even(1361, 5), r.even(583, 5), r.even(386, 5), r.even(14424, 1))
vaccs <- cbind(vaccs, rbind(r.even(1352, 20), r.even(45379, 20), r.even(36858, 10), r.even(23132, 5), r.even(23686, 5), r.even(17056, 5), r.even(6733, 5), r.even(4228, 5), r.even(7299, 5), r.even(395103, 1)))
vaccs <- cbind(vaccs, rbind(r.even(2202, 20), r.even(75231, 20), r.even(60873, 10), r.even(37872, 5), r.even(38787, 5), r.even(27809, 5), r.even(11101, 5), r.even(7071, 5), r.even(13338, 5), r.even(543397, 1)))
vaccs <- cbind(vaccs, rbind(r.even(3563, 20), r.even(124955, 20), r.even(97479, 10), r.even(59808, 5), r.even(60157, 5), r.even(43763, 5), r.even(18352, 5), r.even(13783, 5), r.even(32634, 5), r.even(670040, 1)))
dimnames(vaccs) <- list(ages = 0:80, date = lubridate::as_date(c("20201213","20201220","20201227","20210103")))

vaccs <- vaccs %>% as.tbl_cube(met_name = "jabs") %>% as_tibble() %>% mutate(date = lubridate::as_date(date), Age.Grp = cut(ages, breaks = c(0, 1, 5, 15, 25, 45, 65, 75, Inf), right = FALSE, ordered_result = T)) %>% group_by(Age.Grp, date) %>% summarise(jabs = sum(jabs))

probs <- c(8942, 4969, 13188, 5172, 7416, 7709, 6218)
probs <- rbind(probs, c(70593, 59141, 112177, 78726, 71914, 92577, 75194), c(86555, 80364, 152828, 145645, 106357, 143259, 101894), c(112616, 115341, 211936, 201823, 144815, 195125, 141767))
dimnames(probs) <- list(date = lubridate::as_date(c("20201213","20201220","20201227","20210103")), region = regions)
probs <- apply(probs, 1, function(x) x / sum(x))
probs <- probs %>% as.tbl_cube(met_name = "propn") %>% as_tibble() %>% mutate(date = lubridate::as_date(date))

merged <- vaccs %>% left_join(probs)
merged_wide <- pivot_wider(merged, id_cols = 1:3, names_from = region, values_from = propn)

tmp.samples <- NULL
for(i in 1:nrow(merged_wide))
    tmp.samples <- cbind(tmp.samples, rmultinom(1, merged_wide[i, ]$jabs, c(merged_wide[i, ]$East_of_England, merged_wide[i, ]$London, merged_wide[i, ]$Midlands, merged_wide[i, ]$North_East, merged_wide[i, ]$North_West, merged_wide[i, ]$South_East, merged_wide[i, ]$South_West)))
rownames(tmp.samples) <- colnames(merged_wide)[-(1:3)]
sampled.jabs <- merged_wide[, 1:2]
tmp.samples <- as.data.frame(t(tmp.samples))
sampled.jabs <- bind_cols(sampled.jabs, tmp.samples)

sampled.jabs <- sampled.jabs %>% pivot_longer(cols = -(1:2), names_to = "Region", values_to = "Jabs") %>% right_join(expand.grid(date = lubridate::ymd("20200217"):lubridate::ymd("20210117"), Region = colnames(merged_wide)[-(1:3)], Age.Grp = unique(merged_wide$Age.Grp)) %>% as.data.frame %>% mutate(`date` = lubridate::as_date(`date`))) %>% replace_na(list(Jabs = 0)) %>% arrange(date)

vacc1.files <- paste0(file.path(getwd(), "data", "RTM_format", "NHS", "vaccination", "dummy_dose1_"), regions, ".txt")
names(vacc1.files) <- regions

for(reg in regions){
    region.dat <- pivot_wider(sampled.jabs %>% filter(Region == reg),
                              id_cols = 2,
                              names_from = Age.Grp,
                              values_from = Jabs)
    tmpFile <- vacc1.files[reg]
    dir.create(dirname(tmpFile), recursive = TRUE, showWarnings = FALSE)

    region.dat %>%
        write_tsv(tmpFile,
                  col_names = FALSE)
    }

## ====== DATA ON SECOND VACCINATION

vaccs <- rbind(r.even(0, 20), r.even(0, 20), r.even(0, 10), r.even(0, 5), r.even(0, 5), r.even(0, 5), r.even(0, 5), r.even(0, 5), r.even(0, 5), r.even(0, 1))
vaccs <- cbind(vaccs, rbind(r.even(0, 20), r.even(0, 20), r.even(0, 10), r.even(0, 5), r.even(0, 5), r.even(0, 5), r.even(0, 5), r.even(0, 5), r.even(0, 5), r.even(0, 1)))
vaccs <- cbind(vaccs, rbind(r.even(0, 20), r.even(0, 20), r.even(0, 10), r.even(0, 5), r.even(0, 5), r.even(0, 5), r.even(0, 5), r.even(0, 5), r.even(0, 5), r.even(0, 1)))
vaccs <- cbind(vaccs, rbind(r.even(81, 20), r.even(3700, 20), r.even(3125, 10), r.even(1951, 5), r.even(1876, 5), r.even(1311, 5), r.even(476, 5), r.even(240, 5), r.even(153, 5), r.even(6156, 1)))
dimnames(vaccs) <- list(ages = 0:80, date = lubridate::as_date(c("20201213","20201220","20201227","20210103")))

vaccs <- vaccs %>% as.tbl_cube(met_name = "jabs") %>% as_tibble() %>% mutate(date = lubridate::as_date(date), Age.Grp = cut(ages, breaks = c(0, 1, 5, 15, 25, 45, 65, 75, Inf), right = FALSE, ordered_result = T)) %>% group_by(Age.Grp, date) %>% summarise(jabs = sum(jabs))

probs <- rep(1, 7)
probs <- rbind(probs, rep(1, 7), rep(1,7), c(3710, 1102, 4165, 1102, 4501, 2408, 2041))
dimnames(probs) <- list(date = lubridate::as_date(c("20201213","20201220","20201227","20210103")), region = regions)
probs <- apply(probs, 1, function(x) x / sum(x))
probs <- probs %>% as.tbl_cube(met_name = "propn") %>% as_tibble() %>% mutate(date = lubridate::as_date(date))

merged <- vaccs %>% left_join(probs)
merged_wide <- pivot_wider(merged, id_cols = 1:3, names_from = region, values_from = propn)

tmp.samples <- NULL
for(i in 1:nrow(merged_wide))
    tmp.samples <- cbind(tmp.samples, rmultinom(1, merged_wide[i, ]$jabs, c(merged_wide[i, ]$East_of_England, merged_wide[i, ]$London, merged_wide[i, ]$Midlands, merged_wide[i, ]$North_East, merged_wide[i, ]$North_West, merged_wide[i, ]$South_East, merged_wide[i, ]$South_West)))
rownames(tmp.samples) <- colnames(merged_wide)[-(1:3)]
sampled.jabs <- merged_wide[, 1:2]
tmp.samples <- as.data.frame(t(tmp.samples))
sampled.jabs <- bind_cols(sampled.jabs, tmp.samples)

sampled.jabs <- sampled.jabs %>% pivot_longer(cols = -(1:2), names_to = "Region", values_to = "Jabs") %>% right_join(expand.grid(date = lubridate::ymd("20200217"):lubridate::ymd("20210117"), Region = colnames(merged_wide)[-(1:3)], Age.Grp = unique(merged_wide$Age.Grp)) %>% as.data.frame %>% mutate(`date` = lubridate::as_date(`date`))) %>% replace_na(list(Jabs = 0)) %>% arrange(date)

vacc2.files <- paste0(file.path(getwd(), "data", "RTM_format", "NHS", "vaccination", "dummy_dose2_"), regions, ".txt")
names(vacc2.files) <- regions

for(reg in regions){
    region.dat <- pivot_wider(sampled.jabs %>% filter(Region == reg),
                              id_cols = 2,
                              names_from = Age.Grp,
                              values_from = Jabs)
    tmpFile <- vacc2.files[reg]
    dir.create(dirname(tmpFile), recursive = TRUE, showWarnings = FALSE)

    region.dat %>%
        write_tsv(tmpFile,
                  col_names = FALSE)
    }
