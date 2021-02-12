load("tmp.RData")


readvacdat <- function(reg, ds){
    r <- which(regions == reg)
    fl <- ifelse(dose == 1, vac1.files[r], vacn.files[r])
    read_tsv(fl, col_names = c("sdate", age.labs)) %>% mutate(dose = ds, region = reg) %>% distinct()
}

## Read in current data
x <- data.frame()
for(reg in regions)
    for(ds in 1:2)
        x <- bind_rows(x, readvacdat(reg, ds))

mdate <- max(x$sdate)

newdates <- mdate + 1:42

xadd <- x %>% bind_rows(expand_grid(sdate = newdates, dose = 1:2, region = regions) %>%
                        mutate(`<1yr` = 0,
                               `1-4` = 0,
                               `5-14` = 0,
                               `15-24` = 0,
                               `25-44` = 0,
                               `45-64` = 0,
                               `65-74` = 0,
                               `75+` = 0)
                        )
for(reg in regions){
    r <- which(regions = reg)
    write_tsv(xadd %>% filter(dose == 1, region == reg) %>% arrange(sdate) %>% select(-c(dose, region)),
              vac1.files[r]
