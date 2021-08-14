require(tidyverse)
require(cubelyr)

suppressMessages(extract <- R.utils::extract)

## Get the repo location
proj.dir <- file.path(dirname(dirname(getwd())))

## Load in convolution functions for estimating hospitalisations, ICUs and deaths
source(file.path(proj.dir, "R", "output", "convolution.R"))
## Functions on paramters of gamma distributions
source(file.path(proj.dir, "R", "output", "gamma_fns.R"))

## Other useful functions
drop.from.names <- function(x, value) {
  names(x)[!(names(x) %in% value)]
}

apply.convolution <- function(start, func, over = "date") {
    output.dimnames <- dimnames(start)
    result.length <- length(dimnames(start)[[over]])
    name.over <- dimnames(start)[[over]]
    dims.to.preserve <- drop.from.names(output.dimnames, over)
    result <- start %>%
        apply(dims.to.preserve, conv, b = func)
    names(dimnames(result))[1] <- over
    result <- result %>%
        extract(indices = list(1:result.length), dims = over)
    dimnames(result)[[over]] <- name.over
    return(result %>% aperm(perm = names(dimnames(start))))
}

## Where are our outputs to be found
output.dir <- file.path(proj.dir, "model_runs", "20210806", "Prev459_cm6ons_NHS60cutoff_IFR4bp_18wk2_prev14-0PHE_matrices_20210806_timeuse_household_deaths")
load(file.path(output.dir, "projections_midterm.RData"))

## Estimated and projected incidence stored in `incidence' object

## Initially, not incorporating age
infections <- apply(infections, c("iteration", "date", "region"), sum)

## Delay distribution
warwick.delay <- c(rep(0, 3), 0.101541687310437, 0.100800333047344, 0.0959763158457943, 0.0891720283028308, 0.0814111772059684, 0.0733313959322210, 0.0653516009824231, 0.0577138301669170, 0.0505836728272167, 0.0440378148640618, 0.0381018779042572, 0.0327901892018164, 0.0280853405593643, 0.0239399668017949, 0.0203183211083898, 0.0171850072099235, 0.0144766370432761, 0.0121482487049283, 0.0101690151236128, 0.00848264674110468, 0.00705198507468494, 0.00584947875638519, 0.00483839870685633, 0.00398982407913745, 0.00328286523283823, 0.00269531866872596, 0.00220667028721191, 0.00180248093733307, 0.00147007619932521, 0.00119579517381864)

## CoCIN info - lengths of stay
mixes <- c(0.164, 1 - 0.103 - 0.164, 0.03, 0.103 * (1 - 0.291))
idxs <- table(sample(4, 1000000, replace = TRUE, prob = mixes))
mean.gams <- c(rgamma(idxs[1], shape = 9.8*42, rate = 42),
               rgamma(idxs[2], shape = 8.8*168, rate = 168),
               rgamma(idxs[3], shape = 4.88*15.8, rate = 4.88),
               rgamma(idxs[4], shape = 13.94*17.7, rate = 13.94))
sd.gams <- c(rgamma(idxs[1], shape = 7.4 * 31.39, rate = 31.39),
             rgamma(idxs[2], shape = 8 * 101.06, rate = 101.06),
             rgamma(idxs[3], shape = 10.7*3.24, rate = 3.24),
             rgamma(idxs[4], shape = 9.8 * 9.83, rate = 9.83))
samp <- sapply(1:length(mean.gams), function(x) {
    ab <- gamma.shape.rate.params(c(mean.gams[x], sd.gams[x]^2))
    rgamma(1, shape = ab[1], rate = ab[2])
})
los.dist <- cut(samp, 0:ceiling(max(samp))) %>%
    table() / 1000000


## Beds available - fields icu_prev_acute1 in the NHS SitRep
beds.used <- c(339, 998, 831, 1003, 801, 421, 327)
## beds.used <- c(197, 560, 646, 806, 741, 275, 207)
beds.available <- as_tibble(list(capacity = c(9915, 12760, 15092, 14934, 13331, 10993, 8948), region = dimnames(infections)$region))
beds.total <- beds.used / (beds.used + beds.available$capacity)
beds.day <- lubridate::as_date("20210805")

## Function to translate these infections into hospital occupancy
hosp.occupancy <- function(infections = infections, delay1 = warwick.delay, delay2 = los.dist){
    ## Take a fraction, ihr of your infections
    sev.infects <- infections
    ## Delay over the convolution distribution
    admission <- apply.convolution(sev.infects, delay1)
    ## Apply length of stay - take cumsums to get cumulative numbers departed from hospital
    departure <- apply.convolution(admission, delay2) %>%
        apply(c("iteration", "region"), cumsum)
    ## Similarly we want cumulative admissions
    admission <- apply(admission, c("iteration", "region"), cumsum)
    ## Occupancy is all those who have entered - all those who have left
    admission - departure
}

occupancy <- hosp.occupancy(infections)

## Scale the occupancies to match what we observe on beds.day
dat.idx <- which(lubridate::as_date(as.integer(dimnames(occupancy)[[1]])) == beds.day)

## Get the scalings we need to match current data
scalings <- extract(occupancy, indices = dat.idx, dims = "date", drop = TRUE) %>%
    sweep(2, beds.used, "/")

## Apply these to the occupancies
occupancy <- sweep(occupancy, 2:3, scalings, "/")

## stop()

## Get quantiles and plot
occupancy <- apply(occupancy, c("date", "region"), quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
    as.tbl_cube(met_name = "occupancy") %>%
    as_tibble() %>%
    mutate(Date = lubridate::as_date(date)) %>%
    select(-date) %>%
    pivot_wider(id_cols = c("region", "Date"), names_from = Var1, names_prefix = "q", values_from = occupancy) %>%
    left_join(beds.available)

require(ggplot2)

## First ten days
gg.occupancy <- occupancy %>% filter(Date >= beds.day - 10, Date <= beds.day + 10) %>%
    ggplot(aes(x = Date, y = `q50%`, ymin = `q2.5%`, ymax = `q97.5%`, fill = region)) +
    geom_line() +
    geom_ribbon(alpha = 0.3) +
    geom_ribbon(aes(x=Date,y=`q50%`, ymin = `q25%`, ymax = `q75%`, fill = region), alpha = 0.6) +
    facet_wrap(~region, scales = "free_y") +
    ggtitle("The Next Ten Days") +
    theme(legend.position = "none")
    ## xlim(range()

ggsave("tendays.png", gg.occupancy, width = 14, height = 7)

## Intermediate
gg.intermediate <- occupancy %>% filter(Date >= beds.day - 10, Date <= beds.day + 28) %>%
    ggplot(aes(x = Date, y = `q50%`, ymin = `q2.5%`, ymax = `q97.5%`, fill = region)) +
    geom_line() +
    geom_ribbon(alpha = 0.3) +
    geom_ribbon(aes(x=Date,y=`q50%`, ymin = `q25%`, ymax = `q75%`, fill = region), alpha = 0.6) +
    ## geom_line(aes(x = Date, y = capacity), size = 2, linetype = 2) +
    ## geom_line(aes(x = Date, y = 0.8 * capacity), size = 1.6, linetype = 3) +
    facet_wrap(~region, scales = "free_y") +
    ggtitle("Th Next 28 Days") +
    theme(legend.position = "none")

ggsave("intermediate.png", gg.intermediate, width = 14, height = 7)

## Another fifty days
gg.capacity <- occupancy %>% filter(Date >= beds.day - 10, Date <= beds.day + 78) %>%
    ggplot(aes(x = Date, y = `q50%`, ymin = `q2.5%`, ymax = `q97.5%`, fill = region)) +
    geom_line() +
    geom_ribbon(alpha = 0.3) +
    geom_ribbon(aes(x=Date,y=`q50%`, ymin = `q25%`, ymax = `q75%`, fill = region), alpha = 0.6) +
    ## geom_line(aes(x = Date, y = capacity), size = 2, linetype = 2) +
    ## geom_line(aes(x = Date, y = 0.8 * capacity), size = 1.6, linetype = 3) +
    facet_wrap(~region, scales = "free_y") +
    ggtitle("The Next 78 Days") +
    theme(legend.position = "none")

ggsave("capacity.png", gg.capacity, width = 14, height = 7)

## What's in the bank already

## Set infections to zero beyond 11/10
banked.infections <- infections
banked.infections[, (dat.idx + 1):dim(infections)[2], ] <- 0
banked.occupancy <- banked.infections %>%
    hosp.occupancy() %>%
    sweep(2:3, scalings, "/") %>%
    apply(c("date", "region"), quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
    as.tbl_cube(met_name = "occupancy") %>%
    as_tibble() %>%
    mutate(Date = lubridate::as_date(date)) %>%
    select(-date) %>%
    pivot_wider(id_cols = c("region", "Date"), names_from = Var1, names_prefix = "q", values_from = occupancy) %>%
    left_join(beds.available)

## First ten days
gg.occupancy <- banked.occupancy %>% filter(Date >= beds.day - 10, Date <= beds.day + 10) %>%
    ggplot(aes(x = Date, y = `q50%`, ymin = `q2.5%`, ymax = `q97.5%`, fill = region)) +
    geom_line() +
    geom_ribbon(alpha = 0.3) +
    geom_ribbon(aes(x=Date,y=`q50%`, ymin = `q25%`, ymax = `q75%`, fill = region), alpha = 0.6) +
    facet_wrap(~region, scales = "free_y") +
    ggtitle("The Next Ten Days - of infections already occurred") +
    theme(legend.position = "none")
    ## xlim(range()

ggsave("banked_tendays.png", gg.occupancy, width = 14, height = 7)

## Total lock-down in two weeks?

## Set infections to zero beyond two weeks time
lockdown.infections <- infections
lockdown.infections[, (dat.idx + 14):dim(infections)[2], ] <- 0
lockdown.occupancy <- lockdown.infections %>%
    hosp.occupancy() %>%
    sweep(2:3, scalings, "/") %>%
    apply(c("date", "region"), quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
    as.tbl_cube(met_name = "occupancy") %>%
    as_tibble() %>%
    mutate(Date = lubridate::as_date(date)) %>%
    select(-date) %>%
    pivot_wider(id_cols = c("region", "Date"), names_from = Var1, names_prefix = "q", values_from = occupancy) %>%
    left_join(beds.available)


## Another two weeks
gg.capacity <- lockdown.occupancy %>% filter(Date >= beds.day - 10, Date <= beds.day + 78) %>%
    ggplot(aes(x = Date, y = `q50%`, ymin = `q2.5%`, ymax = `q97.5%`, fill = region)) +
    geom_line() +
    geom_ribbon(alpha = 0.3) +
    geom_ribbon(aes(x=Date,y=`q50%`, ymin = `q25%`, ymax = `q75%`, fill = region), alpha = 0.6) +
    ## geom_line(aes(x = Date, y = capacity), size = 2, linetype = 2) +
    ## geom_line(aes(x = Date, y = 0.8 * capacity), size = 1.6, linetype = 3) +
    facet_wrap(~region, scales = "free_y") +
    ggtitle(glue::glue("Longer-term - if no further infections occur after {format(beds.day + 14, '%d/%m')}")) +
    theme(legend.position = "none")

ggsave("lockdown_capacity.png", gg.capacity, width = 14, height = 7)

## Set infections to zero beyond four weeks time
lockdown2.infections <- infections
lockdown2.infections[, (dat.idx + 28):dim(infections)[2], ] <- 0
lockdown2.occupancy <- lockdown2.infections %>%
    hosp.occupancy() %>%
    sweep(2:3, scalings, "/") %>%
    apply(c("date", "region"), quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
    as.tbl_cube(met_name = "occupancy") %>%
    as_tibble() %>%
    mutate(Date = lubridate::as_date(date)) %>%
    select(-date) %>%
    pivot_wider(id_cols = c("region", "Date"), names_from = Var1, names_prefix = "q", values_from = occupancy) %>%
    left_join(beds.available)


## Another two weeks
gg.capacity <- lockdown2.occupancy %>% filter(Date >= beds.day - 10, Date <= beds.day + 78) %>%
    ggplot(aes(x = Date, y = `q50%`, ymin = `q2.5%`, ymax = `q97.5%`, fill = region)) +
    geom_line() +
    geom_ribbon(alpha = 0.3) +
    geom_ribbon(aes(x=Date,y=`q50%`, ymin = `q25%`, ymax = `q75%`, fill = region), alpha = 0.6) +
    ## geom_line(aes(x = Date, y = capacity), size = 2, linetype = 2) +
    ## geom_line(aes(x = Date, y = 0.8 * capacity), size = 1.6, linetype = 3) +
    facet_wrap(~region, scales = "free_y") +
    ggtitle(glue::glue("Longer-term - if no further infections occur after {format(beds.day + 28, '%d/%m')}")) +
    theme(legend.position = "none")

ggsave("lockdown2_capacity.png", gg.capacity, width = 14, height = 7)


## ## Merge Warwick and PHE data
## warwick.data <- read_csv("Warwick_Output_for_Paul.csv") %>%
##     filter(ValueType == "hospital_prev") %>%
##     mutate(Date = as.Date(paste(Year, Month, Day, sep = "/"), format = "%Y/%m/%d"),
##            region = str_replace_all(Geography, " ", "_"),
##            group = "Warwick") %>%
##     inner_join(beds.available) %>%
##     rename(`q5%` = `Quantile 0.05`,
##            `q95%` = `Quantile 0.95`,
##            `q50%` = Value) %>%
##     select(region, Date, `q5%`, `q50%`, `q95%`, group, capacity) %>%
##     bind_rows(occupancy %>% select(region, Date, `q5%`, `q50%`, `q95%`, capacity) %>% mutate(group = "PHECam")) %>%
##     filter(Date >= beds.day, Date <= beds.day + 78)

## ## Comparison with Warwick's numbers
## gg.compare <- warwick.data %>%
##     ggplot(aes(x = Date, y = `q50%`, ymin = `q5%`, ymax = `q95%`, colour = group, fill = group)) +
##     geom_line() +
##     geom_ribbon(alpha = 0.2) +
##     geom_line(aes(x = Date, y = capacity), size = 2, linetype = 2, colour = "black") +
##     geom_line(aes(x = Date, y = 0.8 * capacity), size = 1.6, linetype = 3, colour = "black") +
##     facet_wrap(~region, scales = "free_y") +
##     ggtitle("Comparison of PHE and Warwick")

## ggsave("comparison.png", gg.compare, width = 14, height = 7)
