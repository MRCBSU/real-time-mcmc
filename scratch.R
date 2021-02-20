for(per in 1:num.bp){

    TAcolname <- paste0("time", per)
    TA <- bind_rows(TA %>% filter(era < (per - 1)) %>% mutate(!!TAcolname := 0),
                    TA %>% filter(era == (per - 1)) %>% mutate(!!TAcolname := time),
                    TA %>% filter(era > (per - 1)) %>% mutate(!!TAcolname := max(TA$time))
                    )
    if(per > 1){ ## Remove some redundant rows
        TAprevcolname <- paste0("time", per - 1)
        TA <- TA %>% filter(!(!!sym(TAprevcolname) == max(TA$time) & !!sym(TAcolname) == 0 & time == 0))
        }
    ## TA <- TA %>%
    ##     mutate(!!TAcolname := ifelse(era <= (per - 1) ? 0 :
    ##                                  ifelse(era == (per - 1) ? time : max(TA$time)))
    ##            )
}

