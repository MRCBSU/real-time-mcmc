for(dt in (d0 + 1):d.end){

    ## Get denominator population sizes
    ijab <- jab.dat %>% left_join(jab.dat %>%
                              filter(dose == "First") %>%
                              group_by(region, age.grp) %>%
                              summarise(sdate = sdate, cumsum.n1 = cumsum(n))) %>%
        left_join(jab.dat %>%
                  filter(dose == "Second") %>%
                  group_by(region, age.grp) %>%
                  summarise(sdate = sdate, cumsum.n2 = cumsum(n)) %>%
                  replace_na(list(cumsum.n2 = 0))) %>%
        mutate(denom.pop1 = pop - cumsum.n1,
               denom.pop2 = cumsum.n1 - cumsum.n2)

    
    ## In reg, how many people are due a second jab.
    capacity <- vacc.guide$pricks[vacc.guide$wc == dt]
    vacs.due <- ijab %>%
        left_join(ijab %>%
                  arrange(sdate) %>% 
                  group_by(region, age.grp, dose) %>%
                  summarise(sdate = sdate, vac2.due = pos.part(dplyr::lag(cumsum.n1, twelve.weeks) - cumsum.n2))
                  ) %>%
        filter(sdate == dt - 1, dose == "Second") %>%
        mutate(n = vac2.due * ifelse(sum(vac2.due) > capacity, capacity / sum(vac2.due), 1),
               sdate = lubridate::as_date(dt)) %>%
        select(sdate, region, age.grp, dose, n, pPfizer, pop)
    
    capacity <- capacity - sum(vacs.due$n)

    ## Which priority groups are now fully immunised subject to uptake.
    tmp2 <-     
        ijab %>% filter(sdate == dt - 1, dose == "First") %>%
        inner_join(xdist, by = c("region", "age.grp")) %>%
        arrange(region, age.grp) %>%
        mutate(exhausted = ((pop * uptake) <= cumsum.n1)) %>%
        mutate(f = f + ifelse(lead(exhausted, default = FALSE), lead(f, default = 0), 0)) %>%
        mutate(f = ifelse(exhausted, 0, f)) %>%
        mutate(vac1due = f * capacity) %>%
        mutate(exceed.f = f * pos.part(1 - ((uptake * pop) - cumsum.n1) / vac1due))

    while(any(tmp2$exceed.f > 0 & !tmp2$exhausted)){
        
        idx <- which(tmp2$exceed.f > 0 & !tmp2$exhausted)
        tmp2$f[idx] <- tmp2$f[idx] - tmp2$exceed.f[idx]
        tmp2$f[idx-1] <- tmp2$f[idx-1] + tmp2$exceed.f[idx]
        tmp2$exhausted[idx] <- TRUE
        tmp2$exhausted[idx-1] <- FALSE
        
        tmp2 <- tmp2 %>%
            mutate(vac1due = f * capacity) %>%
            mutate(exceed.f = f * pos.part(1 - ((uptake * pop) - cumsum.n1) / vac1due))
    }
    
    jab.dat <- jab.dat %>%
        bind_rows(
            vacs.due,
            tmp2 %>%
            mutate(n = vac1due,
                   sdate = lubridate::as_date(dt)) %>%
            select(sdate, region, age.grp, dose, n, pPfizer, pop)
        )
            
            
}
    
