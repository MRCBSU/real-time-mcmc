stopifnot(all(!is.na(beta.rw.vals)))
stopifnot(length(beta.rw.vals) == nbetas*nr)
vals <- beta.rw.vals[static.zero.beta.locs]
if(!flag.earlier_cm & flag.hack_match_june2020_runs) {
    for(i in static.zero.beta.locs[static.zero.beta.locs > 1]) {
        beta.rw.vals[i-1] <- beta.rw.vals[i]
        beta.rw.vals[i] <- 0
    }
}
stopifnot(all(beta.rw.vals[static.zero.beta.locs] == 0))
stopifnot(all(!is.na(beta.rw.props)))
stopifnot(length(beta.rw.props) == nbetas*nr)
if(flag.hack_match_june2020_runs) {
    beta.rw.props[beta.rw.props == 0] <- 0.02
    beta.rw.props[static.zero.beta.locs] <- 0
}
stopifnot(all(beta.rw.props[static.zero.beta.locs] == 0))
stopifnot(all(beta.rw.props[-static.zero.beta.locs] > 0))
stopifnot(length(contact.dist) == length(contact.proposal))
stopifnot(length(contact.reduction) == length(contact.proposal))
stopifnot(contact.reduction[zero.contact.elements] == 0)
stopifnot(contact.proposal[zero.contact.elements] == 0)
stopifnot(contact.proposal[-zero.contact.elements] != 0)
stopifnot(all(!is.na(contact.reduction)))
stopifnot(all(!is.na(contact.proposal)))
stopifnot(all(!is.na(contact.dist)))
