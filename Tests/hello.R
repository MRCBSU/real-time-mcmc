library(Rmpi)

id <- mpi.comm.rank(comm = 0)
np <- mpi.comm.size(comm = 0)
hostname <- mpi.get.processor.name()

cpus_per_task <- Sys.getenv("SLURM_CPUS_PER_TASK")

if(id == 0) {

    first.where.true <- function(vec, predicate) {
    	true.at <- which(predicate(vec))
    	if (length(true.at) == 0) return(NULL)
    	index.to.use <- min(true.at)
    	return(vec[index.to.use])
    }	

    file.length <- nrow(data.table::fread("../data/raw/vaccination/20211012 immunisations SPIM.csv", select = 1, header = T))
    input.colnames <- colnames(data.table::fread("20211003 immunisations SPIM/20211003 immunisations SPIM.csv", nrows = 0))
    nrow_per_process <- floor(file.length / np) 

    possible.col.names <- list(
        age = "age",
        region = c("region_of_residence", "Region_of_Residence"),
        sdate = "vaccination_date",
        type = "product_display_type",
        dose = c("string_dose_number", "dose_number"),
        ltla_code = "ltla_code"
    )

    is.valid.col.name <- function(name) {name %in% input.colnames}
    first.valid.col.name <- function(names) {first.where.true(names, is.valid.col.name)}
    col.names <- lapply(possible.col.names, first.valid.col.name)
    invalid.col.names <- sapply(col.names, is.null)
    if(any(invalid.col.names)) {
        names.invalid.cols <- paste0(names(possible.col.names)[invalid.col.names], collapse = ", ")
        stop(paste("No valid column name(s) for:", names.invalid.cols))
    }

    vacc.date.fmt <-"%d%b%Y"

    vacc.col.args <- list()
    vacc.col.args[[col.names[["type"]]]] <- col_character()
    vacc.col.args[[col.names[["age"]]]] <- col_double()
    vacc.col.args[[col.names[["region"]]]] <- col_character()
    vacc.col.args[[col.names[["dose"]]]] <- col_character()
    vacc.col.args[[col.names[["sdate"]]]] <- col_date(vacc.date.fmt)
    vacc.col.args[[col.names[["ltla_code"]]]] <- col_character()
    vacc.cols <- do.call(cols_only, vacc.col.args)

    get.col.type <- function(x){
        ifelse(identical(x, col_character()), "character",
        ifelse(identical(x, col_double()), "numeric",
        ifelse(identical(x, col_integer()), "integer",
        ifelse(identical(x, col_date(vacc.date.fmt)), "character", "NULL"))))
    }

    fields <- sapply(vacc.col.args[input.colnames], get.col.type) %>%
        unlist(use.names = F) %>%
        as.vector()

    read_in.dat <- list(
        fields = fields,
        file.length = file.length,
        input.colnames = input.colnames,
        nrow_per_process = nrow_per_process
    )

    mpi.bcast.Robj(read_in.dat, rank = 0, comm = 0)

    vacc.df.part <- data.table::fread("20211003 immunisations SPIM/20211003 immunisations SPIM.csv", colClasses = fields, header = F, skip = 1, nrows = file.length - nrow_per_process * (np - 1), nThread = cpus_per_task)
} else {
    read_in.dat <- NULL

    mpi.bcast.Robj(read_in.dat, rank = 0, comm = 0)

    vacc.df.part <- data.table::fread("20211003 immunisations SPIM/20211003 immunisations SPIM.csv", colClasses = fields, header = F, skip = 1, nrows = file.length - nrow_per_process * (np - 1), nThread = cpus_per_task)
}


msg <- sprintf("Hello world from process %03d of %03d, on host %s\n",
    id, np, hostname)
cat(msg)

mpi.barrier(comm = 0)
mpi.finalize()