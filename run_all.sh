#!/usr/bin/env bash
# Usage: ./run_all.sh <deaths_input_file> <Scottish_deaths_input_file> <date_of_run>
set -e

declare -a region
region=($(seq 1 7))

for r in ${region[@]}; do
	Rscript R/data/format_deaths.R $1 $3 $r
	Rscript set_up.R $3 $r
done
Rscript R/data/format_Scottish_deaths.R $2 $3 8
Rscript set_up.R $3 8

for r in ${region[@]}; do
	Rscript run.R $3 $r &
done
Rscript run.R $3 8 &
wait

for r in ${region[@]}; do
	Rscript R/output/build_report.R $3 $r &
done
Rscript R/output/build_report.R $3 8 &
wait

Rscript R/output/combine_region_runs.R $3 1
