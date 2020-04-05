set -e
set -x

declare -a region
region=($(seq 1 7))

for r in ${region[@]}; do
	Rscript R/data/format_deaths.R /data/covid-19/data-raw/manual-downloads/2020-04-04/deaths.csv $r
	Rscript set_up.R $r
done

for r in ${region[@]}; do
	Rscript run.R $r &
done
wait

for r in ${region[@]}; do
	Rscript R/output/build_report.R $r &
done
wait

Rscript R/output/combine_region_runs.R
