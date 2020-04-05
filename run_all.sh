set -e
set -x

declare -a region
region=($(seq 1 7))

for r in ${region[@]}; do
	Rscript R/data/format_deaths.R $1 $r
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
