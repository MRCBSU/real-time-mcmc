set -e

declare -a region
region=($(seq 1 7))

for r in ${region[@]}; do
	Rscript R/data/format_deaths.R $1 $2 $r
	Rscript set_up.R $2 $r
done

for r in ${region[@]}; do
	Rscript run.R $2 $r &
done
wait

for r in ${region[@]}; do
	Rscript R/output/build_report.R $2 $r &
done
wait

Rscript R/output/combine_region_runs.R $2 1
