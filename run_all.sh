set -e

declare -a region
region=($(seq 1 7))

for r in ${region[@]}; do
	Rscript R/data/format_deaths.R ~/Downloads/2020-04-05_deaths.csv $r
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
