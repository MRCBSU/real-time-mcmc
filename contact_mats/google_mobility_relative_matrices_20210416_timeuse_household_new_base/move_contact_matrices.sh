OUTPUT=../../real-time-mcmc/contact_mats/google_mobility_relative_matrices_20210416_timeuse_household
gunzip contact_mat.tar.gz
tar -xvf contact_mat.tar
rm contact_mat.tar

mkdir ${OUTPUT}_old_base
cp * ${OUTPUT}_old_base
cd ${OUTPUT}_old_base
ln -s ../base_matrices/base_matrices_20200529.rds base_matrices.rds

cd -

mkdir ${OUTPUT}_new_base
cp * ${OUTPUT}_new_base
cd ${OUTPUT}_new_base
ln -s ../base_matrices/base_matrices_20201218.rds base_matrices.rds
