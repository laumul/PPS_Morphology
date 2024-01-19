# Laura Mulvey | 10.06.2022

# take the marginal likelihood values produed by the stepping stone and path sampling

model=$1
echo "Sim no.;mk;mk+GV;mk+GVmultistate;mk+G;mk+V;mk+Vmultistate;mk+Gmultistate" > $model"_path_results.csv"

files=$(ls output_$model)

ls output_$model/mk > $model.txt
# loop=$(cat $model.txt | wc -l)

cd output_$model

for i in {1..20}
do
for j in $files
do
sim=$(cat ../$model.txt | sed -n $i'p')
cd $j/$sim

declare "step_$j"=$(sed -n '2p' ss_may_marginal_likelihood.csv)

cd ../..
done
echo "$i;$step_mk;$step_mk_GV;$step_mk_GVmultistate;$step_mk_G;$step_mk_V;$step_mk_Vmultistate;$step_mk_Gmultistate" >> ../$model"_path_results.csv"
done

cd ..
