echo " ----- Getting data from Dalco48 benchmark_rheological_Krieger_Dougherty_v2 case ----- "
cd /home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/examples/benchmark_rheological_Krieger_Dougherty_v2/tmp_trash && \
scp pedro@10.195.70.2:/home/pedro/git/LBDEMcoupling-public/examples/benchmark_rheological_Krieger_Dougherty_v2/outDir/tmp/lattice_average_energy.csv . && \
scp pedro@10.195.70.2:/home/pedro/git/LBDEMcoupling-public/examples/benchmark_rheological_Krieger_Dougherty_v2/outDir/tmp/relative_apparent_viscosity.csv . && \
echo " ----- Data from Dalco48 benchmark_rheological_Krieger_Dougherty_v2 uploaded already -----"

echo " ----- Getting data from Yggdrasil benchmark_rheological_Krieger_Dougherty_v2_D3Q27_without_bulk_viscosity_modification case ----- "
cd /home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/examples/benchmark_rheological_Krieger_Dougherty_v2_D3Q27_without_bulk_viscosity_modification/tmp_trash && \
scp desantan@login1.yggdrasil.hpc.unige.ch:/home/users/d/desantan/scratch/LBDEM/LBDEMcoupling-public/examples/benchmark_rheological_Krieger_Dougherty_v2_D3Q27_without_bulk_viscosity_modification/outDir/tmp/lattice_average_energy.csv . && \
scp desantan@login1.yggdrasil.hpc.unige.ch:/home/users/d/desantan/scratch/LBDEM/LBDEMcoupling-public/examples/benchmark_rheological_Krieger_Dougherty_v2_D3Q27_without_bulk_viscosity_modification/outDir/tmp/relative_apparent_viscosity.csv . && \
echo " ----- Data from Yggdrasil benchmark_rheological_Krieger_Dougherty_v2_D3Q27_without_bulk_viscosity_modification uploaded already -----"

echo " ----- Getting data from Yggdrasil benchmark_rheological_Krieger_Dougherty_v2_D3Q27_with_bulk_viscosity_modification case ----- "
cd /home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/examples/benchmark_rheological_Krieger_Dougherty_v2_D3Q27_with_bulk_viscosity_modification/tmp_trash && \
scp desantan@login1.yggdrasil.hpc.unige.ch:/home/users/d/desantan/scratch/LBDEM/LBDEMcoupling-public/examples/benchmark_rheological_Krieger_Dougherty_v2_D3Q27_with_bulk_viscosity_modification/outDir/tmp/lattice_average_energy.csv . && \
scp desantan@login1.yggdrasil.hpc.unige.ch:/home/users/d/desantan/scratch/LBDEM/LBDEMcoupling-public/examples/benchmark_rheological_Krieger_Dougherty_v2_D3Q27_with_bulk_viscosity_modification/outDir/tmp/relative_apparent_viscosity.csv . && \
echo " ----- Data from Yggdrasil benchmark_rheological_Krieger_Dougherty_v2_D3Q27_with_bulk_viscosity_modification uploaded already -----"


echo " ----- Getting data from Yggdrasil benchmark_rheological_Krieger_Dougherty case ----- "
cd /home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/examples/benchmark_rheological_Krieger_Dougherty/tmp_trash_cluster && \
scp desantan@login1.yggdrasil.hpc.unige.ch:/home/users/d/desantan/scratch/LBDEM/LBDEMcoupling-public/examples/benchmark_rheological_Krieger_Dougherty/outDir/tmp/lattice_average_energy.csv . && \
scp desantan@login1.yggdrasil.hpc.unige.ch:/home/users/d/desantan/scratch/LBDEM/LBDEMcoupling-public/examples/benchmark_rheological_Krieger_Dougherty/outDir/tmp/relative_apparent_viscosity.csv . && \
echo " ----- Data from Yggdrasil benchmark_rheological_Krieger_Dougherty uploaded already -----"

