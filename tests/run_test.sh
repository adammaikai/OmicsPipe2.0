jobname=`echo omics_pipe_test | cut -c1-15`
#test
echo "#PBS -N ${jobname}"
echo '#PBS -l nodes=1:ppn=8'                       
echo '#PBS -l mem=32gb'                            #Memory required for highest mem job
echo '#PBS -l walltime=480:00:00'                  #Walltime 
echo "#PBS -o /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_run"  #Directory for stdout
echo "#PBS -e /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_run"  #Directory for stderr
echo '#PBS -m ea'
echo '#PBS -M kfisch@scripps.edu'                  #Email
echo 
echo 'module load python-addons'                   #Load python
echo 'cd /gpfs/home/kfisch/scripts/omics_pipeline-devel'
echo 'python main_test.py RNAseq /gpfs/home/kfisch/scripts/omics_pipeline-devel/tests/test_params.yaml'                          #Excecute main script
echo 'exit 0'
