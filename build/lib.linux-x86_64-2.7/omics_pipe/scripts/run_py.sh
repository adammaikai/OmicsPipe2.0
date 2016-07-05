jobname=`echo omics_pipe | cut -c1-15`

echo "#PBS -N ${jobname}"
echo '#PBS -l nodes=1:ppn=8'                       
echo '#PBS -l mem=31gb'                            #Memory required for highest mem job
echo '#PBS -l walltime=480:00:00'                  #Walltime 
echo "#PBS -o /gpfs/group/su/kfisch/Bodymap/logs"  #Directory for stdout
echo "#PBS -e /gpfs/group/su/kfisch/Bodymap/logs"  #Directory for stderr
echo '#PBS -m ea'
echo '#PBS -M kfisch@scripps.edu'                  #Email
echo 
echo 'cd /gpfs/group/su/kfisch/OA/src/bodymap/src' #Location of python scripts
echo 'module load python'                   #Load python
echo 'python automate.py'                          #Excecute main script
echo 'exit 0'
