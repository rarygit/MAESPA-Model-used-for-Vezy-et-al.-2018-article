#!/bin/bash

# Shell to use
#$ -S /bin/bash
# Name of the job in SGE
#$ -N my_job
# Name of the queue to use
#$-q cemeb.q
# Maximum hardware time allowed for this job
#$-l h_rt=15:00:00
# Merge output (error and standard with '-j y') into one single file called my_job.out
#$-o my_job.out -j y
# run in the current directory
#$ -cwd

#-----------------Lancement des runs a partir du repertoire dans le scrach/$USER des noeuds :

#Creation dossier temporaire :
#il faut le creer prealablement :
if [ ! -d /export/scrach/$USER ]; then mkdir /export/scrach/$USER/; fi
mkdir -p /export/scrach/$USER/tmp/

#---------Gestion des differentes nomenclatures de fichiers en fonction de $SGE_TASK_ID pour la copie vers les noeuds:

if (( 1<=$SGE_TASK_ID<=9)); then
#Copie des fichiers d entree du /home vers le dossier temporaire du noeud:
cp -r /home/$USER/MAESPA/2-Simulations_FullPeriod2012_Seq/000$SGE_TASK_ID'_'* /export/scrach/$USER/tmp/
fi

if (( 10<=$SGE_TASK_ID<=99)); then
#Copie des fichiers d entree du /home vers le dossier temporaire du noeud:
cp -r /home/$USER/MAESPA/2-Simulations_FullPeriod2012_Seq/00$SGE_TASK_ID'_'* /export/scrach/$USER/tmp/
fi

if (( 100<=$SGE_TASK_ID<=999)); then
#Copie des fichiers d entree du /home vers le dossier temporaire du noeud:
cp -r /home/$USER/MAESPA/2-Simulations_FullPeriod2012_Seq/0$SGE_TASK_ID'_'* /export/scrach/$USER/tmp/
fi


#Copie de l'executable du /home vers le tmp du noeud :
cp -r /home/$USER/MAESPA/2-Simulations_FullPeriod2012_Seq/maespa.out /export/scrach/$USER/tmp/


#------------------------------------------------


#Execution de l'executable dans le dossier temporaire du noeud :
cd /export/scrach/$USER/tmp/ # On se met dans le repertoire tmp du noeud
echo "$SGE_TASK_ID" | ./maespa.out # lancement 

#Compression des resultats dans le fichier $SGE_TASK_ID'_maespa.tar.gz' - Pour le decompresser c'est : tar -xvf $SGE_TASK_ID'_maespa.tar.gz':

tar czf $SGE_TASK_ID'_maespa.tar.gz' $SGE_TASK_ID'_hrflux.dat' $SGE_TASK_ID'_Dayflx.dat' $SGE_TASK_ID'_Maeserr.dat' $SGE_TASK_ID'_watbal.dat' $SGE_TASK_ID'_watbalday.dat' $SGE_TASK_ID'_watlay.dat' $SGE_TASK_ID'_watsoilt.dat' $SGE_TASK_ID'_watupt.dat'

#Copie des résultats compresses dans le home :
cp -rp $SGE_TASK_ID'_maespa.tar.gz' /home/$USER/MAESPA/2-Simulations_FullPeriod2012_Seq/.

#Suppression du dossier temporaire sur les noeuds :
rm -rf /export/scrach/$USER







