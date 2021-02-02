# Make array for batches
list=(6 7 8)
echo ${list[@]}

# Make filepath directory
dir=../Rawdata/P2_4/
echo $dir

# Make path to barcodes file
barcodes=../info/barcodes_lane2.1.tsv
echo $barcodes

for i in ${list[@]}; do
        filepath1=$dir"P2_4_CKDL190143248-1a-8_H723VBBXX_L[$i]_1.fq.gz"
        echo $filepath1
        filepath2=$dir"P2_4_CKDL190143248-1a-8_H723VBBXX_L[$i]_2.fq.gz"
        echo $filepath2
        jobname=proc_radtags_[$i]_2_4
        echo $jobname

        qsub -N $jobname \
                -q general \
                -W group_list=jro0014_lab \
                -W x=FLAGS:ADVRES:jro0014_lab \
                -d /scratch/phyletica/distichus/samples/ \
                 -l walltime=48:00:00 << HereDoc


                        module load gcc/5.3.0
                        module load zlib
                        module load stacks/2.41

                        process_radtags --paired -1 $filepath1 -2 $filepath2 -b $barcodes -o /scratch/phyletica/distichus/samples/ \
                                --renz_1 sbfI --renz_2 mspI --inline_null -c -q -r &> process_radtags.lane2.1.oe 

HereDoc

done
