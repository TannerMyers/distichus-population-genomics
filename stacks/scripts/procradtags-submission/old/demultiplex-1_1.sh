# Make array for batches
list=(6 7 8)
echo ${list[@]}

# Make filepath directory
dir=../Rawdata/P1_1/
echo $dir

# Make path to barcodes file
barcodes=../info/barcodes_lane1_1.tsv
echo $barcodes

for i in ${list[@]}; do
        filepath1=$dir"P1_1_CKDL190143248-1a-1_H723VBBXX_L"$i"_1.fq.gz"
        echo $filepath1
        filepath2=$dir"P1_1_CKDL190143248-1a-1_H723VBBXX_L"$i"_2.fq.gz"
        echo $filepath2
        jobname=proc_radtags_[$i]_1_1
        echo $jobname

        # Make a new new directory with $i in the name
        NEW_DIR="$i"_reads
        echo $NEW_DIR

        if [ ! -d /scratch/phyletica/distichus/samples/$NEW_DIR ]
                then
                        mkdir /scratch/phyletica/distichus/samples/$NEW_DIR
                else 
                        echo "$NEW_DIR already exists"
                #       continue

                cd /scratch/phyletica/distichus/samples/$NEW_DIR

                echo $(pwd)

        qsub -N $jobname \
                -q general \
                -W group_list=jro0014_lab \
                -W x=FLAGS:ADVRES:jro0014_lab \
                -d /scratch/phyletica/distichus/samples/ \
                 -l walltime=48:00:00 << HereDoc

                        module load /home/shared/jro0014_lab/privatemodules/stacks/2.5

                        process_radtags --paired -1 $filepath1 -2 $filepath2 -b $barcodes -o /scratch/phyletica/distichus/samples/ \
                                -o /scratch/phyletica/distichus/samples/$NEW_DIR \
                                --renz_1 sbfI --renz_2 mspI --inline_null -c -q -r &> process_radtags.lane1_1.oe 

HereDoc
                 continue
        fi
        
done
