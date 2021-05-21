# Make array for batches
list=(6 7 8)
echo ${list[@]}

# Make filepath directory
dir=../../Rawdata/P1_1/
echo $dir

# Make path to barcodes file
barcodes=../../info/barcodes_lane1_1.tsv
echo $barcodes

for i in ${list[@]}; do
        
	# sed -i 's/1a-1/1a-*/g' demultiplex*_*.sh 
	# sed -i 's/1_1/*_*/g' demultiplex*_*.sh Substitute specific values in script

	filepath1=$dir"P1_1_CKDL190143248-1a-1_H723VBBXX_L[$i]_1.fq.gz"
        echo $filepath1
        filepath2=$dir"P1_1_CKDL190143248-1a-1_H723VBBXX_L[$i]_2.fq.gz"
        echo $filepath2
        jobname=proc_radtags_[$i]_1_1
        echo $jobname

        qsub -N $jobname \
                -q general \
                -W group_list=jro0014_lab \
                -W x=FLAGS:ADVRES:jro0014_lab \
                -d /scratch/phyletica/distichus/cleaned/cleaned1_1 \
                 -l walltime=48:00:00 << HereDoc


                        module load stacks/2.5

                        process_radtags -P -1 $filepath1 -2 $filepath2 -b $barcodes -o ./ \
                                --renz-1 sbfI --renz-2 mspI --inline_null -c -q -r &> process_radtags.lane1_1.oe 

HereDoc

done
