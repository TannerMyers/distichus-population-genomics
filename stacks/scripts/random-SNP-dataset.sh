[[ -d $out_folder/run_%j ]] || mkdir -p $out_folder/run_%j # the %j is slurm for the job number -- figure out how to do this with PBS

populations \
	--write-random-snp \
	--min-mac 2

cp $out_folder/populations.structure $structure_folder/populations.structure.%j


