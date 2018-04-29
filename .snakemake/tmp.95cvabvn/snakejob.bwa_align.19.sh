#!/bin/bash
# properties = {"rule": "bwa_align", "local": false, "input": ["/home/tamburin/fiona/bacteremia/1.assemble_trimmed/filtered_assemblies/B09.fasta", "/home/tamburin/fiona/bacteremia/1.assemble_trimmed/filtered_assemblies/B09.fasta.amb", "/home/tamburin/fiona/bacteremia/1.assemble_trimmed/filtered_assemblies/B09.fasta.ann", "/home/tamburin/fiona/bacteremia/1.assemble_trimmed/filtered_assemblies/B09.fasta.bwt", "/home/tamburin/fiona/bacteremia/1.assemble_trimmed/filtered_assemblies/B09.fasta.pac", "/home/tamburin/fiona/bacteremia/1.assemble_trimmed/filtered_assemblies/B09.fasta.sa", "input_samples/B09.fq"], "output": ["filtered_bam/B09.filtered.bam"], "wildcards": ["B09"], "params": {}, "log": [], "threads": 8, "resources": {"mem": 32, "time": 6}, "jobid": 19, "cluster": {}}
cd /srv/gsfs0/projects/bhatt/fiona/bacteremia/strainsifter && \
/home/tamburin/tools/miniconda3/bin/python \
-m snakemake filtered_bam/B09.filtered.bam --snakefile /srv/gsfs0/projects/bhatt/fiona/bacteremia/strainsifter/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /srv/gsfs0/projects/bhatt/fiona/bacteremia/strainsifter/.snakemake/tmp.95cvabvn /home/tamburin/fiona/bacteremia/1.assemble_trimmed/filtered_assemblies/B09.fasta /home/tamburin/fiona/bacteremia/1.assemble_trimmed/filtered_assemblies/B09.fasta.amb /home/tamburin/fiona/bacteremia/1.assemble_trimmed/filtered_assemblies/B09.fasta.ann /home/tamburin/fiona/bacteremia/1.assemble_trimmed/filtered_assemblies/B09.fasta.bwt /home/tamburin/fiona/bacteremia/1.assemble_trimmed/filtered_assemblies/B09.fasta.pac /home/tamburin/fiona/bacteremia/1.assemble_trimmed/filtered_assemblies/B09.fasta.sa input_samples/B09.fq --latency-wait 5 \
--benchmark-repeats 1 --attempt 1 \
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --no-hooks --nolock --timestamp  --force-use-threads  --allowed-rules bwa_align  && exit 0 || exit 1

