#!/bin/bash
# properties = {"rule": "pileup", "local": false, "input": ["passed_samples/B09.bam", "/home/tamburin/fiona/bacteremia/1.assemble_trimmed/filtered_assemblies/B09.fasta", "/home/tamburin/fiona/bacteremia/1.assemble_trimmed/filtered_assemblies/B09.fasta.fai"], "output": ["pileup/B09.pileup"], "wildcards": ["B09"], "params": {}, "log": [], "threads": 1, "resources": {"mem": 32, "time": 1}, "jobid": 9, "cluster": {}}
cd /srv/gsfs0/projects/bhatt/fiona/bacteremia/strainsifter && \
/home/tamburin/tools/miniconda3/bin/python \
-m snakemake pileup/B09.pileup --snakefile /srv/gsfs0/projects/bhatt/fiona/bacteremia/strainsifter/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /srv/gsfs0/projects/bhatt/fiona/bacteremia/strainsifter/.snakemake/tmp.95cvabvn passed_samples/B09.bam /home/tamburin/fiona/bacteremia/1.assemble_trimmed/filtered_assemblies/B09.fasta /home/tamburin/fiona/bacteremia/1.assemble_trimmed/filtered_assemblies/B09.fasta.fai --latency-wait 5 \
--benchmark-repeats 1 --attempt 1 \
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --no-hooks --nolock --timestamp  --force-use-threads  --allowed-rules pileup  && exit 0 || exit 1

