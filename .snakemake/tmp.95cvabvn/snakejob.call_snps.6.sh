#!/bin/bash
# properties = {"rule": "call_snps", "local": false, "input": ["pileup/B09.pileup"], "output": ["snp_calls/B09.tsv"], "wildcards": ["B09"], "params": {"min_cvg": 5, "min_freq": 0.8, "min_qual": 20}, "log": [], "threads": 1, "resources": {"mem": 32, "time": 2}, "jobid": 6, "cluster": {}}
cd /srv/gsfs0/projects/bhatt/fiona/bacteremia/strainsifter && \
/home/tamburin/tools/miniconda3/bin/python \
-m snakemake snp_calls/B09.tsv --snakefile /srv/gsfs0/projects/bhatt/fiona/bacteremia/strainsifter/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /srv/gsfs0/projects/bhatt/fiona/bacteremia/strainsifter/.snakemake/tmp.95cvabvn pileup/B09.pileup --latency-wait 5 \
--benchmark-repeats 1 --attempt 1 \
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --no-hooks --nolock --timestamp  --force-use-threads  --allowed-rules call_snps  && exit 0 || exit 1

