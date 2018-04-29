#!/bin/bash
# properties = {"rule": "calc_coverage", "local": false, "input": ["filtered_bam/B09.filtered.bam"], "output": ["coverage/B09.cvg"], "wildcards": ["B09"], "params": {"cvg": 5}, "log": [], "threads": 1, "resources": {"mem": 16, "time": 1}, "jobid": 16, "cluster": {}}
cd /srv/gsfs0/projects/bhatt/fiona/bacteremia/strainsifter && \
/home/tamburin/tools/miniconda3/bin/python \
-m snakemake coverage/B09.cvg --snakefile /srv/gsfs0/projects/bhatt/fiona/bacteremia/strainsifter/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /srv/gsfs0/projects/bhatt/fiona/bacteremia/strainsifter/.snakemake/tmp.95cvabvn filtered_bam/B09.filtered.bam --latency-wait 5 \
--benchmark-repeats 1 --attempt 1 \
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --no-hooks --nolock --timestamp  --force-use-threads  --allowed-rules calc_coverage  && exit 0 || exit 1

