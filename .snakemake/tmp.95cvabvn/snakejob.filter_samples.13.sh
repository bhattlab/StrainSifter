#!/bin/bash
# properties = {"rule": "filter_samples", "local": false, "input": ["coverage/B09.cvg"], "output": ["passed_samples/B09.bam"], "wildcards": ["B09"], "params": {"min_cvg": 5, "min_perc": 0.5}, "log": [], "threads": 1, "resources": {"mem": 1, "time": 1}, "jobid": 13, "cluster": {}}
cd /srv/gsfs0/projects/bhatt/fiona/bacteremia/strainsifter && \
/home/tamburin/tools/miniconda3/bin/python \
-m snakemake passed_samples/B09.bam --snakefile /srv/gsfs0/projects/bhatt/fiona/bacteremia/strainsifter/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /srv/gsfs0/projects/bhatt/fiona/bacteremia/strainsifter/.snakemake/tmp.95cvabvn coverage/B09.cvg --latency-wait 5 \
--benchmark-repeats 1 --attempt 1 \
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --no-hooks --nolock --timestamp  --force-use-threads  --allowed-rules filter_samples  && exit 0 || exit 1

