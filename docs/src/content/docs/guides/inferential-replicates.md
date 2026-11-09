---
title: Inferential replicates
description: Bootstrap replicates of the abundance estimates, written as a Parquet table.
---

`oarfish` has the ability to compute [_inferential replicates_](https://academic.oup.com/nar/article/47/18/e105/5542870) of its quantification estimates. This is performed by bootstrap sampling of the original read mappings, and subsequently performing inference under each resampling. These inferential replicates allow assessing the variance of the point estimate of transcript abundance, and can lead to improved differential analysis at the transcript level, if using a differential testing tool that takes advantage of this information.

The generation of inferential replicates is controlled by the `--num-bootstraps` argument to `oarfish`. The default value is `0`, meaning that no inferential replicates are generated. If you set this to some value greater than `0`, the requested number of inferential replicates will be generated. It is recommended, if generating inferential replicates, to run `oarfish` with multiple threads, since replicate generation is highly-parallelized.

If replicates are generated, they are written to a [`Parquet`](https://parquet.apache.org/) file, starting with the specified output stem and ending with `infreps.pq` (see [Output files](/oarfish/reference/output/)).
