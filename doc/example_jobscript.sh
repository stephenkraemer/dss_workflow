#!/usr/bin/env bash
# properties = {properties}
# Make sure that conda is available, and load the dss_workflow environment
# (which you need to create the first time you use this,
# the workflow does not create the conda environment by itself, i.e.
# no use of the conda directive)
PATH=/icgc/dkfzlsdf/analysis/B080/kraemers/anaconda3/bin:$PATH
source activate dmr_calling
{exec_job}
