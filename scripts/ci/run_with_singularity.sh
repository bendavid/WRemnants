#!/bin/bash
SINGULARITY_BIND="/scratch,/eos,/cvmfs" singularity run /cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/bendavid/cmswmassdocker/wmassdevrolling\:v8 $@
