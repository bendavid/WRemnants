# WRemnants

## Instructions:

Get the code (after forking from the central WMass repository)
```bash
git clone --recursive https://github.com/WMass/WRemnants.git
cd WRemnants/
MY_GIT_USER=$(git config user.github)
git remote set-url origin git@github.com:$MY_GIT_USER/WRemnants.git
git remote add wmass-central git@github.com:WMass/WRemnants.git
git push -f origin main
```

Run the code
```bash
singularity run /scratch/singularity/pythonrootarchdevrolling # needs to be on lxplus8s10
source WRemnants/setup.sh
python WRemnants/scripts/histmakers/mw_with_mu_eta_pt.py
```
