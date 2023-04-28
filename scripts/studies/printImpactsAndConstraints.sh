source local_setup.sh
declare -a combine_folders=("WMass_eta_pt_non-closure-Apr26")
declare -a fitresults_files=("fitresults-w_charge_inclusive" "fitresults-w_plus" "fitresults-w_minus")
#declare -a postfixes=("_A-M-separated" "_A-M-combined" "_A-M-separated-dep-only" "_A-M-combined-dep-only")
declare -a postfixes=("_A-M-combined-one-group-ind")

for combine_folder in "${combine_folders[@]}"; do
    for fitresults_file in "${fitresults_files[@]}"; do
        for postfix in "${postfixes[@]}"; do
            echo ${combine_files_dir}/${combine_folder}/${fitresults_file}${postfix}.root 
            echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
            echo "-----------------------Impacts--------------------------"
            python3 scripts/combine/printImpacts.py -f ${combine_files_dir}/${combine_folder}/${fitresults_file}${postfix}.root | grep '.*nonClosure.*'
            echo "---------------------Constraints------------------------"
            python3 scripts/combine/studies/printConstraints.py -i ${combine_files_dir}/${combine_folder}/${fitresults_file}${postfix}.root
            echo ""
        done
    done
done
