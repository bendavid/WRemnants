#!/bin/bash

all_folders=$(xrdfs root://eosuser.cern.ch// ls -l $1 | grep '^d' | sort -k 5 -r | awk '{print $7}')
newest_folders=$(echo "$all_folders" | head -n $2)
for folder in $all_folders; do
if ! [[ ${newest_folders[@]} =~ ${folder} ]]; then
    xrdfs root://eosuser.cern.ch// ls -l $folder | while read line; do
    if [ "${line:0:1}" == "d" ]; then
        sub_folder=$(echo $line | awk '{print $7}')
        xrdfs root://eosuser.cern.ch// ls -l $sub_folder | while read sub_line; do
        if [ "${sub_line:0:1}" == "d" ]; then
            echo "Not expected to find a folder here"
        else
            file=$(echo $sub_line | awk '{print $7}')
            xrdfs root://eosuser.cern.ch// rm $file
        fi
        done
        xrdfs root://eosuser.cern.ch// rmdir $sub_folder
    else
        file=$(echo $line | awk '{print $7}')
        xrdfs root://eosuser.cern.ch// rm $file
    fi
    done
    xrdfs root://eosuser.cern.ch// rmdir $folder
fi
done