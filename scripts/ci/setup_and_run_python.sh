. ./setup.sh
if [ -e "/home/${USER:0:1}/${USER}/.keytab" ]; then
    echo "Initialize kerberos"
    kinit -kt /home/${USER:0:1}/${USER}/.keytab ${USER}@CERN.CH
fi
python3 $@
