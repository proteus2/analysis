#! /bin/bash

machs=( Meteo Cluster )
hdirs=( /data1/kyh /data/kyh )
        # meteo, cluster

chkerr(){
if [ $1 -ne 0 ] ; then echo "   Failed: $2" ; continue ; fi
}

cd $HOME
cp .cshrc cshrc

imach=0

for hh in ${hdirs[*]}; do

  echo ; echo " Copy to the Machine - ${machs[$imach]}"

  if [ -d "$hh" ] ; then
    cp cshrc $hh/
    chkerr $? "cp cshrc $hh/"
    mv $hh/cshrc $hh/.cshrc
    chkerr $? "mv $hh/cshrc $hh/.cshrc"
  else echo '   Failed: not connected.'
  fi

  imach=$(( $imach + 1 ))

done

rm -f cshrc

exit 0

