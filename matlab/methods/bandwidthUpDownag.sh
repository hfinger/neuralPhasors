#!/bin/bash
#fname=$(cat /dev/urandom | tr -cd 'a-f0-9' | head -c 32)
fname=$(date +%s | sha256sum | base64 | head -c 32)
date1=$(date +"%s%N")
dd if=/dev/urandom of=/net/store/nbp/projects/phasesim/temp/speedtest/randfiles/$fname bs=$1M count=1 >& /dev/null
cp -r /net/store/nbp/projects/phasesim/temp/speedtest/randfiles/$fname /tmp/$fname
rm -r /net/store/nbp/projects/phasesim/temp/speedtest/randfiles/$fname
rm -r /tmp/$fname
date2=$(date +"%s%N")
diff=$(($date2-$date1))
echo $diff
