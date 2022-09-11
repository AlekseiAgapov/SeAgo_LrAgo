#!/usr/bin/bash

./pipeline.sh -p 6 -d ../../syn/first/6_h/no_cipro/
./pipeline.sh -p 6 -d ../../syn/first/6_h/cipro/
./pipeline.sh -p 6 -d ../../syn/first/12_h/no_cipro/
./pipeline.sh -p 6 -d ../../syn/first/12_h/cipro/
echo 'syn is done'
./pipeline.sh -p 6 -d ../../lro/first/6_h/no_cipro/
./pipeline.sh -p 6 -d ../../lro/first/6_h/cipro/
./pipeline.sh -p 6 -d ../../lro/first/12_h/no_cipro/
./pipeline.sh -p 6 -d ../../lro/first/12_h/cipro/
echo 'first lro is over'
./pipeline.sh -p 6 -d ../../lro/second/6_h/no_cipro/
./pipeline.sh -p 6 -d ../../lro/second/6_h/cipro/
./pipeline.sh -p 6 -d ../../lro/second/12_h/no_cipro/
./pipeline.sh -p 6 -d ../../lro/second/12_h/cipro/
echo 'thats it'
