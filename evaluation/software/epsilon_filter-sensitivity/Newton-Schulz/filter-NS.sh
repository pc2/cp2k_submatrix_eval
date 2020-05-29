#!/bin/bash

for i in $(seq -8 1 -3); do
	for j in 1.0 1.7 3.1 5.6; do
	    echo -n "${j}E${i}: "
		(grep "Tr(DK)" ${j}E${i}-NS.out | awk '{print $3}'; grep "SCF   " ${j}E${i}-NS.out | awk '{print $3}'; grep "_fixed_mu" ${j}E${i}-NS.out | awk '{print $7}')|xargs
	done
done
