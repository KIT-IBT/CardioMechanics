#!/bin/bash

# Cellmodels used
EVFILE="OHaraRudy_endo OHaraRudy_mid OHaraRudy_epi CourtemancheEtAl CourtemancheEtAl_APG CourtemancheEtAl_CT"
STRETCH="1.0"
MODE="II"

BCL="0.8"
SIMLEN="800"
TOUTBEGIN="799.2"

i=0
for model in $EVFILE; do
	for stretch in $STRETCH; do
		echo "Simulating $model and stretch = $stretch ..."
		ElphyModelTest -evfile ./${MODE}/${model}.ev -text ${BCL} -textbegin 0.0 -tbegin 0 -tend ${SIMLEN} -stretch ${stretch} -toutbegin ${TOUTBEGIN} -outfile ./${MODE}/${model}_${stretch}.txt -progress #> /dev/null 2>&1
		if [ $? -ne 0 ]
		then
			echo "Something went wrong. Please repeat last step manually and check for errors."
			exit $?
		fi
		echo "Done"
	done
	((i++))
done