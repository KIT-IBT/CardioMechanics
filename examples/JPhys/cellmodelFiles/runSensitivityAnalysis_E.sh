#!/bin/bash

# Cellmodels used
EVFILE="ORd_0_E ORd_1_E ORd_2_E ORd_3_E CRN_0_E CRN_1_E CRN_2_E CRN_3_E"
STRETCH="1.0"
MODE="Sensitivity"

i=0
for model in $EVFILE; do
	for stretch in $STRETCH; do
		echo "Simulating $model and stretch = $stretch ..."
		ElphyModelTest -evfile ./${MODE}/${model}.ev -text 0.8 -textbegin 100.0 -textend 101.0 -tbegin 0 -tend 3.5 -clamp ./${MODE}/stretch_clamp.txt -toutbegin 0.0 -outfile ./${MODE}/${model}_${stretch}.txt -progress -lp #> /dev/null 2>&1
		if [ $? -ne 0 ]
		then
			echo "Something went wrong. Please repeat last step manually and check for errors."
			exit $?
		fi
		echo "Done"
	done
	((i++))
done