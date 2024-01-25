#!/bin/bash

# Cellmodels used
EVFILE="ORd_0 ORd_1 ORd_2 ORd_3 ORd_4 ORd_5 ORd_6 ORd_7 CRN_0 CRN_1 CRN_2 CRN_3 CRN_4 CRN_5 CRN_6 CRN_7"
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