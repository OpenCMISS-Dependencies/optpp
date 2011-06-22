#!/bin/sh

/bin/rm -f failures.log details.log results.dat
touch failures.log details.log results.dat

echo "-------------------------------------------------------" >> details.log;
echo " " >> details.log;
echo "OPT++ xml test results: " `date` >> details.log ;
echo " " >> details.log;
echo "-------------------------------------------------------" >> details.log;

runcommand=$1

for reg_test in tstCG.xml tstFD.xml tstNewton.xml tstQNewton.xml tstPDS.xml tstTRPDS.xml tstNIPS.xml tstHock14.xml tstHock65.xml
do
#
#   First check to see if the program runs at all
#
    if $runcommand ../../bin/opt++.e $reg_test; then
	:
    else
	echo -e "$reg_test crashed" >> failures.log;
    fi
#
#   Now check to see if the test problems passed/failed
#
#    set status = `grep PASSED $reg_test.out >> details.log`
#    set status = `grep FAILED $reg_test.out >> details.log`
#    set status = `grep PASSED $reg_test.out >> results.dat`
#    set status = `grep FAILED $reg_test.out >> results.dat`
done
#
#   If any of the tests resulted in a failure store result in failures log
#   The flag -s checks for existence and nonzero size file
#
#set status = `grep -l FAILED *.out >> failures.log`
#
#  Print out messages depending on whether all tests passed or not
#
if [ -s failures.log ]; then
    echo "----------------------------------------------------------------";
    echo " ";
    echo "  OPT++ xml test results: FAILURES DETECTED";
    echo "  ";
    echo "      " `date` ;
    echo " ";
    echo "The following tests failed:";
    echo `cat failures.log`;
    echo " "
    echo "Please check the *.out files for details";
    echo "---------------------------------------------------------------";

    echo "----------------------------------------------------" >> details.log;
    echo "OPT++ xml test results: FAILURES DETECTED" >> details.log;
    echo "----------------------------------------------------" >> details.log;

else
    echo "----------------------------------------------------------------";
    echo " ";
    echo "  OPT++ xml test results: all tests COMPLETED";
    echo "  Please check the results in tests/xml/*.out";
    echo "  to ensure that they are correct";
    echo "  ";
    echo "      " `date` ;
    echo "----------------------------------------------------------------";

    echo "----------------------------------------------------" >> details.log;
    echo "OPT++ xml test results: all tests COMPLETED" >> details.log;
    echo "----------------------------------------------------" >> details.log;

fi

exit 0
