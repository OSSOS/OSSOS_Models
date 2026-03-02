#! /bin/bash

\rm -f OSSOS-model.out ModelUsed.dat Inner

make

time ./Inner <<EOF
123476790
0
8.0d0
.true.
./OSSOS-model.out
EOF
if [ $? != 0 ]; then
    status=$?
    echo "================================================================="
    echo "=======================  Test FAILED!   ========================="
    echo "================================================================="
    exit ${status}
fi

tail --lines=+3 ModelUsed.dat > zzz1
tail --lines=+3 ModelUsed-check-8.66.dat > zzz2
diff zzz1 zzz2 > /dev/null 2>&1
status=$?
\rm -f zzz1 zzz2
\rm -f OSSOS-model.out ModelUsed.dat
if [ $status != 0 ]; then
    echo "================================================================="
    echo "=======================  Test FAILED!   ========================="
    echo "================================================================="
    exit ${status}
else
    echo "================================================================="
    echo "=====================  Test successful!  ========================"
    echo "================================================================="
fi

exit
