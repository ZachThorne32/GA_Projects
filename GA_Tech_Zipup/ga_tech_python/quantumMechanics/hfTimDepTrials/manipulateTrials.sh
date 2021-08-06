#manipulate AN (duration of trial) and BN (scaler for trials)
for dir in trial*
do
    echo $dir

    cd $dir

    cp ../runhfTD.py .

    let N=`echo $dir | awk -Fl '{print $2}'`
    echo $N
    let AN=$(($N))
    echo $AN
    let BN=$(($N+3))
    if [[ $(($BN%2)) -eq 0 ]]
    then
        echo $BN
    else
        BN=$(($BN+1))
        echo $BN
    fi

    sed -i "s/aaa/$AN/g" runQuantTrials.py
    sed -i "s/bbb/$BN/g" runQuantTrials.py

    cd ..

done
