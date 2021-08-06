#manipulate AN (duration of trial) and BN (scaler for trials)
for dir in trial*
do
    echo $dir

    cd $dir

    cp ../runMonteCarlo.py .

    let N=`echo $dir | awk -Fl '{print $2}'`
    echo $N
    let AN=$N*5000000
    echo $AN
    let BN=$N/2+1
    echo $BN

    sed -i "s/aaa/$AN/g" runMonteCarlo.py
    sed -i "s/bbb/$BN/g" runMonteCarlo.py

    cd ..

done
