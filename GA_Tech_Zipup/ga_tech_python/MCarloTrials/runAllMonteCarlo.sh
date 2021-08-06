for trial in trial*
do
    echo $trial

    cd $trial

    python runMonteCarlo.py

    cd ..

done
