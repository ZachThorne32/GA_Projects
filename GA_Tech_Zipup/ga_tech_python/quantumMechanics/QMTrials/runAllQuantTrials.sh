for trial in trial*
do
    echo $trial

    cd $trial

    python3 runQuantTrials.py

    cd ..

done
