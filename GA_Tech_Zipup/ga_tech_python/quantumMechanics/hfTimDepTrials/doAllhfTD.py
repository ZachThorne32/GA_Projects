for trial in trial*
do
    echo $trial

    cd $trial

    python3 runhfTD.py

    cd ..

done
