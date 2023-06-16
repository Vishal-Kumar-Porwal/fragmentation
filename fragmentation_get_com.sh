#!/bin/bash
rm -r XYZ

mkdir XYZ

echo -n 'Enter the name of the molecule:'
read molecule

echo -n 'Enter the name and path of xyz file:'
read xyz_file                                                          

echo -n 'Enter the name and path of connection file:'
read connection_file

printf "$xyz_file\n$connection_file" | python3  fragmentation.py 

echo 'creating com files'
    
cd XYZ

for filexyz in *.xyz
do
    # 
    line2=$filexyz
    echo $filexyz
    line_1=$(echo $filexyz | rev | cut -c5-| rev)
    echo $line_1
    echo $line2
    line3=$(echo $line2 | rev | cut -c 5)
    echo $line3
    #line_2=$(echo $line_1 | sed 's/ /_/g')
    #echo $line_2
 
    echo  '%nprocshared=8'                                                 >> $molecule"_"$line_1.com
    echo  '%mem=112GB       '                                              >> $molecule"_"$line_1.com
    echo  "%chk=rotor_scan_$molecule"_"$line_1.chk"                          >> $molecule"_"$line_1.com
    echo  '#p opt freq 6-311+g(2d,d,p) b2plypd3'                   >> $molecule"_"$line_1.com
    echo  '               '                                                >> $molecule"_"$line_1.com
    echo  " $molecule bond breaking  :"  $line_1                             >> $molecule"_"$line_1.com
    echo  '               '                                                >> $molecule"_"$line_1.com
    echo  "0  $line3           "                                                >> $molecule"_"$line_1.com
    sed -n '3,$w tempfile' $filexyz
    echo tempfile
    cat tempfile                                                           >> $molecule"_"$line_1.com                                                          
    echo '                '                                                >> $molecule"_"$line_1.com
    #echo $line                                                             >> rotor_scan_$molecule"_"$line_1.com
    #echo '                '                                                >> rotor_scan_$molecule"_"$line_1.com
    #echo '                '                                                >> rotor_scan_$molecule"_"$line_1.com
    


 
done 



