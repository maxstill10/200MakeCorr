#!/bin/bash

FILE="run16_1_runnum.list"

#6sem, Centred, Flatt
mode="6sem"

runnumbers=()

while IFS= read -r line; do
   numbers+=("$line")   
done < "$FILE"

for number in "${numbers[@]}"; do
    mkdir -p /star/data01/pwg/mmorozov/Polarization/200GeV/log/$number
    mkdir -p /star/data01/pwg/mmorozov/Polarization/200GeV/data/$number
    OutFile="submiter_$number.xml"


    cat > "$OutFile" <<EOF
<!-- Task description -->
<job maxFilesPerProcess="5" filesPerHour="8" fileListSyntax="xrootd">
<!-- Decription of the task -->
<shell>singularity exec -e -B /direct -B /star -B /afs -B /gpfs -B /sdcc/lustre02 /cvmfs/star.sdcc.bnl.gov/containers/rhic_sl7.sif</shell>
<command>
    <!-- setenv NODEBUG yes -->

    ln -s /star/u/mmorozov/200MakeCorr/PicoDstAnalyzer.C .
    ln -s /star/u/mmorozov/200MakeCorr/MaxRunner.C .
    ln -s /star/u/mmorozov/200MakeCorr/StRoot .
    ln -s /star/u/mmorozov/200MakeCorr/.sl73_gcc485 .
    ln -s /star/u/mmorozov/200MakeCorr/libStPicoDst.so .

    root -b -l -q MaxRunner.C\(\"\$FILELIST\",\"\$JOBID.root\",\"$mode\",\"$number\",\"run14\"\)
    
    </command>

<!-- Get input files from get_file_lis.pl -->
<input URL="catalog:star.bnl.gov?filetype=daq_reco_picoDst,trgsetupname=AuAu_200_production_2014,production=P16id,runnumber=$number,filename~st_physics,storage!=hpss" preferStorage="local" singleCopy="true" nFiles="all" /> <!-- nFiles="all"  -->
<stdout URL="file:/star/data01/pwg/mmorozov/Polarization/200GeV/log/$number/\$JOBID.out"/>
<stderr URL="file:/star/data01/pwg/mmorozov/Polarization/200GeV/log/$number/\$JOBID.err"/>
<output fromScratch="\$JOBID.root" toURL="file:/star/data01/pwg/mmorozov/Polarization/200GeV/data/$number"/>


  <Generator>
    <Location>/star/data01/pwg/mmorozov/Polarization/200GeV/log/$number</Location>
  </Generator>


</job>

EOF

   star-submit-beta $OutFile
done
