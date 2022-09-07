#!/bin/bash

#on cbsumm10 (Cornell BioHPC server that has BLAST2GO license)
mkdir /workdir/$USER
mkdir /workdir/$USER/20_BLAST2GO
cd /workdir/$USER/20_BLAST2GO

#copy blast2go files
cp /shared_data/blast2go/annotation.prop ./
cp /shared_data/blast2go/go.obo ./

#copy InterProScan output and BLAST output
scp arj66@cbsugaurav.biohpc.cornell.edu:/workdir/arj66/EUPHORBIA_PEPLUS_GENOME/18_InterProScan/ips_output.xml ./
scp arj66@cbsugaurav.biohpc.cornell.edu:/workdir/arj66/EUPHORBIA_PEPLUS_GENOME/19_BLAST_for_BLAST2GO/uniref_blastresults.xml ./

nohup /usr/local/blast2go/blast2go_cli.run -properties annotation.prop -useobo go.obo -loadblast uniref_blastresults.xml -loadips50 ips_output.xml -mapping -annotation -statistics all -saveannot Euphorbia_peplus -savereport Euphorbia_peplus_blast2go -tempfolder ./ >& annotatelogfile &

dos2unix Euphorbia_peplus.annot

cd ..

scp -rp 20_BLAST2GO arj66@cbsugaurav.biohpc.cornell.edu:/workdir/arj66/EUPHORBIA_PEPLUS_GENOME/
