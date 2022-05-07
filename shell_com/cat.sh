cd rawdata

cat */SJT-0710-1*.R1.fastq.gz  > EPS2C-5.R1.fastq.gz
cat */SJT-0710-1*.R2.fastq.gz  > EPS2C-5.R2.fastq.gz

cat */SJT-0710-2*.R1.fastq.gz  > EPS8C-4.R1.fastq.gz
cat */SJT-0710-2*.R2.fastq.gz  > EPS8C-4.R2.fastq.gz

cat */SJT-0710-3*.R1.fastq.gz  > ESC.R1.fastq.gz
cat */SJT-0710-3*.R2.fastq.gz  > ESC.R2.fastq.gz

cat */SJT-0710-4*.R1.fastq.gz  > TdEPS.R1.fastq.gz
cat */SJT-0710-4*.R2.fastq.gz  > TdEPS.R2.fastq.gz

cd ..
mkdir 2nd.trim.map
cd  2nd.trim.map

mkdir EPS2C-5
mkdir EPS8C-4
mkdir ESC
mkdir TdEPS

cp /home1/liziyi/EPS/WGBS/final.trim.map.sh EPS2C-5/
cp /home1/liziyi/EPS/WGBS/final.trim.map.sh EPS8C-4/
cp /home1/liziyi/EPS/WGBS/final.trim.map.sh ESC/
cp /home1/liziyi/EPS/WGBS/final.trim.map.sh TdEPS/



