rsync -auvz skameno@ssh.alma.cl:/home/skameno/public_html/Grid/Stokes/ WORK/
cd WORK
\rm *.Rdata
\rm *.png
\rm *.pdf
\rm *.table
\rm *.CSV
\rm *.csv
\rm *.html
\rm -rf *_files
\rm -rf AMAPOLA
\ls -la *Flux.log | awk '{if ($5 >= 1024) print $9}' > fileList
Rscript ../StokesFluxStat.R fileList
Rscript ../polCalCSV.R
Rscript ../SSOUID.R
Rscript ../ReadAeff.R fileList
Rscript ../StokesText.R
mkdir AMAPOLA
mkdir AMAPOLA/PDF
mkdir AMAPOLA/PNG
mkdir AMAPOLA/Table
cp ../HTML/*.html AMAPOLA
cp -r ../HTML/resources AMAPOLA
cp -r ../HTML/images AMAPOLA
cp Flux.Rdata AMAPOLA
cp AeDF.Rdata AMAPOLA
cp Dterm.Rdata AMAPOLA
cp amapola.txt AMAPOLA
cp PolCalBand*.csv AMAPOLA
cp *.png AMAPOLA/PNG
cp *.pdf AMAPOLA/PDF
mv *.table AMAPOLA/Table
cp *.html AMAPOLA
mv common.flux_files AMAPOLA
mv J*.flux_files AMAPOLA
mv Band*LSTplot_files AMAPOLA
#rsync -a --safe-links --delete AMAPOLA skameno@ssh.alma.cl:/home/skameno/public_html/
rsync -a --safe-links AMAPOLA skameno@ssh.alma.cl:/home/skameno/public_html/
cd ..
