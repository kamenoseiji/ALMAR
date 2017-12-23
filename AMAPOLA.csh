mkdir WORK
cd WORK
Rscript ~/Programs/ALMAR/PlotFlux.R
scp *.html skameno@ssh.alma.cl:/home/skameno/public_html/AMAPOLA/
cd ..
rm -rf WORK
