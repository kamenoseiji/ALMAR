mkdir WORK
cd WORK
scp skameno@ssh.alma.cl:/home/skameno/public_html/Grid/Stokes/\*Flux.log .
ls *Flux.log > fileList
Rscript ../StokesFluxStat.R fileList
scp Flux.Rdata skameno@ssh.alma.cl:/home/skameno/public_html/Grid/Stokes/
Rscript ../PlotFlux.R
scp amapola.txt skameno@ssh.alma.cl:/home/skameno/public_html/AMAPOLA/
scp *.html skameno@ssh.alma.cl:/home/skameno/public_html/AMAPOLA/
cd ..
rm -rf WORK
