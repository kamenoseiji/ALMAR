mkdir WORK
cd WORK
scp skameno@ssh.alma.cl:/home/skameno/public_html/Grid/Stokes/\*Flux.log .
Rscript ../StokesFluxStat.R *Flux.log
scp Flux.Rdata skameno@ssh.alma.cl:/home/skameno/public_html/Grid/Stokes/
Rscript ../PlotFlux.R
scp *.html skameno@ssh.alma.cl:/home/skameno/public_html/AMAPOLA/
cd ..
rm -rf WORK
