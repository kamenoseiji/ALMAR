#mkdir WORK
rsync -auvz skameno@ssh.alma.cl:/home/skameno/public_html/Grid/Stokes/ WORK/
cd WORK
#scp skameno@ssh.alma.cl:/home/skameno/public_html/Grid/Stokes/\*Flux.log .
ls *Flux.log > fileList
Rscript ../StokesFluxStat.R fileList
scp Flux.Rdata skameno@ssh.alma.cl:/home/skameno/public_html/Grid/Stokes/
scp amapola.txt skameno@ssh.alma.cl:/home/skameno/public_html/AMAPOLA/
Rscript ../PlotFlux.R
scp *.html skameno@ssh.alma.cl:/home/skameno/public_html/AMAPOLA/
cd ..
#rm -rf WORK
