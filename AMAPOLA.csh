#mkdir WORK
rsync -auvz skameno@ssh.alma.cl:/home/skameno/public_html/Grid/Stokes/ WORK/
cd WORK
#ls *Flux.log > fileList
ls -la *Flux.log | awk '{if ($5 >= 1024) print $9}' > fileList
Rscript ../StokesFluxStat.R fileList
Rscript ../ReadAeff.R *Flux.log
scp Flux.Rdata skameno@ssh.alma.cl:/home/skameno/public_html/Grid/Stokes/
scp AeDF.Rdata skameno@ssh.alma.cl:/home/skameno/public_html/Grid/Stokes/
scp amapola.txt skameno@ssh.alma.cl:/home/skameno/public_html/AMAPOLA/
Rscript ../PlotFlux.R
scp PolQuery.CSV skameno@ssh.alma.cl:/home/skameno/public_html/AMAPOLA/
scp *.html skameno@ssh.alma.cl:/home/skameno/public_html/AMAPOLA/
cd ..
