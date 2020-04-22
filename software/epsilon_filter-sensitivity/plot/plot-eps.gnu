set terminal epslatex standalone size 4.5in,3.2in color

set output "plot-eps.tex"

col1= 57*65536 + 106*256 + 177 #(57,106,177)
col2=218*65536 + 124*256 +  48 #(218,124,48)
col3= 62*65536 + 150*256 +  81 #(62,150,81)
col4=204*65536 +  37*256 +  41 #(204,37,41)
col5= 83*65536 +  81*256 +  84 #(83,81,84)
col6=107*65536 +  76*256 + 154 #(107,76,154)
col7=146*65536 +  36*256 +  40 #(146,36,40)
col8=148*65536 + 139*256 +  61 #(148,139,61)

set border 31 linewidth 2

set logscale y 10
set xrange [-9:-2]

set xlabel '$\epsilon_\textrm{filter}$'
set ylabel 'Time (s)'

set format x "$10^{%.0f}$"

plot 'SM.dat' using 1:4 title 'Submatrix Method' with lp lw 5 lc rgb col1 pt 6, \
     'NS.dat' using 1:4 title 'Newton-Schulz' with lp lw 5 lc rgb col2 pt 4

set output "bla"
!latex plot-eps.tex
!dvips plot-eps.dvi
!ps2eps -f plot-eps.ps
!rm *.dvi *.aux *.ps *-inc.eps *.tex bla
