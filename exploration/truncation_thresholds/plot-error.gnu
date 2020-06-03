set terminal epslatex standalone size 4.5in,2.5in color
set output "plot-error.tex"

set key vertical maxrows 2 samplen 3

col1= 57*65536 + 106*256 + 177 #(57,106,177)
col2=218*65536 + 124*256 +  48 #(218,124,48)
col3= 62*65536 + 150*256 +  81 #(62,150,81)
col4=204*65536 +  37*256 +  41 #(204,37,41)
col5= 83*65536 +  81*256 +  84 #(83,81,84)
col6=107*65536 +  76*256 + 154 #(107,76,154)
col7=146*65536 +  36*256 +  40 #(146,36,40)
col8=148*65536 + 139*256 +  61 #(148,139,61)

set style line 1 linecolor rgb col1 linewidth 5 dashtype 1 pointtype 6 pointsize default
set style line 2 linecolor rgb col2 linewidth 5 dashtype 1 pointtype 4 pointsize default
set style line 3 linecolor rgb col3 linewidth 5 dashtype 1 pointtype 2 pointsize default
set style line 4 linecolor rgb col4 linewidth 5 dashtype 1 pointtype 8 pointsize default

set border 31 linewidth 2

set xlabel 'Number of atoms'
set ylabel 'Error (meV/atom)'

set logscale y 10
set yrange [0.0001:10]

hartree = 27.2116 # eV

plot 'results.dat' using ($2):(abs($3-$4)/$2*hartree*1000) title '$\epsilon_\textrm{filter}=10^{-4}$' with lp ls 1, \
     'results.dat' using ($2):(abs($3-$5)/$2*hartree*1000) title '$\epsilon_\textrm{filter}=10^{-5}$' with lp ls 2, \
     'results.dat' using ($2):(abs($3-$6)/$2*hartree*1000) title '$\epsilon_\textrm{filter}=10^{-6}$' with lp ls 3, \
     'results.dat' using ($2):(abs($3-$7)/$2*hartree*1000) title '$\epsilon_\textrm{filter}=10^{-7}$' with lp ls 4

set output "bla"
!latex plot-error.tex
!dvips plot-error.dvi
!ps2eps -f plot-error.ps
!rm *.dvi *.aux *.ps *-inc.eps *.tex bla
