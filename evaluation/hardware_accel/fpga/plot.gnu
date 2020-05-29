set terminal epslatex standalone size 4.5in,3.2in color

set key autotitle column nobox samplen 1 noenhanced
unset title
set output "plot.tex"
set multiplot

lw=5

col1= 57*65536 + 106*256 + 177 #(57,106,177)
col2=218*65536 + 124*256 +  48 #(218,124,48)
col3= 62*65536 + 150*256 +  81 #(62,150,81)
col4=204*65536 +  37*256 +  41 #(204,37,41)
col5= 83*65536 +  81*256 +  84 #(83,81,84)
col6=107*65536 +  76*256 + 154 #(107,76,154)
col7=146*65536 +  36*256 +  40 #(146,36,40)
col8=148*65536 + 139*256 +  61 #(148,139,61)

f1="../gpu/half.dat"
f2="../gpu/halfsingle.dat"
f3="../gpu/single.dat"
f3f="single.dat"
f4="../gpu/double.dat"

#set logscale x
#set logscale y
set style fill solid noborder
set tics scale 2

#set xtics ("1" 1,"2" 2,"3" 1e-4,"$10^{-3}$" 1e-3,"$10^{-2}$" 1e-2,"$10^{-1}$" 1e-1,"$10^{-0}$" 1e-0)
#set ytics ("$10^{0}$" 1, "" 2 1,"" 3 1,"" 4 1,"" 5 1,"" 6 1,"" 7 1,"" 8 1,"" 9 1, "$10^{1}$" 10, "" 20 1, "" 30 1, "" 40 1, "" 50 1, "" 60 1, "" 70 1, "" 80 1, "" 90 1,"$10^{2}$" 100,"" 200 1,"" 300 1,"" 400 1,"" 500 1,"" 600 1,"" 700 1,"" 800 1,"" 900 1,"$10^{3}$" 1000)
set xlabel "Sign iteration"
set ylabel "Energy diff. from FP64 [meV/atom]" offset 1
set key font ",12"

set border 31 linewidth 2

#set key tmargin 2.0
#set key at 0.7,750
set key spacing 2.0

EREF=-74.249194445722438
#EREF=0
conv=27.211/(32*3)*1000
plot [4.5:15.5][-4:6] \
f1 u ($4):(conv*($5-EREF)) w lp lc rgb col1 lw lw pt 6 t "GPU FP16", \
f2 u ($4):(conv*($5-EREF)) w lp lc rgb col2 lw lw pt 4 t "GPU FP16'", \
f3 u ($4):(conv*($5-EREF)) w lp lc rgb col3 lw lw pt 2 t "GPU FP32", \
f4 u ($4):(conv*($5-EREF)) w lp lc rgb col4 dashtype "." lw lw pt 8 t "GPU FP64", \
f3f u ($4):(conv*($5-EREF)) w lp lc rgb col5 dashtype "-.-" lw lw pt 2 t "FPGA FP32"

set origin .17, .47
set size .55,.5
#clear
unset key
#unset grid
#unset object
#unset arrow
set ylabel "abs. Energy diff." offset 1
set xlabel ""
EREF=0
conv=1.0
#set xtics .1
#set ytics .5
#set bmargin 1
#set tmargin 1
#set lmargin 3
#set rmargin 1
set logscale y
#set ytics ("$10^{-12}$" 1e-12,"$10^{-10}$" 1e-10,"$10^{-8}$" 1e-8,"$10^{-6}$" 1e-6,"$10^{-4}$" 1e-4,"$10^{-2}$" 1e-2)
set format y "$10^{%T}$"

set logscale y

EREF=-74.249194445722438
#EREF=0
conv=27.211/(32*3)*1000
plot [4.5:15.5][1e-3:6] \
f1 u ($4):(abs(conv*($5-EREF))) w lp lc rgb col1 lw lw pt 6 t "GPU FP16", \
f2 u ($4):(abs(conv*($5-EREF))) w lp lc rgb col2 lw lw pt 4 t "GPU FP16'", \
f3 u ($4):(abs(conv*($5-EREF))) w lp lc rgb col3 lw lw pt 2 t "GPU FP32", \
f4 u ($4):(abs(conv*($5-EREF))) w lp lc rgb col4 lw lw pt 8 t "GPU FP64", \
f3f u ($4):(abs(conv*($5-EREF))) w lp lc rgb col5  lw lw pt 2 t "FPGA FP32"

unset multiplot


set output "bla"
!latex plot.tex
!dvips plot.dvi
!mv plot.ps energy.ps
!rm energy.eps
!ps2eps energy.ps
!rm *.dvi *.aux 





reset
set terminal epslatex standalone size 4.5in,3.2in color

set key autotitle column nobox samplen 1 noenhanced
unset title
set output "plot.tex"

lw=5

f1="../gpu/half.dat"
f2="../gpu/halfsingle.dat"
f3="../gpu/single.dat"
f3f="single.dat"
f4="../gpu/double.dat"

#set logscale x
#set logscale y
set style fill solid noborder
set tics scale 2

set format y "$10^{%T}$"

#set xtics ("1" 1,"2" 2,"3" 1e-4,"$10^{-3}$" 1e-3,"$10^{-2}$" 1e-2,"$10^{-1}$" 1e-1,"$10^{-0}$" 1e-0)
#set ytics ("$10^{0}$" 1, "" 2 1,"" 3 1,"" 4 1,"" 5 1,"" 6 1,"" 7 1,"" 8 1,"" 9 1, "$10^{1}$" 10, "" 20 1, "" 30 1, "" 40 1, "" 50 1, "" 60 1, "" 70 1, "" 80 1, "" 90 1,"$10^{2}$" 100,"" 200 1,"" 300 1,"" 400 1,"" 500 1,"" 600 1,"" 700 1,"" 800 1,"" 900 1,"$10^{3}$" 1000)
set xlabel "Sign iteration"
set ylabel "$||X_k^2-I||_F$"
set key font ",12"

set border 31 linewidth 2

#set key tmargin 2.0
#set key at 0.7,750
set key spacing 2.0
set key right bottom
set logscale y

EREF=0
conv=1.0
plot [2.5:15.5][1e-13:1] \
f1 u ($4):(conv*($6-EREF)) w lp lc rgb col1 lw lw pt 6 t "GPU FP16", \
f2 u ($4):(conv*($6-EREF)) w lp lc rgb col2 lw lw pt 4 t "GPU FP16'", \
f3 u ($4):(conv*($6-EREF)) w lp lc rgb col3 lw lw pt 2 t "GPU FP32", \
f4 u ($4):(conv*($6-EREF)) w lp lc rgb col4 lw lw pt 8 t "GPU FP64", \
f3f u ($4):(conv*($6-EREF)) w lp lc rgb col5 lw lw pt 2 t "FPGA FP32"



set output "bla"
!latex plot.tex
!dvips plot.dvi
!mv plot.ps energy2.ps
!rm energy2.eps
!ps2eps energy2.ps
!rm *.dvi *.aux 
