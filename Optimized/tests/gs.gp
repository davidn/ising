#!/usr/bin/env gnuplot

set term pdf
set output 'gridsize.pdf'
set xlabel '1/number of points'
set ylabel 'critical temperature J/k'
set ytics rotate by 45 add ('Tc(onsager)' 2.269)
set xrange [0:]
fit (y+m*x) 'gridsize' u (1/$1**2):4:5 via y,m
plot '' u (1/$1**2):4:5 notitle w yerrorbars, (y+m*x) notitle
