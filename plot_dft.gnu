reset

set terminal wxt enhanced size 1100,700
set style line 1 lc rgb "orange" lw 2
set style line 2 lc rgb "blue"   lw 2

set grid xtics ytics mxtics mytics
set xlabel "x"
set ylabel "f(x)"

Ns = "10 50 100 1000 10000"

do for [n in Ns] {
    filename = sprintf("dft_graphics_%s.dat", n)
    title = sprintf("Сравнение DFT для N = %s", n)

    set title title

    plot filename using 1:2 with lines ls 1 title "analytic", \
         filename using 1:3 with lines ls 2 title "DFT"

    pause -1   # ← ДЕЛАЕТ ОКНО ИНТЕРАКТИВНЫМ И ЗУМИРУЕМЫМ
}
