#
#



set xrange [0.449:0.450]

# 軸ラベル
set xlabel "t (sec)" font ",16"
set ylabel "y (m)" font ",16"

# 軸の数字のフォントを大きく
set tics font ",14"

# データを凡例なしで描画
# plot 'pos_y.txt' using ($3):($4) w lp notitle
plot 'pos_y.txt' using ($3):($4) with points pointtype 7 pointsize 0.2 notitle

pause -1

# pdf出力

set terminal pdfcairo font ",16" linewidth 2 rounded size 8cm,6cm
set output "plot_closeup.pdf"

replot

unset output
