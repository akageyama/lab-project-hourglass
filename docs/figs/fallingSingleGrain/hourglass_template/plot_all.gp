
set yrange [0:1.2]
set xrange [0:0.7]

# 軸ラベル
set xlabel "t (sec)" font ",16"
set ylabel "y (m)" font ",16"

# 軸の数字のフォントを大きく
set tics font ",14"

# データを凡例なしで描画
plot 'pos_y.txt' using ($3):($4) w l notitle

pause -1

# pdf出力

set terminal pdfcairo font ",16" linewidth 2 rounded size 8cm,6cm
set output "plot_all.pdf"

replot

unset output
