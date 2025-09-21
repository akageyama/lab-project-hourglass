#!/bin/bash

infile="case02.txt"  # Do not remove this file!


#####################################
#  次の二行と 以下の2点を変更する
#     (1) gnuplot script の yrange
#     (2) gnuplot のplot文
#####################################

# N=1  # 平均をとる行数  for case02-N1-rawdata.pdf
# pdffile="case02-N1-rawdata.pdf"

N=48   # 平均をとる行数 for case02-N48-averaged.pdf
pdffile="case02-N48-averaged.pdf"


averageddata=$(mktemp)
workfile=$(mktemp)
plotscript=$(mktemp)


#---------------------------
generate_plotscript()
#---------------------------
{
cat<<EOF
#
set xrange [-1:13]


#####################################
## For case02-N48-averaged.pdf
#####################################
set yrange [0.995:1.005]

#####################################
## For case02-N1-rawdata.pdf
#####################################
#set yrange [0.9:1.1]

# 軸ラベル
set xlabel "t (sec)" font ",16"
set ylabel "weight (kg)" font ",16"

# 軸の数字のフォントを大きく
set tics font ",14"

####################################
# For case02-N48-averaged.pdf
# 1.0の線と理論値とデータをプロット
####################################
plot 1.0 notitle lw 0.5, \
     1.00204 t "theory" lw 0.5, \
     '$averageddata' with points pointtype 7 pointsize 0.2 notitle

#####################################
## For case02-N1-rawdata.pdf
## 1.0の線とデータをプロット
#####################################
#plot 1.0 notitle lw 0.5, \
#     '$averageddata' with points pointtype 7 pointsize 0.2 notitle

pause -1

# pdf出力

set terminal pdfcairo font ",16" linewidth 2 rounded size 8cm,6cm
set output "$pdffile"

replot

unset output
EOF
}


#---------------------------
take_average_in_N_lines()
#---------------------------
{
  awk -v N=$N '
  {
      t_val = $1        # t の値を逐次更新（最後に読んだものが残る）
      y_sum += $2
      count++

      if (count == N) {
          printf "%.8f %.8f\n", t_val, y_sum/N
          y_sum=0; count=0
      }
  }
  END {
      # 余りの行があった場合は平均を出す（このとき t は最後の行の値）
      if (count > 0) {
          printf "%.8f %.8f\n", t_val, y_sum/count
      }
  }' "$workfile"
}




cat $infile | awk '{print $3, $4}' > $workfile

echo "$N 行毎に平均値をとる"
take_average_in_N_lines > $averageddata


echo "gnuplot script作成"
generate_plotscript > $plotscript


echo "plot実行"
gnuplot $plotscript

rm $plotscript $averageddata $workfile
