set pm3d map
set size square

if(!exists("N1")){
    N1 = 10000
}
if(!exists("N2")){
    N2 = 300
}
if(!exists("S")){
    S = 100
}
if(!exists("name")){
    name = "c"
}
if(!exists("outputPath")){
    outputPath = "./analysis/"
}

set palette model RGB
set palette defined ( 0 "black", 1 "white" )
unset key
unset border
unset xtics
unset ytics
unset colorbox

set xrange [0:N]
set yrange [0:N]
set term png font "Helvetica,20" crop

do for [i=50:N2:50]{
    num = sprintf('%09d',i)
    pic_path = outputPath.'/analysis/'.name.'_profile'
    set output pic_path.'/'.name.num.'.png'
    filename = outputPath.'/data/'.name.'-time_'.num.'.dat'
    set title "T = ".num offset 0, -1
    splot filename u 1:2:3 w image palette
}

D = N2 + S - N2 + N2/S * S
do for [i=D:N1:S]{
    num=sprintf('%09d',i)
    pic_path = outputPath.'/analysis/'.name.'_profile'
    set output pic_path.'/'.name.num.'.png'
    filename = outputPath.'/data/'.name.'-time_'.num.'.dat'
    set title "T = ".num offset 0, -1
    splot filename u 1:2:3 w image palette
}
print "Completed printing images"
