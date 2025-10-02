lim = 1000
retardo = 2
f = 0.01

set terminal qt size 800,800 enhanced font 'Verdana,10'
set size square

do for [i=0:lim]{
  titulo = sprintf("paso = %.4d - tiempo = %.5f",i,i*5e-5)  
  set title titulo
  file = sprintf("output/state_%.4d.txt", i)
  plot file every ::1374::7301 u 2:3 w p ps 1 pt 7 lc rgb "blue" not,\
       "" every ::0::1373 u 2:3 w p ps 1 pt 7 lc rgb "black" not
  pause retardo
}

do for [i=0:lim]{
  titulo = sprintf("paso = %.4d - tiempo = %.5f",i,i*5e-5)  
  set title titulo
  file = sprintf("output/state_%.4d.txt", i)
  plot file every ::1374::7301 u 2:3:($4*f):($5*f) w vectors lc rgb "blue" not,\
       "" every ::0::1373 u 2:3 w p ps 1 pt 7 lc rgb "black" not
  pause retardo
}
