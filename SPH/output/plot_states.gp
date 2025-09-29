# ==========================
# Archivo: plot_states.gp
# ==========================

# Número máximo de pasos (ajústalo si es necesario)
lim = 100
f = 0.01

set terminal pngcairo size 800,800 enhanced font 'Verdana,10'
set size square

# ==========================
# Calcular rangos automáticamente
# ==========================
xmin = 1e30
xmax = -1e30
ymin = 1e30
ymax = -1e30

do for [i=0:lim] {
    file = sprintf("state_%.4d.txt", i)
    if (system(sprintf("test -f %s", file)) == 0) {
        stats file using 2:3 nooutput
        if (STATS_min_x < xmin) xmin = STATS_min_x
        if (STATS_max_x > xmax) xmax = STATS_max_x
        if (STATS_min_y < ymin) ymin = STATS_min_y
        if (STATS_max_y > ymax) ymax = STATS_max_y
    }
}

set xrange [xmin:xmax]
set yrange [ymin:ymax]

# ==========================
# Fase 1: posiciones
# ==========================
do for [i=0:lim] {
    file = sprintf("state_%.4d.txt", i)
    if (system(sprintf("test -f %s", file)) == 0) {
        set title sprintf("Paso = %.4d", i)
        set output sprintf("frame1_%.4d.png", i)
        plot file u 2:3 w p ps 1 pt 7 lc rgb "blue" notitle
    }
}

# ==========================
# Fase 2: posiciones + vectores de velocidad
# ==========================
do for [i=0:lim] {
    file = sprintf("state_%.4d.txt", i)
    if (system(sprintf("test -f %s", file)) == 0) {
        set title sprintf("Paso = %.4d", i)
        set output sprintf("frame2_%.4d.png", i)
        plot file u 2:3 w p ps 1 pt 7 lc rgb "blue" notitle, \
             file u 2:3:($4*f):($5*f) w vectors lc rgb "red" notitle
    }
}
