REF='SiO2Spheres_Kruger.SIFlux'
DATA00='/home/homer/work/scuff-em/examples/SiO2Spheres/SiO2Sphere_501.SIFlux'
DATA01='/home/homer/work/scuff-em/examples/SiO2Spheres/SiO2Spheres_501.SIFlux'
DATA02='/home/homer/work/scuff-em/examples/SiO2Spheres/SiO2Sphere_1479.SIFlux'

DATA10='SiO2Sphere_677.SIFlux.JDEPFT'
DATA11='SiO2Spheres_677.SIFlux.JDEPFT'

DATA20='Run2/SiO2Sphere_3298.SIFlux.JDEPFT' 
DATA21='Run2/SiO2Spheres_3298.SIFlux.JDEPFT'

ifeq(x,y,z)=(x==y ? z : 1/0)

set xlabel 'Omega (3e14 rad/sec)'

set logscale y
set format y "%.0e"
set format x "%.1f"

##################################################
##################################################
##################################################
unset logscale x
set yrange [1e-6:]
set xrange [0.01:1]
set terminal x11 1
set title 'Power radiation'
set ylabel 'Radiated Power flux (dimensionless)'
set key at 0.85,1e-3
plot REF    u 2:(abs($3)) t 'Krueger'      w l lw 2            \
    ,DATA00 u 2:(2*$4) t 'SCUFF, N=501'  w p pt 7 ps 1.25          \
    ,DATA10 u 2:(2*abs($4)) t 'BUFF,  N=677' w p pt 7 ps 1.25    \
    ,DATA20 u 2:(2*abs($4)) t 'BUFF,  N=3298' w p pt 7 ps 1.25
#call 'png' 'SiO2Sphere_PowerRadiation.png'

##################################################
##################################################
##################################################
unset logscale x
unset yrange
set autoscale y
set xrange [0.1:1]
set terminal x11 2
set title 'Power transfer'
set ylabel 'Transferred Power flux (dimensionless)'
set key at 0.85,1e-8
plot REF    u 2:4                t 'Krueger'      w l lw 2           \
    ,DATA01  u (ifeq($3,12,$2)):4 t 'SCUFF, N=501' w p pt 7 ps 1.25  \
    ,DATA02  u (ifeq($3,12,$2)):4 t 'SCUFF, N=1479' w p pt 6 ps 1.25 \
    ,DATA11 u (ifeq($3,12,$2)):(ifeq($1,10,abs($4))) t 'BUFF,  N=677' w p pt 9 ps 1.00   \
    ,DATA21 u (ifeq($3,12,$2)):(ifeq($1,10,abs($4))) t 'BUFF,  N=3298' w p pt 8 ps 1.25
#call 'png' 'SiO2Spheres_PowerTransfer.png'

##################################################
##################################################
##################################################
set terminal x11 3
set title 'Force flux (1-->2)'
set key at 0.85,1e-6
set ylabel 'Force density (nanonewtons/watts)'
plot REF   u 2:(abs($5))                t 'Krueger'      w l lw 2           \
    ,DATA01 u (ifeq($3,12,$2)):(abs($5)) t 'SCUFF, N=501' w p pt 7 ps 1.00  \
    ,DATA02 u (ifeq($3,12,$2)):(abs($5)) t 'SCUFF, N=1479' w p pt 6 ps 1.25 \
    ,DATA11 u (ifeq($3,12,$2)):(ifeq($1,10,abs($8))) t 'BUFF, N=677'  w p pt 9 ps 1.00  \
    ,DATA21 u (ifeq($3,12,$2)):(ifeq($1,10,abs($8))) t 'BUFF, N=3298' w p pt 8 ps 1.25 
#call 'png' 'SiO2Spheres_F12.png'

##################################################
##################################################
##################################################
set terminal x11 4
set key at 0.85,1e-6
set title 'Force flux 2-->2'
set ylabel 'Force density (nanonewtons/watts)'
plot REF    u 2:(+$6)                    t 'Krueger (positive)' w l lw 2     \
    ,REF    u 2:(-$6)                    t 'Krueger (negative)' w l lw 2     \
    ,DATA01  u (ifeq($3,22,$2)):(abs($5)) t 'SCUFF, N=501' w p pt 7 ps 1.00  \
    ,DATA02  u (ifeq($3,22,$2)):(abs($5)) t 'SCUFF, N=1479' w p pt 6 ps 1.25 \
    ,DATA11 u (ifeq($3,22,$2)):(ifeq($1,10,abs($8))) t 'BUFF, N=677'  w p pt 9 ps 1.00  \
    ,DATA21 u (ifeq($3,22,$2)):(ifeq($1,10,abs($8))) t 'BUFF, N=3298' w p pt 8 ps 1.25 
#call 'png' 'SiO2Spheres_F22.png'
