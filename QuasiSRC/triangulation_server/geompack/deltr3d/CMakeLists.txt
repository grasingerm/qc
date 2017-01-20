cmake_minimum_required(VERSION 3.1)

add_library(deltr3d SHARED dtris3.f dtriw3.f imptr3.f bnsrt3.f walkt3.f vbfac.f 
nwthou.f nwthin.f nwthed.f swapes.f swapmu.f tetmu.f frstet.f tetlst.f 
order3.f htins.f htdel.f htfun.c htsrc.f availf.f updatf.f ifacty.f imptrf.f 
fndmsw.f swaptf.f itris3.f imptrd.f nwthfc.f)
