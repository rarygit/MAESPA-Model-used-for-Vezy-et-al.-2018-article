


setwd("D:/These/Projects/MAESPA_Validation/Test/make/")
intel= system.time(system2("MAESPAV2.6.exe", input=as.character(1), wait=T))
Intelhrflx= readhrflux("1_hrflux.dat")
gfortran= system.time(system2("maespa.exe", input=as.character(1), wait=T))
GFortranhrflx= readhrflux("1_hrflux.dat")


plot(Intelhrflx$hrPAR[Intelhrflx$Tree==unique(Intelhrflx$Tree)[1]])
lines(GFortranhrflx$hrPAR[GFortranhrflx$Tree==132])
plot(Intelhrflx$TCAN[Intelhrflx$Tree==132])
lines(GFortranhrflx$TCAN[GFortranhrflx$Tree==132])
# Lines have to be exactly overlapping points.



gfortran
intel
# gfortran's compiled MAESPA takes 2x times intel's compiled MAESPA

