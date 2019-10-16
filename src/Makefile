
# for makhmal
MPICC = /usr/local/bin/mpic++ -Wall

CFLAGS = -O3
CC = g++
INC= -I../include
LIBS =  -L../lib -lbabak -lniftiio -lznz -lm -lz -lc -llapack -lblas
LESSLIBS =  -L../lib -lbabak -lniftiio -lznz -lm -lz -lc

all: libbabak.a brainwash 3dwarper acpcdetect applywarp3d scaleImage avgImage one_sample_t ivf unwarp2d cphdr reslice kaiba oi PILtransform atra aa art2fsl fsl2art yuki analyze2nii multiplyTransformation niftiinfo break4D gambling_vancova_setup medianImage threshold ccstats applymask inverseMatrix4 LMdetect

# because of functions in matrixops, linking the libbabak will possibly require lapack libs.

leastSquaresAffineTrans.o: leastSquaresAffineTrans.cpp
	$(CC) -O2 -c $(INC) leastSquaresAffineTrans.cpp

hist1D_plot.o: hist1D_plot.cpp
	$(CC) -O2 -c $(INC) hist1D_plot.cpp

set_dim.o: set_dim.cpp
	$(CC) -O2 -c $(INC) set_dim.cpp

findMSP.o: findMSP.cpp
	$(CC) -O2 -c $(INC) findMSP.cpp

singular_value_decomposition.o: singular_value_decomposition.cpp
	$(CC) -O2 -c $(INC) singular_value_decomposition.cpp

DKI.o: DKI.cpp
	$(CC) -O2 -c $(INC) DKI.cpp

artlib.o: artlib.cpp
	$(CC) -O2 -c $(INC) artlib.cpp

volume.o: volume.cpp
	$(CC) -O2 -c $(INC) volume.cpp

binomial.o: binomial.cpp
	$(CC) -O2 -c $(INC) binomial.cpp

smooth.o: smooth.cpp
	$(CC) -O2 -c $(INC) smooth.cpp

niftiimage.o: niftiimage.cpp
	$(CC) -O2 -c $(INC) niftiimage.cpp

utils.o: utils.cpp
	$(CC) -O2 -c $(INC) utils.cpp

nifti.o: nifti.cpp
	$(CC) -O2 -c $(INC) nifti.cpp

errorMessage.o: errorMessage.cpp
	$(CC) -O2 -c $(INC) errorMessage.cpp

maskOps.o: maskOps.cpp
	$(CC) -O2 -c $(INC) maskOps.cpp

EMFIT.o: EMFIT.cpp
	$(CC) -O2 -c $(INC) EMFIT.cpp

max_cc.o: max_cc.cpp
	$(CC) -O2 -c $(INC) max_cc.cpp

statistics.o: statistics.cpp
	$(CC) -O2 -c $(INC) statistics.cpp

ginverse.o: ginverse.cpp
	$(CC) -O2 -c $(INC) ginverse.cpp

permutation.o: permutation.cpp
	$(CC) -O2 -c $(INC) permutation.cpp

hpsort.o: hpsort.cpp
	$(CC) -O2 -c $(INC) hpsort.cpp

random.o: random.cpp
	$(CC) -O2 -c $(INC) random.cpp

dicomIO.o: dicomIO.cpp
	$(CC) -O2 -c $(INC) dicomIO.cpp

nkiIO.o: nkiIO.cpp
	$(CC) -O2 -c $(INC) nkiIO.cpp

matrixops.o: matrixops.cpp
	$(CC) -O2 -c $(INC) matrixops.cpp

subsets.o: subsets.cpp 
	$(CC) -O2 -c $(INC) subsets.cpp 

analyzeio.o: analyzeio.cpp 
	$(CC) -O2 -c $(INC) analyzeio.cpp 

getoption.o: getoption.cpp 
	$(CC) -O2 -c $(INC) getoption.cpp 

swap.o: swap.cpp 
	$(CC) -O2 -c $(INC) swap.cpp 

fileinfo.o: fileinfo.cpp
	$(CC) -O2 -c $(INC) fileinfo.cpp

histogram.o: histogram.cpp
	$(CC) -O2 -c $(INC) histogram.cpp

resize.o: resize.cpp
	$(CC) -O2 -c $(INC) resize.cpp

gaussian_kernel.o: gaussian_kernel.cpp
	$(CC) -O2 -c $(INC) gaussian_kernel.cpp

convolution.o: convolution.cpp
	$(CC) -O2 -c $(INC) convolution.cpp

resliceImage.o: resliceImage.cpp
	$(CC) -O2 -c $(INC) resliceImage.cpp

matrixCom.o: matrixCom.cpp
	$(CC) -O2 -c $(INC) matrixCom.cpp

registration.o: registration.cpp
	$(CC) -O2 -c $(INC) registration.cpp

legendre.o: legendre.cpp
	$(CC) -O2 -c $(INC) legendre.cpp

cubicspline.o: cubicspline.cpp
	$(CC) -O2 -c $(INC) cubicspline.cpp

medianfilter.o: medianfilter.cpp
	$(CC) -O2 -c $(INC) medianfilter.cpp

directionCode.o: directionCode.cpp
	$(CC) -O2 -c $(INC) directionCode.cpp

checkNiftiFileExtension.o: checkNiftiFileExtension.cpp
	$(CC) -O2 -c $(INC) checkNiftiFileExtension.cpp

getNiftiImageOrientation.o: getNiftiImageOrientation.cpp
	$(CC) -O2 -c $(INC) getNiftiImageOrientation.cpp

isOrientationCodeValid.o: isOrientationCodeValid.cpp
	$(CC) -O2 -c $(INC) isOrientationCodeValid.cpp

PILtransform.o: PILtransform.cpp
	$(CC) -O2 -c $(INC) PILtransform.cpp

reorientVolume.o: reorientVolume.cpp
	$(CC) -O2 -c $(INC) reorientVolume.cpp

rotate.o: rotate.cpp
	$(CC) -O2 -c $(INC) rotate.cpp

setLowHigh.o: setLowHigh.cpp
	$(CC) -O2 -c $(INC) setLowHigh.cpp

compute_cm.o: compute_cm.cpp
	$(CC) -O2 -c $(INC) compute_cm.cpp

standardize.o: standardize.cpp
	$(CC) -O2 -c $(INC) standardize.cpp

sph.o: sph.cpp
	$(CC) -O2 -c $(INC) sph.cpp

landmarks.o: landmarks.cpp
	$(CC) -O2 -c $(INC) landmarks.cpp

getARTHOME.o: getARTHOME.cpp
	$(CC) -O2 -c $(INC) getARTHOME.cpp

update_qsform.o: update_qsform.cpp
	$(CC) -O2 -c $(INC) update_qsform.cpp

libbabak.a: leastSquaresAffineTrans.o hist1D_plot.o set_dim.o findMSP.o singular_value_decomposition.o DKI.o artlib.o volume.o binomial.o smooth.o niftiimage.o utils.o nifti.o errorMessage.o maskOps.o EMFIT.o max_cc.o statistics.o ginverse.o permutation.o hpsort.o random.o dicomIO.o nkiIO.o subsets.o analyzeio.o getoption.o swap.o fileinfo.o histogram.o resize.o gaussian_kernel.o convolution.o resliceImage.o matrixCom.o registration.o legendre.o cubicspline.o medianfilter.o directionCode.o checkNiftiFileExtension.o getNiftiImageOrientation.o isOrientationCodeValid.o PILtransform.o reorientVolume.o rotate.o setLowHigh.o compute_cm.o standardize.o sph.o matrixops.o landmarks.o getARTHOME.o update_qsform.o
	ar -ru ../lib/libbabak.a leastSquaresAffineTrans.o hist1D_plot.o set_dim.o findMSP.o singular_value_decomposition.o DKI.o artlib.o volume.o binomial.o smooth.o niftiimage.o utils.o nifti.o errorMessage.o maskOps.o EMFIT.o max_cc.o statistics.o ginverse.o permutation.o hpsort.o random.o dicomIO.o nkiIO.o subsets.o analyzeio.o getoption.o swap.o fileinfo.o histogram.o resize.o gaussian_kernel.o convolution.o resliceImage.o matrixCom.o registration.o legendre.o cubicspline.o medianfilter.o directionCode.o checkNiftiFileExtension.o getNiftiImageOrientation.o isOrientationCodeValid.o PILtransform.o reorientVolume.o rotate.o setLowHigh.o compute_cm.o standardize.o sph.o matrixops.o landmarks.o getARTHOME.o update_qsform.o

clean: 
	rm -f *.o

# compilation of executables
PILtransform: PILtransform.cxx
	$(CC) $(CFLAGS) -o PILtransform PILtransform.cxx $(INC) $(LIBS)
	mv PILtransform ${ARTHOME}/bin

brainwash: brainwash.cxx
	$(CC) $(CFLAGS) -o brainwash brainwash.cxx $(LIBS) $(INC)
	mv brainwash ${ARTHOME}/bin

3dwarper: 3dwarper.cxx
	$(CC) $(CFLAGS) -o 3dwarper 3dwarper.cxx $(LIBS) $(INC)
	mv 3dwarper ${ARTHOME}/bin

acpcdetect: acpcdetect.cxx
	$(CC) $(CFLAGS) -o acpcdetect acpcdetect.cxx $(LIBS) $(INC)
	mv acpcdetect ${ARTHOME}/bin

applywarp3d: applywarp3d.cxx
	$(CC) $(CFLAGS) -o applywarp3d applywarp3d.cxx $(INC) $(LESSLIBS)
	mv applywarp3d ${ARTHOME}/bin

scaleImage: scaleImage.cxx
	$(CC) $(CFLAGS) -o scaleImage scaleImage.cxx $(INC) $(LIBS)
	mv scaleImage ${ARTHOME}/bin

avgImage: avgImage.cxx
	$(CC) $(CFLAGS) -o avgImage avgImage.cxx $(INC) $(LESSLIBS)
	mv avgImage ${ARTHOME}/bin

one_sample_t: one_sample_t.cxx
	$(CC) $(CFLAGS) -o one_sample_t one_sample_t.cxx $(INC) $(LIBS)
	mv one_sample_t ${ARTHOME}/bin

ivf: ivf.cxx
	$(CC) $(CFLAGS) -o ivf ivf.cxx $(INC) $(LIBS)
	mv ivf ${ARTHOME}/bin

unwarp2d: unwarp2d.cxx
	$(CC) $(CFLAGS) -o unwarp2d unwarp2d.cxx $(INC) $(LIBS)
	mv unwarp2d ${ARTHOME}/bin

cphdr: cphdr.cxx
	$(CC) $(CFLAGS) -o cphdr cphdr.cxx $(INC) $(LESSLIBS)
	mv cphdr ${ARTHOME}/bin

reslice: reslice.cxx
	$(CC) $(CFLAGS) -o reslice reslice.cxx $(INC) $(LIBS)
	mv reslice ${ARTHOME}/bin

applymask: applymask.cxx
	$(CC) $(CFLAGS) -o applymask applymask.cxx $(INC) $(LIBS)
	mv applymask ${ARTHOME}/bin

inverseMatrix4: inverseMatrix4.cxx
	$(CC) $(CFLAGS) -o inverseMatrix4 inverseMatrix4.cxx $(INC) $(LIBS)
	mv inverseMatrix4 ${ARTHOME}/bin

LMdetect: LMdetect.cxx
	$(CC) $(CFLAGS) -o LMdetect LMdetect.cxx $(INC) $(LIBS)
	mv LMdetect ${ARTHOME}/bin

kaiba: kaiba.cxx
	$(CC) $(CFLAGS) -o kaiba kaiba.cxx $(INC) $(LIBS)
	mv kaiba ${ARTHOME}/bin
#Example: $(CC) $(CFLAGS) -o kaiba kaiba.cxx -DBETA_PARAM=0.3 $(INC) $(LIBS)
# If we used the option -DBETA_PARAM=0.3 as in the above example, it will be the same as
# having the following statement in kaiba.cxx
#define BETA_PARAM 0.3

oi: oi.cxx
	$(CC) $(CFLAGS) -o oi oi.cxx $(INC) $(LIBS)
	mv oi ${ARTHOME}/bin

atra: atra.cxx
	$(CC) $(CFLAGS) -o atra atra.cxx $(INC) $(LIBS)
	mv atra ${ARTHOME}/bin

aa: aa.cxx
	$(CC) $(CFLAGS) -o aa aa.cxx $(INC) $(LIBS)
	mv aa ${ARTHOME}/bin

art2fsl: art2fsl.cxx
	$(CC) $(CFLAGS) -o art2fsl art2fsl.cxx $(INC) $(LIBS)
	mv art2fsl ${ARTHOME}/bin

fsl2art: fsl2art.cxx
	$(CC) $(CFLAGS) -o fsl2art fsl2art.cxx $(INC) $(LIBS)
	mv fsl2art ${ARTHOME}/bin

yuki: yuki.cxx
#	$(MPICC) $(CFLAGS) -o yuki yuki.cxx $(INC) $(LIBS)
	$(MPICC) $(CFLAGS) -o yuki yuki.cxx $(INC) -L${ARTHOME}/lib -lbabak -lniftiio -lznz -lm -lz -lc -llapack -lblas
	mv yuki ${ARTHOME}/bin

analyze2nii: analyze2nii.cxx
	$(CC) $(CFLAGS) -o analyze2nii analyze2nii.cxx $(INC) $(LIBS)
	mv analyze2nii ${ARTHOME}/bin

multiplyTransformation: multiplyTransformation.cxx
	$(CC) $(CFLAGS) -o multiplyTransformation multiplyTransformation.cxx $(INC) $(LESSLIBS)
	mv multiplyTransformation ${ARTHOME}/bin

niftiinfo: niftiinfo.cxx
	$(CC) $(CFLAGS) -o niftiinfo niftiinfo.cxx $(INC) $(LESSLIBS)
	mv niftiinfo ${ARTHOME}/bin

break4D: break4D.cxx
	$(CC) $(CFLAGS) -o break4D break4D.cxx $(INC) $(LIBS)
	mv break4D ${ARTHOME}/bin

vancova: vancova.cxx
	$(CC) $(CFLAGS) -o vancova vancova.cxx $(INC) $(LIBS) -ldcdf  -lf2c
	mv vancova ${ARTHOME}/bin

gambling_vancova_setup: gambling_vancova_setup.cxx
	$(CC) $(CFLAGS) -o gambling_vancova_setup gambling_vancova_setup.cxx $(INC) $(LIBS)
	mv gambling_vancova_setup ${ARTHOME}/bin

medianImage: medianImage.cxx
	$(CC) $(CFLAGS) -o medianImage medianImage.cxx $(INC) $(LIBS)
	mv medianImage ${ARTHOME}/bin

threshold: threshold.cxx
	$(CC) $(CFLAGS) -o threshold threshold.cxx $(INC) $(LIBS)
	mv threshold ${ARTHOME}/bin

ccstats: ccstats.cxx
	$(CC) $(CFLAGS) -o ccstats ccstats.cxx $(INC) $(LIBS)
	mv ccstats ${ARTHOME}/bin
