########################################
#                                      #
#        Makefile for OneDimDVR        #
#                                      #
########################################

flag = -O3

analyzers:
	cd Che2wfn     ; make flag=$(flag); cd ..
	cd Che2tran    ; make flag=$(flag); cd ..
	cd density     ; make flag=$(flag); cd ..
	cd population  ; make flag=$(flag); cd ..
	cd transmission; make flag=$(flag); cd ..
	cd Wigner      ; make flag=$(flag); cd ..
	cd expectation ; make flag=$(flag); cd ..

.PHONY: install
install:
	rm -rf OneDimDVR
	mkdir  OneDimDVR
	mv Che2wfn/Che2wfn.exe            OneDimDVR
	mv Che2tran/Che2tran.exe          OneDimDVR
	mv density/density.exe            OneDimDVR
	mv population/population.exe            OneDimDVR
	mv transmission/transmission.exe  OneDimDVR
	mv Wigner/Wigner.exe              OneDimDVR
	mv expectation/expectation.exe    OneDimDVR
	ln -s ../density/animate.py       OneDimDVR/animate-density.py
	ln -s ../Wigner/animate.py        OneDimDVR/animate-Wigner.py
	ln -s ../expectation/visualize.py OneDimDVR/visualize-expectation.py
