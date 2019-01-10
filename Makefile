mex:
	cd src/; \
		matlab -r "build_all; exit" -nodesktop -nosplash; \
		

test:
	cd tests/; \
		matlab -r "test_explogit; test_lcexplogit; exit" -nodesktop -nojvm -nosplash

clean:
	rm -f 	build/*