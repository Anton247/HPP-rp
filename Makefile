videos: run
	make -B density/density.avi >/dev/null
	make -B rest/rest.avi >/dev/null
	make -B move/move.avi >/dev/null 

init:
	make rules.txt
	make config.c

run: hpprp FORCE
	./hpprp

hpprp: hpprp.c config.c collide.c
	g++ -fopenmp -O3 hpprp.c -o hpprp -lrt -Wall 

collide.c: rules.txt system/gen_collide.py
	python3 system/gen_collide.py rules.txt collide.c

rules.txt:
	if [ -f rules.txt ] ; \
		then echo "rules.txt exists, not overriding" ; \
		touch rules.txt ; \
		else echo "rules.txt does not exist, create" ; \
		python3 system/gen_rules.py \
			5-5=0.3333 10-10=0.3333 16-16=0.3333\
			15-15=0.3333 21-21=0.3333 26-26=0.3333 > rules.txt ; \
	fi

config.c:
	if [ -f config.c ] ; \
		then echo "config.c exists, not overriding" ; \
		touch config.c ; \
	else echo "config.c does not exist, create" ; \
		WIDTH=500 HEIGHT=500 ITERS=501 \
		SAVE_PERIOD=50 AVERAGING_RADIUS=1 ENSEMBLE_SIZE=2 \
		INIT1=0.7 INIT2=0.7 INIT4=0.7 INIT8=0.7 \
		INIT16=0.25 INIT32=0 INIT64=0 INIT128=0 \
		SOURCE_WIDTH=20 \
		SRC1=1 SRC2=1 SRC4=1 SRC8=1 \
		SRC16=0.75 SRC32=0 SRC64=0 SRC128=0 \
		python3 system/gen_config_flatwave.py > config.c ; \
	fi

density/density.avi:
	system/image_gen.sh density density/density.avi

rest/rest.avi:
	system/image_gen.sh rest rest/rest.avi

move/move.avi:
	system/image_gen.sh move move/move.avi

clean:
	rm -f density/* rest/* move/* hpprp collide.c

purge:
	make clean
	rm -f config.c rules.txt 

FORCE:
