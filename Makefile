fnpeaks: fnpeaks.c
	gcc -O3 -o $@ $< -lm -funroll-all-loops -ffast-math -msse
