# sodMPI

An older code made in Fortran90 that serves to simulate the evolution of the ground state of the sodium atom in a strong laser field.
The calculation is quite robust so the code was designed to be used on clusters with many active threads.
MPI and openMP were used to optimize for multithreading.

The runtime of the code on 256 threads was about 30 hours per session.
