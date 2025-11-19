#ifndef RANDOM_H

#define RANDOM_h

void innit_genrand64(unsigned long long seed);

void init_by_array64(unsigned long long*, unsigned long long);

unsigned long long genrand64_int64(void);

long long genrand64_int63(void);
/* generates a random number on [0,1]-real-interval */
double genrand64_real1(void);
/* generates a random number on [0,1)-real-interval */
double genrand64_real2(void);
/* generates a random number on (0,1)-real-interval */
double genrand64_real3(void);

#endif
