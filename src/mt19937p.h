/* mt19937p.h
 *
 * This is a version of the Mersenne Twister (aka mt19937) that is designed for
 * parallel use. Usage:
 *
 * struct mt19937p myMT19937;
 * sgenrand(seed,&myMT19937);
 * randomNumber = genrand(&myMT19937);
 *
 * Simply create different mt19937p structures for your different streams.
 * Disclaimer & Warning: This has only been lightly tested and has shown that it
 * generates identical sequences as mt19937.c! Use at your own risk!
 *
 * Massaged into parallel for my Austen McDonald <austen@cc.gatech.edu>
 * http://www.prism.gatech.edu/~gte363v/
 * November, 2000
 */

#ifndef MT19937P
#define MT19937P

#define MATRIX_A 0x9908b0df /* constant vector a */
#define NN 624

struct mt19937p {
	unsigned long mt[NN];
	int mti;
	unsigned long mag01[2];
}; 
typedef struct mt19937p mt19937p;

#ifdef __cplusplus
extern "C" {
#endif

/* initializing the array with a NONZERO seed */
void sgenrand(unsigned long seed, struct mt19937p* config);

double /* generating reals */
/* unsigned long */ /* for integer generation */
genrand(struct mt19937p* config);

unsigned long /* for integer generation */
genrand_int(struct mt19937p* config);

#ifdef __cplusplus
}
#endif
		
#endif
