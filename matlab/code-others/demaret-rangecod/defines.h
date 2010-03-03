#ifndef _DEFINES_H_
#define _DEFINES_H_

#define TRUE  1
#define FALSE 0

#define MAX(a,b)    ( (a) > (b) ? (a) : (b) )
#define MIN(a,b)    ( (b) > (a) ? (a) : (b) )

#define SQUARE(x)   ( (x)*(x) )

#define SIGN(x)     ( ((x) < 0) ? -1 : ((x) > 0) ? 1 : 0 )
#define ABS(x)      ( (x) < 0 ? -(x) : (x) )

#define ODD(x)      ( (x) & 1 == 1 )
#define EVEN(x)     ( (x) & 1 == 0 )

#define MOD2(x)     ( (x) & 0x01 )
#define MOD4(x)     ( (x) & 0x03 )

/* take care with swap: a and b shouldn't be the same variable! */
#define SWAP(a,b)   { (a)=(a)^(b); (b)=(a)^(b); (a)=(a)^(b); }

#define LEN2D(p)    ( SQUARE((p).x) + SQUARE((p).y) )
#define LEN3D(p)    ( SQUARE((p).x) + SQUARE((p).y) + SQUARE((p).z) )

#define DIST2D(p,q) ( SQUARE((p).x-(q).x) + SQUARE((p).y-(q).y) )
#define DIST3D(p,q) ( SQUARE((p).x-(q).x) + SQUARE((p).y-(q).y) + SQUARE((p).z-(q).z) )

#define ADD2D(p,q,r)  { (r).x = (p).x + (q).x; (r).y = (p).y + (q).y; }
#define DIFF2D(p,q,r) { (r).x = (p)->x - (q)->x; (r).y = (p)->y - (q)->y; }

#define DOT2D(p,q)  ( (p).x*(q).x + (p).y*(q).y )

#endif /* of _DEFINES_H_ */
