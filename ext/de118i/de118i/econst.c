/*							econst.c	*/
/*  e type constants used by high precision check routines */

#include "ehead.h"
/* 0.0 */
unsigned short ezero[NE] = {
0, 0000000,0000000,0000000,0000000,0000000,};
extern unsigned short ezero[];
/* 5.0E-1 */
unsigned short ehalf[NE] = {
0, 0000000,0000000,0000000,0100000,0x3ffe,};
extern unsigned short ehalf[];
/* 1.0E0 */
unsigned short eone[NE] = {
0, 0000000,0000000,0000000,0100000,0x3fff,};
extern unsigned short eone[];
/* 2.0E0 */
unsigned short etwo[NE] = {
0, 0000000,0000000,0000000,0100000,0040000,};
extern unsigned short etwo[];
/* 3.2E1 */
unsigned short e32[NE] = {
0, 0000000,0000000,0000000,0100000,0040004,};
extern unsigned short e32[];
/* 6.93147180559945309417232121458176568075500134360255E-1 */
unsigned short elog2[NE] = {
0xc9e4,0x79ab,0150717,0013767,0130562,0x3ffe,};
extern unsigned short elog2[];
/* 1.41421356237309504880168872420969807856967187537695E0 */
unsigned short esqrt2[NE] = {
0x597e,0x6484,0174736,0171463,0132404,0x3fff,};
extern unsigned short esqrt2[];
/* 2/sqrt(PI) =
 * 1.12837916709551257389615890312154517168810125865800E0 */
unsigned short eoneopi[NE] = {
0x71d5,0x688d,0012333,0135202,0110156,0x3fff,};
extern unsigned short eoneopi[];
/* 3.14159265358979323846264338327950288419716939937511E0 */
unsigned short epi[NE] = {
0xc4c6,0xc234,0020550,0155242,0144417,0040000,};
extern unsigned short epi[];
/* 5.7721566490153286060651209008240243104215933593992E-1 */
unsigned short eeul[NE] = {
0xd1be,0xc7a4,0076660,0063743,0111704,0x3ffe,};
extern unsigned short eeul[];
