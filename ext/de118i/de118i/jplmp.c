//------------------------------------------------------------------------------
/* Orbits of the five minor planets as implemented in the
 * JPL DE102 and DE118/DE200 numerical integrations.
 * This program was written with the aid of documentation
 * kindly supplied by E. M. Standish of JPL.
 *
 * Note that in this model the orbits are fixed Keplerian ellipses.
 * Each is specified by a perihelion date and period
 * of revolution, together with a set of Chebyshev
 * expansions for the x, y, and z coordinates covering
 * one period.  The result is in 1950 coordinates, which
 * are then rotated into the J2000 FK5 system.
 */
//------------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>


#ifndef NOINCS

#include "mconf.h"
#include "prec.h"
#include "ssystem.h"
#include "const.h"

#endif

#if LDOUBLE
#if IBMPC
/* Julian date of perihelion
 */
short mpperihdi[] = {
0x410b,0x4a0e,0xbe78,0x94e8,0x4014, XPD
0x7271,0x43bf,0x97f1,0x94e9,0x4014, XPD
0x9082,0x602c,0x1c62,0x94e8,0x4014, XPD
0x90ac,0x93f2,0xd46a,0x94cb,0x4014, XPD
0x5a1d,0x3b64,0x1b1f,0x94c3,0x4014, XPD
};

/* Orbital period, in days
 */
short mpperiodi[] = {
0x7c07,0x187e,0x3c21,0xd212,0x4009, XPD
0xc41e,0x14ce,0xaf3a,0xd255,0x4009, XPD
0x5547,0x4413,0x07dd,0xa5bc,0x4009, XPD
0x443d,0x3ba3,0x9d88,0xa849,0x4009, XPD
0x872b,0xd916,0xf7ce,0xc903,0x4009, XPD
};

/* Chebyshev coefficients for one complete revolution
 */
short chi[] = {
0x2ef3,0x70e8,0xb8e6,0xeada,0xbffd, XPD
0x8620,0x0338,0x48c9,0xcb6b,0x3ffe, XPD
0x7b6b,0x2410,0xf574,0x95a2,0xc000, XPD
0x7684,0x9374,0x50fe,0xe642,0xbffe, XPD
0xf073,0xc30c,0x823d,0xa68d,0x3ffe, XPD
0xc9d5,0xe3fe,0xf69b,0xce0e,0x3ffb, XPD
0xab89,0xcef9,0x1f96,0x97ed,0xbff9, XPD
0xc1a5,0x558e,0x5783,0xb558,0x3ff7, XPD
0x91b7,0xf90b,0xad04,0xdbd4,0xbff7, XPD
0x5ea0,0xef12,0x6321,0xa3bf,0xbff5, XPD
0x30a9,0xe055,0x4d7a,0x87b3,0x3ff3, XPD
0xf94b,0x022b,0x03ca,0x9f56,0xbff1, XPD
0x50e4,0x3cf3,0xf3ed,0xd162,0x3ff1, XPD
0xe899,0x6c6e,0x5e31,0xa98a,0x3fef, XPD
0xfa1d,0x23ec,0x9c7c,0xb2a2,0xbfec, XPD
0x1d2e,0x641a,0x2fd3,0xd453,0x3ffc, XPD
0xd7c0,0xa169,0x104c,0xaa2c,0x3fff, XPD
0xd4e5,0xc6e1,0x322d,0x8748,0x3fff, XPD
0x137a,0x48af,0x16e8,0xc0a0,0xbfff, XPD
0x8dd6,0x7ab6,0x51a9,0x9693,0xbffd, XPD
0x1ae2,0x4178,0x4f29,0xac61,0x3ffc, XPD
0xb606,0x9528,0x2154,0x895a,0x3ff8, XPD
0x1ad3,0xec28,0xb714,0x97b4,0x3ff8, XPD
0x95a0,0x7d95,0x16d9,0xc6be,0x3ff6, XPD
0xd94f,0xc4a7,0x0daf,0x88fc,0xbff6, XPD
0x0fed,0xdd4f,0x90d8,0xf55d,0xbff1, XPD
0xa0e8,0xa08a,0x4435,0x854b,0xbff2, XPD
0x9c20,0x05c8,0xdad6,0xbd4c,0xbff0, XPD
0xd655,0x1ee2,0xa79f,0x8dd4,0x3ff0, XPD
0x13ab,0xc44e,0xb3bc,0xa17f,0x3feb, XPD
0xffc5,0xaff0,0x0443,0xc3dc,0x3ffc, XPD
0x2ee1,0x191e,0x8db2,0xedc3,0x3ffd, XPD
0xfbca,0x5bef,0x0890,0xf995,0x3ffe, XPD
0x4c05,0x8b00,0x5657,0x8691,0xbffe, XPD
0xd402,0x9d7c,0x06fc,0x8ae6,0xbffd, XPD
0x9332,0xe3a5,0x5003,0xf0d9,0x3ffa, XPD
0xb0c5,0x709e,0xb5e2,0xfd66,0x3ff7, XPD
0xfa7a,0x33ee,0x84e5,0xd3f6,0x3ff6, XPD
0x7d23,0x7fc5,0x9030,0xb754,0x3ff6, XPD
0xd544,0x20ba,0xf561,0xbf64,0xbff4, XPD
0xdd09,0x1278,0x75c0,0xe256,0xbff1, XPD
0x8c84,0x654a,0xe7cd,0xba3c,0xbff0, XPD
0xef05,0xd1f3,0xc832,0xae9e,0xbff0, XPD
0x9838,0x76b0,0x527c,0xc62a,0x3fee, XPD
0x6a50,0x9181,0x8e02,0x94f9,0x3feb, XPD
0x8c47,0x7619,0xe834,0xafbb,0x3ffa, XPD
0x4d6a,0xe040,0x3340,0xb118,0x3fff, XPD
0x2c42,0xaea5,0x7368,0xcbf2,0xbfff, XPD
0x26a5,0x107e,0x003e,0xb8f0,0xbfff, XPD
0x78ec,0x227e,0x6b87,0xa7d8,0x3ffd, XPD
0x776c,0x8870,0xce5a,0xc14b,0x3ff9, XPD
0x0ec4,0x64dd,0x82c9,0xb1c2,0x3ffa, XPD
0x5d24,0x1da0,0x578b,0x8cd8,0x3ffa, XPD
0x541b,0x49ff,0x0ad1,0xdb8b,0xbff6, XPD
0xc49c,0xdf10,0xfe16,0x8d68,0x3ff7, XPD
0x1518,0xac23,0x97ee,0x8728,0xbff6, XPD
0xba10,0x5f9b,0x1bea,0xbeec,0xbff4, XPD
0x461d,0xce2a,0x21da,0xa4b0,0xbff2, XPD
0xa48d,0xa0f2,0x32bf,0x9b07,0xbff3, XPD
0xef2e,0x0919,0x1bea,0xc26b,0x3ff0, XPD
0x6970,0xc64d,0x74c0,0x806b,0xbffb, XPD
0x938b,0xe7cf,0x3de2,0xe7bd,0x3ffe, XPD
0xb439,0xdd69,0x6ac6,0x9509,0x4000, XPD
0x10db,0xfa33,0x92f5,0xf200,0xbffe, XPD
0xc3b0,0x23b9,0x6213,0xf54f,0xbffd, XPD
0x2a99,0xc5b8,0xa470,0xfcf0,0x3ff8, XPD
0xb9f6,0xaf91,0x6e9c,0x81e6,0xbffb, XPD
0x9a8b,0xe3e0,0xf81f,0xb84d,0x3ff9, XPD
0xa385,0xdb36,0x0d2d,0xa06f,0x3ff7, XPD
0x6aa0,0x7dd1,0x40e6,0xb90b,0x3ff6, XPD
0x0e5c,0xfb38,0x9a44,0xc589,0x3ff6, XPD
0x20c4,0x1c5d,0x6ed5,0xf9d5,0xbff3, XPD
0xa99d,0xf919,0x0d0d,0xf0b2,0x3ff2, XPD
0xa18a,0x7e4f,0x2bba,0xcadd,0xbff2, XPD
0x055f,0xfbea,0xd6b1,0x8e12,0xbff1, XPD
0xcc99,0xfbeb,0xcf20,0x97cc,0x3ff8, XPD
0x9768,0x3556,0xed1b,0x8f84,0xbffd, XPD
0x3993,0x43ba,0xaeef,0xb02b,0xbffd, XPD
0xc0bc,0xf3b2,0x1196,0x95e0,0x3ffd, XPD
0x0dcc,0x1e08,0x5d61,0x90fc,0x3ffb, XPD
0x6aa4,0x865f,0x305e,0x9ca6,0xbff7, XPD
0x27b1,0x0fe6,0xc833,0x998c,0x3ff8, XPD
0x1e15,0xa40e,0xddad,0xe448,0xbff7, XPD
0x9b53,0xf51a,0x81a6,0xbda4,0xbff4, XPD
0x181c,0xdd8a,0x518d,0xe533,0xbff4, XPD
0x8aca,0xe178,0x60a5,0xe980,0xbff3, XPD
0x8b23,0x944a,0xb479,0x9ab9,0x3ff2, XPD
0x5757,0xee6c,0x2943,0x8e42,0xbff0, XPD
0x16dc,0x5faa,0xd53e,0xfb45,0x3ff0, XPD
0x6400,0x2b18,0x90eb,0xa7f0,0x3fee, XPD
0x555b,0xdcc2,0x4156,0xf43e,0xbffb, XPD
0x3be7,0x40e7,0xd8a6,0xa8d5,0xbfff, XPD
0xe1a1,0x101d,0x0863,0xa7ec,0xbffe, XPD
0x1459,0x02d3,0x061f,0xbe3d,0x3fff, XPD
0x3773,0x59ea,0xaecd,0xb7fc,0x3ffc, XPD
0x95e8,0x2fe3,0xb6c0,0xa17d,0xbffc, XPD
0x5e07,0x5492,0x734b,0xe5ae,0xbff6, XPD
0x5afd,0x3fa2,0x14d1,0xc01a,0xbff8, XPD
0x2114,0xf426,0xc723,0x8892,0xbff6, XPD
0x6732,0x62bc,0x7de4,0x86f5,0x3ff6, XPD
0x9449,0x2d19,0xbf48,0x92b7,0x3ff0, XPD
0x6aa0,0x19fb,0x2e9d,0xcceb,0x3ff2, XPD
0xe484,0xc4ee,0x7659,0x9095,0x3ff0, XPD
0x1dcc,0xf144,0xc00c,0x9247,0xbff0, XPD
0xdfa5,0x5d17,0x2500,0xa9f3,0x3fe9, XPD
0xf8da,0x1811,0x21ee,0xc5af,0xbffd, XPD
0x1f70,0x719b,0xaf10,0x9abc,0x3ffd, XPD
0xcc3d,0x36d5,0x66a8,0x87e9,0xc000, XPD
0xe678,0xc70c,0x5525,0xae5a,0xbffd, XPD
0xfa51,0x8321,0x12d1,0x94ea,0x3ffe, XPD
0x722f,0x798a,0x8c6a,0x9401,0x3ffa, XPD
0xb83c,0x3073,0xf19b,0xb9e5,0xbff8, XPD
0xc320,0xbc1b,0x8ddf,0xb00f,0x3ff6, XPD
0x27d6,0x65dc,0xef54,0xdd13,0xbff7, XPD
0x984e,0xdb57,0xff51,0xf760,0xbff3, XPD
0x31ad,0xaee8,0xc63b,0xed7f,0x3ff1, XPD
0x90b1,0x4db9,0xac2c,0xbbce,0xbff0, XPD
0x7731,0xd1c6,0x7a05,0xea0b,0x3ff1, XPD
0x8527,0xb467,0xbe45,0x8610,0x3fee, XPD
0xd169,0x20b6,0x8e7f,0x898d,0x3feb, XPD
0xd014,0x1c7e,0x8db7,0x8dc9,0xbffc, XPD
0x9b21,0x7f40,0xbee1,0x9631,0x3ffd, XPD
0x1ba5,0x1fb9,0x9120,0xc2f6,0xbffe, XPD
0xd16d,0x1090,0xf860,0xa93b,0xbffd, XPD
0xcfea,0xd7ed,0x7e4d,0xd59d,0x3ffc, XPD
0x55f1,0xd81f,0x3296,0x8fa9,0x3ffa, XPD
0x9e4e,0xf6dc,0x7fdf,0x8555,0xbff7, XPD
0x797e,0x0960,0x5b21,0xaae4,0x3ff6, XPD
0x0359,0x504d,0xf928,0x9e90,0xbff6, XPD
0x4f40,0xf9c3,0xcd59,0xf01d,0xbff3, XPD
0xd618,0x1c5d,0x2494,0xaa58,0x3ff0, XPD
0xb6dd,0xaf8d,0x3132,0xb64b,0xbff0, XPD
0x5030,0x78bf,0xe337,0xa7dd,0x3ff0, XPD
0x3db1,0xb6a3,0x2a13,0x8221,0x3fee, XPD
0x04ff,0x6457,0x45f6,0xc551,0x3fe9, XPD
0x5a94,0x3c56,0x3e1f,0x83f7,0xbffa, XPD
0x0be3,0x034a,0x7513,0xfda3,0x3ffe, XPD
0x28a8,0x2266,0x2da6,0xe1ed,0x3fff, XPD
0x2ad2,0x5864,0xe694,0x84e6,0xbfff, XPD
0xb12f,0x263c,0xdd05,0xbca6,0xbffd, XPD
0x2b50,0xb5a2,0xe4d8,0xac3f,0x3ff9, XPD
0x54e3,0x252a,0xbd8d,0xbcd6,0xbffa, XPD
0xcb52,0x036a,0x4ef8,0xc830,0x3ff9, XPD
0x7f13,0x58e7,0xf385,0x8774,0x3ff7, XPD
0xf419,0xffc4,0x71c3,0xb740,0x3ff6, XPD
0x861a,0xc66e,0xcd61,0x90f0,0x3ff6, XPD
0x43c6,0x2f73,0xe63e,0x8e1c,0xbff4, XPD
0xe9fe,0x75f0,0x3a4f,0x9042,0x3ff2, XPD
0x660d,0xa227,0x10d0,0xce2c,0xbff2, XPD
0xe264,0x95d0,0x9394,0xe3fd,0xbff0, XPD
0xa931,0xb65d,0xa9d9,0xdf89,0xbff9, XPD
0x33c7,0x6312,0x62c2,0x8059,0xbfff, XPD
0x6786,0x1c17,0x64db,0xbf59,0x3fff, XPD
0x8286,0xc022,0x57d4,0x8681,0x3fff, XPD
0x246b,0xc515,0x7f7d,0x9fc7,0xbffd, XPD
0xa745,0xb941,0xda9e,0xae53,0xbff9, XPD
0xa930,0xb04b,0x0c33,0x9ff0,0xbffa, XPD
0x6aac,0xe9ab,0x8d9c,0xca9a,0xbff9, XPD
0x0bd5,0x1fe1,0x85c3,0xe573,0x3ff6, XPD
0xb83e,0x44c4,0x61eb,0xb976,0xbff6, XPD
0x86c5,0xe2ff,0x1366,0xf584,0x3ff5, XPD
0xc16c,0x7a8d,0xc9b9,0x8fd3,0x3ff4, XPD
0xeaea,0xad75,0x5ce4,0xf45c,0x3ff1, XPD
0xf303,0xe966,0xc9ff,0xd0a8,0x3ff2, XPD
0xefe5,0xa5ba,0xec2e,0xc118,0xbff0, XPD
0x695b,0x19a3,0xf63a,0xeed4,0xbff8, XPD
0x3cee,0x8591,0x08c3,0xa023,0xbffd, XPD
0xdcdc,0x164f,0xe7a2,0xcc70,0x3ffe, XPD
0xf945,0x9bd5,0x4d5a,0xa7d1,0x3ffd, XPD
0xd0c9,0x4f37,0x0e9b,0xaab6,0xbffc, XPD
0x5893,0x4535,0xa4e0,0xd980,0xbff7, XPD
0x7ed7,0xfca6,0x618c,0xaae1,0xbff9, XPD
0x9884,0x5e45,0x1ef7,0xfcc8,0xbff7, XPD
0xf482,0x53fa,0x659d,0xf526,0x3ff5, XPD
0x7c45,0xddde,0x2699,0xe765,0xbff4, XPD
0x30b6,0xafc4,0x29a1,0x8328,0x3ff5, XPD
0x9f99,0x5542,0xc9a5,0xb372,0x3ff2, XPD
0x6958,0xfcda,0x30a5,0x828a,0x3ff1, XPD
0xaac5,0x0a8e,0x2609,0x822b,0x3ff1, XPD
0x385d,0x20a1,0x0d7e,0xce4f,0xbfef, XPD
0x67ab,0xfaa3,0x0bf2,0xe2bf,0xbffd, XPD
0x1431,0x3df7,0x25ef,0xa7d1,0x3ffd, XPD
0x2351,0x0015,0x4537,0xae88,0x4000, XPD
0x3f7b,0x0bef,0xb2a4,0xa52c,0xbffd, XPD
0xf810,0x7089,0xa3e6,0xdbb8,0xbffd, XPD
0x061c,0x2fe2,0xdff6,0xf54b,0xbff8, XPD
0x224a,0x8eb7,0xf194,0xe79e,0xbffb, XPD
0x20fc,0xce98,0x453c,0xef1b,0x3ff7, XPD
0x2272,0x99fc,0x2dc4,0x9559,0xbff8, XPD
0x6006,0x9ff4,0x5bc7,0xa2d7,0x3ff6, XPD
0x1ec5,0xab24,0xd820,0xed40,0x3ff6, XPD
0x5259,0xfbe7,0x3c49,0x8e03,0x3ff2, XPD
0x3be4,0xfafe,0xae90,0xc15b,0x3ff5, XPD
0xcd84,0xe436,0xf8c3,0xdf07,0xbff1, XPD
0xc9f0,0xd207,0x4674,0xcd8f,0x3ff2, XPD
0x444a,0x86c7,0x1225,0xd490,0xbffa, XPD
0x3a84,0xd553,0xcb88,0xaeda,0xbfff, XPD
0xdfc0,0x042e,0x6ac5,0xa39d,0x3ffd, XPD
0x2c2a,0x0fdb,0xfa04,0xac19,0x3fff, XPD
0xec6f,0xe120,0x28ac,0xcdfa,0xbffa, XPD
0x0101,0x4363,0x52e6,0xff95,0x3ffa, XPD
0xabc2,0x7c25,0xe8b8,0xd921,0xbff8, XPD
0x48d2,0x83a1,0x4434,0xf922,0xbff9, XPD
0x1020,0xd87a,0x9b17,0x8c01,0xbff5, XPD
0x37ce,0xce84,0x963f,0xa9ab,0xbff8, XPD
0x2487,0x7733,0x9d8e,0xde69,0x3ff3, XPD
0xe9ec,0x0aa9,0xda3d,0x93f7,0xbff4, XPD
0x46c3,0x54d7,0x5d05,0xb543,0x3ff2, XPD
0x8cef,0x7f0d,0x62c6,0xe862,0x3ff3, XPD
0x5180,0x6fc8,0x919a,0xc0b3,0x3fef, XPD
0x4b84,0x8b57,0x0d1c,0xb3cc,0xbffb, XPD
0xb11c,0xf620,0x5ea2,0xd88e,0xbffe, XPD
0x83e9,0x9d21,0xf293,0x8a64,0x3ffe, XPD
0x81d6,0x8300,0x7524,0xd525,0x3ffe, XPD
0xf35a,0x2872,0x040c,0xae3a,0xbffb, XPD
0x71fc,0x2a7c,0xd9fd,0x9e44,0x3ffa, XPD
0xd10e,0xe0ca,0x9042,0xb7a9,0xbff9, XPD
0x72e4,0xc712,0x70d2,0x9a46,0xbff9, XPD
0xc8d7,0x1694,0x8804,0xecd9,0xbff5, XPD
0x28c2,0xaa89,0xa200,0xd222,0xbff7, XPD
0xc2a8,0x29b5,0xe8e7,0xbc20,0x3ff4, XPD
0xdd89,0x6e7a,0xe641,0xb741,0xbff3, XPD
0x89c0,0x4520,0x7873,0x9952,0x3ff3, XPD
0xba02,0x8866,0x3001,0x8fe7,0x3ff3, XPD
0x79a8,0xb088,0x4e04,0xa2ff,0x3ff0, XPD
};

short RK00i[] = {0x9673,0x7f73,0x211a,0xfffb,0x3ffe, XPD};
short RK01i[] = {0x8fe4,0x3565,0xa8be,0xb732,0xbff8, XPD};
short RK02i[] = {0x758a,0x4022,0x473d,0x9f38,0xbff7, XPD};
short RK10i[] = {0xbc29,0x7734,0xa8c3,0xb732,0x3ff8, XPD};
short RK11i[] = {0x65df,0x4b69,0xe72a,0xfffb,0x3ffe, XPD};
short RK12i[] = {0xb6fb,0x11e1,0x30e9,0xe3db,0xbfef, XPD};
short RK20i[] = {0x4b30,0x0e3c,0x4725,0x9f38,0x3ff7, XPD};
short RK21i[] = {0x06d2,0xd97a,0x186c,0xe3ec,0xbfef, XPD};
short RK22i[] = {0x37c9,0x340a,0x39f0,0xffff,0x3ffe, XPD};
short RL00i[] = {0xda14,0x5ac6,0x21c1,0xfffb,0x3ffe, XPD};

short RL01i[] = {0x51b0,0x8730,0x1eee,0xb725,0xbff8, XPD};
short RL02i[] = {0xfe33,0xf924,0x834c,0x9f33,0xbff7, XPD};
short RL10i[] = {0x203d,0xad12,0x1eeb,0xb725,0x3ff8, XPD};
short RL11i[] = {0x6580,0x4bc1,0xe7c5,0xfffb,0x3ffe, XPD};
short RL12i[] = {0xb887,0x8995,0x9199,0xe3d0,0xbfef, XPD};
short RL20i[] = {0x7326,0x18e6,0x835a,0x9f33,0x3ff7, XPD};
short RL21i[] = {0x00a4,0xbe4b,0x6591,0xe3c7,0xbfef, XPD};
short RL22i[] = {0x7bc9,0x0f05,0x39fc,0xffff,0x3ffe, XPD};

#endif
#if MIEEE
long mpperihdi[] = {
0x40140000,0x94e8be78,0x4a0e410b,
0x40140000,0x94e997f1,0x43bf7271,
0x40140000,0x94e81c62,0x602c9082,
0x40140000,0x94cbd46a,0x93f290ac,
0x40140000,0x94c31b1f,0x3b645a1d,
};
long mpperiodi[] = {
0x40090000,0xd2123c21,0x187e7c07,
0x40090000,0xd255af3a,0x14cec41e,
0x40090000,0xa5bc07dd,0x44135547,
0x40090000,0xa8499d88,0x3ba3443d,
0x40090000,0xc903f7ce,0xd916872b,
};
long chi[] = {
0xbffd0000,0xeadab8e6,0x70e82ef3,
0x3ffe0000,0xcb6b48c9,0x03388620,
0xc0000000,0x95a2f574,0x24107b6b,
0xbffe0000,0xe64250fe,0x93747684,
0x3ffe0000,0xa68d823d,0xc30cf073,
0x3ffb0000,0xce0ef69b,0xe3fec9d5,
0xbff90000,0x97ed1f96,0xcef9ab89,
0x3ff70000,0xb5585783,0x558ec1a5,
0xbff70000,0xdbd4ad04,0xf90b91b7,
0xbff50000,0xa3bf6321,0xef125ea0,
0x3ff30000,0x87b34d7a,0xe05530a9,
0xbff10000,0x9f5603ca,0x022bf94b,
0x3ff10000,0xd162f3ed,0x3cf350e4,
0x3fef0000,0xa98a5e31,0x6c6ee899,
0xbfec0000,0xb2a29c7c,0x23ecfa1d,
0x3ffc0000,0xd4532fd3,0x641a1d2e,
0x3fff0000,0xaa2c104c,0xa169d7c0,
0x3fff0000,0x8748322d,0xc6e1d4e5,
0xbfff0000,0xc0a016e8,0x48af137a,
0xbffd0000,0x969351a9,0x7ab68dd6,
0x3ffc0000,0xac614f29,0x41781ae2,
0x3ff80000,0x895a2154,0x9528b606,
0x3ff80000,0x97b4b714,0xec281ad3,
0x3ff60000,0xc6be16d9,0x7d9595a0,
0xbff60000,0x88fc0daf,0xc4a7d94f,
0xbff10000,0xf55d90d8,0xdd4f0fed,
0xbff20000,0x854b4435,0xa08aa0e8,
0xbff00000,0xbd4cdad6,0x05c89c20,
0x3ff00000,0x8dd4a79f,0x1ee2d655,
0x3feb0000,0xa17fb3bc,0xc44e13ab,
0x3ffc0000,0xc3dc0443,0xaff0ffc5,
0x3ffd0000,0xedc38db2,0x191e2ee1,
0x3ffe0000,0xf9950890,0x5beffbca,
0xbffe0000,0x86915657,0x8b004c05,
0xbffd0000,0x8ae606fc,0x9d7cd402,
0x3ffa0000,0xf0d95003,0xe3a59332,
0x3ff70000,0xfd66b5e2,0x709eb0c5,
0x3ff60000,0xd3f684e5,0x33eefa7a,
0x3ff60000,0xb7549030,0x7fc57d23,
0xbff40000,0xbf64f561,0x20bad544,
0xbff10000,0xe25675c0,0x1278dd09,
0xbff00000,0xba3ce7cd,0x654a8c84,
0xbff00000,0xae9ec832,0xd1f3ef05,
0x3fee0000,0xc62a527c,0x76b09838,
0x3feb0000,0x94f98e02,0x91816a50,
0x3ffa0000,0xafbbe834,0x76198c47,
0x3fff0000,0xb1183340,0xe0404d6a,
0xbfff0000,0xcbf27368,0xaea52c42,
0xbfff0000,0xb8f0003e,0x107e26a5,
0x3ffd0000,0xa7d86b87,0x227e78ec,
0x3ff90000,0xc14bce5a,0x8870776c,
0x3ffa0000,0xb1c282c9,0x64dd0ec4,
0x3ffa0000,0x8cd8578b,0x1da05d24,
0xbff60000,0xdb8b0ad1,0x49ff541b,
0x3ff70000,0x8d68fe16,0xdf10c49c,
0xbff60000,0x872897ee,0xac231518,
0xbff40000,0xbeec1bea,0x5f9bba10,
0xbff20000,0xa4b021da,0xce2a461d,
0xbff30000,0x9b0732bf,0xa0f2a48d,
0x3ff00000,0xc26b1bea,0x0919ef2e,
0xbffb0000,0x806b74c0,0xc64d6970,
0x3ffe0000,0xe7bd3de2,0xe7cf938b,
0x40000000,0x95096ac6,0xdd69b439,
0xbffe0000,0xf20092f5,0xfa3310db,
0xbffd0000,0xf54f6213,0x23b9c3b0,
0x3ff80000,0xfcf0a470,0xc5b82a99,
0xbffb0000,0x81e66e9c,0xaf91b9f6,
0x3ff90000,0xb84df81f,0xe3e09a8b,
0x3ff70000,0xa06f0d2d,0xdb36a385,
0x3ff60000,0xb90b40e6,0x7dd16aa0,
0x3ff60000,0xc5899a44,0xfb380e5c,
0xbff30000,0xf9d56ed5,0x1c5d20c4,
0x3ff20000,0xf0b20d0d,0xf919a99d,
0xbff20000,0xcadd2bba,0x7e4fa18a,
0xbff10000,0x8e12d6b1,0xfbea055f,
0x3ff80000,0x97cccf20,0xfbebcc99,
0xbffd0000,0x8f84ed1b,0x35569768,
0xbffd0000,0xb02baeef,0x43ba3993,
0x3ffd0000,0x95e01196,0xf3b2c0bc,
0x3ffb0000,0x90fc5d61,0x1e080dcc,
0xbff70000,0x9ca6305e,0x865f6aa4,
0x3ff80000,0x998cc833,0x0fe627b1,
0xbff70000,0xe448ddad,0xa40e1e15,
0xbff40000,0xbda481a6,0xf51a9b53,
0xbff40000,0xe533518d,0xdd8a181c,
0xbff30000,0xe98060a5,0xe1788aca,
0x3ff20000,0x9ab9b479,0x944a8b23,
0xbff00000,0x8e422943,0xee6c5757,
0x3ff00000,0xfb45d53e,0x5faa16dc,
0x3fee0000,0xa7f090eb,0x2b186400,
0xbffb0000,0xf43e4156,0xdcc2555b,
0xbfff0000,0xa8d5d8a6,0x40e73be7,
0xbffe0000,0xa7ec0863,0x101de1a1,
0x3fff0000,0xbe3d061f,0x02d31459,
0x3ffc0000,0xb7fcaecd,0x59ea3773,
0xbffc0000,0xa17db6c0,0x2fe395e8,
0xbff60000,0xe5ae734b,0x54925e07,
0xbff80000,0xc01a14d1,0x3fa25afd,
0xbff60000,0x8892c723,0xf4262114,
0x3ff60000,0x86f57de4,0x62bc6732,
0x3ff00000,0x92b7bf48,0x2d199449,
0x3ff20000,0xcceb2e9d,0x19fb6aa0,
0x3ff00000,0x90957659,0xc4eee484,
0xbff00000,0x9247c00c,0xf1441dcc,
0x3fe90000,0xa9f32500,0x5d17dfa5,
0xbffd0000,0xc5af21ee,0x1811f8da,
0x3ffd0000,0x9abcaf10,0x719b1f70,
0xc0000000,0x87e966a8,0x36d5cc3d,
0xbffd0000,0xae5a5525,0xc70ce678,
0x3ffe0000,0x94ea12d1,0x8321fa51,
0x3ffa0000,0x94018c6a,0x798a722f,
0xbff80000,0xb9e5f19b,0x3073b83c,
0x3ff60000,0xb00f8ddf,0xbc1bc320,
0xbff70000,0xdd13ef54,0x65dc27d6,
0xbff30000,0xf760ff51,0xdb57984e,
0x3ff10000,0xed7fc63b,0xaee831ad,
0xbff00000,0xbbceac2c,0x4db990b1,
0x3ff10000,0xea0b7a05,0xd1c67731,
0x3fee0000,0x8610be45,0xb4678527,
0x3feb0000,0x898d8e7f,0x20b6d169,
0xbffc0000,0x8dc98db7,0x1c7ed014,
0x3ffd0000,0x9631bee1,0x7f409b21,
0xbffe0000,0xc2f69120,0x1fb91ba5,
0xbffd0000,0xa93bf860,0x1090d16d,
0x3ffc0000,0xd59d7e4d,0xd7edcfea,
0x3ffa0000,0x8fa93296,0xd81f55f1,
0xbff70000,0x85557fdf,0xf6dc9e4e,
0x3ff60000,0xaae45b21,0x0960797e,
0xbff60000,0x9e90f928,0x504d0359,
0xbff30000,0xf01dcd59,0xf9c34f40,
0x3ff00000,0xaa582494,0x1c5dd618,
0xbff00000,0xb64b3132,0xaf8db6dd,
0x3ff00000,0xa7dde337,0x78bf5030,
0x3fee0000,0x82212a13,0xb6a33db1,
0x3fe90000,0xc55145f6,0x645704ff,
0xbffa0000,0x83f73e1f,0x3c565a94,
0x3ffe0000,0xfda37513,0x034a0be3,
0x3fff0000,0xe1ed2da6,0x226628a8,
0xbfff0000,0x84e6e694,0x58642ad2,
0xbffd0000,0xbca6dd05,0x263cb12f,
0x3ff90000,0xac3fe4d8,0xb5a22b50,
0xbffa0000,0xbcd6bd8d,0x252a54e3,
0x3ff90000,0xc8304ef8,0x036acb52,
0x3ff70000,0x8774f385,0x58e77f13,
0x3ff60000,0xb74071c3,0xffc4f419,
0x3ff60000,0x90f0cd61,0xc66e861a,
0xbff40000,0x8e1ce63e,0x2f7343c6,
0x3ff20000,0x90423a4f,0x75f0e9fe,
0xbff20000,0xce2c10d0,0xa227660d,
0xbff00000,0xe3fd9394,0x95d0e264,
0xbff90000,0xdf89a9d9,0xb65da931,
0xbfff0000,0x805962c2,0x631233c7,
0x3fff0000,0xbf5964db,0x1c176786,
0x3fff0000,0x868157d4,0xc0228286,
0xbffd0000,0x9fc77f7d,0xc515246b,
0xbff90000,0xae53da9e,0xb941a745,
0xbffa0000,0x9ff00c33,0xb04ba930,
0xbff90000,0xca9a8d9c,0xe9ab6aac,
0x3ff60000,0xe57385c3,0x1fe10bd5,
0xbff60000,0xb97661eb,0x44c4b83e,
0x3ff50000,0xf5841366,0xe2ff86c5,
0x3ff40000,0x8fd3c9b9,0x7a8dc16c,
0x3ff10000,0xf45c5ce4,0xad75eaea,
0x3ff20000,0xd0a8c9ff,0xe966f303,
0xbff00000,0xc118ec2e,0xa5baefe5,
0xbff80000,0xeed4f63a,0x19a3695b,
0xbffd0000,0xa02308c3,0x85913cee,
0x3ffe0000,0xcc70e7a2,0x164fdcdc,
0x3ffd0000,0xa7d14d5a,0x9bd5f945,
0xbffc0000,0xaab60e9b,0x4f37d0c9,
0xbff70000,0xd980a4e0,0x45355893,
0xbff90000,0xaae1618c,0xfca67ed7,
0xbff70000,0xfcc81ef7,0x5e459884,
0x3ff50000,0xf526659d,0x53faf482,
0xbff40000,0xe7652699,0xddde7c45,
0x3ff50000,0x832829a1,0xafc430b6,
0x3ff20000,0xb372c9a5,0x55429f99,
0x3ff10000,0x828a30a5,0xfcda6958,
0x3ff10000,0x822b2609,0x0a8eaac5,
0xbfef0000,0xce4f0d7e,0x20a1385d,
0xbffd0000,0xe2bf0bf2,0xfaa367ab,
0x3ffd0000,0xa7d125ef,0x3df71431,
0x40000000,0xae884537,0x00152351,
0xbffd0000,0xa52cb2a4,0x0bef3f7b,
0xbffd0000,0xdbb8a3e6,0x7089f810,
0xbff80000,0xf54bdff6,0x2fe2061c,
0xbffb0000,0xe79ef194,0x8eb7224a,
0x3ff70000,0xef1b453c,0xce9820fc,
0xbff80000,0x95592dc4,0x99fc2272,
0x3ff60000,0xa2d75bc7,0x9ff46006,
0x3ff60000,0xed40d820,0xab241ec5,
0x3ff20000,0x8e033c49,0xfbe75259,
0x3ff50000,0xc15bae90,0xfafe3be4,
0xbff10000,0xdf07f8c3,0xe436cd84,
0x3ff20000,0xcd8f4674,0xd207c9f0,
0xbffa0000,0xd4901225,0x86c7444a,
0xbfff0000,0xaedacb88,0xd5533a84,
0x3ffd0000,0xa39d6ac5,0x042edfc0,
0x3fff0000,0xac19fa04,0x0fdb2c2a,
0xbffa0000,0xcdfa28ac,0xe120ec6f,
0x3ffa0000,0xff9552e6,0x43630101,
0xbff80000,0xd921e8b8,0x7c25abc2,
0xbff90000,0xf9224434,0x83a148d2,
0xbff50000,0x8c019b17,0xd87a1020,
0xbff80000,0xa9ab963f,0xce8437ce,
0x3ff30000,0xde699d8e,0x77332487,
0xbff40000,0x93f7da3d,0x0aa9e9ec,
0x3ff20000,0xb5435d05,0x54d746c3,
0x3ff30000,0xe86262c6,0x7f0d8cef,
0x3fef0000,0xc0b3919a,0x6fc85180,
0xbffb0000,0xb3cc0d1c,0x8b574b84,
0xbffe0000,0xd88e5ea2,0xf620b11c,
0x3ffe0000,0x8a64f293,0x9d2183e9,
0x3ffe0000,0xd5257524,0x830081d6,
0xbffb0000,0xae3a040c,0x2872f35a,
0x3ffa0000,0x9e44d9fd,0x2a7c71fc,
0xbff90000,0xb7a99042,0xe0cad10e,
0xbff90000,0x9a4670d2,0xc71272e4,
0xbff50000,0xecd98804,0x1694c8d7,
0xbff70000,0xd222a200,0xaa8928c2,
0x3ff40000,0xbc20e8e7,0x29b5c2a8,
0xbff30000,0xb741e641,0x6e7add89,
0x3ff30000,0x99527873,0x452089c0,
0x3ff30000,0x8fe73001,0x8866ba02,
0x3ff00000,0xa2ff4e04,0xb08879a8,
};
long RK00i[] = {0x3ffe0000,0xfffb211a,0x7f739673};
long RK01i[] = {0xbff80000,0xb732a8be,0x35658fe4};
long RK02i[] = {0xbff70000,0x9f38473d,0x4022758a};
long RK10i[] = {0x3ff80000,0xb732a8c3,0x7734bc29};
long RK11i[] = {0x3ffe0000,0xfffbe72a,0x4b6965df};
long RK12i[] = {0xbfef0000,0xe3db30e9,0x11e1b6fb};
long RK20i[] = {0x3ff70000,0x9f384725,0x0e3c4b30};
long RK21i[] = {0xbfef0000,0xe3ec186c,0xd97a06d2};
long RK22i[] = {0x3ffe0000,0xffff39f0,0x340a37c9};
long RL00i[] = {0x3ffe0000,0xfffb21c1,0x5ac6da14};
long RL01i[] = {0xbff80000,0xb7251eee,0x873051b0};
long RL02i[] = {0xbff70000,0x9f33834c,0xf924fe33};
long RL10i[] = {0x3ff80000,0xb7251eeb,0xad12203d};
long RL11i[] = {0x3ffe0000,0xfffbe7c5,0x4bc16580};
long RL12i[] = {0xbfef0000,0xe3d09199,0x8995b887};
long RL20i[] = {0x3ff70000,0x9f33835a,0x18e67326};
long RL21i[] = {0xbfef0000,0xe3c76591,0xbe4b00a4};
long RL22i[] = {0x3ffe0000,0xffff39fc,0x0f057bc9};
#endif

#if !UNK
#define mpperihd ((long double *)mpperihdi)
#define mpperiod ((long double *)mpperiodi)
#define ch ((long double *)chi)

#define RK00 (*(long double *)RK00i)
#define RK01 (*(long double *)RK01i)
#define RK02 (*(long double *)RK02i)
#define RK10 (*(long double *)RK10i)
#define RK11 (*(long double *)RK11i)
#define RK12 (*(long double *)RK12i)
#define RK20 (*(long double *)RK20i)
#define RK21 (*(long double *)RK21i)
#define RK22 (*(long double *)RK22i)

#define RL00 (*(long double *)RL00i)
#define RL01 (*(long double *)RL01i)
#define RL02 (*(long double *)RL02i)
#define RL10 (*(long double *)RL10i)
#define RL11 (*(long double *)RL11i)
#define RL12 (*(long double *)RL12i)
#define RL20 (*(long double *)RL20i)
#define RL21 (*(long double *)RL21i)
#define RL22 (*(long double *)RL22i)
#endif
#if UNK
/* Julian date of perihelion
 */
long double mpperihd[] = {
2439727.61747L, /* Ceres */
2439781.98561L, /* Pallas */
2439687.09607L, /* Vesta */
2437877.10408L, /* Iris */
2437318.78050L /* Bamberga */
};

/* Orbital period, in days
 */
long double mpperiod[] = {
1680.56984L,
1682.67764L,
1325.87596L,
1346.30048L,
1608.124L
};

/* Chebyshev coefficients for one complete revolution
 */
long double ch[] = {
/* Ceres, X coordinate */
-4.586999684620e-1L,  7.946057787660e-1L, -2.338071216029e0L,
-8.994494077228e-1L,  6.505967522417e-1L,  1.006144777607e-1L,
-1.854568643327e-2L,  5.534212809533e-3L, -6.708702547208e-3L,
-1.249295092718e-3L,  2.588279743209e-4L, -7.597732320643e-5L,
 9.984328008949e-5L,  2.021080252791e-5L, -2.661871914125e-6L,
/* Ceres, Y coordinate */
 2.073485825322e-1L,  1.329469716470e0L,   1.056890747406e0L,
-1.504885543270e0L,  -2.940927047142e-1L,  1.683399552279e-1L,
 8.383305119724e-3L,  9.259394446135e-3L,  3.032570436963e-3L,
-2.090218869642e-3L, -1.169993836602e-4L, -1.271190734330e-4L,
-4.513268808669e-5L,  3.381506984284e-5L,  1.203260096346e-6L,
/* Ceres, Z coordinate */
 1.912689844430e-1L,  4.643825797944e-1L,  9.749303200187e-1L,
-5.256551707962e-1L, -2.712862189644e-1L,  5.880099540648e-2L,
 7.733191309744e-3L,  3.234298176906e-3L,  2.797398760321e-3L,
-7.301115768150e-4L, -1.079262419828e-4L, -4.440257835932e-5L,
-4.163271004842e-5L,  1.181157357331e-5L,  1.109948926800e-6L,

/* Pallas */
 4.290381150293e-2L,  1.383551031759e0L,  -1.593336511721e0L,
-1.444824247651e0L,   3.278230287990e-1L,  2.359571746298e-2L,
 4.339839066079e-2L,  3.438600727700e-2L, -3.349962368965e-3L,
 4.315494601308e-3L, -2.062356083727e-3L, -7.283107794077e-4L,
-1.570587316340e-4L, -2.956926691273e-4L,  4.635295741246e-5L,

-6.270495618154e-2L,  9.052313498467e-1L,  2.328699773983e0L,
-9.453212595492e-1L, -4.791212700237e-1L,  1.543823298115e-2L,
-6.342779555957e-2L,  2.249811612921e-2L,  4.896050868149e-3L,
 2.823546738447e-3L,  3.014183200304e-3L, -4.765200088867e-4L,
 2.295451275798e-4L, -1.934661374033e-4L, -6.774596619107e-5L,

 9.265138887655e-3L, -2.803110213960e-1L, -3.440832774230e-1L,
 2.927251335876e-1L,  7.079384758595e-2L, -4.780553453247e-3L,
 9.371943957602e-3L, -6.966694108343e-3L, -7.234291203889e-4L,
-8.743303801260e-4L, -4.453687185879e-4L,  1.475576496875e-4L,
-3.391705564484e-5L,  5.990810038801e-5L,  1.000998683469e-5L,

/* Vesta */
-1.192593674310e-1L, -1.319026070765e0L,  -6.559453241208e-1L,
 1.486237301962e0L,   1.796748459660e-1L, -1.577061228563e-1L,
-3.504660747540e-3L, -1.172496827075e-2L, -2.083943956635e-3L,
 2.059310183891e-3L,  3.498024045962e-5L,  1.954256820853e-4L,
 3.447147301144e-5L, -3.487593512748e-5L,  3.165561395639e-7L,

-3.861017802374e-1L,  3.022207934926e-1L, -2.123620666763e0L,
-3.405329331031e-1L,  5.816966783051e-1L,  3.613428926437e-2L,
-1.134632677418e-2L,  2.686473977278e-3L, -6.746761188693e-3L,
-4.718377988248e-4L,  1.132484047640e-4L, -4.477675310431e-5L,
 1.116012719484e-4L,  7.990920741835e-6L,  1.024849381860e-6L,

-1.384641783726e-1L,  2.933482790246e-1L, -7.615748122623e-1L,
-3.305356614367e-1L,  2.086086020481e-1L,  3.507346879406e-2L,
-4.069030226609e-3L,  2.607605217269e-3L, -2.419529752218e-3L,
-4.579857152266e-4L,  4.061324790585e-5L, -4.346220957083e-5L,
 4.002255160844e-5L,  7.756325500855e-6L,  3.675324354182e-7L,

/* Iris */
-3.221821085369e-2L,  9.907754107874e-1L,  1.765050607796e0L,
-1.038296530184e0L,  -3.684605664049e-1L,  2.102656075083e-2L,
-4.610322993751e-2L,  2.443709777445e-2L,  4.133814722713e-3L,
 2.796199583826e-3L,  2.211618558954e-3L, -5.421176461335e-4L,
 1.375758192501e-4L, -1.966210666429e-4L, -5.435717896184e-5L,

-2.728732274927e-2L, -1.002727837482e0L,   1.494915587413e0L,
 1.050822237856e0L,  -3.120689240482e-1L, -2.128021907064e-2L,
-3.904728666659e-2L, -2.473189982202e-2L,  3.501148373402e-3L,
-2.829932123194e-3L,  1.873137825390e-3L,  5.486575959091e-4L,
 1.165203058428e-4L,  1.989930460645e-4L, -4.603799673429e-5L,

-1.457714124299e-2L, -3.127672899889e-1L,  7.985977907868e-1L,
 3.277687237841e-1L, -1.667101175588e-1L, -6.637650018584e-3L,
-2.085942318799e-2L, -7.714285965209e-3L,  1.870345977899e-3L,
-8.827023325163e-4L,  1.000647623044e-3L,  1.711353200638e-4L,
 6.224622956062e-5L,  6.206920105118e-5L, -2.459392543221e-5L,

/* Bamberga */
-4.428638204675e-1L,  3.277675489944e-1L,  2.727067283355e0L,
-3.226066422703e-1L, -4.291430681990e-1L, -1.497170325524e-2L,
-1.130961297475e-1L,  7.296952063929e-3L, -9.115500159538e-3L,
 2.484760195898e-3L,  3.620198026209e-3L,  1.354338061210e-4L,
 1.475205478064e-3L, -1.063495348004e-4L,  1.960369766505e-4L,

-5.189520920287e-2L, -1.366052095235e0L,   3.195603719233e-1L,
 1.344542743667e0L,  -5.028739823144e-2L,  6.239826567272e-2L,
-1.325271345735e-2L, -3.041184731784e-2L, -1.068163091916e-3L,
-1.035585091379e-2L,  4.242182929449e-4L, -5.644537880111e-4L,
 1.728659993507e-4L,  4.432379144528e-4L,  2.297180181498e-5L,

-8.779154055939e-2L, -8.459223888259e-1L,  5.406028375993e-1L,
 8.326028074398e-1L, -8.507159387687e-2L,  3.863988067556e-2L,
-2.241972137474e-2L, -1.883241687594e-2L, -1.807020047676e-3L,
-6.412819960435e-3L,  7.176534798334e-4L, -3.495357888632e-4L,
 2.924387939936e-4L,  2.744733357682e-4L,  3.886158089997e-5
};

long double RK00 =  0.9999256791774783L;
long double RK01 = -0.0111815116768724L;
long double RK02 = -0.0048590038154553L;
long double RK10 =  0.0111815116959975L;
long double RK11 =  0.9999374845751042L;
long double RK12 = -0.0000271625775175L;
long double RK20 =  0.0048590037714450L;
long double RK21 = -0.0000271704492210L;
long double RK22 =  0.9999881946023742L;

long double RL00 =  0.9999257180268403L;
long double RL01 = -0.0111782838886141L;
long double RL02 = -0.0048584357372842L;
long double RL10 =  0.0111782838782385L;
long double RL11 =  0.9999375206641666L;
long double RL12 = -0.0000271576311202L;
long double RL20 =  0.0048584357611567L;
long double RL21 = -0.0000271533600777L;
long double RL22 =  0.9999881973626738L;
#endif /* UNK */
#endif /* LDOUBLE */

#if !(LDOUBLE)
#if UNK
/* Julian date of perihelion
 */
DOUBLE mpperihd[] = {
2439727.61747, /* Ceres */
2439781.98561, /* Pallas */
2439687.09607, /* Vesta */
2437877.10408, /* Iris */
2437318.78050 /* Bamberga */
};

/* Orbital period, in days
 */
DOUBLE mpperiod[] = {
1680.56984,
1682.67764,
1325.87596,
1346.30048,
1608.124
};

/* Chebyshev coefficients for one complete revolution
 */
DOUBLE ch[] = {
/* Ceres, X coordinate */
-4.586999684620e-1,  7.946057787660e-1, -2.338071216029e0,
-8.994494077228e-1,  6.505967522417e-1,  1.006144777607e-1,
-1.854568643327e-2,  5.534212809533e-3, -6.708702547208e-3,
-1.249295092718e-3,  2.588279743209e-4, -7.597732320643e-5,
 9.984328008949e-5,  2.021080252791e-5, -2.661871914125e-6,
/* Ceres, Y coordinate */
 2.073485825322e-1,  1.329469716470e0,   1.056890747406e0,
-1.504885543270e0,  -2.940927047142e-1,  1.683399552279e-1,
 8.383305119724e-3,  9.259394446135e-3,  3.032570436963e-3,
-2.090218869642e-3, -1.169993836602e-4, -1.271190734330e-4,
-4.513268808669e-5,  3.381506984284e-5,  1.203260096346e-6,
/* Ceres, Z coordinate */
 1.912689844430e-1,  4.643825797944e-1,  9.749303200187e-1,
-5.256551707962e-1, -2.712862189644e-1,  5.880099540648e-2,
 7.733191309744e-3,  3.234298176906e-3,  2.797398760321e-3,
-7.301115768150e-4, -1.079262419828e-4, -4.440257835932e-5,
-4.163271004842e-5,  1.181157357331e-5,  1.109948926800e-6,

/* Pallas */
 4.290381150293e-2,  1.383551031759e0,  -1.593336511721e0,
-1.444824247651e0,   3.278230287990e-1,  2.359571746298e-2,
 4.339839066079e-2,  3.438600727700e-2, -3.349962368965e-3,
 4.315494601308e-3, -2.062356083727e-3, -7.283107794077e-4,
-1.570587316340e-4, -2.956926691273e-4,  4.635295741246e-5,

-6.270495618154e-2,  9.052313498467e-1,  2.328699773983e0,
-9.453212595492e-1, -4.791212700237e-1,  1.543823298115e-2,
-6.342779555957e-2,  2.249811612921e-2,  4.896050868149e-3,
 2.823546738447e-3,  3.014183200304e-3, -4.765200088867e-4,
 2.295451275798e-4, -1.934661374033e-4, -6.774596619107e-5,

 9.265138887655e-3, -2.803110213960e-1, -3.440832774230e-1,
 2.927251335876e-1,  7.079384758595e-2, -4.780553453247e-3,
 9.371943957602e-3, -6.966694108343e-3, -7.234291203889e-4,
-8.743303801260e-4, -4.453687185879e-4,  1.475576496875e-4,
-3.391705564484e-5,  5.990810038801e-5,  1.000998683469e-5,

/* Vesta */
-1.192593674310e-1, -1.319026070765e0,  -6.559453241208e-1,
 1.486237301962e0,   1.796748459660e-1, -1.577061228563e-1,
-3.504660747540e-3, -1.172496827075e-2, -2.083943956635e-3,
 2.059310183891e-3,  3.498024045962e-5,  1.954256820853e-4,
 3.447147301144e-5, -3.487593512748e-5,  3.165561395639e-7,

-3.861017802374e-1,  3.022207934926e-1, -2.123620666763e0,
-3.405329331031e-1,  5.816966783051e-1,  3.613428926437e-2,
-1.134632677418e-2,  2.686473977278e-3, -6.746761188693e-3,
-4.718377988248e-4,  1.132484047640e-4, -4.477675310431e-5,
 1.116012719484e-4,  7.990920741835e-6,  1.024849381860e-6,

-1.384641783726e-1,  2.933482790246e-1, -7.615748122623e-1,
-3.305356614367e-1,  2.086086020481e-1,  3.507346879406e-2,
-4.069030226609e-3,  2.607605217269e-3, -2.419529752218e-3,
-4.579857152266e-4,  4.061324790585e-5, -4.346220957083e-5,
 4.002255160844e-5,  7.756325500855e-6,  3.675324354182e-7,

/* Iris */
-3.221821085369e-2,  9.907754107874e-1,  1.765050607796e0,
-1.038296530184e0,  -3.684605664049e-1,  2.102656075083e-2,
-4.610322993751e-2,  2.443709777445e-2,  4.133814722713e-3,
 2.796199583826e-3,  2.211618558954e-3, -5.421176461335e-4,
 1.375758192501e-4, -1.966210666429e-4, -5.435717896184e-5,

-2.728732274927e-2, -1.002727837482e0,   1.494915587413e0,
 1.050822237856e0,  -3.120689240482e-1, -2.128021907064e-2,
-3.904728666659e-2, -2.473189982202e-2,  3.501148373402e-3,
-2.829932123194e-3,  1.873137825390e-3,  5.486575959091e-4,
 1.165203058428e-4,  1.989930460645e-4, -4.603799673429e-5,

-1.457714124299e-2, -3.127672899889e-1,  7.985977907868e-1,
 3.277687237841e-1, -1.667101175588e-1, -6.637650018584e-3,
-2.085942318799e-2, -7.714285965209e-3,  1.870345977899e-3,
-8.827023325163e-4,  1.000647623044e-3,  1.711353200638e-4,
 6.224622956062e-5,  6.206920105118e-5, -2.459392543221e-5,

/* Bamberga */
-4.428638204675e-1,  3.277675489944e-1,  2.727067283355e0,
-3.226066422703e-1, -4.291430681990e-1, -1.497170325524e-2,
-1.130961297475e-1,  7.296952063929e-3, -9.115500159538e-3,
 2.484760195898e-3,  3.620198026209e-3,  1.354338061210e-4,
 1.475205478064e-3, -1.063495348004e-4,  1.960369766505e-4,

-5.189520920287e-2, -1.366052095235e0,   3.195603719233e-1,
 1.344542743667e0,  -5.028739823144e-2,  6.239826567272e-2,
-1.325271345735e-2, -3.041184731784e-2, -1.068163091916e-3,
-1.035585091379e-2,  4.242182929449e-4, -5.644537880111e-4,
 1.728659993507e-4,  4.432379144528e-4,  2.297180181498e-5,

-8.779154055939e-2, -8.459223888259e-1,  5.406028375993e-1,
 8.326028074398e-1, -8.507159387687e-2,  3.863988067556e-2,
-2.241972137474e-2, -1.883241687594e-2, -1.807020047676e-3,
-6.412819960435e-3,  7.176534798334e-4, -3.495357888632e-4,
 2.924387939936e-4,  2.744733357682e-4,  3.886158089997e-5
};

DOUBLE RK00 =  0.9999256791774783;
DOUBLE RK01 = -0.0111815116768724;
DOUBLE RK02 = -0.0048590038154553;
DOUBLE RK10 =  0.0111815116959975;
DOUBLE RK11 =  0.9999374845751042;
DOUBLE RK12 = -0.0000271625775175;
DOUBLE RK20 =  0.0048590037714450;
DOUBLE RK21 = -0.0000271704492210;
DOUBLE RK22 =  0.9999881946023742;

DOUBLE RL00 =  0.9999257180268403;
DOUBLE RL01 = -0.0111782838886141;
DOUBLE RL02 = -0.0048584357372842;
DOUBLE RL10 =  0.0111782838782385;
DOUBLE RL11 =  0.9999375206641666;
DOUBLE RL12 = -0.0000271576311202;
DOUBLE RL20 =  0.0048584357611567;
DOUBLE RL21 = -0.0000271533600777;
DOUBLE RL22 =  0.9999881973626738;
#endif /* UNK */
#endif /* not LDOUBLE */

//------------------------------------------------------------------------------

#define NCOF 15
DOUBLE  pcofs[NCOF];
DOUBLE vcofs[NCOF];

DOUBLE FLOOR(), SQRT();

//------------------------------------------------------------------------------
aroids (JD, y)

  DOUBLE   JD;
  register DOUBLE y[];
{

  DOUBLE t, psum, vsum, per;
  DOUBLE v[3], p[3];
  DOUBLE *pans;
  int i, j, k, iob;

  //fprintf (stderr, "aroids .. \n");


  pans = &y[6*IAROIDS]; /* output pointer */

  for( iob = 0; iob < 5; iob++ )
  {
    per = mpperiod[iob];

    t = (JD - mpperihd[iob])/per;
    t -= FLOOR(t);
    t = Two * t - One;

    /* Chebyshev expansion
      */
    pcofs[0] = One;
    pcofs[1] = t;
    vcofs[0] = Zero;
    vcofs[1] = One;
    t *= Two;
    for( i=2; i<NCOF; i++ )
    {
      pcofs[i] = t * pcofs[i-1]  -  pcofs[i-2];
      vcofs[i] = Two * pcofs[i-1] + t * vcofs[i-1] - vcofs[i-2];
    }

    j = 45 * iob;

    for( k=0; k<3; k++ )
    {
      psum = Zero;
      vsum = Zero;
      for( i=NCOF-1; i>=0; i-- )
      {
        psum += pcofs[i] * ch[j+i];
        vsum += vcofs[i] * ch[j+i];
      }
      p[k] = psum;
      v[k] = vsum / ( Half * per );
      j += NCOF;
    }

#if DE200
    r118_200(p);
    r118_200(v);
#endif

    for( j=0; j<3; j++ )
    {
      /* Convert from heliocentric to barycentric coordinates */
      v[j] += y[(6*ISUN)+2*j];
      p[j] += y[(6*ISUN+1)+2*j];

      *pans++ = v[j];
      *pans++ = p[j];
    }
  }

}
//------------------------------------------------------------------------------

/* Rotate from DE118 (FK4) to DE200 (FK5)
 */
//------------------------------------------------------------------------------
r118_200( vec )

  DOUBLE vec[];
{
  DOUBLE x, y, z;

  x = RK00*vec[0] + RK01*vec[1] + RK02*vec[2];
  y = RK10*vec[0] + RK11*vec[1] + RK12*vec[2];
  z = RK20*vec[0] + RK21*vec[1] + RK22*vec[2];

  vec[0] = x;
  vec[1] = y;
  vec[2] = z;

}
//------------------------------------------------------------------------------

/* Rotate DE200 to DE118 */
r200_118( vec )
  DOUBLE vec[];
{
  DOUBLE x, y, z;

  x = RK00*vec[0] + RK10*vec[1] + RK20*vec[2];
  y = RK01*vec[0] + RK11*vec[1] + RK21*vec[2];
  z = RK02*vec[0] + RK12*vec[1] + RK22*vec[2];

  vec[0] = x;
  vec[1] = y;
  vec[2] = z;

}
//------------------------------------------------------------------------------
/* Rotate DE102 to DE200 */
//------------------------------------------------------------------------------
r102_200( vec )
  DOUBLE vec[];
{
  DOUBLE x, y, z;

  x = RL00*vec[0] + RL01*vec[1] + RL02*vec[2];
  y = RL10*vec[0] + RL11*vec[1] + RL12*vec[2];
  z = RL20*vec[0] + RL21*vec[1] + RL22*vec[2];

  vec[0] = x;
  vec[1] = y;
  vec[2] = z;
}
//------------------------------------------------------------------------------
/* Rotate DE200 to DE102 */
//------------------------------------------------------------------------------
r200_102( vec )
  DOUBLE vec[];
{
  DOUBLE x, y, z;

  x = RL00*vec[0] + RL10*vec[1] + RL20*vec[2];
  y = RL01*vec[0] + RL11*vec[1] + RL21*vec[2];
  z = RL02*vec[0] + RL12*vec[1] + RL22*vec[2];

  vec[0] = x;
  vec[1] = y;
  vec[2] = z;

}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
