//------------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>

#include "mconf.h"
#include "prec.h"
#include "ssystem.h"

//------------------------------------------------------------------------------

#if LDOUBLE
static union
{
  long double ld;
  unsigned short s[8];
}  funion;
#else
static union
{
  double d;
  unsigned short s[4];
}  funion;
#endif

//------------------------------------------------------------------------------
long double as_to_ld (char *s) 
{

  long double x;

#if LDOUBLE
  //printf ("LDOUBLE \n");
  asctoe64 (s, funion.s);
  x = funion.ld;
#else
  //printf ("DOUBLE \n");
  asctoe53 (s, funion.s);
  x = funion.d;
#endif

  return (x);
}
//------------------------------------------------------------------------------

/* rdnums.c
 */

static FILE *in;
char *fgets();

struct number
	{
	DOUBLE *val;
	int n;
	};


extern DOUBLE JD0, JD1, C, EMRAT, KG, AU, RADS, RADM, RADE, AE, AM;
extern DOUBLE K2M, LOVENO, PHASE, PSLP1, PSLPI, PSLPIA, PSLPIB;
extern DOUBLE JDEPOCH, LBET, LGAM, LGAMBET;
extern DOUBLE CMRSQ, CMR2, BMR2, AMR2;

extern DOUBLE yn0[], yn1[], GMs[], Je[], Jm[], Cnm[], Snm[]; // ���������?
              // � ��� ���� ������ �� ���������� ??? 

//------------------------------------------------------------------------------

struct number nums[] = {

#ifdef _0

  {&JD0,     1},               // Julian date of initial state vector
  {yn0,      6*(NBODY+FMASS)}, // ���� �� ������ Lunar librations (Euler angles)


  {&JD1,     1},               // Julian date of test state vector
  {yn1,      6*(NBODY+FMASS)},


  {&C,       1}, // Speed of light, au/d ;2.99792458e5  CLIGHT
  {&EMRAT,   1}, // Earth's mass divided by Moon's mass 
  {&KG,      1}, // Gaussian gravitational constant

  {GMs, NTOTAL}, // masses (GM values) of all objects:
                 // 
                 // 0.0,  /*place holder for libration */
                 // Mercury, Venus, ... Bamberga 

  {&AU,      1}, // Astronomical unit in meters
  {&RADS,    1},
  {&RADM,    1},
  {&RADE,    1},
  {&AE,      1},
  {Je,       3},
  {&AM,      1},
  {Jm,       3},

  {Cnm,      8},
  {Snm,      8},

  {&K2M,     1},
  {&LOVENO,  1},
  {&PHASE,   1},
  {&PSLP1,   1},
  {&PSLPI,   1},
  {&PSLPIA,  1},
  {&PSLPIB,  1},
  {&JDEPOCH, 1},
  {&LBET,    1},
  {&LGAM,    1},
  {&LGAMBET, 1},
  {&CMRSQ,   1},
  {&CMR2,    1},
  {&BMR2,    1},
  {&AMR2,    1},

#endif

  {NULL, 0}
};

//------------------------------------------------------------------------------
set_planets_0 ()
{


  // ; Initial solar system state vector

  JD0 = as_to_ld ("2440400.50");  // Julian date of initial state vector


  // ; yn0[]

  yn0[ 0] = as_to_ld ("1.00015500951082782e-4");  // phidot Lunar librations (Euler angles) 
  yn0[ 1] = as_to_ld ("6.03343019729665792e-3");  // phi
  yn0[ 2] = as_to_ld ("1.49039876740924137e-5");  // thetadot
  yn0[ 3] = as_to_ld ("3.82393007236771923e-1");  // theta
  yn0[ 4] = as_to_ld ("-1.19434270765326686e-4"); // psidot
  yn0[ 5] = as_to_ld ("1.28120436555489181e0");   // psi

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~```


  yn0[ 6] = as_to_ld ("3.7085069798210382186e-003"); /* Mercury */
  yn0[ 7] = as_to_ld ("3.6030663368975339466e-001"); // x
  yn0[ 8] = as_to_ld ("2.4854958767430945324e-002"); // ydot
  yn0[ 9] = as_to_ld ("-9.4812876741771684223e-002"); // y
  yn0[10] = as_to_ld ("1.2929109014677844626e-002"); // zdot
  yn0[11] = as_to_ld ("-8.7466840117233140436e-002"); // z

  yn0[12] = as_to_ld ("1.1156645711264016669e-002"); // /* Venus */
  yn0[13] = as_to_ld ("6.0786466491731583464e-001"); 
  yn0[14] = as_to_ld ("1.5494075513638794325e-002"); 
  yn0[15] = as_to_ld ("-3.5518362463675619232e-001"); 
  yn0[16] = as_to_ld ("6.2773904546696609267e-003"); 
  yn0[17] = as_to_ld ("-1.9824142909855515515e-001"); 

  yn0[18] = as_to_ld ("1.6833251020051496668e-002"); /* EMB */
  yn0[19] = as_to_ld ("1.0820747754938311664e-001"); 
  yn0[20] = as_to_ld ("1.5602036176919105255e-003"); 
  yn0[21] = as_to_ld ("-9.2711110739430602933e-001"); 
  yn0[22] = as_to_ld ("6.7646174015847273137e-004"); 
  yn0[23] = as_to_ld ("-4.0209347855944090112e-001"); 

  yn0[24] = as_to_ld ("1.4481919298277924969e-002"); /* Mars */
  yn0[25] = as_to_ld ("-1.2796408611369531836e-001"); 
  yn0[26] = as_to_ld ("8.0528538390447499843e-005"); 
  yn0[27] = as_to_ld ("-1.3262618005333617013e+000"); 
  yn0[28] = as_to_ld ("-3.5188931029397090065e-004"); 
  yn0[29] = as_to_ld ("-6.0530808652523961512e-001"); 

  yn0[30] = as_to_ld ("1.0053452569924098185e-003"); /* Jupiter */
  yn0[31] = as_to_ld ("-5.3896824544609061333e+000"); 
  yn0[32] = as_to_ld ("-6.5298425191689416643e-003"); 
  yn0[33] = as_to_ld ("-7.7026549518616593034e-001"); 
  yn0[34] = as_to_ld ("-2.8258787532429609536e-003"); 
  yn0[35] = as_to_ld ("-1.9866431165907522014e-001"); 

  yn0[36] = as_to_ld ("-3.1594662504012930114e-003"); /* Saturn */
  yn0[37] = as_to_ld ("7.9527768530257360864e+000"); 
  yn0[38] = as_to_ld ("4.3714634278372622354e-003"); 
  yn0[39] = as_to_ld ("4.5078822184006553686e+000"); 
  yn0[40] = as_to_ld ("1.9441395169137103763e-003"); 
  yn0[41] = as_to_ld ("1.5201955253183338898e+000"); 

  yn0[42] = as_to_ld ("1.7108310564806817248e-004"); /* Uranus */
  yn0[43] = as_to_ld ("-1.8278236586353147533e+001"); 
  yn0[44] = as_to_ld ("-3.7646704682815900043e-003"); 
  yn0[45] = as_to_ld ("-9.5764572881482056433e-001"); 
  yn0[46] = as_to_ld ("-1.6519678610257000136e-003"); 
  yn0[47] = as_to_ld ("-1.6132190397271035415e-001"); 

  yn0[48] = as_to_ld ("2.6225242764289213785e-003"); /* Neptune */
  yn0[49] = as_to_ld ("-1.6367191358770888335e+001"); 
  yn0[50] = as_to_ld ("-1.5277473123858904045e-003"); 
  yn0[51] = as_to_ld ("-2.3760896725373076342e+001"); 
  yn0[52] = as_to_ld ("-6.9183197562182804864e-004"); 
  yn0[53] = as_to_ld ("-9.3213866179497290101e+000"); 

  yn0[54] = as_to_ld ("2.8177758090360373050e-004"); /* Pluto */
  yn0[55] = as_to_ld ("-3.0447680255169362534e+001"); 
  yn0[56] = as_to_ld ("-3.1469590804946202045e-003"); 
  yn0[57] = as_to_ld ("-5.3177934960261367037e-001"); 
  yn0[58] = as_to_ld ("-1.0794238049289112837e-003"); 
  yn0[59] = as_to_ld ("9.0596584886274922101e+000"); 

  yn0[60] = as_to_ld ("5.98752118417335956e-4");   /* MOON */
  yn0[61] = as_to_ld ("-8.35703164195952601e-4"); 
  yn0[62] = as_to_ld ("-1.74153713527242722e-4"); 
  yn0[63] = as_to_ld ("-1.98543915768166071e-3"); 
  yn0[64] = as_to_ld ("-8.84771962437116281e-5"); 
  yn0[65] = as_to_ld ("-1.08326877048661754e-3"); 

  yn0[66] = as_to_ld ("-2.8369446340813151639e-007");  /* Sun */
  yn0[67] = as_to_ld ("4.5144118714356666407e-003"); 
  yn0[68] = as_to_ld ("5.1811944086463255444e-006"); 
  yn0[69] = as_to_ld ("7.2282841152065867346e-004"); 
  yn0[70] = as_to_ld ("2.2306588340621263489e-006"); 
  yn0[71] = as_to_ld ("2.4659100492567986271e-004"); /* Asteroid states follow */


}
//------------------------------------------------------------------------------
set_planets_1 ()
{

  // ; Test state vector

  JD1 = as_to_ld ("2440800.50");  // Julian date of test state vector


  // ;  yn1[]

  yn1[ 0] = as_to_ld ("1.19798509852348462e-4");  /* Librations */
  yn1[ 1] = as_to_ld (" 3.26749396278042651e-2"); 
  yn1[ 2] = as_to_ld ("-6.85863913398886621e-5"); 
  yn1[ 3] = as_to_ld (" 3.85255316277713118e-1"); 
  yn1[ 4] = as_to_ld ("-9.87921606265917336e-5"); 
  yn1[ 5] = as_to_ld (" 1.25594900035445747e0"); 

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~```


  yn1[ 6] = as_to_ld (" 1.1640059959070938016e-002");   /* Mercury */
  yn1[ 7] = as_to_ld ("-3.4350297408563828055e-001"); 
  yn1[ 8] = as_to_ld ("-1.7964875644231269124e-002"); 
  yn1[ 9] = as_to_ld ("-2.5118910315168408963e-001"); 
  yn1[10] = as_to_ld ("-1.0817311005210173814e-002"); 
  yn1[11] = as_to_ld ("-9.9170891227838064646e-002"); 

  yn1[12] = as_to_ld (" 1.8597282203783061216e-002");  /* Venus */
  yn1[13] = as_to_ld ("-2.6914630075961837156e-001"); 
  yn1[14] = as_to_ld ("-6.5864799228651037249e-003"); 
  yn1[15] = as_to_ld ("-6.1655661757510340864e-001"); 
  yn1[16] = as_to_ld ("-4.1428517490915759141e-003"); 
  yn1[17] = as_to_ld ("-2.6070731202257053265e-001"); 

  yn1[18] = as_to_ld (" 1.3083581197550196016e-002");  /* EMB */
  yn1[19] = as_to_ld (" 6.4260215465460413745e-001"); 
  yn1[20] = as_to_ld (" 9.8846107560980890857e-003"); 
  yn1[21] = as_to_ld ("-7.2083597294516633464e-001"); 
  yn1[22] = as_to_ld (" 4.2862797207098424858e-003"); 
  yn1[23] = as_to_ld ("-3.1263839868006851566e-001"); 

  yn1[24] = as_to_ld ("-1.0289202815795734684e-002");  /* Mars */
  yn1[25] = as_to_ld ("-1.0394941428468679346e+000"); 
  yn1[26] = as_to_ld ("-7.0723334989418045152e-003"); 
  yn1[27] = as_to_ld (" 1.1522403836312279294e+000"); 
  yn1[28] = as_to_ld ("-2.9703138309275424941e-003"); 
  yn1[29] = as_to_ld (" 5.5677401157133424433e-001"); 

  yn1[30] = as_to_ld (" 4.6228653666202124358e-003");   /* Jupiter */
  yn1[31] = as_to_ld ("-4.2373615722355993645e+000"); 
  yn1[32] = as_to_ld ("-5.0555918692615562352e-003"); 
  yn1[33] = as_to_ld ("-3.1464733091277224306e+000"); 
  yn1[34] = as_to_ld ("-2.2818455722493645242e-003"); 
  yn1[35] = as_to_ld ("-1.2462241659031373286e+000"); 

  yn1[36] = as_to_ld ("-4.2556740044874090445e-003");   /* Saturn */
  yn1[37] = as_to_ld (" 6.4633979534971069954e+000"); 
  yn1[38] = as_to_ld (" 3.5633047172379600052e-003"); 
  yn1[39] = as_to_ld (" 6.1037564142437140792e+000"); 
  yn1[40] = as_to_ld (" 1.6573768541155527860e-003"); 
  yn1[41] = as_to_ld (" 2.2444618449141959815e+000"); 

  yn1[42] = as_to_ld (" 5.2274776790945661154e-004");  /* Uranus */
  yn1[43] = as_to_ld ("-1.8139340904007687965e+001"); 
  yn1[44] = as_to_ld ("-3.7317457890160985748e-003"); 
  yn1[45] = as_to_ld ("-2.4578888288247007605e+000"); 
  yn1[46] = as_to_ld ("-1.6425078479167742140e-003"); 
  yn1[47] = as_to_ld ("-8.2063943185823257564e-001");
 
  yn1[48] = as_to_ld (" 2.6898588576179062356e-003");  /* Neptune */
  yn1[49] = as_to_ld ("-1.5304564931392374665e+001"); 
  yn1[50] = as_to_ld ("-1.4254385368896515149e-003"); 
  yn1[51] = as_to_ld ("-2.4351618762495032811e+001"); 
  yn1[52] = as_to_ld ("-6.5161737589378315099e-004"); 
  yn1[53] = as_to_ld ("-9.5901150170962550989e+000");
 
  yn1[54] = as_to_ld (" 3.9518201228036738355e-004");  /* Pluto */
  yn1[55] = as_to_ld ("-3.0312344985766555864e+001"); 
  yn1[56] = as_to_ld ("-3.1426251869952901448e-003"); 
  yn1[57] = as_to_ld ("-1.7898553889579568107e+000"); 
  yn1[58] = as_to_ld ("-1.1124350707669052540e-003"); 
  yn1[59] = as_to_ld (" 8.6212535443336042115e+000"); 

  yn1[60] = as_to_ld ("-4.58408893102436927e-4");  /* MOON */
  yn1[61] = as_to_ld ("-1.60130208907259920e-3"); 
  yn1[62] = as_to_ld ("-2.78697065392994499e-4"); 
  yn1[63] = as_to_ld (" 1.95720839700677492e-3"); 
  yn1[64] = as_to_ld ("-1.70697259101564644e-4"); 
  yn1[65] = as_to_ld (" 9.68483037769266028e-4"); 

  yn1[66] = as_to_ld ("-3.4432210186490844452e-006");  /* Sun */
  yn1[67] = as_to_ld (" 3.7798301065137354516e-003"); 
  yn1[68] = as_to_ld ("  4.0360524289263751236e-006"); 
  yn1[69] = as_to_ld (" 2.6305582306528893588e-003"); 
  yn1[70] = as_to_ld (" 1.8100237472517860139e-006"); 
  yn1[71] = as_to_ld (" 1.0818800652186113487e-003"); 


  return;
}
//------------------------------------------------------------------------------
void
set_const_by_hand ()
{


  printf ( "\n" );
  printf ( "Reading in physical parameters...\n" );
  printf ( "\n" );


  set_planets_0 ();

  set_planets_1 ();


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  C     = as_to_ld ("173.144632720536344565")  ; // Speed of light, au/d ;2.99792458e5 CLIGHT
  EMRAT = as_to_ld ("8.13005869999999999376E1"); // Earth's mass divided by Moon's mass 
  KG    = as_to_ld ("0.01720209895")           ; // Gaussian gravitational constant 

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  GMs[ 0] = as_to_ld ("0.0");      // place holder for libration 
  GMs[ 1] = as_to_ld ("4.91254745145081175785E-11"); // Mercury 

  GMs[ 2] = as_to_ld ("7.24345620963276523095E-10"); // Venus 
  GMs[ 3] = as_to_ld ("8.88769273403302327042E-10"); // Earth 

  GMs[ 4] = as_to_ld ("9.54952894222405763492E-11"); // Mars 
  GMs[ 5] = as_to_ld ("2.82534210344592625472E-7");  // Jupiter 
  GMs[ 6] = as_to_ld ("8.45946850483065929285E-8");  // Saturn 
  GMs[ 7] = as_to_ld ("1.28881623813803488851E-8");  // Uranus 
  GMs[ 8] = as_to_ld ("1.53211248128427618918E-8");  // Neptune 
  GMs[ 9] = as_to_ld ("2.27624775186369921644E-12"); // Pluto 
  GMs[10] = as_to_ld ("1.09318924524284471369E-11"); // Moon 
  GMs[11] = as_to_ld ("2.95912208285591102582E-4");  // Sun 
  GMs[12] = as_to_ld ("1.746e-13");                  // Ceres 
  GMs[13] = as_to_ld ("3.2e-14");                    // Pallas 
  GMs[14] = as_to_ld ("4.08e-14");                   // Vesta 
  GMs[15] = as_to_ld ("1.6e-15");                    // Iris 
  GMs[16] = as_to_ld ("2.6e-15");                    // Bamberga 


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  AU    =  as_to_ld ("1.4959787066e11");  /* Astronomical unit in meters */
  RADS  =  as_to_ld ("6.96e5");  /* Radius of the Sun, kilometers */
  RADM  =  as_to_ld ("1738.0");  /* Radius of the Moon (AM is the value actually used) */
  RADE  =  as_to_ld ("6378.14");  /* Radius of the Earth (AE is the value actually used) */
  AE    =  as_to_ld ("4.26352325064817808471E-5");  /* equatorial radius of Earth, in au */

  Je[0] =  as_to_ld ("1.08263e-3");  /* J2E */
  Je[1] =  as_to_ld ("-2.54e-6");    /* J3E */
  Je[2] =  as_to_ld ("-1.61e-6");    /* J4E */

  AM = as_to_ld ("1.16178124180819699985E-5");  /* equatorial radius of Moon, in au */

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Jm[0] = as_to_ld ("2.02150907893e-4");   /* J2M */
  Jm[1] = as_to_ld ("1.21260448837e-5");   /* J3M */
  Jm[2] = as_to_ld ("-1.45383007072e-7");  /* J4M */

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Cnm[0] = as_to_ld ("2.230351309e-5");     /* C22 */
  Cnm[1] = as_to_ld ("3.07082741328e-5");   /* C31 */
  Cnm[2] = as_to_ld ("4.88840471683e-6");   /* C32 */
  Cnm[3] = as_to_ld ("1.43603108489e-6");   /* C33 */
  Cnm[4] = as_to_ld ("-7.17780149806e-6");  /* C41 */
  Cnm[5] = as_to_ld ("-1.43951838385e-6");  /* C42 */
  Cnm[6] = as_to_ld ("-8.54788154819e-8");  /* C43 */
  Cnm[7] = as_to_ld ("-1.5490389313e-7");   /* C44 */
  

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Snm[0] = as_to_ld ("0.0");               // S22 
  Snm[1] = as_to_ld ("5.61066891941e-6");  // S31 
  Snm[2] = as_to_ld ("1.68743052295e-6");  // S32 
  Snm[3] = as_to_ld ("-3.343544677e-7");   // S33 
  
  Snm[4] = as_to_ld ("2.94743374914e-6");  // S41 
  Snm[5] = as_to_ld ("-2.8843721272e-6");  // S42 
  Snm[6] = as_to_ld ("-7.88967312839e-7"); // S43 
  Snm[7] = as_to_ld (" 5.6404155572e-8");  // S44 


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  K2M    = as_to_ld ("2.22160557154e-2"); 
  LOVENO = as_to_ld ("0.30"); 
  PHASE  = as_to_ld ("4.0700012e-2"); 

  //; synodic period of the moon

  PSLP1  = as_to_ld ("2.29971502189818919e-1");        // DE118 initial psidot: 
  PSLPI  = as_to_ld ("2.732158222801637949485e1");     // = 2 pi / PSLP1 = period in days 
  PSLPIA = as_to_ld ("2.73215820789337158203125e1");   // PSLPA + PSLPB = PSLPI 
  PSLPIB = as_to_ld ("1.490826636745416323896386e-7"); //
  JDEPOCH= as_to_ld ("2440400.5");      //
  LBET = as_to_ld ("6.31686773468e-4"); //  (C-A)/B 
  LGAM = as_to_ld ("2.28022183594e-4"); //  (B-A)/C 

  //; LGAMBET = (LGAMd - LBETd)/(1.0-LBETd*LGAMd);
  //; CMR2 = CMRSQ*AMd*AMd;
  //; BMR2 = CMRSQ*AMd*AMd*(1.0 + LGAMd)/(1.0 + LBETd);
  //; AMR2 = CMRSQ*AMd*AMd*(1.0 - LBETd * LGAMd)/(1.0 + LBETd);

  LGAMBET = as_to_ld ("-4.03664648017289733947E-4");      //
  CMRSQ   = as_to_ld ("0.3906895261319410091");           //
  CMR2    = as_to_ld ("5.2732758299330413853212309E-11"); //
  BMR2    = as_to_ld ("5.2711485389894113965222358E-11"); //
  AMR2    = as_to_ld ("5.2699461151199766018387453E-11"); //

  return;
}
//------------------------------------------------------------------------------
rdnums ()
{

  char str[128];

#if LDOUBLE
  long double x;
  long double *px;
#else
  double x;
  double *px;
#endif

  struct number *pnum;
  int i, n;

  //printf ( "\n" );
  //printf ( "Reading in physical parameters...\n" );
  //printf ( "\n" );

  in = fopen ( "aconst.h", "r" );
  if( in == NULL )
  {
    printf( "can't find aconst.h\n" );
    exit (0);
  }

  pnum = nums; // ��������� �� ������ ������� �������� �������� ������

inloop:

  px = pnum->val; // ��������� �� ������, ���� ����� ��������� ������
  n  = pnum->n;   // ������� ������ ����� ���������

  fprintf (stderr, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ n = %d \n", n);

  if (n <= 0)
    goto done;

  for (i = 0; i < n; i++ )  // ���������, ��������� �����
  {

  nxtline:
    if (fgets (str, 80, in) == NULL)
    {
      printf ("aconst.h file i/o error\n");
      exit(0);
    }

    if( str[0] == ';' ) // ��� �����������
      goto nxtline;

    x = as_to_ld (str); 


    // �������� ��������
    // 
    fprintf (stderr, "%f ", (double) x);


    if( *px != x )
    {

#if LDOUBLE
      funion.ld = *px;
      e64toasc( funion.s, str, 25 );
      printf( "%s changed to ", str );
      funion.ld = x;
      e64toasc( funion.s, str, 25 );
      printf( "%s\n", str );
#else
      funion.d = *px;
      e53toasc( funion.s, str, 18 );
      printf( "%s changed to ", str );
      funion.d = x;
      e53toasc( funion.s, str, 18 );
      printf( "%s\n", str );
#endif

      *px = x;
    }

    ++px; // 
  }

  ++pnum; // ������ �� ����. ������� � ��������

  fprintf (stderr, "\n");

  goto inloop;

done:

  fclose (in);

  
  // ������ �����������:
  //
  //set_const_by_hand ();

  return;
}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
