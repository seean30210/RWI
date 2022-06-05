#define  PHYSICS                        MHD
#define  DIMENSIONS                     2
#define  GEOMETRY                       POLAR
#define  BODY_FORCE                     POTENTIAL
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            4

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   EIGHT_WAVES
#define  BACKGROUND_FIELD               NO
#define  AMBIPOLAR_DIFFUSION            NO
#define  RESISTIVITY                    NO
#define  HALL_MHD                       NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  Mstar                          0
#define  Mdisk                          1
#define  Mplanet                        2
#define  Viscosity                      3

/* [Beg] user-defined constants (do not change this line) */

#define  LIMITER                        VANLEER_LIM
#define  UNIT_LENGTH                    (25.0*CONST_au)
#define  UNIT_DENSITY                   (0.000949065*2.2*CONST_Msun/(UNIT_LENGTH*UNIT_LENGTH*UNIT_LENGTH))
#define  UNIT_VELOCITY                  sqrt(CONST_G*g_inputParam[Mstar]*CONST_Msun/UNIT_LENGTH)
#define  UNIT_MAGFIELD                  (sqrt(4.*CONST_PI*UNIT_DENSITY*

/* [End] user-defined constants (do not change this line) */
