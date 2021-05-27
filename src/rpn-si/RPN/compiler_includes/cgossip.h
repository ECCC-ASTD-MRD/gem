/*
 * Generated automatically by fh2h.pl
 * !!! DO NOT EDIT !!!
 * Edit the original fortran header file instead
 * or fix fh2h.pl if there is a translation bug.
 */


#ifndef FH2H_FGOSSIP_H
#define FH2H_FGOSSIP_H


#ifdef __c
extern "C" {
#endif


#ifndef IMPLICIT
#define IMPLICIT  /* Only to point out implicit types */
#endif

/*------ fortran header (without commons and data statements) ----------*/

#define INIT_ERROR (-1)
#define SERVER_ERROR (-2)
#define CONNECTION_ERROR (-3)
#define READ_ERROR (-4)
#define WRITE_ERROR (-5)
#define READ_TIMEOUT (-6)
#define WRITE_TIMEOUT (-7)
#define READ_TYPE_ERROR (-8)
#define WRITE_TYPE_ERROR (-9)
#define DATA_LENGTH_ERROR (-10)
#define SEND_COMMAND_ERROR (-11)

/*------ common blocks -------------------------------------------------*/

/*------ data statements -----------------------------------------------*/

#ifndef NO_STATIC_DATA


#endif  /* #ifndef NO_STATIC_DATA */

/*------ end of fortran header -----------------------------------------*/

#ifdef __c
}
#endif


#endif  /* #ifndef FH2H_FGOSSIP_H */
