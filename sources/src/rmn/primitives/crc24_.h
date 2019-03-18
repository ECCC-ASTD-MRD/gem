/**
 * \file crc24_.h
 * Functions and types for CRC checks.
 *
 * Generated on Fri Feb  8 10:07:03 2013,
 * by pycrc v0.8, http://www.tty1.net/pycrc/
 * using the configuration:
 *    Width        = 24
 *    Poly         = 0x864cfb
 *    XorIn        = 0xb704ce
 *    ReflectIn    = False
 *    XorOut       = 0x000000
 *    ReflectOut   = False
 *    Algorithm    = table-driven
 *****************************************************************************/
#ifndef __CRC24__H__
#define __CRC24__H__

#include <stdlib.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif


/**
 * The definition of the used algorithm.
 *****************************************************************************/
#define CRC_ALGO_TABLE_DRIVEN 1


/**
 * The type of the CRC values.
 *
 * This type must be big enough to contain at least 24 bits.
 *****************************************************************************/
typedef uint32_t crc24_t;


/**
 * Calculate the initial crc value.
 *
 * \return     The initial crc value.
 *****************************************************************************/
static crc24_t crc24_init(void)
{
    return 0xb704ce;
}


/**
 * Update the crc value with new data.
 *
 * \param crc      The current crc value.
 * \param data     Pointer to a buffer of \a data_len bytes.
 * \param data_len Number of bytes in the \a data buffer.
 * \return         The updated crc value.
 *****************************************************************************/
crc24_t crc24_update(crc24_t crc, const unsigned char *data, size_t data_len);
crc24_t crc24_update_le(crc24_t crc, const unsigned char *data, size_t data_len, int mask);

/**
 * Calculate the final crc value.
 *
 * \param crc  The current crc value.
 * \return     The final crc value.
 *****************************************************************************/
static crc24_t crc24_finalize(crc24_t crc)
{
    return crc ^ 0x000000;
}


#ifdef __cplusplus
}           /* closing brace for extern "C" */
#endif

#endif      /* __CRC24__H__ */
