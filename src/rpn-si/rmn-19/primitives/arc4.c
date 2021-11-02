/* 
 * Copyright (c) 2006-2007, Christophe Devine
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer
 *       in the documentation and/or other materials provided with the
 *       distribution.
 *     * Neither the name of the XySSL nor the names of its contributors
 *       may be used to endorse or promote products derived from this
 *       software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
/*
 *  The ARCFOUR algorithm was publicly disclosed on 94/09.
 *
 *  http://groups.google.com/group/sci.crypt/msg/10a300c9d21afca0
 */

#define XYSSL_ARC4_C

#if defined(XYSSL_ARC4_C)

#include "arc4.h"

/*
 * ARC4 key schedule
 */
void arc4_setup( arc4_context *ctx, unsigned char *key, int keylen )
{
    int i, j, k, a;
    unsigned char *m;
    unsigned char scrap[1024];

    ctx->x = 0;
    ctx->y = 0;
    m = ctx->m;

    for( i = 0; i < 256; i++ )
        m[i] = (unsigned char) i;

    j = k = 0;

    for( i = 0; i < 256; i++, k++ )
    {
        if( k >= keylen ) k = 0;

        a = m[i];
        j = ( j + a + key[k] ) & 0xFF;
        m[i] = m[j];
        m[j] = (unsigned char) a;
    }
#ifndef ORIGINAL_TEST
    arc4_crypt(ctx,scrap,1024);  /* waste 1024 key bytes to reduce vulnerability */
#endif
}

/*
 * ARC4 cipher function
 */
void arc4_crypt( arc4_context *ctx, unsigned char *buf, int buflen )
{
    int i, x, y, a, b;
    unsigned char *m;

    x = ctx->x;
    y = ctx->y;
    m = ctx->m;

    for( i = 0; i < buflen; i++ )
    {
        x = ( x + 1 ) & 0xFF; a = m[x];
        y = ( y + a ) & 0xFF; b = m[y];

        m[x] = (unsigned char) b;
        m[y] = (unsigned char) a;

        buf[i] = (unsigned char)
            ( buf[i] ^ m[(unsigned char)( a + b )] );
    }

    ctx->x = x;
    ctx->y = y;
}

#if defined(SELF_TEST)

#include <string.h>
#include <stdio.h>

/*
 * ARC4 tests vectors as posted by Eric Rescorla in sep. 1994:
 * used only with -DORIGINAL_TEST (otherwise we perform a crypt/decrypt test)
 *
 * http://groups.google.com/group/comp.security.misc/msg/10a300c9d21afca0
 */
static const unsigned char arc4_test_key[3][8] =
{
    { 0x01, 0x23, 0x45, 0x67, 0x89, 0xAB, 0xCD, 0xEF },
    { 0x01, 0x23, 0x45, 0x67, 0x89, 0xAB, 0xCD, 0xEF },
    { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 }
};

static const unsigned char arc4_test_pt[3][8] =
{
    { 0x01, 0x23, 0x45, 0x67, 0x89, 0xAB, 0xCD, 0xEF },
    { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 },
    { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 }
};

static const unsigned char arc4_test_ct[3][8] =
{
    { 0x75, 0xB7, 0x87, 0x80, 0x99, 0xE0, 0xC5, 0x96 },
    { 0x74, 0x94, 0xC2, 0xE7, 0x10, 0x4B, 0x08, 0x79 },
    { 0xDE, 0x18, 0x89, 0x41, 0xA3, 0x37, 0x5D, 0x3A }
};

/*
 * Checkup routine
 */
int arc4_self_test( int verbose )
{
    int i;
    unsigned char buf[8];
    unsigned char buf2[8];
    arc4_context ctx;
    arc4_context ctx2;

    for( i = 0; i < 3; i++ )
    {
        if( verbose != 0 )
            printf( "  ARC4 test #%d: \n", i + 1 );

        memcpy( buf, arc4_test_pt[i], 8 );
        arc4_setup( &ctx, (unsigned char *) arc4_test_key[i], 8 );
        arc4_crypt( &ctx, buf, 8 );

        memcpy( buf2,buf,8); /* now we decrypt and expect to get back original byte stream */
        arc4_setup( &ctx2, (unsigned char *) arc4_test_key[i], 8 );
        arc4_crypt( &ctx2, buf2, 8 );

        if( memcmp( buf2, arc4_test_pt[i], 8 ) != 0 )
        {
            if( verbose != 0 )
                printf( "decryption after encryption failed\n" );
        }
	else
	{
            if( verbose != 0 )
                printf( "decryption after encryption OK\n" );
	}
#ifdef ORIGINAL_TEST
        if( memcmp( buf, arc4_test_ct[i], 8 ) != 0 )
        {
            if( verbose != 0 )
                printf( "encryption failed\n" );
        }
	else
	{
            if( verbose != 0 )
                printf( "encryption OK\n" );
	}

        if( verbose != 0 )
            printf( "passed\n" );
#endif
    }
    if( verbose != 0 )
        printf( "\n" );

    return( 0 );
}

main()
{
arc4_self_test(1);
}
#endif

#endif
