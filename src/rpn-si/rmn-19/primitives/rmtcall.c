/* RMNLIB - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
 *                          Environnement Canada
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#if !defined (_FLOAT1)
#include <rpnmacros.h>
/*                                                                   */
/*  RMTCALL: effectuer un appel a une fonction FORTRAN dont on a     */
/*           obtenu l'adresse au moyen de LOCF en passant une        */
/*           liste d'adresses pour les parametres (40 adresses       */
/*           DANS LE FORMAT UTILISE PAR LOCF).                       */
/*           ex : STATUT = RMTCALL(ADRESSE,LISTE)                    */
/*                                                                   */
/*           INTEGER LISTE(40), IADRF, STATUS                        */
/*           EXTERNAL AMOI                                           */
/*           INTEGER AMOI, RMTCALL                                   */
/*             .......                                               */
/*           LISTE(1) = LOCF(.....)                                  */
/*             .......                                               */
/*           LISTE(40) = LOCF(.....)                                 */
/*           IADRF = LOCF(AMOI)                                      */
/*           STATUT = RMTCALL(IADRF,LISTE)                           */
/*                                                                   */
/*  notes: toutes ces fonctions travaillent avec des MOTS.           */
/*         Pour FORTRAN, un INTEGER ou un REAL occupent un mot.      */
/*         NE PAS UTILISER pour une variable de type CHARACTER.      */
wordint f77name(rmtcall)(unsigned long long *entry_in,unsigned long long *args_in)
{
  typedef wordint *W_ptr;
  W_ptr args[41];
  int i;
  union {
    long long ptr_sub;
    wordint (* entry)();
    } callee;
  union {
    long long ptr;
    wordint * entry;
    } arg_temp;
 
  callee.ptr_sub = *entry_in;  
  for (i=0; i<41; i++) {
    arg_temp.ptr = args_in[i];
    args[i] = arg_temp.entry;
    }

    
   /* printca(args[0]); */
   return((*callee.entry)(
                 args[ 0],args[ 1],args[ 2],args[ 3],args[ 4],
                 args[ 5],args[ 6],args[ 7],args[ 8],args[ 9],
                 args[10],args[11],args[12],args[13],args[14],
                 args[15],args[16],args[17],args[18],args[19],
                 args[20],args[21],args[22],args[23],args[24],
                 args[25],args[26],args[27],args[28],args[29],
                 args[30],args[31],args[32],args[33],args[34],
                 args[35],args[36],args[37],args[38],args[39],
                 args[40]));


}
#endif
