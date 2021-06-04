#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <fcntl.h>

#include <rpnmacros.h>

/* ENTRY = ID + TEXT + NewLine
   ID = 8 chars, TEXT = 55 chars */

#define MAXTEXT   90
#define MAXID      9
#define ENTRYSIZE 128

#define READ_PTR 00L
#define WRIT_PTR 20L
#define LIMT_PTR 40L
#define DATA_POS 60L

/* description of the format of the log file :

byte #      line #     contents
00            1         read pointer
20            2         write pointer
40            3         entry size
60            4         message # 1
60+ENTRYSIZE  5         message # 2
 and so on 
*/
