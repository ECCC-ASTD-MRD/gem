#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

int mkstemp_(char *template) {
   int fd = mkstemp(template);

   // Close the file for now since it will be opened again later
   close(fd);
   return fd;
}
