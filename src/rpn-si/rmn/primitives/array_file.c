#include <rpnmacros.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <malloc.h>
/*
 contents of test file test_array_file.f90 follow

 program test_array_file
 integer array(100)
 call array_from_file(array,size(array),'test_file') ! read test_file
 print *,array(1),array(2)
 call array_to_file(array,size(array),'test_file2')  ! write test_file2
 stop
 end
*/
/*
  FORTRAN callable subroutine array_from_file, read a file into a FORTRAN integer array

  call array_from_file(array,nw,file_name)
  integer array(nw)
  character *(*) file_name

  upon return
  array(1) contains the index of the last integer used in array
  array(2) contains the number of characters read
  array(3:array(1)) contain the data from the file

  if the file is too big to fit into array, array(1) = nw-2 upon return
*/
void f77name(array_from_file)(ftnword *array, ftnword *nw, char *file_name, F2Cl nc_file_name)
{
  char *buffer;
  int fd;

  buffer=(char *)malloc(nc_file_name+1);        /* allocate buffer for file name */
  strncpy(buffer,file_name,nc_file_name);       /* copy file name into buffer */
  buffer[nc_file_name]='\0';
  nc_file_name--;                               /* and eliminate trailing blanks */
  while ( buffer[nc_file_name]==' ' && nc_file_name>0 ) {
    buffer[nc_file_name]='\0';
    nc_file_name--;
    }
  fd=open(buffer,O_RDONLY);                     /* open file read only */
  array[1]=read(fd,(char *)(array+2),(*nw-2)*sizeof(ftnword));  /* read data from file up to array capacity */
  array[0]=2 + ((array[1]+sizeof(ftnword)-1))/sizeof(ftnword);  /* set las index used in array */
  close(fd);
}

/*
  FORTRAN callable subroutine array_to_file, write an "file in FORTRAN integer array" to a file
  this file should have been previously read by array_from_file

  call array_to_file(array,nw,file_name)
  integer array(nw)
  character *(*) file_name

  array(1) contains the index of the last integer used in array
  array(2) contains the number of characters to write into file
  array(3:array(1)) contain the data to be written into the file

  if the file is too big to fit into array, array(1) = nw-2
*/
void f77name(array_to_file)(ftnword *array, ftnword *nw, char *file_name, F2Cl nc_file_name)
{
  char *buffer;
  int fd;
  int nc;
  off_t len;

  buffer=(char *)malloc(nc_file_name+1);        /* allocate buffer for file name */
  strncpy(buffer,file_name,nc_file_name);       /* copy file name into buffer */
  buffer[nc_file_name]='\0';
  nc_file_name--;                               /* and eliminate trailing blanks */
  while ( buffer[nc_file_name]==' ' && nc_file_name>0 ) {
    buffer[nc_file_name]='\0';
    nc_file_name--;
    }
  fd=open(buffer,O_CREAT+O_RDWR,0777);          /* open file for writing */
  len=write(fd,(char *)(array+2),array[1]);     /* write data into file */
  ftruncate(fd,len);                            /* make sure to truncate after write */
  close(fd);                                    /* close file */
}
