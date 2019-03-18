#include<stdio.h>


read_buffer( int , int *, int );
write_buffer(int , int *, int );


main( int argc, char *argv[])
{
   int i, fd;
   int *Buf_in, *Buf_out;
   int *buffer;
   int nb_elmt = 100;

   if(argc != 2)
      {
       fprintf(stderr,"commande file_name\n");
       exit(1);
       }
   
   if((fd = open(argv[1],O_RDWR | O_CREAT,0744)) == -1)
      {
      fprintf(stderr,"Can't open %s\n",argv[1]);
      exit(1);
      }

   if((buffer = malloc(1024 * 2)) == NULL)
      {
       fprintf(stderr,"can't allocate memory for buffer\n");
       exit(1);
       }
   
   Buf_in = read_buffer(fd, buffer, nb_elmt);

   for(i=0;i<100;i++)
      fprintf(stdout," VAL = %d\n",buffer[i]);

   Buf_out = write_buffer(fd,buffer, nb_elmt);

   }
   
