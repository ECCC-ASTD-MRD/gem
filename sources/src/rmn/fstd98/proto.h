#if !defined(F2Cl)
#define F2Cl int
#endif
int fnom_index(int iun);
int error_msg(char *function_name, int errcode, int errlevel);
int file_index(int iun);
ftnword f77name(xdfopn)(ftnword *fiun,char *mode,ftnword_2 *pri,ftnword *fnpri,
             ftnword_2 *aux,ftnword *fnaux,char *appl,F2Cl l1,F2Cl l2);
int c_xdfopn(int iun,char *mode,word_2 *pri,int npri,
             word_2 *aux,int naux,char *appl);
ftnword f77name(xdfcls)(ftnword *fiun);
int c_xdfcls(int iun);
ftnword f77name(xdfsta)(ftnword *fiun,ftnword *stat,ftnword *fnstat,
			ftnword_2 *pri,ftnword *fnpri,
			ftnword_2 *aux,ftnword *fnaux,
			char *vers,char *appl,F2Cl l1,F2Cl l2);
int c_xdfsta(int iun,word *stat,int nstat,
                    word_2 *pri,int npri,word_2 *aux,int naux,
                    char *vers,char *appl);
ftnword f77name(xdfimp)(ftnword *fiun,ftnword *stat,ftnword *fnstat,
                    ftnword_2 *pri,ftnword_2 *aux,
                    char *vers,char *appl,F2Cl l1,F2Cl l2);
int c_xdfimp(int iun,word *stat,int nstat,word_2 *pri,word_2 *aux,
                    char *vers,char *appl);
ftnword f77name(xdfini)(ftnword *fiun,word *buf,ftnword *fidtyp,
			ftnword *keys,ftnword *fnkeys,
			ftnword *info,ftnword *fninfo);
int c_xdfini(int iun,buffer_interface_ptr buf,int idtyp,
             word *keys,int nkeys,word *info,int ninfo);

void build_burp_prim_keys(burp_record *brpk, word *keys,
                                 burp_record *mask, word *mskkeys,
				 int index, int mode);
void build_burp_info_keys(word *buf, word *keys, int index, int mode);
void build_fstd_info_keys(word *buf, word *keys, int index, int mode);
void build_fstd_prim_keys(word *buf, word *keys, word *mask, word *mskkeys,
				int index, int mode);
ftnword f77name(xdfadd)(word *buf,word *donnees,
                        ftnword *fnelm, ftnword *fnbits, ftnword *fdatyp);
int c_xdfadd(word *buffer, word *donnees, int nelm, int nbits, int datyp);
ftnword f77name(xdfprm)(ftnword *fhandle,ftnword *addr,ftnword *lng,
                        ftnword *idtyp,ftnword *primk,ftnword *fnprim);
int c_xdfprm(int handle,int *addr,int *lng,int *idtyp,
	     word *primk,int nprim);
ftnword f77name(xdfhdr)(word *buf,ftnword *addr,ftnword *lng,
                        ftnword *idtyp,ftnword *primk,ftnword *fnprim,
			ftnword *info, ftnword *fninfo);
int c_xdfhdr(buffer_interface_ptr buf ,int *addr,int *lng,int *idtyp,
	     word *primk,int nprim,word *info,int ninfo);
ftnword f77name(xdfloc)(ftnword *fiun, ftnword *fhandle, ftnword *primk,
			ftnword *fnprim);
int c_xdfloc(int iun, int handle, word *primk,int nprim);
int c_xdfloc2(int iun, int handle, word *primk,int nprim, word *mskkeys);
ftnword f77name(xdfget)(ftnword *fhandle, word *buf);
int c_xdfget(int handle, buffer_interface_ptr buf);
ftnword f77name(xdfput)(ftnword *fiun, ftnword *fhandle,
			word *buf);
int c_xdfput(int iun, int handle, buffer_interface_ptr buf);
ftnword f77name(xdfopt)(char *foptname, char *foptc, ftnword *foptv,
			F2Cl l1, F2Cl l2);
int c_xdfopt(char *optname, char *optc, int optv);
ftnword f77name(xdfgop)(char *foptname, char *foptc, ftnword *foptv,
			F2Cl l1, F2Cl l2);
int c_xdfgop(char *optname, char *optc, int *optv);
ftnword f77name(xdfins)(word *buf,word *donnees,
                        ftnword *fbitpos, ftnword *fnelm,
			ftnword *fnbits, ftnword *fdatyp);
int c_xdfins(word *buffer, word *donnees, int bitpos,
             int nelm, int nbits, int datyp);
ftnword f77name(xdfxtr)(word *buf,word *donnees,
                        ftnword *fbitpos, ftnword *fnelm,
			ftnword *fnbits, ftnword *fdatyp);
int c_xdfxtr(word *buffer, word *donnees, int bitpos,
             int nelm, int nbits, int datyp);
ftnword f77name(xdfrep)(word *buf,word *donnees,
                        ftnword *fbitpos, ftnword *fnelm,
			ftnword *fnbits, ftnword *fdatyp);
int c_xdfrep(word *buffer, word *donnees, int bitpos,
	     int nelm, int nbits, int datyp);
ftnword f77name(xdfcut)(word *buf,
			ftnword *fbitpos, ftnword *fnelm,
			ftnword *fnbits, ftnword *fdatyp);
int c_xdfcut(void *buffer, int bitpos, int nelm, int nbits, int datyp);
ftnword f77name(xdfupd)(ftnword *fiun,word *buf,
			ftnword *fidtyp,
			ftnword *keys,ftnword *fnkeys,
			ftnword *info,ftnword *fninfo);
int c_xdfupd(int iun,buffer_interface_ptr buf,int idtyp,
             word *keys,int nkeys,word *info,int ninfo);
ftnword f77name(xdfuse)(ftnword *fsrc_unit, ftnword *fdest_unit);
int c_xdfuse(int src_unit, int dest_unit);
ftnword f77name(xdfcle)(char *fkeyname,ftnword *fbit1,ftnword *flkey,
			ftnword *ftkey,ftnword *fdesc1,ftnword *fdesc2,F2Cl l1);
int c_xdfcle(char *keyname,int bit1,int lkey,int tkey,int *desc1,int *desc2);
int c_qdfmsig(int iun, char* newappl);
ftnword f77name(qdfmsig)(ftnword *fiun,char *appl,F2Cl l1);
ftnword f77name(mrbdel)(word *buf, ftnword *f_number);
int c_mrbdel(void *buffer, int number);
ftnword f77name(xdflnk)(ftnword *liste, ftnword *fn);
int c_xdflnk(word *liste, int n);
int c_xdfunl(word *liste, int n);
ftnword f77name(fstvoi)(ftnword *f_iun,char *options,F2Cl l1);
int c_fstvoi(int iun,char *options);
ftnword f77name(fstouv)(ftnword *f_iun, char *options, F2Cl l1);
int c_fstouv(int iun, char *options);
ftnword f77name(secateur)(char *filename, ftnword *f_where, F2Cl l1);
void c_fst_env_var(char *cles, int index, char *content);
