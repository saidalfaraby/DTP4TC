#ifndef BISON_MDPPARSER_TAB_H
# define BISON_MDPPARSER_TAB_H

#ifndef YYSTYPE
typedef union {
  double val;
  int ival;
  DdNode * dnode;
  DdNode ** p_dnode;
  double *darray;
  char *nme;
} yystype;
# define YYSTYPE yystype
# define YYSTYPE_IS_TRIVIAL 1
#endif
# define	NAME	257
# define	REAL	258
# define	INTGR	259
# define	OPP	260
# define	CLP	261
# define	VARIABLES	262
# define	DISCOUNT	263
# define	TOLERANCE	264
# define	REWARD	265
# define	VALUE	266
# define	ACTION	267
# define	ENDACTION	268
# define	OBSERVATIONS	269
# define	OBSERVE	270
# define	ENDOBSERVE	271
# define	BELIEF	272
# define	ENDBELIEF	273
# define	DISJ	274
# define	CONJ	275
# define	OSB	276
# define	CSB	277
# define	HORIZON	278
# define	VAL	279
# define	COST	280
# define	PRIME	281
# define	UNNORM	282
# define	STARTDD	283
# define	STARTOBSDD	284
# define	ENDDD	285


extern YYSTYPE yylval;

#endif /* not BISON_MDPPARSER_TAB_H */
