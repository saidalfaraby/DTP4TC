/* A Bison parser, made from mdpparser.y
   by GNU bison 1.35.  */

#define YYBISON 1  /* Identify Bison output.  */

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

#line 4 "mdpparser.y"

  #include "POMDP.h"

  //#define YYSTYPE double
  #define MAXDDS 100000
  //#define REDUCE_PRECISION 1
  //#define DONOTNORMALIZE 1
  //#define PARSERDEBUG 1
  int yylex( void );
  int yyerror( const char *s);
  void error(const char *s, const char *v);
  int curr_primed_ovar, curr_ovar, curr_oval, curr_obs_oval;
  DdNode *dds[MAXDDS];
  char *ddnames[MAXDDS];
  DdNode **goodState;
  DdNode **goodStatep;
  double cpt[MAXVALS];
  double pvsum;
  int valindex, level, pvcount, conj, pvindex, numdds, startndds;
  int oldconj[64], conjlevel;
  int numobsparams;
  double obsparams[128];
  bool doing_reward, doing_dd, doing_obs, doing_action, unnormalized;
  int levelIndex[MAXVARS];
  // tells if a level is primed variable
  bool levelPrime[MAXVARS];
  int branchCount[MAXVARS];
  extern POMDP *__theMDP;
  extern POMDP *__thePOMDP;
  map<const char*, int, ltstr> ddindices;

#line 35 "mdpparser.y"
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
#ifndef YYDEBUG
# define YYDEBUG 0
#endif



#define	YYFINAL		147
#define	YYFLAG		-32768
#define	YYNTBASE	32

/* YYTRANSLATE(YYLEX) -- Bison token number corresponding to YYLEX. */
#define YYTRANSLATE(x) ((unsigned)(x) <= 285 ? yytranslate[x] : 93)

/* YYTRANSLATE[YYLEX] -- Bison token number corresponding to YYLEX. */
static const char yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     3,     4,     5,
       6,     7,     8,     9,    10,    11,    12,    13,    14,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31
};

#if YYDEBUG
static const short yyprhs[] =
{
       0,     0,     2,     4,     5,     8,    11,    13,    22,    27,
      29,    32,    37,    39,    42,    45,    46,    47,    53,    55,
      58,    63,    64,    65,    70,    73,    75,    76,    80,    81,
      83,    86,    88,    90,    92,    93,    98,   101,   105,   112,
     115,   116,   118,   119,   120,   124,   125,   126,   127,   133,
     136,   138,   139,   143,   146,   148,   149,   153,   154,   160,
     161,   167,   168,   174,   176,   180,   181,   187,   189,   192,
     194,   196,   199,   202,   204,   207,   209,   210,   216,   219,
     222,   226,   230,   232,   233,   237,   240,   241,   245,   247,
     249,   252
};
static const short yyrhs[] =
{
      36,     0,    33,     0,     0,    34,    35,     0,    35,    54,
       0,    54,     0,    37,    42,    46,    51,    52,    86,    88,
      90,     0,     6,     8,    38,     7,     0,    39,     0,    38,
      39,     0,     6,    40,    41,     7,     0,     3,     0,    41,
       3,     0,     3,     3,     0,     0,     0,     6,    15,    43,
      44,     7,     0,    45,     0,    44,    45,     0,     6,     3,
       4,     7,     0,     0,     0,    18,    47,    48,    19,     0,
      48,    49,     0,    49,     0,     0,     3,    50,    71,     0,
       0,    28,     0,    52,    53,     0,    53,     0,    57,     0,
      54,     0,     0,    56,    55,    71,    31,     0,    29,     3,
       0,    30,     3,     3,     0,    58,    59,    68,    62,    60,
      14,     0,    13,     3,     0,     0,     4,     0,     0,     0,
      26,    61,    71,     0,     0,     0,     0,    16,    63,    65,
      64,    17,     0,    65,    66,     0,    66,     0,     0,     3,
      67,    71,     0,    68,    69,     0,    69,     0,     0,     3,
      70,    71,     0,     0,    22,    21,    72,    79,    23,     0,
       0,    22,    20,    73,    79,    23,     0,     0,     6,    78,
      74,    80,     7,     0,    84,     0,     6,    83,     7,     0,
       0,     6,    16,    75,    76,     7,     0,    77,     0,    77,
      76,     0,     4,     0,     3,     0,     3,    27,     0,    79,
      71,     0,    71,     0,    80,    81,     0,    81,     0,     0,
       6,     3,    82,    71,     7,     0,    83,    85,     0,    85,
      85,     0,     6,     4,     7,     0,     6,     3,     7,     0,
       4,     0,     0,    11,    87,    71,     0,     9,     4,     0,
       0,     9,    89,    71,     0,    91,     0,    92,     0,    10,
       4,     0,    24,     4,     0
};

#endif

#if YYDEBUG
/* YYRLINE[YYN] -- source line where rule number YYN was defined. */
static const short yyrline[] =
{
       0,    66,    67,    69,    69,    91,    92,    95,   100,   140,
     141,   143,   149,   157,   160,   165,   165,   165,   182,   183,
     185,   195,   196,   196,   208,   209,   212,   212,   241,   243,
     247,   248,   250,   251,   253,   253,   271,   275,   287,   296,
     306,   309,   314,   320,   320,   344,   344,   344,   344,   352,
     353,   355,   355,   383,   384,   386,   386,   472,   472,   490,
     490,   507,   507,   534,   535,   559,   559,   596,   597,   599,
     604,   612,   620,   643,   650,   656,   668,   668,   700,   706,
     715,   730,   751,   787,   787,   805,   811,   811,   818,   819,
     821,   827
};
#endif


#if (YYDEBUG) || defined YYERROR_VERBOSE

/* YYTNAME[TOKEN_NUM] -- String name of the token TOKEN_NUM. */
static const char *const yytname[] =
{
  "$", "error", "$undefined.", "NAME", "REAL", "INTGR", "OPP", "CLP", 
  "VARIABLES", "DISCOUNT", "TOLERANCE", "REWARD", "VALUE", "ACTION", 
  "ENDACTION", "OBSERVATIONS", "OBSERVE", "ENDOBSERVE", "BELIEF", 
  "ENDBELIEF", "DISJ", "CONJ", "OSB", "CSB", "HORIZON", "VAL", "COST", 
  "PRIME", "UNNORM", "STARTDD", "STARTOBSDD", "ENDDD", "input", 
  "alphaorbelief", "@1", "ablist", "mdp", "varilist", "varlist", "vardec", 
  "varname", "vallist", "obslist", "@2", "oblist", "obdec", "ibelief", 
  "@3", "belieflist", "belief", "@4", "unn", "actslist", "actionordd", 
  "thedd", "@5", "ddname", "action", "actionname", "oldactioncost", 
  "actioncost", "@6", "observation", "@7", "@8", "obsfunlist", "obsfun", 
  "@9", "acttreelist", "acttree", "@10", "theadd", "@11", "@12", "@13", 
  "@14", "obsparamlist", "obsparam", "rootnodename", "con_dis_add", 
  "currentadd", "subadd", "@15", "primeadd", "constnode", "constadd", 
  "reward", "@16", "disc", "@17", "tolhor", "tol", "hor", 0
};
#endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives. */
static const short yyr1[] =
{
       0,    32,    32,    34,    33,    35,    35,    36,    37,    38,
      38,    39,    40,    41,    41,    42,    43,    42,    44,    44,
      45,    46,    47,    46,    48,    48,    50,    49,    51,    51,
      52,    52,    53,    53,    55,    54,    56,    56,    57,    58,
      59,    59,    60,    61,    60,    62,    63,    64,    62,    65,
      65,    67,    66,    68,    68,    70,    69,    72,    71,    73,
      71,    74,    71,    71,    71,    75,    71,    76,    76,    77,
      78,    78,    79,    79,    80,    80,    82,    81,    83,    83,
      84,    84,    85,    87,    86,    88,    89,    88,    90,    90,
      91,    92
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN. */
static const short yyr2[] =
{
       0,     1,     1,     0,     2,     2,     1,     8,     4,     1,
       2,     4,     1,     2,     2,     0,     0,     5,     1,     2,
       4,     0,     0,     4,     2,     1,     0,     3,     0,     1,
       2,     1,     1,     1,     0,     4,     2,     3,     6,     2,
       0,     1,     0,     0,     3,     0,     0,     0,     5,     2,
       1,     0,     3,     2,     1,     0,     3,     0,     5,     0,
       5,     0,     5,     1,     3,     0,     5,     1,     2,     1,
       1,     2,     2,     1,     2,     1,     0,     5,     2,     2,
       3,     3,     1,     0,     3,     2,     0,     3,     1,     1,
       2,     2
};

/* YYDEFACT[S] -- default rule to reduce with in state S when YYTABLE
   doesn't specify something else to do.  Zero means the default is an
   error. */
static const short yydefact[] =
{
       3,     0,     2,     0,     1,    15,     0,     0,     0,     4,
       6,    34,     0,    21,     0,     0,     9,    36,     0,     5,
       0,    16,    22,    28,    12,     0,     8,    10,    37,     0,
       0,     0,    63,     0,     0,    29,     0,     0,     0,    70,
      82,    65,    61,     0,     0,    59,    57,    35,     0,     0,
      18,    26,     0,    25,     0,     0,    31,    33,    32,    40,
      14,    13,    11,    81,    71,    80,     0,     0,    82,    64,
      78,    79,     0,     0,     0,    17,    19,     0,    23,    24,
      39,    83,    30,     0,    41,     0,    69,     0,    67,     0,
       0,    75,    73,     0,     0,     0,    27,     0,    86,     0,
      55,    45,    54,    66,    68,    76,    62,    74,    60,    72,
      58,    20,    84,    85,     0,     0,     0,     7,    88,    89,
       0,    46,    42,    53,     0,    87,    90,    91,    56,     0,
      43,     0,     0,    51,    47,    50,     0,    38,    77,     0,
       0,    49,    44,    52,    48,     0,     0,     0
};

static const short yydefgoto[] =
{
     145,     2,     3,     9,     4,     5,    15,    16,    25,    38,
      13,    33,    49,    50,    23,    34,    52,    53,    77,    36,
      55,    56,    57,    20,    11,    58,    59,    85,   131,   136,
     122,   129,   140,   134,   135,   139,   101,   102,   120,    92,
      73,    72,    67,    66,    87,    88,    42,    93,    90,    91,
     124,    43,    32,    44,    83,    97,    99,   114,   117,   118,
     119
};

static const short yypact[] =
{
      -1,     1,-32768,     9,-32768,    24,    38,    47,    48,     9,
  -32768,-32768,    -8,    35,    49,    34,-32768,-32768,    51,-32768,
       6,-32768,-32768,    27,-32768,    53,-32768,-32768,-32768,    29,
      15,    28,-32768,    52,    57,-32768,   -11,    58,     8,    -6,
      55,-32768,-32768,    30,    59,-32768,-32768,-32768,    61,    36,
  -32768,-32768,    10,-32768,    62,    -7,-32768,-32768,-32768,    63,
  -32768,-32768,-32768,-32768,-32768,-32768,    64,    60,-32768,-32768,
  -32768,-32768,     6,     6,    65,-32768,-32768,     6,-32768,-32768,
  -32768,-32768,-32768,    66,-32768,    67,-32768,    69,    64,    68,
      40,-32768,-32768,     2,     4,    71,-32768,     6,    75,     7,
  -32768,     0,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,
  -32768,-32768,-32768,-32768,     6,    76,    77,-32768,-32768,-32768,
       6,-32768,    46,-32768,     6,-32768,-32768,-32768,-32768,    79,
  -32768,    70,    78,-32768,    79,-32768,     6,-32768,-32768,     6,
      72,-32768,-32768,-32768,-32768,    83,    86,-32768
};

static const short yypgoto[] =
{
  -32768,-32768,-32768,-32768,-32768,-32768,-32768,    73,-32768,-32768,
  -32768,-32768,-32768,    41,-32768,-32768,-32768,    39,-32768,-32768,
  -32768,    32,    11,-32768,-32768,-32768,-32768,-32768,-32768,-32768,
  -32768,-32768,-32768,-32768,   -42,-32768,-32768,    -5,-32768,   -20,
  -32768,-32768,-32768,-32768,    13,-32768,-32768,    20,-32768,    12,
  -32768,-32768,-32768,     5,-32768,-32768,-32768,-32768,-32768,-32768,
  -32768
};


#define	YYLAST		119


static const short yytable[] =
{
      31,    63,    54,   100,    81,     1,    54,    21,    29,     6,
      29,    61,    29,    51,    10,    62,   121,   115,     7,     8,
      19,    64,     7,     8,    30,   108,    30,   110,    30,    78,
      12,   116,    39,    40,    68,    45,    46,    69,     7,     8,
      14,    26,    48,    75,    14,    41,    89,   106,    70,    71,
      17,    18,    24,    22,    28,    35,    37,    96,    48,    47,
      51,    60,    65,    68,    74,    80,    89,    84,    86,    95,
     100,   105,   130,   109,   109,    98,   103,   112,   111,   113,
     126,   127,   133,   146,   137,   138,   147,    82,    27,   144,
      76,    79,   141,    94,   125,     0,   123,     0,     0,     0,
     128,   104,   107,     0,   132,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   142,     0,     0,   143
};

static const short yycheck[] =
{
      20,     7,    13,     3,    11,     6,    13,    15,     6,     8,
       6,     3,     6,     3,     3,     7,    16,    10,    29,    30,
       9,    27,    29,    30,    22,    23,    22,    23,    22,    19,
       6,    24,     3,     4,     4,    20,    21,     7,    29,    30,
       6,     7,     6,     7,     6,    16,     6,     7,    43,    44,
       3,     3,     3,    18,     3,    28,     3,    77,     6,    31,
       3,     3,     7,     4,     3,     3,     6,     4,     4,     4,
       3,     3,    26,    93,    94,     9,     7,    97,     7,     4,
       4,     4,     3,     0,    14,     7,     0,    55,    15,    17,
      49,    52,   134,    73,   114,    -1,   101,    -1,    -1,    -1,
     120,    88,    90,    -1,   124,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   136,    -1,    -1,   139
};
/* -*-C-*-  Note some compilers choke on comments on `#line' lines.  */
#line 3 "/usr/share/bison/bison.simple"

/* Skeleton output parser for bison,

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002 Free Software
   Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.  */

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

/* This is the parser code that is written into each bison parser when
   the %semantic_parser declaration is not specified in the grammar.
   It was written by Richard Stallman by simplifying the hairy parser
   used when %semantic_parser is specified.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

#if ! defined (yyoverflow) || defined (YYERROR_VERBOSE)

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# if YYSTACK_USE_ALLOCA
#  define YYSTACK_ALLOC alloca
# else
#  ifndef YYSTACK_USE_ALLOCA
#   if defined (alloca) || defined (_ALLOCA_H)
#    define YYSTACK_ALLOC alloca
#   else
#    ifdef __GNUC__
#     define YYSTACK_ALLOC __builtin_alloca
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning. */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
# else
#  if defined (__STDC__) || defined (__cplusplus)
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   define YYSIZE_T size_t
#  endif
#  define YYSTACK_ALLOC malloc
#  define YYSTACK_FREE free
# endif
#endif /* ! defined (yyoverflow) || defined (YYERROR_VERBOSE) */


#if (! defined (yyoverflow) \
     && (! defined (__cplusplus) \
	 || (YYLTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  short yyss;
  YYSTYPE yyvs;
# if YYLSP_NEEDED
  YYLTYPE yyls;
# endif
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAX (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# if YYLSP_NEEDED
#  define YYSTACK_BYTES(N) \
     ((N) * (sizeof (short) + sizeof (YYSTYPE) + sizeof (YYLTYPE))	\
      + 2 * YYSTACK_GAP_MAX)
# else
#  define YYSTACK_BYTES(N) \
     ((N) * (sizeof (short) + sizeof (YYSTYPE))				\
      + YYSTACK_GAP_MAX)
# endif

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  register YYSIZE_T yyi;		\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (0)
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAX;	\
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (0)

#endif


#if ! defined (YYSIZE_T) && defined (__SIZE_TYPE__)
# define YYSIZE_T __SIZE_TYPE__
#endif
#if ! defined (YYSIZE_T) && defined (size_t)
# define YYSIZE_T size_t
#endif
#if ! defined (YYSIZE_T)
# if defined (__STDC__) || defined (__cplusplus)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# endif
#endif
#if ! defined (YYSIZE_T)
# define YYSIZE_T unsigned int
#endif

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		-2
#define YYEOF		0
#define YYACCEPT	goto yyacceptlab
#define YYABORT 	goto yyabortlab
#define YYERROR		goto yyerrlab1
/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */
#define YYFAIL		goto yyerrlab
#define YYRECOVERING()  (!!yyerrstatus)
#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yychar1 = YYTRANSLATE (yychar);				\
      YYPOPSTACK;						\
      goto yybackup;						\
    }								\
  else								\
    { 								\
      yyerror ("syntax error: cannot back up");			\
      YYERROR;							\
    }								\
while (0)

#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Compute the default location (before the actions
   are run).

   When YYLLOC_DEFAULT is run, CURRENT is set the location of the
   first token.  By default, to implement support for ranges, extend
   its range to the last symbol.  */

#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)       	\
   Current.last_line   = Rhs[N].last_line;	\
   Current.last_column = Rhs[N].last_column;
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#if YYPURE
# if YYLSP_NEEDED
#  ifdef YYLEX_PARAM
#   define YYLEX		yylex (&yylval, &yylloc, YYLEX_PARAM)
#  else
#   define YYLEX		yylex (&yylval, &yylloc)
#  endif
# else /* !YYLSP_NEEDED */
#  ifdef YYLEX_PARAM
#   define YYLEX		yylex (&yylval, YYLEX_PARAM)
#  else
#   define YYLEX		yylex (&yylval)
#  endif
# endif /* !YYLSP_NEEDED */
#else /* !YYPURE */
# define YYLEX			yylex ()
#endif /* !YYPURE */


/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (0)
/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
#endif /* !YYDEBUG */

/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   SIZE_MAX < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#if YYMAXDEPTH == 0
# undef YYMAXDEPTH
#endif

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif

#ifdef YYERROR_VERBOSE

# ifndef yystrlen
#  if defined (__GLIBC__) && defined (_STRING_H)
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
#   if defined (__STDC__) || defined (__cplusplus)
yystrlen (const char *yystr)
#   else
yystrlen (yystr)
     const char *yystr;
#   endif
{
  register const char *yys = yystr;

  while (*yys++ != '\0')
    continue;

  return yys - yystr - 1;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined (__GLIBC__) && defined (_STRING_H) && defined (_GNU_SOURCE)
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
#   if defined (__STDC__) || defined (__cplusplus)
yystpcpy (char *yydest, const char *yysrc)
#   else
yystpcpy (yydest, yysrc)
     char *yydest;
     const char *yysrc;
#   endif
{
  register char *yyd = yydest;
  register const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif
#endif

#line 315 "/usr/share/bison/bison.simple"


/* The user can define YYPARSE_PARAM as the name of an argument to be passed
   into yyparse.  The argument should have type void *.
   It should actually point to an object.
   Grammar actions can access the variable by casting it
   to the proper pointer type.  */

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
#  define YYPARSE_PARAM_ARG void *YYPARSE_PARAM
#  define YYPARSE_PARAM_DECL
# else
#  define YYPARSE_PARAM_ARG YYPARSE_PARAM
#  define YYPARSE_PARAM_DECL void *YYPARSE_PARAM;
# endif
#else /* !YYPARSE_PARAM */
# define YYPARSE_PARAM_ARG
# define YYPARSE_PARAM_DECL
#endif /* !YYPARSE_PARAM */

/* Prevent warning if -Wstrict-prototypes.  */
#ifdef __GNUC__
# ifdef YYPARSE_PARAM
int yyparse (void *);
# else
int yyparse (void);
# endif
#endif

/* YY_DECL_VARIABLES -- depending whether we use a pure parser,
   variables are global, or local to YYPARSE.  */

#define YY_DECL_NON_LSP_VARIABLES			\
/* The lookahead symbol.  */				\
int yychar;						\
							\
/* The semantic value of the lookahead symbol. */	\
YYSTYPE yylval;						\
							\
/* Number of parse errors so far.  */			\
int yynerrs;

#if YYLSP_NEEDED
# define YY_DECL_VARIABLES			\
YY_DECL_NON_LSP_VARIABLES			\
						\
/* Location data for the lookahead symbol.  */	\
YYLTYPE yylloc;
#else
# define YY_DECL_VARIABLES			\
YY_DECL_NON_LSP_VARIABLES
#endif


/* If nonreentrant, generate the variables here. */

#if !YYPURE
YY_DECL_VARIABLES
#endif  /* !YYPURE */

int
yyparse (YYPARSE_PARAM_ARG)
     YYPARSE_PARAM_DECL
{
  /* If reentrant, generate the variables here. */
#if YYPURE
  YY_DECL_VARIABLES
#endif  /* !YYPURE */

  register int yystate;
  register int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Lookahead token as an internal (translated) token number.  */
  int yychar1 = 0;

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack. */
  short	yyssa[YYINITDEPTH];
  short *yyss = yyssa;
  register short *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  register YYSTYPE *yyvsp;

#if YYLSP_NEEDED
  /* The location stack.  */
  YYLTYPE yylsa[YYINITDEPTH];
  YYLTYPE *yyls = yylsa;
  YYLTYPE *yylsp;
#endif

#if YYLSP_NEEDED
# define YYPOPSTACK   (yyvsp--, yyssp--, yylsp--)
#else
# define YYPOPSTACK   (yyvsp--, yyssp--)
#endif

  YYSIZE_T yystacksize = YYINITDEPTH;


  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;
#if YYLSP_NEEDED
  YYLTYPE yyloc;
#endif

  /* When reducing, the number of symbols on the RHS of the reduced
     rule. */
  int yylen;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;
#if YYLSP_NEEDED
  yylsp = yyls;
#endif
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed. so pushing a state here evens the stacks.
     */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyssp >= yyss + yystacksize - 1)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack. Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	short *yyss1 = yyss;

	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  */
# if YYLSP_NEEDED
	YYLTYPE *yyls1 = yyls;
	/* This used to be a conditional around just the two extra args,
	   but that might be undefined if yyoverflow is a macro.  */
	yyoverflow ("parser stack overflow",
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),
		    &yyls1, yysize * sizeof (*yylsp),
		    &yystacksize);
	yyls = yyls1;
# else
	yyoverflow ("parser stack overflow",
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),
		    &yystacksize);
# endif
	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyoverflowlab;
# else
      /* Extend the stack our own way.  */
      if (yystacksize >= YYMAXDEPTH)
	goto yyoverflowlab;
      yystacksize *= 2;
      if (yystacksize > YYMAXDEPTH)
	yystacksize = YYMAXDEPTH;

      {
	short *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyoverflowlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);
# if YYLSP_NEEDED
	YYSTACK_RELOCATE (yyls);
# endif
# undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;
#if YYLSP_NEEDED
      yylsp = yyls + yysize - 1;
#endif

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyssp >= yyss + yystacksize - 1)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;


/*-----------.
| yybackup.  |
`-----------*/
yybackup:

/* Do appropriate processing given the current state.  */
/* Read a lookahead token if we need one and don't already have one.  */
/* yyresume: */

  /* First try to decide what to do without reference to lookahead token.  */

  yyn = yypact[yystate];
  if (yyn == YYFLAG)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* yychar is either YYEMPTY or YYEOF
     or a valid token in external form.  */

  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  /* Convert token to internal form (in yychar1) for indexing tables with */

  if (yychar <= 0)		/* This means end of input. */
    {
      yychar1 = 0;
      yychar = YYEOF;		/* Don't call YYLEX any more */

      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yychar1 = YYTRANSLATE (yychar);

#if YYDEBUG
     /* We have to keep this `#if YYDEBUG', since we use variables
	which are defined only if `YYDEBUG' is set.  */
      if (yydebug)
	{
	  YYFPRINTF (stderr, "Next token is %d (%s",
		     yychar, yytname[yychar1]);
	  /* Give the individual parser a way to print the precise
	     meaning of a token, for further debugging info.  */
# ifdef YYPRINT
	  YYPRINT (stderr, yychar, yylval);
# endif
	  YYFPRINTF (stderr, ")\n");
	}
#endif
    }

  yyn += yychar1;
  if (yyn < 0 || yyn > YYLAST || yycheck[yyn] != yychar1)
    goto yydefault;

  yyn = yytable[yyn];

  /* yyn is what to do for this token type in this state.
     Negative => reduce, -yyn is rule number.
     Positive => shift, yyn is new state.
       New state is final state => don't bother to shift,
       just return success.
     0, or most negative number => error.  */

  if (yyn < 0)
    {
      if (yyn == YYFLAG)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }
  else if (yyn == 0)
    goto yyerrlab;

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Shift the lookahead token.  */
  YYDPRINTF ((stderr, "Shifting token %d (%s), ",
	      yychar, yytname[yychar1]));

  /* Discard the token being shifted unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  *++yyvsp = yylval;
#if YYLSP_NEEDED
  *++yylsp = yylloc;
#endif

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  yystate = yyn;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to the semantic value of
     the lookahead token.  This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];

#if YYLSP_NEEDED
  /* Similarly for the default location.  Let the user run additional
     commands if for instance locations are ranges.  */
  yyloc = yylsp[1-yylen];
  YYLLOC_DEFAULT (yyloc, (yylsp - yylen), yylen);
#endif

#if YYDEBUG
  /* We have to keep this `#if YYDEBUG', since we use variables which
     are defined only if `YYDEBUG' is set.  */
  if (yydebug)
    {
      int yyi;

      YYFPRINTF (stderr, "Reducing via rule %d (line %d), ",
		 yyn, yyrline[yyn]);

      /* Print the symbols being reduced, and their result.  */
      for (yyi = yyprhs[yyn]; yyrhs[yyi] > 0; yyi++)
	YYFPRINTF (stderr, "%s ", yytname[yyrhs[yyi]]);
      YYFPRINTF (stderr, " -> %s\n", yytname[yyr1[yyn]]);
    }
#endif

  switch (yyn) {

case 3:
#line 69 "mdpparser.y"
{
  startndds = numdds;
  doing_reward = true;
;
    break;}
case 4:
#line 72 "mdpparser.y"
{
  // look through the list of dds and find those whose names start with OPTALPHA - these are the alpha vectors
  // or those that start wtih BELIEFSAMPLE - this is a belief sample
  for (int i=startndds; i<numdds; i++) {
    if (strncmp(ddnames[i],"OPTALPHA",8) == 0) {
      //dds[i] is an alpha vector
      __thePOMDP->alphas[__thePOMDP->numalphas] = dds[i];
      Cudd_Ref(__thePOMDP->alphas[__thePOMDP->numalphas]);
      __thePOMDP->numalphas++;
    }
    if (strncmp(ddnames[i],"BELIEFSAMPLE",12) == 0) {
      //dds[i] is an belief sample
      __thePOMDP->ufbeliefs[__thePOMDP->numbeliefs] = dds[i];
      Cudd_Ref(__thePOMDP->ufbeliefs[__thePOMDP->numbeliefs]);
      __thePOMDP->numbeliefs++;
    }
  }
  doing_reward = false;
;
    break;}
case 7:
#line 96 "mdpparser.y"
{
  ddindices.clear();
;
    break;}
case 8:
#line 101 "mdpparser.y"
{
  // now we have the variables, so allocate for the NewPrime Diagrams
  // These are diagrams with a single multi-valued primed variable in them,
  // which is translated into log2(nvals) binary variables
  // so we only need numorigvars NewPrime diagrams for each action
  int i;
  __theMDP->NewPrime = (DdNode ***)malloc(MAXACT*(sizeof(DdNode **)));
  for(i=0;i<MAXACT;i++)
    __theMDP->NewPrime[i] = (DdNode **)malloc(__theMDP->numorigvars*(sizeof(DdNode*)));
  // build the good-state ADDs
  buildGoodStateADDs(__theMDP->orig_vars,__theMDP->vars,__theMDP->prime_vars,__theMDP->numorigvars);
  doing_reward = false;
  doing_action = false;
  doing_obs = false;
  doing_dd = false;
  pvindex = 0;
  pvsum = 0.0;

  numdds = 0;
  // we also want to create a bunch of 'generic' dds like ones for variables that
  // stay the same, which will be named SAME<varname>
  for (i=0; i<__theMDP->numorigvars; i++) {
#ifdef PARSERDEBUG
    fprintf(stderr,"************************** about to parse the dd for %d *****\n",i);
#endif
    ddnames[numdds] = new char[128];
    sprintf(ddnames[numdds],"SAME%s",__theMDP->orig_vars[i].name);
    dds[numdds] = __theMDP->buildSameDD(i);
    Cudd_Ref(dds[numdds]);
    // store in hash table
    ddindices[ddnames[numdds]] = numdds;
#ifdef PARSERDEBUG
    fprintf(stderr,"************************** parsed the dd %s at %d*****\n",ddnames[numdds],ddindices[ddnames[numdds]]);
    Cudd_PrintDebug(gbm,dds[numdds],4,100);
#endif
    numdds++;
  }
;
    break;}
case 11:
#line 144 "mdpparser.y"
{
  __theMDP->newADDVar();
  __theMDP->numorigvars++;
;
    break;}
case 12:
#line 150 "mdpparser.y"
{ 
  __theMDP->orig_vars[__theMDP->numorigvars].name = strdup(yyvsp[0].nme);
#ifdef PARSERDEBUG
  fprintf(stderr,"name of %d orig_var %s\n",__theMDP->numorigvars,__theMDP->orig_vars[__theMDP->numorigvars].name);
#endif
;
    break;}
case 13:
#line 158 "mdpparser.y"
{ 
  addVar(__theMDP->orig_vars,__theMDP->numorigvars,yyvsp[0].nme);
;
    break;}
case 14:
#line 160 "mdpparser.y"
{ 
  addVar(__theMDP->orig_vars,__theMDP->numorigvars,yyvsp[-1].nme);
  addVar(__theMDP->orig_vars,__theMDP->numorigvars,yyvsp[0].nme);
;
    break;}
case 16:
#line 166 "mdpparser.y"
{
  __theMDP->numorigobs = 0;
;
    break;}
case 17:
#line 168 "mdpparser.y"
{
  // now we have the observations, so allocate for the ObsFun Diagrams
  // These are diagrams with a single multi-valued primed variable in them,
  // which is translated into log2(nvals) binary variables
  // so we only need numorigvars NewPrime diagrams for each action
  __theMDP->totalObsFun = new DdNode **[MAXACT];
  for (int i=0; i<MAXACT; i++) 
    __theMDP->totalObsFun[i] = new DdNode *[__theMDP->numorigobs];
#ifdef PARSERDEBUG
  fprintf(stderr,"allocating for %d original observations in MDP\n", __theMDP->numorigobs);
#endif
  // also allocate for the initBeliefState
  __theMDP->initBeliefState = new DdNode *[__theMDP->numorigvars];
;
    break;}
case 20:
#line 186 "mdpparser.y"
{
  __theMDP->orig_obs[__theMDP->numorigobs].name = strdup(yyvsp[-2].nme);
  __theMDP->orig_obs[__theMDP->numorigobs].type = ((int) floor(yyvsp[-1].val));
#ifdef PARSERDEBUG
  fprintf(stderr,"adding observation %d which is %s of type %d ",__theMDP->numorigobs,__theMDP->orig_obs[__theMDP->numorigobs].name, __theMDP->orig_obs[__theMDP->numorigobs].type);
#endif
  __theMDP->numorigobs++;
;
    break;}
case 22:
#line 197 "mdpparser.y"
{
  // initialize initBelief
  __theMDP->initBelief = One;
  Cudd_Ref(__theMDP->initBelief);
#ifdef PARSERDEBUG
  fprintf(stderr,"initialized initBelief to One\n");
#endif
;
    break;}
case 23:
#line 204 "mdpparser.y"
{
  pvcount = 0;
;
    break;}
case 26:
#line 213 "mdpparser.y"
{
  // find this variable in orig_vars 
  if ((curr_primed_ovar = findOVar(__theMDP->orig_vars,__theMDP->numorigvars,yyvsp[0].nme)) < 0) {
    fprintf(stderr,"could not find original variable %s\n",yyvsp[0].nme);
    exit(0);
  }
#ifdef PARSERDEBUG
  fprintf(stderr,"found belief over variable %s\n",__theMDP->orig_vars[curr_primed_ovar].name);
#endif
  level = 0;
;
    break;}
case 27:
#line 223 "mdpparser.y"
{
#ifdef PARSERDEBUG
  fprintf(stderr,"parsed belief over variable %s\n",__theMDP->orig_vars[curr_primed_ovar].name);
  Cudd_PrintDebug(gbm,yyvsp[0].dnode,4,100);
#endif
  __theMDP->initBeliefState[curr_primed_ovar]  = yyvsp[0].dnode;
  Cudd_Ref(__theMDP->initBeliefState[curr_primed_ovar]);
  DdNode *temp1 =  Cudd_addApply(gbm,Cudd_addTimes,__theMDP->initBelief,yyvsp[0].dnode);
  Cudd_Ref(temp1);
  Cudd_RecursiveDeref(gbm,yyvsp[0].dnode);
  Cudd_RecursiveDeref(gbm,__theMDP->initBelief);
  __theMDP->initBelief = temp1;
#ifdef PARSERDEBUG
  fprintf(stderr,"Initial Belief is currently\n");
  Cudd_PrintDebug(gbm,__theMDP->initBelief,4,100);
#endif
;
    break;}
case 28:
#line 241 "mdpparser.y"
{
  unnormalized = false;
;
    break;}
case 29:
#line 243 "mdpparser.y"
{
  unnormalized = true;
;
    break;}
case 34:
#line 254 "mdpparser.y"
{
  conjlevel = 0;
;
    break;}
case 35:
#line 256 "mdpparser.y"
{
  dds[numdds] = yyvsp[-1].dnode;
  Cudd_Ref(dds[numdds]);
  Cudd_RecursiveDeref(gbm,yyvsp[-1].dnode);
  // store in hash table
  ddindices[ddnames[numdds]] = numdds;
#ifdef PARSERDEBUG
    fprintf(stderr,"************************** parsed the dd %s at %d*****\n",ddnames[numdds],ddindices[ddnames[numdds]]);
    Cudd_PrintDebug(gbm,dds[numdds],4,100);
#endif
  numdds++;
  doing_dd = false;
  doing_obs = false;
;
    break;}
case 36:
#line 272 "mdpparser.y"
{
  ddnames[numdds] = strdup(yyvsp[0].nme);
  doing_dd = true;
;
    break;}
case 37:
#line 276 "mdpparser.y"
{
  ddnames[numdds] = strdup(yyvsp[-1].nme);
  // second name is observation variable this is  - find it in the list
  if ((curr_primed_ovar = findOVar(__theMDP->orig_obs, __theMDP->numorigobs, yyvsp[0].nme)) < 0) {
    fprintf(stderr,"could not find original observation %s\n",yyvsp[0].nme);
    exit(0);
  }
  doing_dd = true;
  doing_obs = true;
;
    break;}
case 38:
#line 288 "mdpparser.y"
{  
  if (pvcount != __theMDP->numorigvars)
    error("missing primed variable in action",__theMDP->actionlist[__theMDP->numactions].name);
  pvcount = 0;
  __theMDP->numactions++;
  doing_action = false;
;
    break;}
case 39:
#line 297 "mdpparser.y"
{
  /* name of action is $2, increment action counter */
#ifdef PARSERDEBUG
  fprintf(stderr,"------------- action name %s\n",yyvsp[0].nme); 
#endif
  __theMDP->actionlist[__theMDP->numactions].name = strdup(yyvsp[0].nme);
  doing_action = true;
;
    break;}
case 40:
#line 306 "mdpparser.y"
{
  // the action cost is zero
  setActionCost(__theMDP->actionCost,__theMDP->numactions,0.0);
;
    break;}
case 41:
#line 309 "mdpparser.y"
{
  setActionCost(__theMDP->actionCost,__theMDP->numactions,yyvsp[0].val);
;
    break;}
case 42:
#line 314 "mdpparser.y"
{
  // the action cost is zero
  //fprintf(stderr,"action cost 0\n");
  setActionCost(__theMDP->actionCost, __theMDP->numactions,0.0);
  __theMDP->actionCostNoDummy[__theMDP->numactions] =   __theMDP->actionCost[__theMDP->numactions]; 
  Cudd_Ref(__theMDP->actionCostNoDummy[__theMDP->numactions]);
;
    break;}
case 43:
#line 320 "mdpparser.y"
{ pvindex = 0; level = 0; doing_reward = true; ;
    break;}
case 44:
#line 320 "mdpparser.y"
{
  // assign the action cost
  __theMDP->actionCost[__theMDP->numactions] = yyvsp[0].dnode;
  Cudd_Ref(__theMDP->actionCost[__theMDP->numactions]);
  Cudd_RecursiveDeref(gbm,yyvsp[0].dnode);

  __theMDP->actionCostNoDummy[__theMDP->numactions] =   __theMDP->actionCost[__theMDP->numactions]; 
  Cudd_Ref(__theMDP->actionCostNoDummy[__theMDP->numactions]);

  // TEMPORARILY REMOVE THIS FOR PRINTING OUT
  removeAllDummys(__theMDP->actionCost+__theMDP->numactions,__theMDP->numorigvars);

  // multiply by -1 (since its a cost)
  MixGauss value(-1.0);
  DdNode *temp = Cudd_addConst(gbm,&value);
  Cudd_Ref(temp);
  DdNode *temp2 = Cudd_addApply(gbm,Cudd_addTimes,__theMDP->actionCost[__theMDP->numactions],temp);
  Cudd_Ref(temp2);
  Cudd_RecursiveDeref(gbm,temp);
  Cudd_RecursiveDeref(gbm,__theMDP->actionCost[__theMDP->numactions]);
  __theMDP->actionCost[__theMDP->numactions] = temp2;
  doing_reward = false;
;
    break;}
case 46:
#line 345 "mdpparser.y"
{
  doing_obs = true;
  curr_primed_ovar = 0;
;
    break;}
case 47:
#line 348 "mdpparser.y"
{
  doing_obs = false;
;
    break;}
case 51:
#line 356 "mdpparser.y"
{
  // find this variable in orig_obs
#ifdef PARSERDEBUG
  fprintf(stderr,"looking for %s in orig observation %d\n",yyvsp[0].nme,curr_primed_ovar);
#endif
  if ((curr_primed_ovar = findOVar(__theMDP->orig_obs, __theMDP->numorigobs, yyvsp[0].nme)) < 0) {
    fprintf(stderr,"could not find original observation %s\n",yyvsp[0].nme);
    exit(0);
  }
#ifdef PARSERDEBUG
  fprintf(stderr,"doing observation %d which is %s\n",curr_primed_ovar,__theMDP->orig_obs[curr_primed_ovar].name);
#endif
;
    break;}
case 52:
#line 368 "mdpparser.y"
{
  // observation function in theadd ($3)
  // for the curr_primed_ovar value of the original observation
  // build it at that value of orig var and add it to the running sum
#ifdef PARSERDEBUG
  fprintf(stderr,"got the add for %s\n",__theMDP->orig_obs[curr_primed_ovar].name);
  pdd(yyvsp[0].dnode);
#endif
  __theMDP->totalObsFun[__theMDP->numactions][curr_primed_ovar] = yyvsp[0].dnode;
  Cudd_Ref(__theMDP->totalObsFun[__theMDP->numactions][curr_primed_ovar]);
  //#ifdef PARSERDEBUG
  fprintf(stderr,"observation function %d for observation %d is \n",__theMDP->numactions,curr_primed_ovar);
  pdd(__theMDP->totalObsFun[__theMDP->numactions][curr_primed_ovar]);
  //#endif  
;
    break;}
case 55:
#line 387 "mdpparser.y"
{
  /* primed variable is $1 theadd is its action tree*/
#ifdef PARSERDEBUG
   fprintf(stderr,"primed var %s\n",yyvsp[0].nme);
#endif
  // find this variable in orig_vars 
  if ((curr_primed_ovar = findOVar(__theMDP->orig_vars,__theMDP->numorigvars,yyvsp[0].nme)) < 0) {
    fprintf(stderr,"could not find original variable %s\n",yyvsp[0].nme);
    exit(0);
  }
  level = 0;
  conjlevel = 0;
  // the corresponding binary variables are vars[i] and prime_vars[i]
  // where i:orig_vars[curr_primed_ovar].var1index...orig_vars[curr_primed_ovar].var1index+orig_vars[curr_primed_ovar].nbvars-1
;
    break;}
case 56:
#line 401 "mdpparser.y"
{
  // the new prime diagram is that returned in theadd
  // its $3 because the above action is also counted
  DdNode *tmp2;
  DdNode *temp = sumOutPrime(yyvsp[0].dnode,__theMDP->orig_vars+curr_primed_ovar,__theMDP->prime_vars);
  Cudd_Ref(temp);
  if (unnormalized) {
#ifndef DONOTNORMALIZE
    // renormalize - this can add in Dummy states!
#ifdef PARSERDEBUG
    fprintf(stderr,"dividing this:\n");
    pdd(yyvsp[0].dnode);
    fprintf(stderr,"by this (previous summed over %d):\n",curr_primed_ovar);
    pdd(temp);
#endif    
    tmp2 = Cudd_addApply(gbm,Cudd_addDivide,yyvsp[0].dnode,temp);

    Cudd_Ref(tmp2);

    Cudd_RecursiveDeref(gbm,temp);
    Cudd_RecursiveDeref(gbm,yyvsp[0].dnode);
    temp = tmp2;
    // remove dummy states
    //removeAllDummys(&temp,__theMDP->numorigvars);
    //removeAllDummysp(&temp,__theMDP->numorigvars);
    __theMDP->NewPrime[__theMDP->numactions][curr_primed_ovar] = temp;
    Cudd_Ref(__theMDP->NewPrime[__theMDP->numactions][curr_primed_ovar]);
#else
    __theMDP->NewPrime[__theMDP->numactions][curr_primed_ovar] = yyvsp[0].dnode;
    Cudd_Ref(__theMDP->NewPrime[__theMDP->numactions][curr_primed_ovar]);
    Cudd_RecursiveDeref(gbm,yyvsp[0].dnode);
#endif

  } else {
    // check if normalized
    // if normalized - temp should be exactly 1 only for non-dummy values
    // so, remove dummys, take 1-temp and then remove the dummys again - 
    // the result should be Zero
    removeAllDummys(&temp,__theMDP->numorigvars);
    removeAllDummysp(&temp,__theMDP->numorigvars);
    tmp2 = Cudd_addApply(gbm,Cudd_addMinus,One,temp);
    Cudd_Ref(tmp2);
    Cudd_RecursiveDeref(gbm,temp);
    temp = tmp2;
    removeAllDummys(&temp,__theMDP->numorigvars);
    removeAllDummysp(&temp,__theMDP->numorigvars);
    if (temp != Zero) {
      fprintf(stderr,"CPT for primed variable %s is not normalized:\n",yyvsp[-2].nme);
      Cudd_PrintMinterm(gbm,yyvsp[0].dnode);
      fprintf(stderr,"curr_primed_ovar: %d\n summed diagram:\n",curr_primed_ovar);
      Cudd_PrintDebug(gbm,temp,4,100);
      exit(-1);
    } else {
      __theMDP->NewPrime[__theMDP->numactions][curr_primed_ovar] = yyvsp[0].dnode;
      Cudd_Ref(__theMDP->NewPrime[__theMDP->numactions][curr_primed_ovar]);
      Cudd_RecursiveDeref(gbm,yyvsp[0].dnode);
    }
  } 
  Cudd_RecursiveDeref(gbm,temp);
  removeDummy(&(__theMDP->NewPrime[__theMDP->numactions][curr_primed_ovar]),curr_primed_ovar);


#ifdef PARSERDEBUG
  fprintf(stderr,"the add for %d action %d ovar is \n",__theMDP->numactions, curr_primed_ovar);
  Cudd_PrintDebug(gbm,__theMDP->NewPrime[__theMDP->numactions][curr_primed_ovar],4,100);
#endif
  pvcount++;
;
    break;}
case 57:
#line 473 "mdpparser.y"
{
  oldconj[conjlevel] = conj;
#ifdef PARSERDEBUG
  fprintf(stderr,"saving oldconj %d at level %d before setting conj to 1\n",oldconj[conjlevel],conjlevel);
#endif
  conj = 1;
  conjlevel++;
;
    break;}
case 58:
#line 480 "mdpparser.y"
{
  yyval.dnode = yyvsp[-1].dnode;
  Cudd_Ref(yyval.dnode);
  Cudd_RecursiveDeref(gbm,yyvsp[-1].dnode);

#ifdef PARSERDEBUG
  fprintf(stderr,"resetting conj to oldconj at level %d : %d\n",conjlevel-1,oldconj[conjlevel-1]);
#endif
  conjlevel--;
  conj = oldconj[conjlevel];
;
    break;}
case 59:
#line 491 "mdpparser.y"
{
  oldconj[conjlevel] = conj;
#ifdef PARSERDEBUG
  fprintf(stderr,"saving oldconj %d at level %d before setting conj to 1\n",oldconj[conjlevel],conjlevel);
#endif
  conj = 0;
  conjlevel++;
;
    break;}
case 60:
#line 498 "mdpparser.y"
{
  yyval.dnode = yyvsp[-1].dnode;
  Cudd_Ref(yyval.dnode);
  Cudd_RecursiveDeref(gbm,yyvsp[-1].dnode);
#ifdef PARSERDEBUG
  fprintf(stderr,"resetting conj to oldconj at level %d : %d\n",conjlevel-1,oldconj[conjlevel-1]);
#endif
  conjlevel--;
  conj = oldconj[conjlevel];
;
    break;}
case 61:
#line 508 "mdpparser.y"
{
  /* current root node is $2 */
#ifdef PARSERDEBUG
  fprintf(stderr,"********************************  root node %s\n",yyvsp[0].nme);
#endif
  if ((curr_ovar = findOVar(__theMDP->orig_vars,__theMDP->numorigvars,yyvsp[0].nme)) < 0) 
    error("could not find variable",yyvsp[0].nme);
  levelIndex[level] = curr_ovar;
  branchCount[level] = 0;
;
    break;}
case 62:
#line 518 "mdpparser.y"
{
  yyval.dnode = yyvsp[-1].dnode;
  Cudd_Ref(yyval.dnode);
  Cudd_RecursiveDeref(gbm,yyvsp[-1].dnode);
#ifdef PARSERDEBUG
  fprintf(stderr,"***************************** the add is\n");
  Cudd_PrintDebug(gbm,yyval.dnode,4,100);
#endif
  if ((curr_ovar = findOVar(__theMDP->orig_vars,__theMDP->numorigvars,yyvsp[-3].nme)) < 0) 
    error("could not find variable",yyvsp[-3].nme);
  if (branchCount[level] != __theMDP->orig_vars[curr_ovar].nvals) {
    fprintf(stderr,"%d branches missing from variable %s primed variable %s in action %s\n",
	    __theMDP->orig_vars[curr_ovar].nvals-branchCount[level],
	    __theMDP->orig_vars[curr_ovar].name,__theMDP->orig_vars[curr_primed_ovar].name,__theMDP->actionlist[__theMDP->numactions].name);
    exit(-1);
  }
;
    break;}
case 64:
#line 536 "mdpparser.y"
{

  // check to make sure it matches the primed variables number
  if (!doing_reward) {
    if (doing_obs) {
      fprintf(stderr,"forgot keyword 'observe' in observation parameter list of observation %s\n", __theMDP->orig_obs[curr_primed_ovar].valname[curr_obs_oval]);
      exit(0);
    } else {
      if (pvindex != __theMDP->orig_vars[curr_primed_ovar].nvals) 
	error("Some cpt values missing from primed variable",__theMDP->orig_vars[curr_primed_ovar].name);
    }
    yyval.dnode = yyvsp[-1].dnode;
  } else if (doing_reward && pvindex != 1) {
    error("only one reward value per state","please!");
  } else {
    yyval.dnode = yyvsp[-1].dnode;
  }
  Cudd_Ref(yyval.dnode);
  Cudd_RecursiveDeref(gbm,yyvsp[-1].dnode);
  
  // reset for next time
  pvindex = 0;
  pvsum = 0.0;
;
    break;}
case 65:
#line 560 "mdpparser.y"
{
  numobsparams = 0;
;
    break;}
case 66:
#line 562 "mdpparser.y"
{
#ifdef PARSERDEBUG
  fprintf(stderr,"got observation parameters for observation %s\n", __theMDP->orig_obs[curr_primed_ovar].valname[curr_obs_oval]);
  for (int k=0; k<numobsparams; k++)
    fprintf(stderr,"%f ",obsparams[k]);
#endif
  // now, construct the MixGauss with this set of parameters
  int thetype = __theMDP->orig_obs[curr_primed_ovar].type;
  // thetype describes the kind of MixGauss we're processing here:
  // 0: fvdim = 0 is a discrete multinomial density in the mixweights - numobsparams is nmix and can be anything > 2
  // N>0: fvdim = N is a N-D covariance Gaussian mixture
  //       now the numparams should be either fvdim (type>0, in which case there is one mixture component)
  //       or M*(fvdim+1) for M mixture components
  if (thetype == 0 && numobsparams < 2) {
    fprintf(stderr,"makes no sense - type 0 observation function but < 2 weights... what's the point?\n");
    exit(0);
  } else if (thetype>=1 && 
	     (numobsparams < thetype+thetype*(thetype+1)/2 || 
	      (numobsparams > thetype+thetype*(thetype+1)/2 &&
	       numobsparams%(1+thetype+thetype*(thetype+1)/2) != 0))) {
    fprintf(stderr,"makes no sense - type %d observation function but %d  weights\n",thetype,numobsparams);
    exit(0);
  } 
  MixGauss value(thetype,numobsparams,obsparams);
  yyval.dnode = Cudd_addConst(gbm,&value);
  Cudd_Ref(yyval.dnode);
  //delete value;
#ifdef PARSERDEBUG
  fprintf(stderr,"parsed osbservation function for %s\n", __theMDP->orig_obs[curr_primed_ovar].name);
  pdd(yyval.dnode);
#endif
;
    break;}
case 69:
#line 600 "mdpparser.y"
{
  obsparams[numobsparams++] = yyvsp[0].val;
;
    break;}
case 70:
#line 605 "mdpparser.y"
{
#ifdef PARSERDEBUG
  fprintf(stderr,"parsing root node name %s\n",yyvsp[0].nme);
#endif
  yyval.nme = yyvsp[0].nme;
  levelPrime[level] = false;
;
    break;}
case 71:
#line 612 "mdpparser.y"
{
#ifdef PARSERDEBUG
  fprintf(stderr,"saw a prime on variable %s\n",yyvsp[-1].nme);
#endif
  yyval.nme = yyvsp[-1].nme;
  levelPrime[level] = true;
;
    break;}
case 72:
#line 621 "mdpparser.y"
{
  int i;
  DdNode *temp;
  if (conj == 1) {
#ifdef PARSERDEBUG
    fprintf(stderr,"multiplying these two\n");
#endif    
    temp = Cudd_addApply(gbm,Cudd_addTimes,yyvsp[-1].dnode,yyvsp[0].dnode);
  } else {
#ifdef PARSERDEBUG
    fprintf(stderr,"adding these two\n");
#endif    
    temp = Cudd_addApply(gbm,Cudd_addPlus,yyvsp[-1].dnode,yyvsp[0].dnode);
  }
#ifdef PARSERDEBUG
  Cudd_PrintDebug(gbm,yyvsp[-1].dnode,1,2);
  Cudd_PrintDebug(gbm,yyvsp[0].dnode,1,2);
#endif
  Cudd_Ref(temp);
  Cudd_RecursiveDeref(gbm,yyvsp[-1].dnode);
  Cudd_RecursiveDeref(gbm,yyvsp[0].dnode);
  yyval.dnode = temp;
;
    break;}
case 73:
#line 644 "mdpparser.y"
{
  yyval.dnode = yyvsp[0].dnode;
  Cudd_Ref(yyval.dnode);
  Cudd_RecursiveDeref(gbm,yyvsp[0].dnode);
;
    break;}
case 74:
#line 651 "mdpparser.y"
{
  yyval.dnode = Cudd_addApply(gbm,Cudd_addPlus,yyvsp[-1].dnode,yyvsp[0].dnode);
  Cudd_Ref(yyval.dnode);
  Cudd_RecursiveDeref(gbm,yyvsp[-1].dnode);
  Cudd_RecursiveDeref(gbm,yyvsp[0].dnode);
;
    break;}
case 75:
#line 657 "mdpparser.y"
{
  yyval.dnode = yyvsp[0].dnode;
  Cudd_Ref(yyval.dnode);
  Cudd_RecursiveDeref(gbm,yyvsp[0].dnode);
;
    break;}
case 76:
#line 669 "mdpparser.y"
{
#ifdef PARSERDEBUG
  fprintf(stderr,"branch %s level increasing to %d\n",yyvsp[0].nme,level+1);
#endif
  level++;
;
    break;}
case 77:
#line 674 "mdpparser.y"
{
  level--;
#ifdef PARSERDEBUG
  fprintf(stderr,"parsed branch %s\n",yyvsp[-3].nme); 
  Cudd_PrintDebug(gbm,yyvsp[-1].dnode,4,100);
#endif
  if ((curr_oval = findOVal(__theMDP->orig_vars, levelIndex[level],yyvsp[-3].nme)) < 0) {
    fprintf(stderr,"could not find value %s in variable %s action %s\n",
	    yyvsp[-3].nme,__theMDP->orig_vars[levelIndex[level]].name,__theMDP->actionlist[__theMDP->numactions].name);
    exit(-1);
  }
  // need to check if curr_ovar is a primed variable here
  if (!levelPrime[level] && !doing_obs)
    yyval.dnode = buildCubeCPT(__theMDP->vars,curr_oval,__theMDP->orig_vars[levelIndex[level]],yyvsp[-1].dnode);
  else
    yyval.dnode = buildCubeCPT(__theMDP->prime_vars,curr_oval,__theMDP->orig_vars[levelIndex[level]],yyvsp[-1].dnode);
  Cudd_RecursiveDeref(gbm,yyvsp[-1].dnode);
#ifdef PARSERDEBUG
  fprintf(stderr,"branch CPT is\n");
  Cudd_PrintDebug(gbm,yyval.dnode,4,100);
#endif
  branchCount[level]++;
;
    break;}
case 78:
#line 701 "mdpparser.y"
{
  yyval.dnode = Cudd_addApply(gbm,Cudd_addPlus,yyvsp[-1].dnode,yyvsp[0].dnode);
  Cudd_Ref(yyval.dnode);
  Cudd_RecursiveDeref(gbm,yyvsp[-1].dnode);
  Cudd_RecursiveDeref(gbm,yyvsp[0].dnode);
;
    break;}
case 79:
#line 707 "mdpparser.y"
{
  yyval.dnode = Cudd_addApply(gbm,Cudd_addPlus,yyvsp[-1].dnode,yyvsp[0].dnode);
  Cudd_Ref(yyval.dnode);
  Cudd_RecursiveDeref(gbm,yyvsp[-1].dnode);
  Cudd_RecursiveDeref(gbm,yyvsp[0].dnode);
;
    break;}
case 80:
#line 716 "mdpparser.y"
{
  //MixGauss *value = new MixGauss();
#ifdef REDUCE_PRECISION
  //reduce the precision of the inputs 
  MixGauss value(((double) ((int) (yyvsp[-1].val*10)))/10.0);
  //(*value).set(((double) ((int) ($2*10)))/10.0);
#else
  MixGauss value(yyvsp[-1].val);
  //(*value).set($2);
#endif
  //fprintf(stderr,"setting %f to %f\n",$2,((double) ((int) ($2*100)))/100.0);
  yyval.dnode = Cudd_addConst(gbm,&value);
  Cudd_Ref(yyval.dnode);
  //delete value;
;
    break;}
case 81:
#line 731 "mdpparser.y"
{
  // search for NAME in list of defined DDs
  // by looking in hash table
#ifdef PARSERDEBUG
  fprintf(stderr,"********************** looking for dd %s\n",yyvsp[-1].nme);
#endif
  map<const char*, int, ltstr>::iterator ddit = ddindices.find(yyvsp[-1].nme);
  if (ddit==ddindices.end()) {
    fprintf(stderr,"could not find defined DD %s in list\n",yyvsp[-1].nme);
    exit(-1);
  } else {
#ifdef PARSERDEBUG
    fprintf(stderr,"************************** Found defined DD %s in list at %d\n",yyvsp[-1].nme,(*ddit).second);
    Cudd_PrintDebug(gbm,dds[(*ddit).second],1,2);
#endif
    yyval.dnode = dds[(*ddit).second];
    Cudd_Ref(yyval.dnode);
  }
;
    break;}
case 82:
#line 752 "mdpparser.y"
{
  //MixGauss *value = new MixGauss();
#ifdef REDUCE_PRECISION
  //reduce the precision of the inputs 
  MixGauss value(((double) ((int) (yyvsp[0].val*10)))/10.0);
  //(*value).set(((double) ((int) ($1*10)))/10.0);
#else
  MixGauss value(yyvsp[0].val);
  //(*value).set($1);
#endif
  //fprintf(stderr,"setting %f to %f\n",$1,((double) ((int) ($1*100)))/100.0);

  DdNode *temp = Cudd_addConst(gbm,&value);
  Cudd_Ref(temp);
  //delete value;
  if (!doing_reward) {
#ifdef PARSERDEBUG
    fprintf(stderr,"found a constant %f at variable %s\n",yyvsp[0].val,__theMDP->orig_vars[curr_primed_ovar].name);
#endif
    // buildCubeCPT does its own reffing
    if (!doing_obs) {
      yyval.dnode = buildCubeCPT(__theMDP->prime_vars,pvindex,__theMDP->orig_vars[curr_primed_ovar],temp);
    } else {
      fprintf(stderr," what the fuck? Should not be happening in parser!\n");
    }
    Cudd_RecursiveDeref(gbm,temp);
  } else {
    //fprintf(stderr," found a constant %f pvindex %d\n",$1,pvindex);
    yyval.dnode = temp;
  }    
  pvindex++;
  pvsum += yyvsp[0].val;
;
    break;}
case 83:
#line 788 "mdpparser.y"
{ 
#ifdef PARSERDEBUG
  fprintf(stderr, "doing reward\n");
#endif
  pvindex = 0; level = 0; doing_reward = true; 
;
    break;}
case 84:
#line 793 "mdpparser.y"
{
  __theMDP->RewardD = yyvsp[0].dnode;
  Cudd_Ref(__theMDP->RewardD);
  Cudd_RecursiveDeref(gbm,yyvsp[0].dnode);

  __theMDP->RewardDNoDummy =   __theMDP->RewardD; 
  Cudd_Ref(__theMDP->RewardDNoDummy);

  removeAllDummys(&(__theMDP->RewardD),__theMDP->numorigvars);
;
    break;}
case 85:
#line 806 "mdpparser.y"
{ 
  MixGauss discountMixGauss;
  discountMixGauss.set(yyvsp[0].val);
  __theMDP->discount = Cudd_addConst(gbm,&discountMixGauss);
  Cudd_Ref(__theMDP->discount);
;
    break;}
case 86:
#line 811 "mdpparser.y"
{
  doing_reward = true;
;
    break;}
case 87:
#line 813 "mdpparser.y"
{
  __theMDP->discount = yyvsp[0].dnode;
  Cudd_Ref(__theMDP->discount);
;
    break;}
case 90:
#line 822 "mdpparser.y"
{ 
  tolerance = yyvsp[0].val; 
  __theMDP->horizon = -1.0;
;
    break;}
case 91:
#line 828 "mdpparser.y"
{
  tolerance = 0.0;
  __theMDP->horizon = int(yyvsp[0].val);
;
    break;}
}

#line 705 "/usr/share/bison/bison.simple"


  yyvsp -= yylen;
  yyssp -= yylen;
#if YYLSP_NEEDED
  yylsp -= yylen;
#endif

#if YYDEBUG
  if (yydebug)
    {
      short *yyssp1 = yyss - 1;
      YYFPRINTF (stderr, "state stack now");
      while (yyssp1 != yyssp)
	YYFPRINTF (stderr, " %d", *++yyssp1);
      YYFPRINTF (stderr, "\n");
    }
#endif

  *++yyvsp = yyval;
#if YYLSP_NEEDED
  *++yylsp = yyloc;
#endif

  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTBASE] + *yyssp;
  if (yystate >= 0 && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTBASE];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;

#ifdef YYERROR_VERBOSE
      yyn = yypact[yystate];

      if (yyn > YYFLAG && yyn < YYLAST)
	{
	  YYSIZE_T yysize = 0;
	  char *yymsg;
	  int yyx, yycount;

	  yycount = 0;
	  /* Start YYX at -YYN if negative to avoid negative indexes in
	     YYCHECK.  */
	  for (yyx = yyn < 0 ? -yyn : 0;
	       yyx < (int) (sizeof (yytname) / sizeof (char *)); yyx++)
	    if (yycheck[yyx + yyn] == yyx)
	      yysize += yystrlen (yytname[yyx]) + 15, yycount++;
	  yysize += yystrlen ("parse error, unexpected ") + 1;
	  yysize += yystrlen (yytname[YYTRANSLATE (yychar)]);
	  yymsg = (char *) YYSTACK_ALLOC (yysize);
	  if (yymsg != 0)
	    {
	      char *yyp = yystpcpy (yymsg, "parse error, unexpected ");
	      yyp = yystpcpy (yyp, yytname[YYTRANSLATE (yychar)]);

	      if (yycount < 5)
		{
		  yycount = 0;
		  for (yyx = yyn < 0 ? -yyn : 0;
		       yyx < (int) (sizeof (yytname) / sizeof (char *));
		       yyx++)
		    if (yycheck[yyx + yyn] == yyx)
		      {
			const char *yyq = ! yycount ? ", expecting " : " or ";
			yyp = yystpcpy (yyp, yyq);
			yyp = yystpcpy (yyp, yytname[yyx]);
			yycount++;
		      }
		}
	      yyerror (yymsg);
	      YYSTACK_FREE (yymsg);
	    }
	  else
	    yyerror ("parse error; also virtual memory exhausted");
	}
      else
#endif /* defined (YYERROR_VERBOSE) */
	yyerror ("parse error");
    }
  goto yyerrlab1;


/*--------------------------------------------------.
| yyerrlab1 -- error raised explicitly by an action |
`--------------------------------------------------*/
yyerrlab1:
  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      /* return failure if at end of input */
      if (yychar == YYEOF)
	YYABORT;
      YYDPRINTF ((stderr, "Discarding token %d (%s).\n",
		  yychar, yytname[yychar1]));
      yychar = YYEMPTY;
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */

  yyerrstatus = 3;		/* Each real token shifted decrements this */

  goto yyerrhandle;


/*-------------------------------------------------------------------.
| yyerrdefault -- current state does not do anything special for the |
| error token.                                                       |
`-------------------------------------------------------------------*/
yyerrdefault:
#if 0
  /* This is wrong; only states that explicitly want error tokens
     should shift them.  */

  /* If its default is to accept any token, ok.  Otherwise pop it.  */
  yyn = yydefact[yystate];
  if (yyn)
    goto yydefault;
#endif


/*---------------------------------------------------------------.
| yyerrpop -- pop the current state because it cannot handle the |
| error token                                                    |
`---------------------------------------------------------------*/
yyerrpop:
  if (yyssp == yyss)
    YYABORT;
  yyvsp--;
  yystate = *--yyssp;
#if YYLSP_NEEDED
  yylsp--;
#endif

#if YYDEBUG
  if (yydebug)
    {
      short *yyssp1 = yyss - 1;
      YYFPRINTF (stderr, "Error: state stack now");
      while (yyssp1 != yyssp)
	YYFPRINTF (stderr, " %d", *++yyssp1);
      YYFPRINTF (stderr, "\n");
    }
#endif

/*--------------.
| yyerrhandle.  |
`--------------*/
yyerrhandle:
  yyn = yypact[yystate];
  if (yyn == YYFLAG)
    goto yyerrdefault;

  yyn += YYTERROR;
  if (yyn < 0 || yyn > YYLAST || yycheck[yyn] != YYTERROR)
    goto yyerrdefault;

  yyn = yytable[yyn];
  if (yyn < 0)
    {
      if (yyn == YYFLAG)
	goto yyerrpop;
      yyn = -yyn;
      goto yyreduce;
    }
  else if (yyn == 0)
    goto yyerrpop;

  if (yyn == YYFINAL)
    YYACCEPT;

  YYDPRINTF ((stderr, "Shifting error token, "));

  *++yyvsp = yylval;
#if YYLSP_NEEDED
  *++yylsp = yylloc;
#endif

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

/*---------------------------------------------.
| yyoverflowab -- parser overflow comes here.  |
`---------------------------------------------*/
yyoverflowlab:
  yyerror ("parser stack overflow");
  yyresult = 2;
  /* Fall through.  */

yyreturn:
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
  return yyresult;
}
#line 833 "mdpparser.y"


#include "lex.yy.c"
// function to add in zeros for all the dummy variables
// do this by mutliplying each NewPrime diagram by 
// a 'goodstate' ADD for each variable. The goodstate ADD
// for an original (mv) variable is the sum of all cubes over
// binary variables corresponding to valid original values 
// and 0 to all others (which correspond to dummy values)
void removeDummy(DdNode **np, int cpovar) {
   DdNode *tmp,*tmp2;
  tmp = Cudd_addApply(gbm,Cudd_addTimes,goodState[cpovar],np[0]);
  Cudd_Ref(tmp);
  Cudd_RecursiveDeref(gbm,np[0]);
  np[0] = tmp;
}
// same for primed variables
void removeDummyp(DdNode **np, int cpovar) {
   DdNode *tmp,*tmp2;
  tmp = Cudd_addApply(gbm,Cudd_addTimes,goodStatep[cpovar],np[0]);
  Cudd_Ref(tmp);
  Cudd_RecursiveDeref(gbm,np[0]);
  np[0] = tmp;
}
void removeAllDummys(DdNode **np, int nv) {
  int v;
  for (v=0; v<nv; v++) 
    removeDummy(np,v);
}
void removeAllDummysp(DdNode **np, int nv) {
  int v;
  for (v=0; v<nv; v++) 
    removeDummyp(np,v);
}
// add a 1 to np for each dummy state
void addDummyStates(DdNode **np, int nv) {
  int v;
  DdNode *tmp,*tmp2;
  for (v=0; v<nv; v++) {
    tmp = Cudd_addApply(gbm,Cudd_addMinus,One,goodState[v]);
    Cudd_Ref(tmp);
    tmp2 = Cudd_addApply(gbm,Cudd_addPlus,np[0],tmp);
    Cudd_Ref(tmp2);
    Cudd_RecursiveDeref(gbm,tmp);
    Cudd_RecursiveDeref(gbm,np[0]);
    np[0] = tmp2;
  }
    
}
// function to build the goodState ADDs once at the start
void buildGoodStateADDs(onum *ov, rnum *vars, rnum *pvars, int nov) {
  int a,v,va,pv,i,j;
  goodState = new DdNode*[nov];
  goodStatep = new DdNode*[nov];
  DdNode *tmp,*tmp2;
  // construct the goodState ADDs for each varable
  for (v=0; v < nov; v++) {
    goodState[v] = Zero;
    Cudd_Ref(goodState[v]);
    for (va=0; va< ov[v].nvals; va++) {
      tmp = buildCubeCPT(vars,va,ov[v],One);
      tmp2 = Cudd_addApply(gbm,Cudd_addPlus,tmp,goodState[v]);
      Cudd_Ref(tmp2);
      Cudd_RecursiveDeref(gbm,tmp);
      Cudd_RecursiveDeref(gbm,goodState[v]);
      goodState[v] = tmp2;
    }
    // build prime good state ADDs too
    goodStatep[v] = Zero;
    Cudd_Ref(goodStatep[v]);
    for (va=0; va< ov[v].nvals; va++) {
      tmp = buildCubeCPT(pvars,va,ov[v],One);
      tmp2 = Cudd_addApply(gbm,Cudd_addPlus,tmp,goodStatep[v]);
      Cudd_Ref(tmp2);
      Cudd_RecursiveDeref(gbm,tmp);
      Cudd_RecursiveDeref(gbm,goodStatep[v]);
      goodStatep[v] = tmp2;
    }
    //fprintf(stderr,"goodState for %d is \n",v);
    //Cudd_PrintDebug(gbm,goodState[v],4,100);
  }
}

void setActionCost(DdNode **ac, int nac, double val) {
  //  fprintf(stderr,"Action cost is %f for action %d\n",val,nac);
  //MixGauss *value = new MixGauss();
  //  (*value).set((double)(-1.0*fabs(val)));
  MixGauss value((double)(-1.0*fabs(val)));
  ac[nac] = Cudd_addConst(gbm,&value);
  Cudd_Ref(ac[nac]);
  //delete value;
}
// builds an add  rooted at a cube over the 
// variables in v corresponding to ovar's coval branch
// with add at that leaf
DdNode * buildCubeCPT(rnum *v, int coval, onum ovar, DdNode *add) {
  // the returned  is a DdNode* referring to the add at this branch
  //build the cube for this branch
  int nbv = ovar.nbvars;
  int nbvals = int(pow(2.0,nbv));
  int *phase = new int[nbv];
  int i;
  for (i=0; i<nbv; i++)
    phase[i] = 0;
  //phase should be nbvals-coval in binary
  int tmp = nbvals-coval-1;
  i=nbv-1;

  while (tmp > 0 && i >=0) {
    phase[i--] = tmp%2;
    tmp = tmp/2;
  }
  /*
  fprintf(stderr,"variable %s coval %d nbv: %d phase: ",ovar.name,tmp,nbv);
  for (i=0; i<nbv; i++)
    fprintf(stderr,"%d ",phase[i]);
  fprintf(stderr,"\n");
  */
  DdNode **arrayofvars = new DdNode*[nbv];
  for (i=0; i<nbv; i++) {
    arrayofvars[i] = v[ovar.var1index+i].add_var;
    Cudd_Ref(arrayofvars[i]);
  }
  DdNode *branch = Cudd_addComputeCube(gbm,arrayofvars,phase,nbv);
  Cudd_Ref(branch);
  for (i=0; i<nbv; i++) 
    Cudd_RecursiveDeref(gbm,arrayofvars[i]);

  //fprintf(stderr,"built branch:\n");
  //Cudd_PrintDebug(gbm,branch,4,100);

  DdNode *cubecpt = Cudd_addApply(gbm,Cudd_addTimes,branch,add);
  Cudd_Ref(cubecpt);
  Cudd_RecursiveDeref(gbm,branch);

  //fprintf(stderr,"built cube:\n");
  //Cudd_PrintDebug(gbm,cubecpt,4,100);

  delete [] phase;
  delete [] arrayofvars;
  return cubecpt;
}


int findOVal(onum *ov, int ovar, const char *findval) {
  int oval = 0;
  while (oval < ov[ovar].nvals && strcmp(ov[ovar].valname[oval],findval) != 0)
    oval++;
  if (oval < ov[ovar].nvals)
    return oval;
  else 
    return -1;
}


// returns the index of orignal variable findname, or -1 if not found
int findOVar(onum *ov, int nov, const char *findname) {
  int ovar = 0;
  while (ovar < nov && strcmp(ov[ovar].name,findname) != 0)
    ovar++;
  if (ovar < nov)
    return ovar;
  else 
    return -1;
}
// adds a new multi-valued variable to a list of such variables
void addVar(onum *ovar, int vnum, char * varn) {
  ovar[vnum].valname[ovar[vnum].nvals] = strdup(varn);
#ifdef PARSERDEBUG
  fprintf(stderr,"adding %d value (%s) to var %d\n",ovar[vnum].nvals,ovar[vnum].valname[ovar[vnum].nvals],vnum);
#endif
  ovar[vnum].nvals++;
}
//renormalizes a cube of primed variables
// result is not reffed
// does not deref cube
DdNode * renormalizeCube(DdNode *cube, double pvsum) {
  // renormalize the sum - divide by pvsum
  //MixGauss *value = new MixGauss(1.0/pvsum);
  MixGauss value(1.0/pvsum);
  DdNode *temp = Cudd_addConst(gbm,&value);
  Cudd_Ref(temp);
  DdNode *temp2 = Cudd_addApply(gbm,Cudd_addTimes,cube,temp);
  Cudd_RecursiveDeref(gbm,temp);
  //delete value;
  return temp2;
}

int yyerror (const char *s)  /* Called by yyparse on error */
{
  fprintf (stderr,"%s\n", s);
  exit(-1);
}
void error(const char *s, const char *v) {
  fprintf(stderr, "%s %s\n",s,v);
  exit(-1);
}

