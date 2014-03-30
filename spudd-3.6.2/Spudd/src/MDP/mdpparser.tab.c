/* A Bison parser, made by GNU Bison 2.7.12-4996.  */

/* Bison implementation for Yacc-like parsers in C
   
      Copyright (C) 1984, 1989-1990, 2000-2013 Free Software Foundation, Inc.
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.
   
   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.7.12-4996"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* Copy the first part of user declarations.  */
/* Line 371 of yacc.c  */
#line 4 "mdpparser.y"


  //#define USE_POMDP
#ifndef USE_POMDP
  #include "MDP.h"
#else
  #include "POMDP.h"
#endif



  //#define YYSTYPE double
  #define MAXDDS 100000
  //#define REDUCE_PRECISION 1000
  //#define PARSERDEBUG 1
  //#define PARSERDEBUG2 1
  int yylex( void );
  int yyerror( const char *s);
  void error(const char *s, const char *v);
  int curr_primed_ovar, curr_ovar, curr_oval, renormvar;
  DdNode *dds[MAXDDS];
  char *ddnames[MAXDDS];
  DdNode **goodState;
  DdNode **goodStatep;
  double cpt[MAXVALS];
  double pvsum;
  int valindex, level, pvcount, conj, pvindex, numdds;
  int oldconj[64], conjlevel;
  bool doing_reward, doing_dd, doing_obs, doing_action, unnormalized;
  int levelIndex[MAXVARS];
  // tells if a level is primed variable
  bool levelPrime[MAXVARS];
  int branchCount[MAXVARS];
#ifndef USE_POMDP
  extern MDP *__theMDP;
#else
  extern POMDP *__theMDP;
#endif
  map<const char*, int, ltstr> ddindices;

/* Line 371 of yacc.c  */
#line 109 "mdpparser.tab.c"

# ifndef YY_NULL
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULL nullptr
#  else
#   define YY_NULL 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* In a future release of Bison, this section will be replaced
   by #include "mdpparser.tab.h".  */
#ifndef YY_YY_MDPPARSER_TAB_H_INCLUDED
# define YY_YY_MDPPARSER_TAB_H_INCLUDED
/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     NAME = 258,
     REAL = 259,
     INTGR = 260,
     OPP = 261,
     CLP = 262,
     VARIABLES = 263,
     DISCOUNT = 264,
     TOLERANCE = 265,
     REWARD = 266,
     VALUE = 267,
     ACTION = 268,
     ENDACTION = 269,
     OBSERVATIONS = 270,
     OBSERVE = 271,
     ENDOBSERVE = 272,
     BELIEF = 273,
     ENDBELIEF = 274,
     DISJ = 275,
     CONJ = 276,
     RENORM = 277,
     OSB = 278,
     CSB = 279,
     HORIZON = 280,
     VAL = 281,
     COST = 282,
     PRIME = 283,
     UNNORM = 284,
     STARTDD = 285,
     ENDDD = 286
   };
#endif


#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{
/* Line 387 of yacc.c  */
#line 44 "mdpparser.y"

  double val;
  int ival;
  DdNode * dnode;
  DdNode ** p_dnode;
  double *darray;
  char *nme;


/* Line 387 of yacc.c  */
#line 193 "mdpparser.tab.c"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif

extern YYSTYPE yylval;

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */

#endif /* !YY_YY_MDPPARSER_TAB_H_INCLUDED  */

/* Copy the second part of user declarations.  */

/* Line 390 of yacc.c  */
#line 221 "mdpparser.tab.c"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef __attribute__
/* This feature is available in gcc versions 2.5 and later.  */
# if (! defined __GNUC__ || __GNUC__ < 2 \
      || (__GNUC__ == 2 && __GNUC_MINOR__ < 5))
#  define __attribute__(Spec) /* empty */
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif


/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(N) (N)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int yyi)
#else
static int
YYID (yyi)
    int yyi;
#endif
{
  return yyi;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)				\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack_alloc, Stack, yysize);			\
	Stack = &yyptr->Stack_alloc;					\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (YYID (0))
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  6
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   97

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  32
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  58
/* YYNRULES -- Number of rules.  */
#define YYNRULES  85
/* YYNRULES -- Number of states.  */
#define YYNSTATES  141

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   286

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
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
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint8 yyprhs[] =
{
       0,     0,     3,     5,    14,    19,    21,    24,    29,    31,
      34,    37,    38,    43,    45,    48,    53,    56,    59,    61,
      62,    63,    68,    71,    73,    74,    78,    79,    81,    84,
      86,    88,    90,    91,    96,    99,   106,   109,   110,   112,
     113,   114,   118,   119,   120,   121,   127,   130,   132,   133,
     137,   140,   142,   143,   147,   148,   154,   155,   161,   162,
     170,   171,   177,   179,   183,   185,   188,   191,   193,   196,
     198,   199,   205,   208,   211,   215,   219,   221,   222,   226,
     229,   230,   234,   236,   238,   241
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      33,     0,    -1,    34,    -1,    35,    40,    45,    50,    51,
      83,    85,    87,    -1,     6,     8,    36,     7,    -1,    37,
      -1,    36,    37,    -1,     6,    38,    39,     7,    -1,     3,
      -1,    39,     3,    -1,     3,     3,    -1,    -1,     6,    15,
      41,     7,    -1,    42,    -1,    41,    42,    -1,     6,    44,
      43,     7,    -1,    43,     3,    -1,     3,     3,    -1,     3,
      -1,    -1,    -1,    18,    46,    47,    19,    -1,    47,    48,
      -1,    48,    -1,    -1,     3,    49,    70,    -1,    -1,    29,
      -1,    51,    52,    -1,    52,    -1,    56,    -1,    53,    -1,
      -1,    55,    54,    70,    31,    -1,    30,     3,    -1,    57,
      58,    67,    61,    59,    14,    -1,    13,     3,    -1,    -1,
       4,    -1,    -1,    -1,    27,    60,    70,    -1,    -1,    -1,
      -1,    16,    62,    64,    63,    17,    -1,    64,    65,    -1,
      65,    -1,    -1,     3,    66,    70,    -1,    67,    68,    -1,
      68,    -1,    -1,     3,    69,    70,    -1,    -1,    23,    21,
      71,    76,    24,    -1,    -1,    23,    20,    72,    76,    24,
      -1,    -1,    23,    22,     3,    28,    73,    70,    24,    -1,
      -1,     6,    75,    74,    77,     7,    -1,    81,    -1,     6,
      80,     7,    -1,     3,    -1,     3,    28,    -1,    76,    70,
      -1,    70,    -1,    77,    78,    -1,    78,    -1,    -1,     6,
       3,    79,    70,     7,    -1,    80,    82,    -1,    82,    82,
      -1,     6,     4,     7,    -1,     6,     3,     7,    -1,     4,
      -1,    -1,    11,    84,    70,    -1,     9,     4,    -1,    -1,
       9,    86,    70,    -1,    88,    -1,    89,    -1,    10,     4,
      -1,    25,     4,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,    77,    77,    80,    85,   138,   138,   141,   147,   155,
     157,   162,   163,   185,   185,   188,   195,   199,   206,   214,
     216,   216,   230,   230,   234,   234,   268,   270,   275,   275,
     278,   278,   281,   281,   300,   306,   315,   324,   327,   332,
     338,   338,   365,   366,   368,   366,   373,   373,   376,   376,
     392,   392,   395,   395,   507,   507,   525,   525,   542,   542,
     567,   567,   596,   597,   625,   632,   641,   666,   673,   679,
     691,   691,   727,   733,   741,   753,   776,   814,   814,   831,
     836,   836,   844,   844,   847,   853
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "NAME", "REAL", "INTGR", "OPP", "CLP",
  "VARIABLES", "DISCOUNT", "TOLERANCE", "REWARD", "VALUE", "ACTION",
  "ENDACTION", "OBSERVATIONS", "OBSERVE", "ENDOBSERVE", "BELIEF",
  "ENDBELIEF", "DISJ", "CONJ", "RENORM", "OSB", "CSB", "HORIZON", "VAL",
  "COST", "PRIME", "UNNORM", "STARTDD", "ENDDD", "$accept", "input", "mdp",
  "varilist", "varlist", "vardec", "varname", "vallist", "obslist",
  "oblist", "obdec", "obsvallist", "obsname", "ibelief", "$@1",
  "belieflist", "belief", "$@2", "unn", "actslist", "actionordd", "thedd",
  "$@3", "ddname", "action", "actionname", "oldactioncost", "actioncost",
  "$@4", "observation", "$@5", "$@6", "obsfunlist", "obsfun", "$@7",
  "acttreelist", "acttree", "$@8", "theadd", "$@9", "$@10", "$@11", "$@12",
  "rootnodename", "con_dis_add", "currentadd", "subadd", "$@13",
  "primeadd", "constnode", "constadd", "reward", "$@14", "disc", "$@15",
  "tolhor", "tol", "hor", YY_NULL
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    32,    33,    34,    35,    36,    36,    37,    38,    39,
      39,    40,    40,    41,    41,    42,    43,    43,    44,    45,
      46,    45,    47,    47,    49,    48,    50,    50,    51,    51,
      52,    52,    54,    53,    55,    56,    57,    58,    58,    59,
      60,    59,    61,    62,    63,    61,    64,    64,    66,    65,
      67,    67,    69,    68,    71,    70,    72,    70,    73,    70,
      74,    70,    70,    70,    75,    75,    76,    76,    77,    77,
      79,    78,    80,    80,    81,    81,    82,    84,    83,    85,
      86,    85,    87,    87,    88,    89
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     8,     4,     1,     2,     4,     1,     2,
       2,     0,     4,     1,     2,     4,     2,     2,     1,     0,
       0,     4,     2,     1,     0,     3,     0,     1,     2,     1,
       1,     1,     0,     4,     2,     6,     2,     0,     1,     0,
       0,     3,     0,     0,     0,     5,     2,     1,     0,     3,
       2,     1,     0,     3,     0,     5,     0,     5,     0,     7,
       0,     5,     1,     3,     1,     2,     2,     1,     2,     1,
       0,     5,     2,     2,     3,     3,     1,     0,     3,     2,
       0,     3,     1,     1,     2,     2
};

/* YYDEFACT[STATE-NAME] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,     0,     2,    11,     0,     1,     0,    19,     0,
       0,     5,     0,    20,    26,     8,     0,     4,     6,     0,
       0,    13,     0,    27,     0,     0,     0,    18,     0,    12,
      14,    24,     0,    23,     0,     0,     0,    29,    31,    32,
      30,    37,    10,     9,     7,     0,     0,     0,    21,    22,
      36,    34,    77,    28,     0,     0,    38,     0,    17,    16,
      15,     0,     0,    25,    62,     0,    80,     0,     0,    52,
      42,    51,    64,    76,    60,     0,     0,    56,    54,     0,
      78,    79,     0,     0,     0,     3,    82,    83,    33,     0,
      43,    39,    50,    75,    65,    74,     0,    76,    63,    72,
      73,     0,     0,     0,    81,    84,    85,    53,     0,    40,
       0,     0,     0,    69,    67,     0,     0,    58,    48,    44,
      47,     0,    35,    70,    61,    68,    57,    66,    55,     0,
       0,     0,    46,    41,     0,     0,    49,    45,     0,    59,
      71
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     2,     3,     4,    10,    11,    16,    26,     8,    20,
      21,    46,    28,    14,    22,    32,    33,    47,    24,    36,
      37,    38,    55,    39,    40,    41,    57,   110,   121,    91,
     108,   131,   119,   120,   130,    70,    71,    89,   114,   102,
     101,   129,    96,    74,   115,   112,   113,   134,    75,    64,
      76,    54,    65,    67,    82,    85,    86,    87
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -48
static const yytype_int8 yypact[] =
{
      -1,    24,    17,   -48,     1,    28,   -48,    36,    23,    49,
      33,   -48,    47,   -48,    25,   -48,    52,   -48,   -48,    53,
      37,   -48,    54,   -48,   -10,    55,     7,   -48,    56,   -48,
     -48,   -48,    10,   -48,    57,    58,    -9,   -48,   -48,   -48,
     -48,    59,   -48,   -48,   -48,    61,    12,     5,   -48,   -48,
     -48,   -48,   -48,   -48,    62,     5,   -48,    63,   -48,   -48,
     -48,    42,    16,   -48,   -48,     5,    66,     6,    31,   -48,
       9,   -48,    -6,    60,   -48,    26,    68,   -48,   -48,    70,
     -48,   -48,     5,    71,    72,   -48,   -48,   -48,   -48,     5,
     -48,    38,   -48,   -48,   -48,   -48,    73,   -48,   -48,   -48,
     -48,     5,     5,    50,   -48,   -48,   -48,   -48,    74,   -48,
      67,    77,    41,   -48,   -48,     0,     3,   -48,   -48,    74,
     -48,     5,   -48,   -48,   -48,   -48,   -48,   -48,   -48,     5,
       5,    69,   -48,   -48,     5,    64,   -48,   -48,    78,   -48,
     -48
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -48,   -48,   -48,   -48,   -48,    79,   -48,   -48,   -48,   -48,
      75,   -48,   -48,   -48,   -48,   -48,    65,   -48,   -48,   -48,
      48,   -48,   -48,   -48,   -48,   -48,   -48,   -48,   -48,   -48,
     -48,   -48,   -48,   -29,   -48,   -48,    21,   -48,   -47,   -48,
     -48,   -48,   -48,   -48,    -8,   -48,   -20,   -48,   -48,   -48,
     -26,   -48,   -48,   -48,   -48,   -48,   -48,   -48
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const yytype_uint8 yytable[] =
{
      63,    93,    52,    34,    34,     1,    61,     7,    68,    61,
      43,    61,    69,    31,    44,    59,    83,     6,    80,    60,
      35,    35,    94,    62,   126,    90,    62,   128,    62,    48,
      97,    84,     5,    98,     9,   104,    77,    78,    79,     9,
      17,    13,   107,    19,    29,    72,    73,   111,   124,    99,
     100,    12,    15,    19,    23,    25,    27,    31,    42,    45,
      50,    51,    88,    56,    58,   109,    69,    95,   127,   127,
      81,    66,    97,   103,   133,   105,   106,   118,   117,   111,
     123,   122,   135,   136,    53,   140,   137,   138,   139,    18,
     132,    92,   125,     0,   116,    30,     0,    49
};

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-48)))

#define yytable_value_is_error(Yytable_value) \
  YYID (0)

static const yytype_int16 yycheck[] =
{
      47,     7,    11,    13,    13,     6,     6,     6,    55,     6,
       3,     6,     3,     3,     7,     3,    10,     0,    65,     7,
      30,    30,    28,    23,    24,    16,    23,    24,    23,    19,
       4,    25,     8,     7,     6,    82,    20,    21,    22,     6,
       7,    18,    89,     6,     7,     3,     4,     6,     7,    75,
      76,    15,     3,     6,    29,     3,     3,     3,     3,     3,
       3,     3,    31,     4,     3,    27,     3,     7,   115,   116,
       4,     9,     4,     3,   121,     4,     4,     3,    28,     6,
       3,    14,   129,   130,    36,     7,    17,   134,    24,    10,
     119,    70,   112,    -1,   102,    20,    -1,    32
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     6,    33,    34,    35,     8,     0,     6,    40,     6,
      36,    37,    15,    18,    45,     3,    38,     7,    37,     6,
      41,    42,    46,    29,    50,     3,    39,     3,    44,     7,
      42,     3,    47,    48,    13,    30,    51,    52,    53,    55,
      56,    57,     3,     3,     7,     3,    43,    49,    19,    48,
       3,     3,    11,    52,    83,    54,     4,    58,     3,     3,
       7,     6,    23,    70,    81,    84,     9,    85,    70,     3,
      67,    68,     3,     4,    75,    80,    82,    20,    21,    22,
      70,     4,    86,    10,    25,    87,    88,    89,    31,    69,
      16,    61,    68,     7,    28,     7,    74,     4,     7,    82,
      82,    72,    71,     3,    70,     4,     4,    70,    62,    27,
      59,     6,    77,    78,    70,    76,    76,    28,     3,    64,
      65,    60,    14,     3,     7,    78,    24,    70,    24,    73,
      66,    63,    65,    70,    79,    70,    70,    17,    70,    24,
       7
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  However,
   YYFAIL appears to be in use.  Nevertheless, it is formally deprecated
   in Bison 2.4.2's NEWS entry, where a plan to phase it out is
   discussed.  */

#define YYFAIL		goto yyerrlab
#if defined YYFAIL
  /* This is here to suppress warnings from the GCC cpp's
     -Wunused-macros.  Normally we don't worry about that warning, but
     some users do, and we want to make it easy for users to remove
     YYFAIL uses, which will produce warnings from Bison 2.5.  */
#endif

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))

/* Error token number */
#define YYTERROR	1
#define YYERRCODE	256


/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */
#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

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
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  YYUSE (yytype);
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
#else
static void
yy_stack_print (yybottom, yytop)
    yytype_int16 *yybottom;
    yytype_int16 *yytop;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULL, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULL;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - Assume YYFAIL is not used.  It's too flawed to consider.  See
       <http://lists.gnu.org/archive/html/bison-patches/2009-12/msg00024.html>
       for details.  YYERROR is fine as it does not invoke this
       function.
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULL, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YYUSE (yytype);
}




/* The lookahead symbol.  */
int yychar;


#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval YY_INITIAL_VALUE(yyval_default);

/* Number of syntax errors so far.  */
int yynerrs;


/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       `yyss': related to states.
       `yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;

	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),
		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss_alloc, yyss);
	YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

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

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 3:
/* Line 1787 of yacc.c  */
#line 80 "mdpparser.y"
    {
  ddindices.clear();
}
    break;

  case 4:
/* Line 1787 of yacc.c  */
#line 85 "mdpparser.y"
    {
  // now we have the variables, so allocate for the NewPrime Diagrams
  // These are diagrams with a single multi-valued primed variable in them,
  // which is translated into log2(nvals) binary variables
  // so we only need numorigvars NewPrime diagrams for each action
  int i,j;
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
#ifdef PARSERDEBUG2
    Cudd_PrintDebug(gbm,dds[numdds],4,100);
#endif
#endif
    numdds++;

    // and a bunch of dds - one for each value of each variable
    for (j=0; j<__theMDP->orig_vars[i].nvals; j++) {
      ddnames[numdds] = new char[128];
      sprintf(ddnames[numdds],"%s%s",__theMDP->orig_vars[i].name,__theMDP->orig_vars[i].valname[j]);
      dds[numdds] = __theMDP->buildOneCubeOrig(i,j,true);
      Cudd_Ref(dds[numdds]);
      ddindices[ddnames[numdds]] = numdds;
      numdds++;
    }
  }
}
    break;

  case 7:
/* Line 1787 of yacc.c  */
#line 141 "mdpparser.y"
    {
  __theMDP->newADDVar();
  __theMDP->numorigvars++;
}
    break;

  case 8:
/* Line 1787 of yacc.c  */
#line 147 "mdpparser.y"
    { 
  __theMDP->orig_vars[__theMDP->numorigvars].name = strdup((yyvsp[(1) - (1)].nme));
#ifdef PARSERDEBUG
  fprintf(stderr,"name of %d orig_var %s\n",__theMDP->numorigvars,__theMDP->orig_vars[__theMDP->numorigvars].name);
#endif
}
    break;

  case 9:
/* Line 1787 of yacc.c  */
#line 155 "mdpparser.y"
    { 
  addVar(__theMDP->orig_vars,__theMDP->numorigvars,(yyvsp[(2) - (2)].nme));
}
    break;

  case 10:
/* Line 1787 of yacc.c  */
#line 157 "mdpparser.y"
    { 
  addVar(__theMDP->orig_vars,__theMDP->numorigvars,(yyvsp[(1) - (2)].nme));
  addVar(__theMDP->orig_vars,__theMDP->numorigvars,(yyvsp[(2) - (2)].nme));
}
    break;

  case 12:
/* Line 1787 of yacc.c  */
#line 163 "mdpparser.y"
    {
#ifdef USE_POMDP
  // now we have the observations, so allocate for the ObsFun Diagrams
  // These are diagrams with a single multi-valued primed variable in them,
  // which is translated into log2(nvals) binary variables
  // so we only need numorigvars NewPrime diagrams for each action
  int i;
  __theMDP->ObsFun = (DdNode ***)malloc(MAXACT*(sizeof(DdNode **)));
#ifdef PARSERDEBUG
  fprintf(stderr,"allocating for %d original observations in MDP\n", __theMDP->numorigobs);
#endif
  for(i=0;i<MAXACT;i++)
    __theMDP->ObsFun[i] = (DdNode **)malloc(__theMDP->numorigobs*(sizeof(DdNode*)));
  __theMDP->computeTotalNumObs();
  // also allocate for the beliefSTate
  __theMDP->beliefState = (DdNode **) malloc(__theMDP->numorigvars*sizeof(DdNode *));
  // build the good-state ADDs
  // do we have to do this for observations?
  //buildGoodStateADDs(__theMDP->orig_vars,__theMDP->vars,__theMDP->prime_vars,__theMDP->numorigvars);
#endif  
}
    break;

  case 15:
/* Line 1787 of yacc.c  */
#line 188 "mdpparser.y"
    {
#ifdef USE_POMDP
  __theMDP->newADDObs();
  __theMDP->numorigobs++;
#endif
}
    break;

  case 16:
/* Line 1787 of yacc.c  */
#line 195 "mdpparser.y"
    {
#ifdef USE_POMDP
  addVar(__theMDP->orig_obs,__theMDP->numorigobs,(yyvsp[(2) - (2)].nme));
#endif
}
    break;

  case 17:
/* Line 1787 of yacc.c  */
#line 199 "mdpparser.y"
    { 
#ifdef USE_POMDP
  addVar(__theMDP->orig_obs,__theMDP->numorigobs,(yyvsp[(1) - (2)].nme));
  addVar(__theMDP->orig_obs,__theMDP->numorigobs,(yyvsp[(2) - (2)].nme));
#endif
}
    break;

  case 18:
/* Line 1787 of yacc.c  */
#line 206 "mdpparser.y"
    {
#ifdef USE_POMDP
  __theMDP->orig_obs[__theMDP->numorigobs].name = strdup((yyvsp[(1) - (1)].nme));
#ifdef PARSERDEBUG
  fprintf(stderr,"name of %d orig_obs %s\n",__theMDP->numorigobs,__theMDP->orig_obs[__theMDP->numorigobs].name);
#endif
#endif
}
    break;

  case 20:
/* Line 1787 of yacc.c  */
#line 216 "mdpparser.y"
    {
  // initialize initBelief
#ifdef USE_POMDP
  __theMDP->initBelief = One;
  Cudd_Ref(__theMDP->initBelief);
#ifdef PARSERDEBUG
  fprintf(stderr,"initialized initBelief to One\n");
#endif
#endif
}
    break;

  case 21:
/* Line 1787 of yacc.c  */
#line 225 "mdpparser.y"
    {
  pvcount = 0;
}
    break;

  case 24:
/* Line 1787 of yacc.c  */
#line 234 "mdpparser.y"
    {
  // find this variable in orig_vars 
  if ((curr_primed_ovar = findOVar(__theMDP->orig_vars,__theMDP->numorigvars,(yyvsp[(1) - (1)].nme))) < 0) {
    fprintf(stderr,"could not find original variable %s\n",(yyvsp[(1) - (1)].nme));
    exit(0);
  }
#ifdef PARSERDEBUG
  fprintf(stderr,"found belief over variable %s\n",__theMDP->orig_vars[curr_primed_ovar].name);
#endif
  level = 0;
}
    break;

  case 25:
/* Line 1787 of yacc.c  */
#line 244 "mdpparser.y"
    {
#ifdef USE_POMDP
#ifdef PARSERDEBUG
  fprintf(stderr,"parsed belief over variable %s\n",__theMDP->orig_vars[curr_primed_ovar].name);
#ifdef PARSERDEBUG2
  Cudd_PrintDebug(gbm,(yyvsp[(3) - (3)].dnode),4,100);
#endif
#endif
  __theMDP->beliefState[curr_primed_ovar]  = (yyvsp[(3) - (3)].dnode);
  Cudd_Ref(__theMDP->beliefState[curr_primed_ovar]);
  DdNode *temp1 =  Cudd_addApply(gbm,Cudd_addTimes,__theMDP->initBelief,(yyvsp[(3) - (3)].dnode));
  Cudd_Ref(temp1);
  Cudd_RecursiveDeref(gbm,(yyvsp[(3) - (3)].dnode));
  Cudd_RecursiveDeref(gbm,__theMDP->initBelief);
  __theMDP->initBelief = temp1;
#ifdef PARSERDEBUG
  fprintf(stderr,"Initial Belief is currently\n");
#ifdef PARSERDEBUG2
  Cudd_PrintDebug(gbm,__theMDP->initBelief,4,100);
#endif
#endif
#endif
}
    break;

  case 26:
/* Line 1787 of yacc.c  */
#line 268 "mdpparser.y"
    {
  unnormalized = false;
}
    break;

  case 27:
/* Line 1787 of yacc.c  */
#line 270 "mdpparser.y"
    {
  unnormalized = true;
}
    break;

  case 32:
/* Line 1787 of yacc.c  */
#line 281 "mdpparser.y"
    {
  conjlevel = 0;
}
    break;

  case 33:
/* Line 1787 of yacc.c  */
#line 283 "mdpparser.y"
    {
  dds[numdds] = (yyvsp[(3) - (4)].dnode);
  Cudd_Ref(dds[numdds]);
  Cudd_RecursiveDeref(gbm,(yyvsp[(3) - (4)].dnode));
  // store in hash table
  ddindices[ddnames[numdds]] = numdds;
#ifdef PARSERDEBUG
    fprintf(stderr,"************************** parsed the dd %s at %d*****\n",ddnames[numdds],ddindices[ddnames[numdds]]);
#ifdef PARSERDEBUG2
    Cudd_PrintDebug(gbm,dds[numdds],4,100);
#endif
#endif
  numdds++;
  doing_dd = false;
}
    break;

  case 34:
/* Line 1787 of yacc.c  */
#line 300 "mdpparser.y"
    {
  ddnames[numdds] = strdup((yyvsp[(2) - (2)].nme));
  doing_dd = true;
}
    break;

  case 35:
/* Line 1787 of yacc.c  */
#line 306 "mdpparser.y"
    {  
  if (pvcount != __theMDP->numorigvars)
    error("missing primed variable in action",__theMDP->actionlist[__theMDP->numactions].name);
  pvcount = 0;
  __theMDP->numactions++;
  doing_action = false;
}
    break;

  case 36:
/* Line 1787 of yacc.c  */
#line 315 "mdpparser.y"
    {
  /* name of action is $2, increment action counter */
#ifdef PARSERDEBUG
  fprintf(stderr,"------------- action name %s\n",(yyvsp[(2) - (2)].nme)); 
#endif
  __theMDP->actionlist[__theMDP->numactions].name = strdup((yyvsp[(2) - (2)].nme));
  doing_action = true;
}
    break;

  case 37:
/* Line 1787 of yacc.c  */
#line 324 "mdpparser.y"
    {
  // the action cost is zero
  setActionCost(__theMDP->actionCost,__theMDP->numactions,0.0);
}
    break;

  case 38:
/* Line 1787 of yacc.c  */
#line 327 "mdpparser.y"
    {
  setActionCost(__theMDP->actionCost,__theMDP->numactions,(yyvsp[(1) - (1)].val));
}
    break;

  case 39:
/* Line 1787 of yacc.c  */
#line 332 "mdpparser.y"
    {
  // the action cost is zero
  //fprintf(stderr,"action cost 0\n");
  setActionCost(__theMDP->actionCost, __theMDP->numactions,0.0);
  __theMDP->actionCostNoDummy[__theMDP->numactions] =   __theMDP->actionCost[__theMDP->numactions]; 
  Cudd_Ref(__theMDP->actionCostNoDummy[__theMDP->numactions]);
}
    break;

  case 40:
/* Line 1787 of yacc.c  */
#line 338 "mdpparser.y"
    { pvindex = 0; level = 0; doing_reward = true; }
    break;

  case 41:
/* Line 1787 of yacc.c  */
#line 338 "mdpparser.y"
    {
  // assign the action cost
  __theMDP->actionCost[__theMDP->numactions] = (yyvsp[(3) - (3)].dnode);
  Cudd_Ref(__theMDP->actionCost[__theMDP->numactions]);
  Cudd_RecursiveDeref(gbm,(yyvsp[(3) - (3)].dnode));

  __theMDP->actionCostNoDummy[__theMDP->numactions] =   __theMDP->actionCost[__theMDP->numactions]; 
  Cudd_Ref(__theMDP->actionCostNoDummy[__theMDP->numactions]);

  // TEMPORARILY REMOVE THIS FOR PRINTING OUT
#ifdef COSTSPRIMED
  removeAllDummysp(__theMDP->actionCost+__theMDP->numactions,__theMDP->numorigvars);
#else
  removeAllDummys(__theMDP->actionCost+__theMDP->numactions,__theMDP->numorigvars);
#endif
  // multiply by -1 (since its a cost)
  Pair *value = new Pair(-1.0);
  DdNode *temp = Cudd_addConst(gbm,value);
  Cudd_Ref(temp);
  DdNode *temp2 = Cudd_addApply(gbm,Cudd_addTimes,__theMDP->actionCost[__theMDP->numactions],temp);
  Cudd_Ref(temp2);
  Cudd_RecursiveDeref(gbm,temp);
  Cudd_RecursiveDeref(gbm,__theMDP->actionCost[__theMDP->numactions]);
  __theMDP->actionCost[__theMDP->numactions] = temp2;
  doing_reward = false;
}
    break;

  case 43:
/* Line 1787 of yacc.c  */
#line 366 "mdpparser.y"
    {
  doing_obs = true;
}
    break;

  case 44:
/* Line 1787 of yacc.c  */
#line 368 "mdpparser.y"
    {
  doing_obs = false;
}
    break;

  case 48:
/* Line 1787 of yacc.c  */
#line 376 "mdpparser.y"
    {
#ifdef USE_POMDP
  // find this variable in orig_obs
  if ((curr_primed_ovar = findOVar(__theMDP->orig_obs,__theMDP->numorigobs,(yyvsp[(1) - (1)].nme))) < 0) {
    fprintf(stderr,"could not find original observation %s\n",(yyvsp[(1) - (1)].nme));
    exit(0);
  }
#endif
}
    break;

  case 49:
/* Line 1787 of yacc.c  */
#line 384 "mdpparser.y"
    {
#ifdef USE_POMDP
  // observation function in theadd ($3)
  __theMDP->ObsFun[__theMDP->numactions][curr_primed_ovar] = (yyvsp[(3) - (3)].dnode);
  Cudd_Ref(__theMDP->ObsFun[__theMDP->numactions][curr_primed_ovar]);
#endif
}
    break;

  case 52:
/* Line 1787 of yacc.c  */
#line 395 "mdpparser.y"
    {
  /* primed variable is $1 theadd is its action tree*/
#ifdef PARSERDEBUG
   fprintf(stderr,"primed var %s\n",(yyvsp[(1) - (1)].nme));
#endif
  // find this variable in orig_vars 
  if ((curr_primed_ovar = findOVar(__theMDP->orig_vars,__theMDP->numorigvars,(yyvsp[(1) - (1)].nme))) < 0) {
    fprintf(stderr,"could not find original variable %s\n",(yyvsp[(1) - (1)].nme));
    exit(0);
  }
  level = 0;
  conjlevel = 0;
  // the corresponding binary variables are vars[i] and prime_vars[i]
  // where i:orig_vars[curr_primed_ovar].var1index...orig_vars[curr_primed_ovar].var1index+orig_vars[curr_primed_ovar].nbvars-1
}
    break;

  case 53:
/* Line 1787 of yacc.c  */
#line 409 "mdpparser.y"
    {
  // the new prime diagram is that returned in theadd
  // its $3 because the above action is also counted
  removeDummyp(&((yyvsp[(3) - (3)].dnode)),curr_primed_ovar);
  DdNode *tmp2;
  DdNode *temp = sumOutPrime((yyvsp[(3) - (3)].dnode),__theMDP->orig_vars+curr_primed_ovar,__theMDP->prime_vars);
  Cudd_Ref(temp);


  if (unnormalized) {
    // renormalize - this can add in Dummy states!
#ifdef PARSERDEBUG2
    fprintf(stderr,"dividing this:\n");
    pdd((yyvsp[(3) - (3)].dnode));
    fprintf(stderr,"by this (previous summed over %d):\n",curr_primed_ovar);
    pdd(temp);
#endif    
    tmp2 = Cudd_addApply(gbm,Cudd_addDivide,(yyvsp[(3) - (3)].dnode),temp);

    Cudd_Ref(tmp2);

    Cudd_RecursiveDeref(gbm,temp);
    Cudd_RecursiveDeref(gbm,(yyvsp[(3) - (3)].dnode));
    temp = tmp2;
    // remove dummy states
    //removeAllDummys(&temp,__theMDP->numorigvars);
    //removeAllDummysp(&temp,__theMDP->numorigvars);


#ifdef REDUCE_PRECISION
    // reduce the precision of the final diagram
    Pair *prec = new Pair();
    (*prec).set(REDUCE_PRECISION);
    DdNode *predd = Cudd_addConst(gbm,prec);
    Cudd_Ref(predd);
    tmp2 = Cudd_addApply(gbm,reducePrecision,temp,predd);
    Cudd_Ref(tmp2);
    Cudd_RecursiveDeref(gbm,temp);
    Cudd_RecursiveDeref(gbm,predd);
    temp = tmp2;
    // and normalize again
    tmp2 = sumOutPrime(temp,__theMDP->orig_vars+curr_primed_ovar,__theMDP->prime_vars);
    Cudd_Ref(tmp2);
    tmp1 = Cudd_addApply(gbm,Cudd_addDivide,temp,tmp2);
    Cudd_Ref(tmp1);
    Cudd_RecursiveDeref(gbm,temp);
    Cudd_RecursiveDeref(gbm,tmp2);
    temp = tmp1;
#endif



    __theMDP->NewPrime[__theMDP->numactions][curr_primed_ovar] = temp;
    Cudd_Ref(__theMDP->NewPrime[__theMDP->numactions][curr_primed_ovar]);
    
    

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
      fprintf(stderr,"CPT for primed variable %s is not normalized:\n",(yyvsp[(1) - (3)].nme));
      Cudd_PrintMinterm(gbm,(yyvsp[(3) - (3)].dnode));
      fprintf(stderr,"curr_primed_ovar: %d\n summed diagram:\n",curr_primed_ovar);
      Cudd_PrintDebug(gbm,temp,4,100);
      exit(-1);
    } else {
      __theMDP->NewPrime[__theMDP->numactions][curr_primed_ovar] = (yyvsp[(3) - (3)].dnode);
      Cudd_Ref(__theMDP->NewPrime[__theMDP->numactions][curr_primed_ovar]);
      Cudd_RecursiveDeref(gbm,(yyvsp[(3) - (3)].dnode));
    }
  } 
  Cudd_RecursiveDeref(gbm,temp);
  // removeDummy(&(__theMDP->NewPrime[__theMDP->numactions][curr_primed_ovar]),curr_primed_ovar);


#ifdef PARSERDEBUG
  fprintf(stderr,"the add for %d action %d ovar is \n",__theMDP->numactions, curr_primed_ovar);
#ifdef PARSERDEBUG2
  Cudd_PrintDebug(gbm,__theMDP->NewPrime[__theMDP->numactions][curr_primed_ovar],4,100);
#endif
#endif
  pvcount++;
}
    break;

  case 54:
/* Line 1787 of yacc.c  */
#line 507 "mdpparser.y"
    {
  oldconj[conjlevel] = conj;
#ifdef PARSERDEBUG
  fprintf(stderr,"saving oldconj %d at level %d before setting conj to 1\n",oldconj[conjlevel],conjlevel);
#endif
  conj = 1;
  conjlevel++;
}
    break;

  case 55:
/* Line 1787 of yacc.c  */
#line 514 "mdpparser.y"
    {
  (yyval.dnode) = (yyvsp[(4) - (5)].dnode);
  Cudd_Ref((yyval.dnode));
  Cudd_RecursiveDeref(gbm,(yyvsp[(4) - (5)].dnode));

#ifdef PARSERDEBUG
  fprintf(stderr,"resetting conj to oldconj at level %d : %d\n",conjlevel-1,oldconj[conjlevel-1]);
#endif
  conjlevel--;
  conj = oldconj[conjlevel];
}
    break;

  case 56:
/* Line 1787 of yacc.c  */
#line 525 "mdpparser.y"
    {
  oldconj[conjlevel] = conj;
#ifdef PARSERDEBUG
  fprintf(stderr,"saving oldconj %d at level %d before setting conj to 1\n",oldconj[conjlevel],conjlevel);
#endif
  conj = 0;
  conjlevel++;
}
    break;

  case 57:
/* Line 1787 of yacc.c  */
#line 532 "mdpparser.y"
    {
  (yyval.dnode) = (yyvsp[(4) - (5)].dnode);
  Cudd_Ref((yyval.dnode));
  Cudd_RecursiveDeref(gbm,(yyvsp[(4) - (5)].dnode));
#ifdef PARSERDEBUG
  fprintf(stderr,"resetting conj to oldconj at level %d : %d\n",conjlevel-1,oldconj[conjlevel-1]);
#endif
  conjlevel--;
  conj = oldconj[conjlevel];
}
    break;

  case 58:
/* Line 1787 of yacc.c  */
#line 542 "mdpparser.y"
    {
  oldconj[conjlevel] = conj;
#ifdef PARSERDEBUG
  fprintf(stderr,"saving oldconj %d at level %d before setting conj to 1\n",oldconj[conjlevel],conjlevel);
#endif
  conj = 2;
  conjlevel++;
}
    break;

  case 59:
/* Line 1787 of yacc.c  */
#line 549 "mdpparser.y"
    {
  //renormalize $6 over variable $3
  if ((renormvar = findOVar(__theMDP->orig_vars,__theMDP->numorigvars,(yyvsp[(3) - (7)].nme))) < 0) {
    fprintf(stderr,"could not find original variable %s\n",(yyvsp[(3) - (7)].nme));
    exit(0);
  }
  DdNode *temp = NULL;
  __theMDP->normalizeFunction((yyvsp[(6) - (7)].dnode),&temp,renormvar);
  (yyval.dnode) = temp;
  Cudd_Ref((yyval.dnode));
  Cudd_RecursiveDeref(gbm,temp);
  Cudd_RecursiveDeref(gbm,(yyvsp[(6) - (7)].dnode));
#ifdef PARSERDEBUG
  fprintf(stderr,"resetting conj to oldconj at level %d : %d\n",conjlevel-1,oldconj[conjlevel-1]);
#endif
  conjlevel--;
  conj = oldconj[conjlevel];
}
    break;

  case 60:
/* Line 1787 of yacc.c  */
#line 567 "mdpparser.y"
    {
  /* current root node is $2 */
#ifdef PARSERDEBUG
  fprintf(stderr,"********************************  root node %s\n",(yyvsp[(2) - (2)].nme));
#endif
  if ((curr_ovar = findOVar(__theMDP->orig_vars,__theMDP->numorigvars,(yyvsp[(2) - (2)].nme))) < 0) 
    error("could not find variable",(yyvsp[(2) - (2)].nme));
  levelIndex[level] = curr_ovar;
  branchCount[level] = 0;
}
    break;

  case 61:
/* Line 1787 of yacc.c  */
#line 577 "mdpparser.y"
    {
  (yyval.dnode) = (yyvsp[(4) - (5)].dnode);
  Cudd_Ref((yyval.dnode));
  Cudd_RecursiveDeref(gbm,(yyvsp[(4) - (5)].dnode));
#ifdef PARSERDEBUG
  fprintf(stderr,"***************************** the add is\n");
#ifdef PARSERDEBUG2
  Cudd_PrintDebug(gbm,(yyval.dnode),4,100);
#endif
#endif
  if ((curr_ovar = findOVar(__theMDP->orig_vars,__theMDP->numorigvars,(yyvsp[(2) - (5)].nme))) < 0) 
    error("could not find variable",(yyvsp[(2) - (5)].nme));
  if (branchCount[level] != __theMDP->orig_vars[curr_ovar].nvals) {
    fprintf(stderr,"%d branches missing from variable %s primed variable %s in action %s\n",
	    __theMDP->orig_vars[curr_ovar].nvals-branchCount[level],
	    __theMDP->orig_vars[curr_ovar].name,__theMDP->orig_vars[curr_primed_ovar].name,__theMDP->actionlist[__theMDP->numactions].name);
    exit(-1);
  }
}
    break;

  case 63:
/* Line 1787 of yacc.c  */
#line 597 "mdpparser.y"
    {

  // check to make sure it matches the primed variables number
  if (!doing_reward) {
    if (doing_obs) {
#ifdef USE_POMDP
      if (pvindex != __theMDP->orig_obs[curr_primed_ovar].nvals) 
	error("Some cpt values missing from observation",__theMDP->orig_obs[curr_primed_ovar].name);
#endif
    } else {
      if (pvindex != __theMDP->orig_vars[curr_primed_ovar].nvals) 
	error("Some cpt values missing from primed variable",__theMDP->orig_vars[curr_primed_ovar].name);
    }
    (yyval.dnode) = (yyvsp[(2) - (3)].dnode);
  } else if (doing_reward && pvindex != 1) {
    error("only one reward value per state","please!");
  } else {
    (yyval.dnode) = (yyvsp[(2) - (3)].dnode);
  }
  Cudd_Ref((yyval.dnode));
  Cudd_RecursiveDeref(gbm,(yyvsp[(2) - (3)].dnode));

  // reset for next time
  pvindex = 0;
  pvsum = 0.0;
}
    break;

  case 64:
/* Line 1787 of yacc.c  */
#line 625 "mdpparser.y"
    {
#ifdef PARSERDEBUG
  fprintf(stderr,"parsing root node name %s\n",(yyvsp[(1) - (1)].nme));
#endif
  (yyval.nme) = (yyvsp[(1) - (1)].nme);
  levelPrime[level] = false;
}
    break;

  case 65:
/* Line 1787 of yacc.c  */
#line 632 "mdpparser.y"
    {
#ifdef PARSERDEBUG
  fprintf(stderr,"saw a prime on variable %s\n",(yyvsp[(1) - (2)].nme));
#endif
  (yyval.nme) = (yyvsp[(1) - (2)].nme);
  levelPrime[level] = true;
}
    break;

  case 66:
/* Line 1787 of yacc.c  */
#line 641 "mdpparser.y"
    {
  int i;
  DdNode *temp;
  if (conj == 1) {
#ifdef PARSERDEBUG
    fprintf(stderr,"multiplying these two\n");
#endif    
    temp = Cudd_addApply(gbm,Cudd_addTimes,(yyvsp[(1) - (2)].dnode),(yyvsp[(2) - (2)].dnode));
  } else {
#ifdef PARSERDEBUG
    fprintf(stderr,"adding these two\n");
#endif    
    temp = Cudd_addApply(gbm,Cudd_addPlus,(yyvsp[(1) - (2)].dnode),(yyvsp[(2) - (2)].dnode));
  }
#ifdef PARSERDEBUG
#ifdef PARSERDEBUG2
  Cudd_PrintDebug(gbm,(yyvsp[(1) - (2)].dnode),1,2);
  Cudd_PrintDebug(gbm,(yyvsp[(2) - (2)].dnode),1,2);
#endif
#endif
  Cudd_Ref(temp);
  Cudd_RecursiveDeref(gbm,(yyvsp[(1) - (2)].dnode));
  Cudd_RecursiveDeref(gbm,(yyvsp[(2) - (2)].dnode));
  (yyval.dnode) = temp;
}
    break;

  case 67:
/* Line 1787 of yacc.c  */
#line 666 "mdpparser.y"
    {
  (yyval.dnode) = (yyvsp[(1) - (1)].dnode);
  Cudd_Ref((yyval.dnode));
  Cudd_RecursiveDeref(gbm,(yyvsp[(1) - (1)].dnode));
}
    break;

  case 68:
/* Line 1787 of yacc.c  */
#line 673 "mdpparser.y"
    {
  (yyval.dnode) = Cudd_addApply(gbm,Cudd_addPlus,(yyvsp[(1) - (2)].dnode),(yyvsp[(2) - (2)].dnode));
  Cudd_Ref((yyval.dnode));
  Cudd_RecursiveDeref(gbm,(yyvsp[(1) - (2)].dnode));
  Cudd_RecursiveDeref(gbm,(yyvsp[(2) - (2)].dnode));
}
    break;

  case 69:
/* Line 1787 of yacc.c  */
#line 679 "mdpparser.y"
    {
  (yyval.dnode) = (yyvsp[(1) - (1)].dnode);
  Cudd_Ref((yyval.dnode));
  Cudd_RecursiveDeref(gbm,(yyvsp[(1) - (1)].dnode));
}
    break;

  case 70:
/* Line 1787 of yacc.c  */
#line 691 "mdpparser.y"
    {
#ifdef PARSERDEBUG
  fprintf(stderr,"branch %s level increasing to %d\n",(yyvsp[(2) - (2)].nme),level+1);
#endif
  level++;
}
    break;

  case 71:
/* Line 1787 of yacc.c  */
#line 696 "mdpparser.y"
    {
  level--;
#ifdef PARSERDEBUG
  fprintf(stderr,"parsed branch %s\n",(yyvsp[(2) - (5)].nme)); 
#ifdef PARSERDEBUG2
  Cudd_PrintDebug(gbm,(yyvsp[(4) - (5)].dnode),4,100);
#endif
#endif
  if ((curr_oval = findOVal(__theMDP->orig_vars, levelIndex[level],(yyvsp[(2) - (5)].nme))) < 0) {
    fprintf(stderr,"could not find value %s in variable %s action %s\n",
	    (yyvsp[(2) - (5)].nme),__theMDP->orig_vars[levelIndex[level]].name,__theMDP->actionlist[__theMDP->numactions].name);
    exit(-1);
  }
  // need to check if curr_ovar is a primed variable here
  if (!levelPrime[level] && !doing_obs)
    (yyval.dnode) = buildCubeCPT(__theMDP->vars,curr_oval,__theMDP->orig_vars[levelIndex[level]],(yyvsp[(4) - (5)].dnode));
  else
    (yyval.dnode) = buildCubeCPT(__theMDP->prime_vars,curr_oval,__theMDP->orig_vars[levelIndex[level]],(yyvsp[(4) - (5)].dnode));
  Cudd_RecursiveDeref(gbm,(yyvsp[(4) - (5)].dnode));
#ifdef PARSERDEBUG
  fprintf(stderr,"branch CPT is\n");
#ifdef PARSERDEBUG2
  Cudd_PrintDebug(gbm,(yyval.dnode),4,100);
#endif
#endif
  branchCount[level]++;
}
    break;

  case 72:
/* Line 1787 of yacc.c  */
#line 727 "mdpparser.y"
    {
  (yyval.dnode) = Cudd_addApply(gbm,Cudd_addPlus,(yyvsp[(1) - (2)].dnode),(yyvsp[(2) - (2)].dnode));
  Cudd_Ref((yyval.dnode));
  Cudd_RecursiveDeref(gbm,(yyvsp[(1) - (2)].dnode));
  Cudd_RecursiveDeref(gbm,(yyvsp[(2) - (2)].dnode));
}
    break;

  case 73:
/* Line 1787 of yacc.c  */
#line 733 "mdpparser.y"
    {
  (yyval.dnode) = Cudd_addApply(gbm,Cudd_addPlus,(yyvsp[(1) - (2)].dnode),(yyvsp[(2) - (2)].dnode));
  Cudd_Ref((yyval.dnode));
  Cudd_RecursiveDeref(gbm,(yyvsp[(1) - (2)].dnode));
  Cudd_RecursiveDeref(gbm,(yyvsp[(2) - (2)].dnode));
}
    break;

  case 74:
/* Line 1787 of yacc.c  */
#line 741 "mdpparser.y"
    {
  Pair *value = new Pair();
#ifdef REDUCE_PRECISION
  //reduce the precision of the inputs 
  (*value).set(((double) ((int) ((yyvsp[(2) - (3)].val)*REDUCE_PRECISION)))/(1.0*REDUCE_PRECISION));
#else
  (*value).set((yyvsp[(2) - (3)].val));
#endif
  //fprintf(stderr,"setting %f to %f\n",$2,((double) ((int) ($2*100)))/100.0);
  (yyval.dnode) = Cudd_addConst(gbm,value);
  Cudd_Ref((yyval.dnode));
}
    break;

  case 75:
/* Line 1787 of yacc.c  */
#line 753 "mdpparser.y"
    {
  // search for NAME in list of defined DDs
  // by looking in hash table
#ifdef PARSERDEBUG
  fprintf(stderr,"********************** looking for dd %s\n",(yyvsp[(2) - (3)].nme));
#endif
  map<const char*, int, ltstr>::iterator ddit = ddindices.find((yyvsp[(2) - (3)].nme));
  if (ddit==ddindices.end()) {
    fprintf(stderr,"could not find defined DD %s in list\n",(yyvsp[(2) - (3)].nme));
    exit(-1);
  } else {
#ifdef PARSERDEBUG
    fprintf(stderr,"************************** Found defined DD %s in list at %d\n",(yyvsp[(2) - (3)].nme),(*ddit).second);
#ifdef PARSERDEBUG2
    Cudd_PrintDebug(gbm,dds[(*ddit).second],1,2);
#endif
#endif
    (yyval.dnode) = dds[(*ddit).second];
    Cudd_Ref((yyval.dnode));
  }
}
    break;

  case 76:
/* Line 1787 of yacc.c  */
#line 776 "mdpparser.y"
    {

  Pair *value = new Pair();
#ifdef REDUCE_PRECISION
  //reduce the precision of the inputs 
  (*value).set(((double) ((int) ((yyvsp[(1) - (1)].val)*REDUCE_PRECISION)))/(1.0*REDUCE_PRECISION));
#else
  (*value).set((yyvsp[(1) - (1)].val));
#endif
  //fprintf(stderr,"setting %f to %f\n",$1,((double) ((int) ($1*100)))/100.0);


  DdNode *temp = Cudd_addConst(gbm,value);
  delete value;
  Cudd_Ref(temp);
  if (!doing_reward) {
#ifdef PARSERDEBUG
    fprintf(stderr,"found a constant %f at variable %s\n",(yyvsp[(1) - (1)].val),__theMDP->orig_vars[curr_primed_ovar].name);
#endif
    // buildCubeCPT does its own reffing
    if (doing_obs) {
#ifdef USE_POMDP
      (yyval.dnode) = buildCubeCPT(__theMDP->obs,pvindex,__theMDP->orig_obs[curr_primed_ovar],temp);
#endif
    } else {
      (yyval.dnode) = buildCubeCPT(__theMDP->prime_vars,pvindex,__theMDP->orig_vars[curr_primed_ovar],temp);
    }
    Cudd_RecursiveDeref(gbm,temp);
  } else {
    //fprintf(stderr," found a constant %f pvindex %d\n",$1,pvindex);
    (yyval.dnode) = temp;
  }    
  pvindex++;
  pvsum += (yyvsp[(1) - (1)].val);
}
    break;

  case 77:
/* Line 1787 of yacc.c  */
#line 814 "mdpparser.y"
    { 
#ifdef PARSERDEBUG
  fprintf(stderr, "doing reward\n");
#endif
  pvindex = 0; level = 0; doing_reward = true; 
}
    break;

  case 78:
/* Line 1787 of yacc.c  */
#line 819 "mdpparser.y"
    {
  __theMDP->RewardD = (yyvsp[(3) - (3)].dnode);
  Cudd_Ref(__theMDP->RewardD);
  Cudd_RecursiveDeref(gbm,(yyvsp[(3) - (3)].dnode));

  __theMDP->RewardDNoDummy =   __theMDP->RewardD; 
  Cudd_Ref(__theMDP->RewardDNoDummy);

  removeAllDummys(&(__theMDP->RewardD),__theMDP->numorigvars);
}
    break;

  case 79:
/* Line 1787 of yacc.c  */
#line 831 "mdpparser.y"
    { 
  Pair discountPair;
  discountPair.set((yyvsp[(2) - (2)].val));
  __theMDP->discount = Cudd_addConst(gbm,&discountPair);
  Cudd_Ref(__theMDP->discount);
}
    break;

  case 80:
/* Line 1787 of yacc.c  */
#line 836 "mdpparser.y"
    {
  doing_reward = true;
}
    break;

  case 81:
/* Line 1787 of yacc.c  */
#line 838 "mdpparser.y"
    {
  __theMDP->discount = (yyvsp[(3) - (3)].dnode);
  Cudd_Ref(__theMDP->discount);
}
    break;

  case 84:
/* Line 1787 of yacc.c  */
#line 847 "mdpparser.y"
    { 
  tolerance = (yyvsp[(2) - (2)].val); 
  __theMDP->horizon = -1.0;
}
    break;

  case 85:
/* Line 1787 of yacc.c  */
#line 853 "mdpparser.y"
    {
  tolerance = 0.0;
  __theMDP->horizon = int((yyvsp[(2) - (2)].val));
}
    break;


/* Line 1787 of yacc.c  */
#line 2552 "mdpparser.tab.c"
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

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

#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}


/* Line 2050 of yacc.c  */
#line 858 "mdpparser.y"


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
#ifdef PARSERDEBUG2
    //Cudd_PrintDebug(gbm,goodState[v],4,100);
#endif

  }
}

void setActionCost(DdNode **ac, int nac, double val) {
  //  fprintf(stderr,"Action cost is %f for action %d\n",val,nac);
  Pair *value = new Pair();
  (*value).set((double)(-1.0*fabs(val)));
  ac[nac] = Cudd_addConst(gbm,value);
  Cudd_Ref(ac[nac]);
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
#ifdef PARSERDEBUG2
  //Cudd_PrintDebug(gbm,branch,4,100);
#endif

  DdNode *cubecpt = Cudd_addApply(gbm,Cudd_addTimes,branch,add);
  Cudd_Ref(cubecpt);
  Cudd_RecursiveDeref(gbm,branch);

  //fprintf(stderr,"built cube:\n");
#ifdef PARSERDEBUG2
  //Cudd_PrintDebug(gbm,cubecpt,4,100);
#endif

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
  Pair *value = new Pair(1.0/pvsum);
  DdNode *temp = Cudd_addConst(gbm,value);
  Cudd_Ref(temp);
  DdNode *temp2 = Cudd_addApply(gbm,Cudd_addTimes,cube,temp);
  Cudd_RecursiveDeref(gbm,temp);
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

