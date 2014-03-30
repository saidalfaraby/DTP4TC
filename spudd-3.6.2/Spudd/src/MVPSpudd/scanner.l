%{
#include <stdlib.h>
#include <string.h>
#include "parser.tab.h"

#define DEBUG 0

%}

%option noyywrap
%option yylineno

ws  ([ \t\n])+
letter [A-Za-z_]
digit  [0-9]
id     {letter}({letter}|{digit})*
number (\-)?{digit}*(\.{digit}+)?
intgr  {digit}+
op     \(
cp     \)
osb    \[
csb    \]
disj   \+
conj   \*

%%
"//"[^}\n]*\n          /* eat up one line comments */
{ws}      ;
debugon   {/*yydebug = 1;*/}
debugoff  {/*yydebug = 0;*/}
variables {if (DEBUG) ECHO; return(VARIABLES);}
reward    {if (DEBUG) ECHO; return(REWARD);}
value     {if (DEBUG) ECHO; return(VALUE);}
discount  {if (DEBUG) ECHO; return(DISCOUNT);}
horizon   {if (DEBUG) ECHO; return(HORIZON);}
tolerance {if (DEBUG) ECHO; return(TOLERANCE);}
action    {if (DEBUG) ECHO; return(ACTION);}
endaction {if (DEBUG) ECHO; return(ENDACTION);}
{op}      {if (DEBUG) ECHO; return(OPP);}
{cp}      {if (DEBUG) ECHO; return(CLP);}
{osb}     {if (DEBUG) ECHO; return(OSB);}
{csb}     {if (DEBUG) ECHO; return(CSB);}
{id}      {if (DEBUG) ECHO; yylval.nme = strdup(yytext); return(NAME);}
{number}  {if (DEBUG) ECHO; yylval.val = atof(yytext); return(REAL);}
{intgr}   {if (DEBUG) ECHO; yylval.ival = atoi(yytext); return(INTGR);}
{disj}    {if (DEBUG) ECHO; return(DISJ);}
{conj}    {if (DEBUG) ECHO; return(CONJ);}
%%

