%{
#include <stdlib.h>
#include <string.h>
#include "mdpparser.tab.h"

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
renorm \#
prime  \'

%%
"//"[^}\n]*\n          /* eat up one line comments */
{ws}      ;
debugon   {/*yydebug = 1;*/}
debugoff  {/*yydebug = 0;*/}
variables {if (DEBUG) ECHO; return(VARIABLES);}
observations {if (DEBUG) ECHO; return(OBSERVATIONS);}
observe   {if (DEBUG) ECHO; return(OBSERVE);}
endobserve {if (DEBUG) ECHO; return(ENDOBSERVE);}
belief   {if (DEBUG) ECHO; return(BELIEF);}
endbelief   {if (DEBUG) ECHO; return(ENDBELIEF);}
reward    {if (DEBUG) ECHO; return(REWARD);}
cost      {if (DEBUG) ECHO; return(COST); }
value     {if (DEBUG) ECHO; return(VALUE);}
discount  {if (DEBUG) ECHO; return(DISCOUNT);}
horizon   {if (DEBUG) ECHO; return(HORIZON);}
tolerance {if (DEBUG) ECHO; return(TOLERANCE);}
action    {if (DEBUG) ECHO; return(ACTION);}
endaction {if (DEBUG) ECHO; return(ENDACTION);}
dd        {if (DEBUG) ECHO; return(STARTDD);}
enddd     {if (DEBUG) ECHO; return(ENDDD);}
unnormalized {if (DEBUG) ECHO; return(UNNORM);}
unnormalised {if (DEBUG) ECHO; return(UNNORM);}
{op}      {if (DEBUG) ECHO; return(OPP);}
{cp}      {if (DEBUG) ECHO; return(CLP);}
{osb}     {if (DEBUG) ECHO; return(OSB);}
{csb}     {if (DEBUG) ECHO; return(CSB);}
{id}      {if (DEBUG) ECHO; yylval.nme = strdup(yytext); return(NAME);}
{number}  {if (DEBUG) ECHO; yylval.val = atof(yytext); return(REAL);}
{intgr}   {if (DEBUG) ECHO; yylval.ival = atoi(yytext); return(INTGR);}
{disj}    {if (DEBUG) ECHO; return(DISJ);}
{conj}    {if (DEBUG) ECHO; return(CONJ);}
{renorm}    {if (DEBUG) ECHO; return(RENORM);}
{prime}     {if (DEBUG) ECHO; return(PRIME);}
%%

