\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\usepackage{enpStyle}
\ProvidesClass{afit-etd}[2011/05/11 version 1.0 Air Force Institute of
Technology Electronic Thesis and Dissertation Class]

% This is a full rewrite of the AFIT class file to clean up the code and 
% remove some inconsistencies caused by references to older packages.  The
% entire class is now contained within this one file.  The format was made
% to be consistent with the AFIT Style Guide with only minor differences as
% were present in the existing class file.  This version is fully compatible
% with pdfLaTeX
%
% According to the documentation, the previous templates were maintained by 
% Mark Roth et al, Joel Young, and Maj Todd Hale.



% define some macros that allow the class to pass options to the underlying
% article.cls without causing erroneous error messages about unused options
\newcommand* \@leqno{\@gobble} 
\newcommand* \@fleqn{\@gobble} 
\newcommand* \@draft{\@gobble} 
\newcommand* \@fontsize{12pt}

\newif \if@pdftex            
\newif \if@nonumbers 
\newif \if@draft             
\newif \if@sstables 
\newif \if@cmfont
\newif \if@flyleafframe

\DeclareOption  {draft}     {\renewcommand*\@draft{draft}\@drafttrue}
\DeclareOption  {leqno}     {\renewcommand*\@fleqn{fleqn}}
\DeclareOption  {fleqn}     {\renewcommand*\@leqno{leqno}}
\DeclareOption  {11pt}      {\renewcommand*\@fontsize{11pt}}
\DeclareOption  {12pt}      {\renewcommand*\@fontsize{12pt}}
\DeclareOption  {cmfont}    {\@cmfonttrue}
\DeclareOption  {pdftex}    {\@pdftextrue}
\DeclareOption  {nonumbers} {\@nonumberstrue}
\DeclareOption  {sstables}  {\@sstablestrue}

\ExecuteOptions {12pt}
\ProcessOptions

\LoadClass[\@draft,\@fleqn,\@leqno,\@fontsize]{article}

\RequirePackage {lastpage} % count the number of pages
\RequirePackage {url}      % intelligent line breaks for URLs and DOIs
\RequirePackage {graphicx} % fundamental package for including graphics
\urlstyle{same}

\AtEndOfClass{
  \RequirePackage[bookmarks, breaklinks, colorlinks,
                  linkcolor=black, urlcolor=black, citecolor=black]{hyperref}
}
\RequirePackage{amsmath,amsthm} % AMS math/thm must be defined prior to txfonts
\if@cmfont\else\RequirePackage{txfonts}\fi % Times fonts versus Computer Modern

% provide the means to switch back and forth between single and double spacing, 
% adapted from the previous OU LaTeX class.
\newcommand* \singlespace{\renewcommand*\baselinestretch{1.0}\large\normalsize}
\newcommand* \@dsp{\renewcommand*\baselinestretch{1.655}\large\normalsize}

\newlength \@lineskip
\newlength \@fracskip
\newlength \@doubleskip
\setlength \@lineskip   {1.000\baselineskip}
\setlength \@fracskip   {0.655\baselineskip}
\setlength \@doubleskip {1.655\baselineskip}

\newcommand* \draftmode{\@drafttrue}
\AtBeginDocument{
  \if@draft  
    \singlespace
  \else
    \@dsp
  \fi
}

\newcounter{@appcount}
\AtEndDocument{
  \if@appendix
    \ifnum \c@chapter=\c@@appcount \else 
    \ClassWarningNoLine{afit-etd}{appendix count has changed. 
                            Rerun to get labels right}
    \fi
    \immediate\write\@auxout{\string\setcounter{@appcount}{\number\c@chapter}}
  \else
    \immediate\write\@auxout{\string\setcounter{@appcount}{0}}
  \fi    
}

\setlength\paperheight    {11.0in}
\setlength\paperwidth     { 8.5in}
                          
\setlength\hoffset        {0in}
\setlength\voffset        {0in}
\setlength\marginparwidth {0in}     % no marginal notes
\setlength\marginparsep   {0in}
\setlength\marginparpush  {0in}

\setlength\oddsidemargin  {0.25in}   % (in addition to the 1in default)
\setlength\evensidemargin {0.25in}   % (in addition to the 1in default)
\setlength\topmargin      {0in}      % (in addition to the 1in default)
\setlength\headheight     {0in} 
\setlength\headsep        {0in} 

\setlength  \footskip{.5in} % set page number .75" from bottom
                          
\setlength\textwidth {\paperwidth}
\setlength\textheight{\paperheight}

\addtolength \textwidth {-2.5in}    % 1.25" for left and right margins
\addtolength \textheight{-2.25in}   % 1" top and 1.25" bottom margins
\addtolength \textheight{-\headheight}
\addtolength \textheight{-\headsep}

\renewcommand \topfraction      {1} % allow floats to cover full page
\renewcommand \textfraction     {0}
\renewcommand \floatpagefraction{.5}

\def\ps@headings{%
  \def\@oddhead {}
  \def\@evenhead{}
  \def\@oddfoot {\hfil\thepage\hfil}
  \def\@evenfoot{}
  \let\@mkboth\markboth
}

% MACROS FOR THE PREFATORY PAGES

\newcommand*\graduation[2]{
  \begingroup
    \def\month{#1}
    \def\monthMar{March}
    \def\monthJun{June}
    \def\monthSep{September}
    \def\monthDec{December}
    \ifx\month\monthMar\endgroup\newcommand*\@gradmonth{#1}
    \else\ifx\month\monthJun\endgroup\newcommand*\@gradmonth{#1}
    \else\ifx\month\monthSep\endgroup\newcommand*\@gradmonth{#1}
    \else\ifx\month\monthDec\endgroup\newcommand*\@gradmonth{#1}
    \else
      \ClassWarning{afit-etd}{graduation month not recognized (expecting March,
        June, September, or December)}
    \fi\fi\fi\fi
  \newcommand*\@gradyear{#2}   
}

\global\let\@coauthor\@empty
\global\let\@signedDate\@empty

\newcommand*\department[1]  {\newcommand*\@department{#1}}
\newcommand*\school[1]      {\newcommand*\@school{#1}}

\newcommand*\dean[1]          {\newcommand*\@dean{#1}}
\newcommand*\deanSignedDate[1]{\newcommand*\@deanSignedDate{#1}}

\newcommand\committee[1]    {\newcommand\@committee{#1}}
\newcommand*\signedDate[1]  {\newcommand*\@signedDate{#1}}

\renewcommand\title[1]      {\renewcommand\@title{#1}}

\renewcommand*\author[1]    {\renewcommand*\@author{#1}}
\newcommand*\rank[1]        {\newcommand*\@rank{#1}}
\newcommand*\previousdegrees[1]{\newcommand\@previousdegrees{#1}}

\newcommand*\coauthor[1]    {\renewcommand*\@coauthor{#1}}
\newcommand*\corank[1]      {\newcommand*\@corank{#1}}
\newcommand*\copreviousdegrees[1]{\newcommand\@copreviousdegrees{#1}}

\newcommand*\@thesistype{THESIS}
\newcommand*\@PhDdegree{Doctor of Philosophy}
\newcommand*\degree[1]{
  \newcommand*\@degree{#1}%
  \ifx\@degree\@PhDdegree
    \newcommand*\@doctype{DISSERTATION}
  \else
    \newcommand*\@doctype{THESIS}
  \fi
}

\newif \if@isreport
\newcommand*\report[1][REPORT]{\newcommand*\@doctype{#1}\@isreporttrue}

\newcommand*\designator[1]{\newcommand*\@designator{#1}}
\newcommand*\@distribution{DISTRIBUTION STATEMENT A\\\singlespace{APPROVED FOR PUBLIC RELEASE; DISTRIBUTION IS UNLIMITED}}
\newcommand*\distribution[1]{\renewcommand\@distribution{#1}}

\newcommand*\flyleafframe{\@flyleafframetrue} 

% PREFATORY PAGES

\newcommand{\@centeredpage}[1]{  
  \newpage
  \vbox to \textheight{
    \singlespace
    \vfil
    \noindent #1 \par
    \vfil
  }
}

\newcommand\makePrefatoryPages{
  
  % error checking
  \begingroup
    \newif\if@ferr
    \@ferrfalse
    \@ifundefined{@gradmonth}{
      \@ferrtrue\ClassWarningNoLine{afit-etd}{missing graduation month}}{}
    \@ifundefined{@department}{
      \@ferrtrue\ClassWarningNoLine{afit-etd}{missing department name}}{}
    \@ifundefined{@author}{
      \@ferrtrue\ClassWarningNoLine{afit-etd}{missing author's name}}{}
    \@ifundefined{@title}{
      \@ferrtrue\ClassWarningNoLine{afit-etd}{missing document's title}}{}
    \@ifundefined{@dean}{
      \@ferrtrue\ClassWarningNoLine{afit-etd}{missing dean's information}}{}
    \if@ferr 
      \ClassError{afit-etd}{missing required prefatory definitions}
                 {All prefatory macros must be defined in the preamble}
    \fi
  \endgroup

  \if@pdftex 
    \hypersetup{pdftitle={\@title},pdfauthor={\@author}} 
    % The following line accomodates transparency in PDF files
    \pdfpageattr {/Group << /S /Transparency /I true /CS /DeviceRGB>>}
  \fi
  
  \pagenumbering{roman}

  % print the prefatory pages
  \@cover
  \if@isreport\else
    \@disclaimer 
    \@titlepage      
    % \@approvalpage 
	\includepdf[pages={1}]{SignaturePage.pdf}
    \@abstract  
    %\@abstract       
    \@dedication     
    \@acknowledgments
  \fi
  \tableofcontents 
  \@listoffigures  
  \@listoftables   
  \@listofsymbols   
  \@listofabbreviations   
  \@makelist

  \newpage
  \pagenumbering{arabic}

  % free up some memory
  \global\let\@pdftex          \@empty
  \global\let\@centeredpage    \@empty
  \global\let\@titlepage       \@empty
  \global\let\@disclaimer      \@empty
  \global\let\@approved        \@empty
  \global\let\@approvalpage    \@empty
  \global\let\@abstract        \@empty
  \global\let\@dedication      \@empty
  \global\let\@acknowledgments \@empty
  \global\let\@listoftables    \@empty
  \global\let\@listoffigures   \@empty
  \global\let\disclaimer       \relax
  \global\let\abstract         \relax
  \global\let\dedication       \relax
  \global\let\acknowledgments  \relax
  \global\let\acknowledgements \relax
  \global\let\tableofcontents  \relax
  \global\let\graduation       \relax
  \global\let\department       \relax
  \global\let\author           \relax
  \global\let\title            \relax
  \global\let\dean             \relax    
  \global\let\@gradyear        \@empty
  \global\let\@gradmonth       \@empty
  \global\let\@department      \@empty
  \global\let\@author          \@empty
  \global\let\@dean            \@empty    
  \global\let\@signwidth       \@empty
}
  
\newcommand\@cover{
  \thispagestyle{empty}
  \begingroup
  %\singlespace
    \if@flyleafframe\else\setlength\fboxrule{0pt}\fi
    \begin{center}
      \addvspace{1in}
      \includegraphics[width=3in]{afitlogo}\par
      \vfill
      \framebox[4in]{
        \vbox to 2in{
          \vspace{\fboxsep}
          {\large\expandafter\uppercase\expandafter{\bf{\@title}}\par}
          \vfill
          \@doctype
          \vfill
          \ifx\@coauthor\@empty
          \@author, \@rank
          \else
          \begin{minipage}[t]{2in}
            \begin{center}
              \@author\\
              \@rank
            \end{center}
          \end{minipage}
          \begin{minipage}[t]{2in}
            \begin{center}
              \@coauthor\\
              \@corank
            \end{center}
          \end{minipage}
          \fi\\[1\@lineskip minus 1\@fracskip]
          \@designator
          \vspace{\fboxsep}
        } % vbox
      }\\ % framebox
      \vfill
      {\bfseries
        \singlespace
        {\sffamily DEPARTMENT OF THE AIR FORCE\\AIR UNIVERSITY}\\[\@lineskip]
        {\itshape\Large AIR FORCE INSTITUTE OF TECHNOLOGY}\par
        \vspace{-.6\@lineskip}
        \rule{4.8in}{3pt}\par
        \vspace{-.7\@lineskip}
        \rule{4.8in}{1pt}\\[.5\@lineskip]
        {\sffamily Wright-Patterson Air Force Base, Ohio}\\[\@doubleskip]
      }
      \@distribution
    \end{center}
  \endgroup
  \newpage\setcounter{page}{1}\thispagestyle{empty}
}

\newcommand\@disclaimer{ 
  \@centeredpage{The views expressed in this
    \expandafter\lowercase\expandafter{\@doctype} are those of the
    author and do not reflect the official policy or position of the
    United States Air Force, the Department of Defense, or the United
    States Government. This material is declared a work of the U.S. Government and is not
    subject to copyright protection in the United States” } 
}
\newcommand\disclaimer[1]{\renewcommand\@disclaimer{\@centeredpage{#1}}}


\newcommand\@titlepage{
  \newpage
  \thispagestyle{empty}
  \noindent\@designator\\[1.5\@doubleskip]
  \begin{center}
    {\large\expandafter\uppercase\expandafter{\@title}\\[\@doubleskip]}
    \@doctype\\[1.5\@doubleskip]
    Presented to the Faculty\\
    \ifx\@doctype\@thesistype\@department\\\fi
    \@school\\
    Air Force Institute of Technology\\
    Air University\\
    Air Education and Training Command\\
    in Partial Fulfillment of the Requirements for the\\
    Degree of \@degree \\[1.5\@doubleskip]
    \ifx\@coauthor\@empty
      \@author, \@previousdegrees\\
      \@rank\\
    \else
      \begin{minipage}[t]{2in}
        \begin{center}
          \singlespace
          \@author, \@previousdegrees\\
          \@rank\\
        \end{center}
      \end{minipage}
      \hfill
      \begin{minipage}[t]{2in}
        \begin{center}
          \singlespace
          \@coauthor, \@copreviousdegrees\\
          \@corank\\
        \end{center}
      \end{minipage}\\[\@lineskip] % accounts for second line in minipage
    \fi
    \vspace{\@doubleskip}
    \@gradmonth\ \@gradyear\\[\@doubleskip]
    \vfill
    \@distribution
  \end{center}
  \setcounter{footnote}{0}
}


% The following loop mechanism for parsing the list of committee members was
% taken from the previous AFIT class with only minor changes.

% Comma delimited list of committee members is the parameter.
% e.g.
%\committee{Iam N. Charge, PhD\\Committee Chair,
%           Maj Painin Deass, PhD\\Committee Member,
%           Howcomb I. Last, PhD\\Committee Member}
%
%prints out signaturelines in aproval page
\newcommand\sign@turelines{
  \@ifundefined{@committee}{
    \@warning{Missing \string\committee\space command!}%
    \def\@comem{MEMBER HERE\\TITLE HERE}%
    \sign@line\sign@line\sign@line
  }%
  {\@for\@comem:=\@committee\do{\sign@line}}%
}
% format of signature line...
\newcommand*\@sigline[2]{
  \makebox[3.25in][l]{#1}\hspace{1in}\makebox[1.25in]{#2}\\
}
\newcommand\sign@line{\par
  \begingroup
    \vspace{2\@lineskip}
    \ifx\@signedDate\@empty
      \vspace{\@lineskip}
    \else
      \@sigline{// signed //}{\@signedDate}
    \fi
      \@sigline{\hrulefill}{\hrulefill}
      \@sigline{\@comem}{Date}
  \endgroup
}

\newcommand{\@approvalpage}{
  \newpage
  \thispagestyle{empty}
  \begingroup
    %\singlespace
    \noindent\@designator\\[2\@doubleskip]
    \begin{center}
      {\large\expandafter\uppercase\expandafter{\@title}\\[2\@doubleskip]}
      \ifx\@coauthor\@empty
        \@author, \@previousdegrees\\
        \@rank
      \else
      \begin{minipage}[t]{2in}
        \begin{center}
          \@author, \@previousdegrees\\
          \@rank
        \end{center}
      \end{minipage}
      \hfill
      \begin{minipage}[t]{2in}
        \begin{center}
          \@coauthor, \@copreviousdegrees\\
          \@corank
        \end{center}
      \end{minipage}
      \fi
    \end{center}
    \vfill
    \noindent Approved:
    \sign@turelines
    \addvspace{\@doubleskip}
    \ifx\@doctype\@thesistype\else
      \noindent Accepted:\par
      \vspace{2\@lineskip}
      \ifx\@deanSignedDate\@empty
        \vspace{\@lineskip}
      \else
        \@sigline{// signed //}{\@deanSignedDate}
      \fi
      \@sigline{\hrulefill}{\hrulefill}
      \@sigline{\@dean}{Date}
      Dean, \@school\par
    \fi
    \vfill
  \endgroup
}

%==================================================================
%																			
%    		ABSTRACT Environment              
%																			
%==================================================================
%\newenvironment{abstract}{
%    \newpage
%    \thispagestyle{plain}
%    %\addcontentsline{toc}{chapter}{Abstract}
%    %\addtocontents{toc}{\protect\addvspace{10\p@}}
%    \noindent\@designator\\[1.5\@doubleskip]
%    %\noindent \@docdesignator               
%    \begin{center}
%	\large\bfseries Abstract
%    \end{center}
%    \vspace{2em}
%    } %{\newpage}

\renewcommand\abstract[1]{
\newcommand\@abstract{\@schaptertoc{Abstract} #1}

}

%\newcommand\@abstract{\@schaptertoc{Abstract} #1}}

% \newcommand\dedication[1]{
% \newcommand\@dedication{\@schaptertoc{Dedication} #1}}

\newcommand\dedication[1]{
   \newcommand\@dedication{
        % \newpage
	    \begin{center}
			% \large\bf Dedication
			\refstepcounter{schapter}
			\addcontentsline{toc}{chapter}{Dedication}
			\@centeredpage{\textit{#1}}
		\end{center}
     }
}

\let\@acknowledgments\relax
\newcommand\acknowledgments[1]{
  \@ifundefined{@acknowledgments}
    {\newcommand\@acknowledgments{\@schaptertoc{Acknowledgments}#1\par\addvspace{\@doubleskip}\hspace{.55\linewidth}\@author}}
    {\ClassWarning{afit-etd}{ignoring subsequent \protect\acknowledgments}}
}
\newcommand\acknowledgements[1]{
  \@ifundefined{@acknowledgments}
    {\newcommand\@acknowledgments{\@schaptertoc{Acknowledgements}#1\par\addvspace{\@doubleskip}\hspace{.55\linewidth}\@author}}
    {\ClassWarning{afit-etd}{ignoring subsequent \protect\acknowledgements}}
}

\newcommand\vita[1]{\@schaptertoc{Vita}#1}


\renewcommand\tableofcontents{
  \clearpage
  \if@pdftex \pdfbookmark[0]{Table of Contents}{TOC} \fi
  \chapter*{Table of Contents}
  \begingroup 
    \singlespace
    \hfill Page\\
    \@starttoc{toc}
  \endgroup
}

\newcounter{listctr}
\newcounter{listcnt}
\renewcommand\thelistcnt{\@Alph\c@listcnt}
\renewcommand\thelistctr{\@Alph\c@listctr}

\newcommand\@makelist{
  \ifnum \c@listcnt>\c@listctr
    \refstepcounter{listctr}
    \csname @listof\thelistctr \endcsname
    \@makelist
  \fi
}
  
\newcommand\addlistof[2]{
  \refstepcounter{listcnt}
  \expandafter\gdef\csname @listof\thelistcnt\endcsname{\@addlistof{#1}{\relax}{}{}{#2}}
}

%% \@addlistof{Figures}{Figure}{lof}{\relax}
%\newcommand\@addlistof[4]{
%  \@schaptertoc{List of #1} 
%  \begingroup
%    %\singlespace
%    \ifx #4\relax\else\renewcommand*\l@figure{\@dottedtocline{1}{0em}{#4}}\fi
%    \ifx #2\relax
%      #4
%    \else
%      \setcounter{tocdepth}{1}
%      \noindent #2 \hfill Page\\
%      \@starttoc{#3}
%    \fi
%  \endgroup
%}

% \@addlistof{Figures}{Figure}{lof}{\relax}

\newcommand\@addlistof[4]{
  \@schaptertoc{List of #1} 
  \begingroup
    %\singlespace
    \ifx #4\relax\else\renewcommand*\l@figure{\@dottedtocline{1}{0em}{#4}}\fi
    %\ifx #4\relax\else\renewcommand*\l@figure{\@dottedtocline{1}{0.25em}{#4}}\fi
    \ifx #2\relax
      #4
    \else
      \setcounter{tocdepth}{1}
      \noindent #2 \hfill Page\\
      \begin{singlespace}\@starttoc{#3}\end{singlespace}
    \fi
  \endgroup
}

\newcommand\@listoffigures{
  \@addlistof{Figures}{Figure}{lof}{\relax}
   %\listoffigures
}

\newcommand\@listoftables{
  \@addlistof{Tables}{Table}{lot}{\relax}
}

\newcommand\listsymbolname{List of Symbols}

\newcommand\@listofsymbols{\relax}
\newcommand\listofsymbols[1][\relax]{\renewcommand\@listofsymbols{\@addlistof{Symbols}{Symbol}{los}{#1}}}
\newcommand\addsymbol[2]{#2\addcontentsline{los}{figure}{\numberline{#2}{#1}}}

\newcommand\@listofabbreviations{\relax}
\newcommand\listabbrevname{List of Abbreviations}
\newcommand\listofabbreviations[1][\relax]{\renewcommand\@listofabbreviations{\@addlistof{Abbreviations}{Abbreviation}{loa}{#1}}}
\newcommand\addabbrev[2]{#2\addcontentsline{loa}{figure}{\numberline{#2}{#1}}}

\newcommand\notables {\renewcommand\@listoftables{}}
\newcommand\nofigures{\renewcommand\@listoffigures{}}

\newcommand\references[1]{
  \clearpage
  \refstepcounter{schapter}
  \addcontentsline{toc}{chapter}{\bibname}
  \begingroup
    \singlespace 
    #1
  \endgroup
}

\newcommand*{\bibname}{Bibliography}
\renewenvironment{thebibliography}[1]{ 
  \@schaptertoc{\bibname}%
  \list{ \@biblabel{\@arabic\c@enumiv}}%
       { \settowidth\labelwidth{\@biblabel{#1}}%
         \leftmargin\labelwidth
         \advance\leftmargin\labelsep
         \usecounter{enumiv}%
         \let\p@enumiv\@empty
         \renewcommand\theenumiv{\@arabic\c@enumiv}
       }%
   \sloppy
   \clubpenalty4000
   \@clubpenalty \clubpenalty
   \widowpenalty4000
   \sfcode`\.\@m
   \singlespace
}
{\def\@noitemerr
  {\@latex@warning{Empty `thebibliography' environment}}%
  \endlist
}

\newif\if@appendix 
\renewcommand\appendix{
  \@appendixtrue
  \setcounter{chapter}{0}
  \renewcommand\thechapter {\@Alph\c@chapter}
  \ifnum \c@@appcount>1
    \newcommand*\@appname{Appendix\space\thechapter}
  \else
    \newcommand*\@appname{Appendix}
  \fi
  \immediate\write\@auxout{\string\@writefile{toc}{%
      \string\setcounter{tocdepth}{0}}} % limit toc to Appendix title only
}

\newcounter{chapter}
\newcounter{schapter}
\renewcommand\thesection  {\thechapter.\@arabic\c@section}

\newcommand\chapter{
  \clearpage
   \global\@topnum\z@
   \@afterindenttrue
   \secdef\@chapter\@schapter
}

\def\@schapter#1{ % used by some bibliography packages
  \@makechapterhead{#1}
  \vspace{\@fracskip}
}

\newcommand\@schaptertoc[2][schapter]{
  \clearpage
  \global\@topnum\z@
  \@afterindenttrue
  \refstepcounter{#1}
  \addcontentsline{toc}{chapter}{#2}
  \@makechapterhead{#2}
}

\if@nonumbers\setcounter{secnumdepth}{0}\fi
\def\@chapter[#1]#2{
  \refstepcounter{chapter}
  \if@appendix
    \typeout{Appendix \thechapter.}
    \addcontentsline{toc}{chapter}{\@appname:\space#1}
    \@makechapterhead{\@appname.\space#2}
  \else
    \typeout{Chapter\space\thechapter.}
    \addcontentsline{toc}{chapter}{\protect\numberline{\thechapter.}#1}
    \@makechapterhead{\thechapter.\quad#2}
  \fi
}

\def\@makechapterhead#1{
  \begingroup
    \ifnum\value{chapter}=1
      \if@appendix\else
        \begin{center}
          \rule{0pt}{3\@lineskip}
          \large\@title
          \vspace{\@doubleskip}
        \end{center}
      \fi
    \else
      \vspace{.25in}
    \fi
    \centering \bfseries #1 \par \nobreak
    \vspace{-\@fracskip}
    \addvspace{2\@lineskip}
  \endgroup
  \@afterheading
}

\renewcommand\section{\@startsection {section}{1}{\z@}%
  {2.1ex \@plus -1ex \@minus -.2ex}{1.5ex \@plus.2ex}%
  {\normalfont\normalsize\bfseries}}
\renewcommand\subsection{\@startsection{subsection}{2}{\parindent}%
  {2.1ex\@plus -1ex \@minus -.2ex}{-2ex}%
  {\normalfont\normalsize\itshape}}
\renewcommand\subsubsection{\@startsection{subsubsection}{3}{2\parindent}%
  {2.1ex\@plus -1ex \@minus -.2ex}{-2ex}%
  {\normalfont\normalsize\itshape}}
\renewcommand\paragraph{\@startsection{paragraph}{4}{\z@}%
  {0ex\@plus -1ex \@minus -.2ex}{-1em}{\normalfont\normalsize\bfseries}}
\renewcommand\subparagraph{\@startsection{subparagraph}{5}{\parindent}%
  {0ex \@plus1ex \@minus .2ex}{-1em}{\normalfont\normalsize\bfseries}}

% The following is extracted directly from the LaTeX source code so that a 
% minor modification could be made to add a period for subsections and deeper
% sections.  This is a simpler change than was previously implemented.
\def\@sect#1#2#3#4#5#6[#7]#8{%
  \ifnum #2>\c@secnumdepth
  \let\@svsec\@empty
  \else
  \refstepcounter{#1}%
  \protected@edef\@svsec{\@seccntformat{#1}\relax}%
  \fi
  \@tempskipa #5\relax
  \ifdim \@tempskipa>\z@
  \begingroup
  #6{%
    \@hangfrom{\hskip #3\relax\@svsec}%
    \interlinepenalty \@M #8\@@par}%
  \endgroup
  \csname #1mark\endcsname{#7}%
  \addcontentsline{toc}{#1}{%
    \ifnum #2>\c@secnumdepth \else
    \protect\numberline{\csname the#1\endcsname}%
    \fi
    #7}%
  \else
  \def\@svsechd{%
    #6{\hskip #3\relax
      \@svsec #8\ifnum #2>1.\fi}% add period after subsections and deeper
    \csname #1mark\endcsname{#7}%
    \addcontentsline{toc}{#1}{%
      \ifnum #2>\c@secnumdepth \else
      \protect\numberline{\csname the#1\endcsname}%
      \fi
      #7}}%
  \fi
  \@xsect{#5}
}
\def\@ssect#1#2#3#4#5{%
  \@tempskipa #3\relax
  \ifdim \@tempskipa>\z@
  \begingroup
  #4{%
    \@hangfrom{\hskip #1}%
    \interlinepenalty \@M #5\@@par}%
  \endgroup
  \else
  \def\@svsechd{#4{\hskip #1\relax #5\ifnum #2>1.\fi}}% add period
  \fi
  \@xsect{#3}
}

\if@nonumbers
  \renewcommand* \thechapter      {\@Roman\c@chapter}
  \renewcommand* \theequation     {\@arabic\c@equation}
  \renewcommand* \thefigure       {\@arabic\c@figure}
  \renewcommand* \thetable        {\@arabic\c@table}
  \renewcommand* \l@figure        {\@dottedtocline{1} {0em}   {2.3em}}
  %\renewcommand* \l@figure        {\@dottedtocline{1} {0em}   {4em}}
  \renewcommand* \l@table         {\@dottedtocline{1} {0em}   {2.3em}}
  \renewcommand* \l@section       {\@dottedtocline{1} {1.5em} {0em}}
  \renewcommand* \l@subsection    {\@dottedtocline{2} {3.0em} {0em}}
  \renewcommand* \l@subsubsection {\@dottedtocline{3} {4.5em} {0em}}
  \renewcommand* \l@paragraph     {\@dottedtocline{4} {6.0em} {0em}}
  \renewcommand* \l@subparagraph  {\@dottedtocline{5} {7.5em} {0em}}
\else
  \@addtoreset {equation} {chapter}
  \@addtoreset {figure}   {chapter}
  \@addtoreset {table}    {chapter}
  \@addtoreset {section}  {chapter}
  \renewcommand* \theequation     {\thechapter.\@arabic\c@equation}
  \renewcommand* \thefigure       {\thechapter.\@arabic\c@figure}
  \renewcommand* \thetable        {\thechapter.\@arabic\c@table}
  \renewcommand* \l@figure        {\@dottedtocline{1} {0em}   {2.3em}}
  \renewcommand* \l@table         {\@dottedtocline{1} {0em}   {2.3em}}
  \renewcommand* \l@section       {\@dottedtocline{1} {1.5em} {2.3em}}
  \renewcommand* \l@subsection    {\@dottedtocline{2} {3.8em} {3.2em}}
  \renewcommand* \l@subsubsection {\@dottedtocline{3} {7.0em} {4.1em}}
  \renewcommand* \l@paragraph     {\@dottedtocline{4} {10em}  {5em}}
  \renewcommand* \l@subparagraph  {\@dottedtocline{5} {12em}  {6em}}
\fi

\setcounter{tocdepth}{3}

\def\fps@figure{tbh}
\def\fnum@figure{Figure\nobreakspace\thefigure}
\renewenvironment{figure}
               {\vspace\@fracskip \@float{figure}}
               {\end@float}

\def\fps@table{tbh}
\def\fnum@table{Table\nobreakspace\thetable}
\renewenvironment{table}
               { \if@sstables \singlespace \fi
                 \vspace\@fracskip \@float{table}
               }
               { \vspace{-\@lineskip}
                 \end@float
               }
\renewenvironment{table*}
               {\if@sstables \singlespace \fi
                \@dblfloat{table}}
               {\end@dblfloat}

\setlength\floatsep      {2\@doubleskip \@plus 1\p@ \@minus 1\p@}
\setlength\textfloatsep  {1\@doubleskip \@plus 1\p@ \@minus 1\p@}
\setlength\intextsep     {1\@doubleskip \@plus 1\p@ \@minus 1\p@}
\setlength\@fptop        {0\p@} % align float pages with the top of the page
\addtolength\textfloatsep\@lineskip

\setlength\abovecaptionskip{1ex}
\setlength\belowcaptionskip{1ex}
\long\def\@makecaption#1#2{%
  \expandafter\ifx#1\fnum@table 
  \else 
    \vskip \abovecaptionskip 
  \fi
  \sbox\@tempboxa{#1. #2}%
  \ifdim \wd\@tempboxa >\hsize
    \begingroup
	  %\singlespace \bf{#1. #2}\par
      \singlespace \rm{#1. #2}\par			% Changed this from \bf to \rm so that all captions are roman, not bold
    \endgroup
  \else
    \global \@minipagefalse
    \begingroup
      \singlespace \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
    \endgroup
  \fi
  \expandafter\ifx#1\fnum@table 
   \vskip \belowcaptionskip
  \fi
}

\newcommand\l@chapter[2]{%
  \addpenalty{-\@highpenalty}%
  \setlength\@tempdima{1.5em}
  \begingroup
    \leavevmode
    \@dottedtocline{0}{0em}{1.5em}{#1}{#2}\par
    \nobreak
  \endgroup
}

\renewcommand\@makefntext[1]{%
  \noindent
  \hb@xt@2em{\hss\@makefnmark}
  \singlespace\footnotesize #1
}

\pagestyle{headings}
\raggedright
\setlength\parindent{2em} % must follow the \raggedright or will be reset to 0
\pretolerance=2000  
\clubpenalty =10000  
\widowpenalty=10000 
