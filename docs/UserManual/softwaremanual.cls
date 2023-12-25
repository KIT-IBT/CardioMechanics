\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{softwaremanual}

% load the base class
\LoadClass[11pt,a4paper]{article}

\RequirePackage[english]{babel}
\RequirePackage[a4paper,top=3.5cm,bottom=3.5cm,left=3cm,right=3cm]{geometry}
\RequirePackage{calc}
\RequirePackage[T1]{fontenc}
\RequirePackage{titlesec, blindtext, color}
\RequirePackage{siunitx}
\sisetup{output-exponent-marker=\ensuremath{\mathrm{e}}}
\RequirePackage{graphicx}
\graphicspath{{Graphics/}}
\RequirePackage{pdfsync}
\RequirePackage{float}
\RequirePackage{microtype}
\RequirePackage{amsmath}
\RequirePackage{dsfont}
\RequirePackage{cases}
\DeclareMathAlphabet{\mathcal}{OMS}{cmsy}{m}{n}
\RequirePackage{booktabs}
\RequirePackage{algorithm}
\RequirePackage{algpseudocode}
\RequirePackage{physics}
\RequirePackage{amssymb}
\RequirePackage[numbers,sort&compress]{natbib}
\RequirePackage[font=footnotesize,labelfont=bf]{caption}

\RequirePackage[pdftex,hidelinks]{hyperref}
\hypersetup{colorlinks=true,linkcolor=blue}

% Colored tables
\RequirePackage{colortbl}	
\definecolor{lightGray}{rgb}{0.95,0.95,0.95}
\definecolor{lightBlue}{rgb}{0.76,0.92,1}
\definecolor{lightGreen}{cmyk}{0.8,0,1,0}
\definecolor{lightYellow}{rgb}{0.99,0.99,0.81}
\definecolor{darkWhite}{rgb}{0.99,0.99,0.95}
\definecolor{gray}{rgb}{0.4,0.4,0.4}
\definecolor{darkblue}{rgb}{0.0,0.0,0.6}
\definecolor{cyan}{rgb}{0.0,0.6,0.6}

\RequirePackage{listings}
\lstdefinestyle{mystyle}{
    backgroundcolor=\color{lightGray},
    basicstyle=\ttfamily\footnotesize,
    commentstyle=\color{lightGreen},
    keywordstyle=\color{darkblue},
    stringstyle=\color{magenta},
    breakatwhitespace=false,         
    breaklines=true,                 
    captionpos=b,                    
    keepspaces=true,                 
    numbers=left,                    
    numbersep=5pt,                  
    showspaces=false,                
    showstringspaces=false,
    showtabs=false,                  
    tabsize=2
}
\lstset{style=mystyle}


% extra front matter commands
\DeclareRobustCommand*{\softwarelogo}[1]{\gdef\@logocmd{#1}}
\DeclareRobustCommand*{\version}[1]{\gdef\@ver{Version #1}}

% default values
\gdef\@logocmd{\@latex@warning{No logo command specified for title page. Ignoring.}}
\gdef\@ver{\@latex@warning{No version number specified for title page. Ignoring.}}

% modify the stock \maketitle command
\renewcommand\maketitle{
\begin{titlepage}
\centering
\@logocmd\par
{\huge\bfseries\@title\par}
\bigskip
{\Large\@ver\par}
\vfill
{\Large \@author\par}
\@date
\end{titlepage}
}

% font size of equation number
\makeatletter
\def\maketag@@@#1{\hbox{\m@th\normalfont\footnotesize#1}}
\makeatother

% space after Figure/Table
\setlength{\textfloatsep}{1.75\baselineskip plus 1pt minus 1pt}
\setlength{\floatsep}{1.75\baselineskip plus 1pt minus 1pt}			% Set space between different float figures
\setlength{\intextsep}{1.75\baselineskip plus 1pt minus 1pt}		% Set space around floats in text

% space between Figure caption and Figure
\captionsetup[figure]{skip=1.1\baselineskip}					

% space between Table caption and table
\captionsetup[table]{skip=1.1\baselineskip}

% space between equation and text
\setlength{\belowdisplayskip}{2.0\baselineskip plus 0pt minus 0pt}
\setlength{\abovedisplayskip}{2.0\baselineskip plus 0pt minus 0pt}

% space around headings
\titlespacing*{\section}{0ex}{3.5ex plus 0ex minus 0ex}{3.5ex plus 0ex minus 0ex}
\titlespacing*{\subsection}{0ex}{3.25ex plus 0ex minus 0ex}{3.25ex plus 0ex minus 0ex}
\titlespacing*{\subsubsection}{0ex}{3.25ex plus 0ex minus 0ex}{3.25ex plus 0ex minus 0ex}

% table font size
\AtBeginEnvironment{tabular}{\footnotesize}
