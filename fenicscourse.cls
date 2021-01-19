% Set up class
\ProvidesClass{fenics}[2011/05/19 FEniCS Course]
\NeedsTeXFormat{LaTeX2e}
\LoadClass[ps2pdf,serif]{beamer}

% Last changed: 2016-12-07

% Packages
\RequirePackage{etex}
\RequirePackage{graphicx}
\RequirePackage{psfrag}
\RequirePackage{amsfonts}
\RequirePackage{amssymb}
\RequirePackage{amsmath}
\RequirePackage{amsthm}
\RequirePackage{amscd}
\RequirePackage{stmaryrd}
\RequirePackage{fancyvrb}
\RequirePackage{color}
\RequirePackage{multimedia}
\RequirePackage{chemarr}
\RequirePackage{cancel}
\RequirePackage[english]{babel}
\RequirePackage[utf8]{inputenc}
\RequirePackage{codehighlight}
\RequirePackage{tikz}
\RequirePackage{pstricks}
\RequirePackage{algorithmic}
\RequirePackage{import}
\RequirePackage{pifont}

\usetikzlibrary{arrows,shapes,backgrounds,decorations,mindmap}

% Colors
\definecolor{fenicsweb}{rgb}{0.60674,0.11719,0.15625}
\definecolor{fenicsred}{rgb}{0.8,0.2,0.0}
\definecolor{grey}{rgb}{0.5,0.5,0.5}

% All sorts of things pertinent to styling beamer
\setbeamertemplate{navigation symbols}{}
\setbeamertemplate{footline}[frame number]
\usecolortheme[rgb={0.8,0.2,0.0}]{structure}
\usefonttheme{structurebold}
\useinnertheme{circles}
\setbeamercolor{frametitle}{fg=fenicsred}
\mode<presentation> {\setbeamercovered{transparent}}

% Environments
\DefineVerbatimEnvironment{code}{Verbatim}{frame=single,rulecolor=\color{blue}}
\DefineVerbatimEnvironment{macrocode}{Verbatim}{commandchars=\\\{\}}

% Title page
\newcommand{\fenicstitle}[5]{
\title{\textcolor{fenicsred}{#1}}
\author{#2}
\institute{\normalsize #3}
\date{\small #4 \\ \vspace{0.5cm} #5}}
\newcommand{\fenicssection}[1]{
\begin{frame}
\vspace{0.5cm}
\begin{center}
\Large\color{fenicsred}{#1}\normalsize
\end{center}
\end{frame}}

% Lecture title page
\newcommand{\fenicslecture}[2]{%
\begin{frame}
\rput{90}(0,1){\includegraphics[width=\textwidth]{png/fenics_mesh.png}} \\
\vspace{1cm}
{\Huge \bf FEniCS Course} \\
\bigskip
{\LARGE #1} \\
\bigskip
\emph{Contributors} \\
#2

\vfill\hfill
\includegraphics[width=0.55\textwidth]{png/fenics_logo_text.png}
%\rput[bl]{90}(7.0,1.0){\includegraphics[width=0.55\textwidth]{png/fenics_logo_text.png}}
\end{frame}
}

% Overview title page
\newcommand{\fenicslectureoverview}[2]{%
\begin{frame}
\rput[tr]{90}(9.5,1){\includegraphics[width=\textwidth]{png/fenics_mesh.png}} \\
\vspace{1cm}
{\Huge \bf FEniCS Course} \\
\bigskip
{\LARGE Overview and Introduction} \\
\vspace{2cm}
\begin{columns}[b]
\begin{column}{0.5\textwidth}
\emph{Lecturer} \\
#1 \\
#2
\end{column}
\begin{column}{0.5\textwidth}
\begin{flushright}
  \includegraphics[width=\textwidth]{png/fenics_logo_text.png}
\end{flushright}
\end{column}
\end{columns}
\end{frame}
}

% Math macros
\newcommand{\dt}{\, \mathrm{d}t}
\newcommand{\dx}{\, \mathrm{d}x}
\newcommand{\dX}{\, \mathrm{d}X}
\newcommand{\ds}{\, \mathrm{d}s}
\newcommand{\dS}{\, \mathrm{d}S}
\newcommand{\ddx}[1]{\mathrm{d}x_{#1}}
\newcommand{\R}{\mathbb{R}}
\newcommand{\foralls}{\forall \;}
\newcommand{\goal}{\mathcal{M}}
\newcommand{\bdry}{\partial T}
\newcommand{\Dt}{\Delta t}
\newcommand{\inner}[2]{\langle #1, #2 \rangle}
\DeclareMathOperator{\Div}{\mathrm{div}}
\DeclareMathOperator{\Curl}{\mathrm{curl}}
\DeclareMathOperator{\Grad}{\mathrm{grad}}
\DeclareMathOperator{\Trace}{\mathrm{tr}}

% Notation macros
\newcommand{\emp}[1]{\texttt{#1}}
\newcommand{\nedelec}{N\'ed\'elec}
\newcommand{\babuska}{Babu\v{s}ka}
\newcommand{\colemph}[1]{{\color{fenicsred}#1}}
\newcommand{\btVFill}{\vskip0pt plus 1filll}

% Remove footnote rule
\renewcommand{\footnoterule}{}

% Other macros
\newcommand{\reference}[3]%
{\let\thefootnote\relax\footnotetext{\hspace{-1.2cm}\tiny #1, \emph{#2} (#3)}}
\renewcommand{\note}[1]%
{\let\thefootnote\relax\footnotetext{\hspace{-1.2cm}\tiny #1}}
\newcommand{\tab}{\hspace*{2em}}

% Useful for two page printing
%\usepackage{pgfpages}
%\pgfpagesuselayout{2 on 1}[a4paper,border shrink=5mm]
