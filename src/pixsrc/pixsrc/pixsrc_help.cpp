#ifndef PS_SIT
#ifdef INTEGER_TYPE_LL
typedef long long PS_SIT;
typedef unsigned long long PS_unsignedSIT;
#endif
#ifdef INTEGER_TYPE_L
typedef long PS_SIT;
typedef unsigned long PS_unsignedSIT;
#endif
#ifdef INTEGER_TYPE_I
typedef int PS_SIT;
typedef unsigned int PS_unsignedSIT;
#endif
#endif

#ifndef PS_FPT
#ifdef SINGLE_PRECISION
typedef float PS_FPT;
#endif
#ifdef DOUBLE_PRECISION
typedef double PS_FPT;
#endif
#endif



#include "pixsrc_help.hpp"
#include "pixsrc_init.hpp"
#include "pixsrc_memory_templates.cpp"
#include "pixsrc_printer_templates.cpp"

#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <algorithm>
#include <unistd.h>

//#define CLR_MAG "\x1b[1;35;40m"
//#define CLR_GRE "\x1b[1;32;40m"
//#define CLR_RED "\x1b[1;31;40m"
//#define CLR_YEL "\x1b[1;33;40m"
//#define CLR_0   "\x1b[0m"
#define CLR_MAG ""
#define CLR_GRE ""
#define CLR_RED ""
#define CLR_YEL ""
#define CLR_0   ""

// copied from http://blog.minhazulhaque.com/2012/09/find-and-replace-all-occurrences-in-cpp-string.html.html
void ps_find_and_replace(string& source, string const& find, string const& replace)
{
    for(string::size_type i = 0; (i = source.find(find, i)) != string::npos;)
    {
        source.replace(i, find.length(), replace);
        i += replace.length();
    }
}

void convert2latex (string& str)
{
    ps_find_and_replace (str, "_", "\\_");
    ps_find_and_replace (str, "<", "\\textless{}");
    ps_find_and_replace (str, ">", "\\textgreater{}");
    ps_find_and_replace (str, "chi^2", "$\\chi^{2}$");
}

void pixsrc_help::general ()
{
    std::cout << std::endl <<
        "pixsrc :: Welcome!"                                                      << std::endl <<
        "pixsrc :: Help for using pixsrc:"                                        << std::endl <<
        "pixsrc :: "                                                              << std::endl <<
        "pixsrc :: All of the following commands are to be run from the"          << std::endl <<
        "pixsrc :: lensmodel prompt."                                             << std::endl <<
        "pixsrc :: "                                                              << std::endl <<
        "pixsrc :: For BRIEF information about input parameters, try"             << std::endl <<
        "pixsrc :: pixsrc help parms"                                             << std::endl <<
        "pixsrc :: "                                                              << std::endl <<
        "pixsrc :: For DETAILED information about input parameters, try"          << std::endl <<
        "pixsrc :: pixsrc help parms more"                                        << std::endl <<
        "pixsrc :: or see the manual."                                                              << std::endl <<
        "pixsrc :: "                                                              << std::endl <<
        "pixsrc :: For information about INPUT files, try"                        << std::endl <<
        "pixsrc :: pixsrc help input"                                             << std::endl <<
        "pixsrc :: or see the manual."                                                              << std::endl <<
        "pixsrc :: "                                                              << std::endl <<
        "pixsrc :: For information about OUTPUT files, try"                       << std::endl <<
        "pixsrc :: pixsrc help output"                                            << std::endl <<
        "pixsrc :: or see the manual."                                                              << std::endl <<
        "pixsrc :: "                                                              << std::endl <<
        "pixsrc :: To create a TEMPLATE for quick-starting, try"                   << std::endl <<
        "pixsrc :: pixsrc help template"                                         << std::endl <<
        "pixsrc :: Directories \"pixsrc\\_in\" and \"pixsrc\\_out\" will be created." << std::endl <<
        "pixsrc :: Several other files will also be created."                      << std::endl <<
        "pixsrc :: You can also modify one of the examples to get started quickly."                                                              << std::endl <<
        "pixsrc :: "                                                              << std::endl <<
        "pixsrc :: To create an up-to-date pixsrc manual, try"                    << std::endl <<
        "pixsrc :: pixsrc help manual"                                            << std::endl <<
        "pixsrc :: The manual will be in LaTeX file format. pixsrc will try compile it into a PDF using pdflatex."    << std::endl <<
        "pixsrc :: The LaTeX file's and other files' names have the format \"pixsrc\\_manual*\"."                       << std::endl <<
        "pixsrc :: "                                                              << std::endl <<
        std::endl;
}

void pixsrc_help::parms ()
{
    PS_SIT size, maxnamesize=0, maxvalsize=0;
    ps_parms_struct *pstruct;
    INIT parmscreator (NULL, &pstruct, &size);

    // get max lengths of strings
    for (PS_SIT i=0; i<size; ++i)
    {
        if ((PS_SIT)pstruct[i].sname.length()>maxnamesize)
            maxnamesize = (PS_SIT)pstruct[i].sname.length();
        if ((PS_SIT)pstruct[i].svalue.length()>maxvalsize)
            maxvalsize = (PS_SIT)pstruct[i].svalue.length();
    }

    // pad strings
    maxnamesize +=1;
    maxvalsize  +=1;
    string whitespace (maxnamesize+maxvalsize, ' ');

    for(PS_SIT i=0; i<size; ++i)
    {
        while ((PS_SIT)pstruct[i].sname.length()<maxnamesize)
            pstruct[i].sname += " ";
        while ((PS_SIT)pstruct[i].svalue.length()<maxvalsize)
            pstruct[i].svalue += " ";
    }

    // print parms
    std::cout << std::endl;
    for(PS_SIT k=0; k<ps_parm_cat_numcategories; ++k)
    {
        std::cout << whitespace;
        for (PS_SIT i=0; i<(PS_SIT)pstruct[i].scategory[k].length()+8; ++i) std::cout << "#";
        std::cout << std::endl;
        std::cout << whitespace << "### " << pstruct[0].scategory[k] << " ###" << std::endl;
        std::cout << whitespace;
        for (PS_SIT i=0; i<(PS_SIT)pstruct[i].scategory[k].length()+8; ++i) std::cout << "#";
        std::cout << std::endl << std::endl;
        for (PS_SIT i=0; i<size; ++i)
        {
            if (k==pstruct[i].category)
            {
                std::cout << pstruct[i].sname << pstruct[i].svalue
                          << "# " << pstruct[i].qdescr[0];
                for (PS_SIT j=1; j<(PS_SIT)pstruct[i].qdescr.size(); ++j)
                    std::cout << std::endl << whitespace << "# " << pstruct[i].qdescr[j];
                std::cout << std::endl;
            }
        }
        std::cout << std::endl;
    }

    delete pstruct->pindex;
    delete [] pstruct;
}

void pixsrc_help::parms_more ()
{
    PS_SIT size;
    ps_parms_struct *pstruct;
    INIT parmscreator (NULL, &pstruct, &size);

    string whitespace ("  ");
    string whitespace2 ("   ");

    // print parms
    std::cout << std::endl;
    for(PS_SIT k=0; k<ps_parm_cat_numcategories; ++k)
    {
        std::cout << whitespace;
        for (PS_SIT i=0; i<(PS_SIT)pstruct[i].scategory[k].length()+8; ++i) std::cout << "#";
        std::cout << std::endl;
        std::cout << whitespace << "### " << pstruct[0].scategory[k] << " ###" << std::endl;
        std::cout << whitespace;
        for (PS_SIT i=0; i<(PS_SIT)pstruct[i].scategory[k].length()+8; ++i) std::cout << "#";
        std::cout << std::endl << std::endl;
        for (PS_SIT i=0; i<size; ++i)
        {
            if (k==pstruct[i].category)
            {
                std::cout << pstruct[i].sname << std::endl;
                for (PS_SIT j=0; j<(PS_SIT)pstruct[i].qdescr.size(); ++j)
                    std::cout << whitespace << pstruct[i].qdescr[j] << std::endl;
                for (PS_SIT j=0; j<(PS_SIT)pstruct[i].notes.size(); ++j)
                    std::cout << whitespace << pstruct[i].notes[j] << std::endl;
                std::cout << whitespace << "Possible inputs.." << std::endl;
                for (PS_SIT j=0; j<(PS_SIT)pstruct[i].entries.size(); ++j)
                {
                    for (PS_SIT m=0; m<(PS_SIT)pstruct[i].entries[j][0].size(); ++m)
                        std::cout << whitespace2 << pstruct[i].entries[j][0][m] << std::endl;
                    for (PS_SIT m=0; m<(PS_SIT)pstruct[i].entries[j][1].size(); ++m)
                        std::cout << whitespace2 << pstruct[i].entries[j][1][m] << std::endl;
                }
                std::cout << "default value: " << pstruct[i].svalue << std::endl << std::endl;
            }
        }
        std::cout << std::endl;
    }

    delete pstruct->pindex;
    delete [] pstruct;
}

void pixsrc_help::manual ()
{
    std::cout << std::endl;
    PRINTER print2screen("pixsrc","Creating LaTeX file.",NULL);

    // make pixsrc out directory
    string cmd;
    cmd = "mkdir -p " + string(CONSTANT dir_out);
    if (system(cmd.c_str()))
        PRINTER printerror("pixsrc","failed to run: " + cmd, NULL);
    cmd = "rm -f " + string(CONSTANT dir_out) + "pixsrc_manual.*";
    if (system(cmd.c_str()))
        PRINTER printwarning("pixsrc","failed to run: " + cmd, NULL);

    // write PNG files 
    string pngfn;
    PS_SIT p;
    unsigned char *q;
    
    pngfn = string(CONSTANT dir_out) + "pixsrc_manual_simple.png";
    p = 32294;
    MEMORY ps_malloc (&q,p);
#include "pixsrc_help_simple.cpp"
    PRINTER printbinary <unsigned char> (-1, "", "", 0, pngfn, q, p, 1);
    MEMORY ps_free (q);

    pngfn = string(CONSTANT dir_out) + "pixsrc_manual_typical.png";
    p = 39750;
    MEMORY ps_malloc (&q,p);
#include "pixsrc_help_typical.cpp"
    PRINTER printbinary <unsigned char> (-1, "", "", 0, pngfn, q, p, 1);
    MEMORY ps_free (q);

    pngfn = string(CONSTANT dir_out) + "pixsrc_manual_uvplane.png";
    p = 39548;
    MEMORY ps_malloc (&q,p);
#include "pixsrc_help_uvplane.cpp"
    PRINTER printbinary <unsigned char> (-1, "", "", 0, pngfn, q, p, 1);
    MEMORY ps_free (q);

    pngfn = string(CONSTANT dir_out) + "pixsrc_manual_uvres.png";
    p = 7421;
    MEMORY ps_malloc (&q,p);
#include "pixsrc_help_uvres.cpp"
    PRINTER printbinary <unsigned char> (-1, "", "", 0, pngfn, q, p, 1);
    MEMORY ps_free (q);

    pngfn = string(CONSTANT dir_out) + "pixsrc_manual_perturbation.png";
    p = 23262;
    MEMORY ps_malloc (&q,p);
#include "pixsrc_help_perturbation.cpp"
    PRINTER printbinary <unsigned char> (-1, "", "", 0, pngfn, q, p, 1);
    MEMORY ps_free (q);

    pngfn = string(CONSTANT dir_out) + "pixsrc_manual_pertkappa.png";
    p = 6288;
    MEMORY ps_malloc (&q,p);
#include "pixsrc_help_pertkappa.cpp"
    PRINTER printbinary <unsigned char> (-1, "", "", 0, pngfn, q, p, 1);
    MEMORY ps_free (q);

    // done writing PNG files


    vector<string> outvec;
    string filename = string (CONSTANT dir_out) + "pixsrc_manual.tex";

    outvec.push_back ("\\documentclass{article}");
    outvec.push_back ("");
    outvec.push_back ("\\usepackage[top=1in, bottom=1.25in, left=1.25in, right=1.25in]{geometry}");
    outvec.push_back ("\\usepackage{hyperref}");
    outvec.push_back ("\\usepackage{float}");
    outvec.push_back ("\\usepackage{graphicx}");
    outvec.push_back ("");
    outvec.push_back ("\\begin{document}");
    outvec.push_back ("");
    outvec.push_back ("\\title{pixsrc user's manual}");
    outvec.push_back ("\\author{Amitpal S. Tagore}");
    outvec.push_back ("");
    outvec.push_back ("\\maketitle");
    outvec.push_back ("\\tableofcontents");
    outvec.push_back ("\\newpage");
    outvec.push_back ("");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("\\section{What is pixsrc?}");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("");
    outvec.push_back ("pixsrc is software for de-lensing objects that have been strongly gravitationally lensed.");
    outvec.push_back ("Among other things, users of pixsrc can:");
    outvec.push_back ("\\begin{itemize}");
    outvec.push_back ("\\item Reconstruct a pixelated source on two different irregular grids.");
    outvec.push_back ("\\item Use a combination of functional forms, like a S\\'{e}rsic profile, to model the source.");
    outvec.push_back ("\\item Expand the source's surface brightness analytically using shapelets.");
    outvec.push_back ("\\item Place a variety of priors on the source in a Bayesian framework.");
    outvec.push_back ("\\item Use a number of pre-source-reconstruction tests to quickly reject bad lens models.");
    outvec.push_back ("\\item Use CUDA to enable GPU computing.");
    outvec.push_back ("\\item Add lens potential perturbations.");
    outvec.push_back ("\\item Model interferometric data in the uv-plane.");
    outvec.push_back ("\\item Compute and/or take into account interpolation errors.");
    outvec.push_back ("\\item Simultaneously analyze multiple data sets of multiple lensed sources or images in different photometric/radio bands.");
    outvec.push_back ("\\end{itemize}");
    outvec.push_back ("");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("\\section{Installation}");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("");
    outvec.push_back ("Setup environment to recognize shared dependencies and run install script.");
    outvec.push_back ("");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("\\section{Quick start}");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("");
    outvec.push_back ("To get started quickly, the easiest way is to modify one of the pixsrc examples.");
    outvec.push_back ("Or, you can have pixsrc create a blank template.");
    outvec.push_back ("From the Unix prompt, type");
    outvec.push_back ("");
    outvec.push_back ("unixprompt\\$ lensmodel pixsrc help template");
    outvec.push_back ("");
    outvec.push_back ("This will create the ``pixsrc\\_in'' directory and several input files, which you can then modify.");
    outvec.push_back ("The basename will be ``bname'' and the imagename will be ``iname''.");
    outvec.push_back ("");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("\\section{Input files}");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("");

    vector<string> inputvec = HELP input (0);
    outvec.insert (outvec.end(), inputvec.begin(), inputvec.end());

    outvec.push_back ("");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("\\section{Parameters file}");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("");

    PS_SIT size;
    ps_parms_struct *pstruct;
    INIT parmscreator (NULL, &pstruct, &size);
    for(PS_SIT k=0; k<ps_parm_cat_numcategories; ++k)
    {
	std::transform(pstruct[0].scategory[k].begin()+1, pstruct[0].scategory[k].end(), 
		       pstruct[0].scategory[k].begin()+1, tolower);
        outvec.push_back ("\\subsection{"+pstruct[0].scategory[k]+"}");
        outvec.push_back ("\\begin{description}");
        for (PS_SIT i=0; i<size; ++i)
        {
            if (k==pstruct[i].category)
            {
                outvec.push_back ("\\item \\textbf{"+pstruct[i].sname+"} \\hfill \\\\");
                for (PS_SIT j=0; j<(PS_SIT)pstruct[i].qdescr.size(); ++j)
                    outvec.push_back (pstruct[i].qdescr[j]);
                outvec.push_back ("");
                for (PS_SIT j=0; j<(PS_SIT)pstruct[i].notes.size(); ++j)
                    outvec.push_back (pstruct[i].notes[j]);
                outvec.push_back("");
                outvec.push_back ("\\medskip");
                for (PS_SIT j=0; j<(PS_SIT)pstruct[i].entries.size(); ++j)
                {
                    for (PS_SIT m=0; m<(PS_SIT)pstruct[i].entries[j][0].size(); ++m)
                        outvec.push_back (pstruct[i].entries[j][0][m]);
                    outvec.push_back("");
                    for (PS_SIT m=0; m<(PS_SIT)pstruct[i].entries[j][1].size(); ++m)
                        outvec.push_back (pstruct[i].entries[j][1][m]);
                    outvec.push_back("");
                    outvec.push_back ("\\smallskip");
                }
                outvec.push_back("");
                outvec.push_back ("\\medskip");
                outvec.push_back ("default: "+pstruct[i].svalue);
            }
        }
        outvec.push_back ("\\end{description}");
    }
    delete pstruct->pindex;
    delete [] pstruct;

    outvec.push_back ("");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("\\section{Output files}");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("");

    vector<string> outputvec = HELP output (0);
    outvec.insert (outvec.end(), outputvec.begin(), outputvec.end());


    outvec.push_back ("");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("\\section{Examples}");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("");
    outvec.push_back ("There are several pixsrc examples provided that highlight difference features.");
    outvec.push_back ("In the directory named \"pixsrc\\_examples\", there are several subdirectories.");
    outvec.push_back ("Each of these contains a lensmodel script that ends in \".in\", which can be used to run a pixsrc example.");
    outvec.push_back ("The true lens model for each case is the same: two singular isothermal spheres (with different Einstein radii), separated by approximately 0.54 arcseconds.");
    outvec.push_back ("Snapshots of some of the more important outputs from pixsrc are are shown below.");
    outvec.push_back ("");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("\\subsection{Simple example}");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("");
    outvec.push_back ("In this case, the lens model is fixed to the true value, and a source reconstruction is done using the default values.");
    outvec.push_back ("The only required inputs to pixsrc, given in the pixsrc\\_in/simple.parameters file, describe the world coordinate system, the point spread function, and the noise level in the data.");
    outvec.push_back ("");
    outvec.push_back ("\\begin{figure}[H]");
    outvec.push_back ("\\centering");
    outvec.push_back ("\\includegraphics[width=0.95\\textwidth]{pixsrc_manual_simple}");
    outvec.push_back ("\\caption{Top row, left to right: data, model, model before PSF convolution, residuals. Bottom row, left to right: reconstructed source, source-plane noise, source-plane signal-to-noise, PSF.}");
    outvec.push_back ("\\end{figure}");
    outvec.push_back ("");
    outvec.push_back ("");
    outvec.push_back ("");
    outvec.push_back ("");
    outvec.push_back ("");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("\\subsection{Typical (multiple datasets) example}");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("");
    outvec.push_back ("In this example, two objects are modelled simultaneously using the same (true) lens model.");
    outvec.push_back ("Unlike the simpler example, some of the default pixsrc settings in the parameters file have been changed.");
    outvec.push_back ("The object in the top row of Fig.~\ref{fig:typical} is modelled using a parametric source model (an elliptical S\\'{e}rsic profile), while the objetc in the bottom row is modelled using grid-free shapelets (a basis-set expansion of the source's surface brightness).");
    outvec.push_back ("");
    outvec.push_back ("\\begin{figure}[H]");
    outvec.push_back ("\\centering");
    outvec.push_back ("\\includegraphics[width=0.95\\textwidth]{pixsrc_manual_typical}");
    outvec.push_back ("\\caption{Top row, left to right: data, model, model before PSF convolution, residuals, reconstructed source. Bottom row, left to right: data, model, model before PSF convolution, residuals, reconstructed source.}");
    outvec.push_back ("\\label{fig:typical}");
    outvec.push_back ("\\end{figure}");
    outvec.push_back ("");
    outvec.push_back ("");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("\\subsection{Lens optimization example}");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("");
    outvec.push_back ("This example shows how to efficiently minimize the lens model parameters.");
    outvec.push_back ("In an outer loop, the lens model parameters are varied, while in the inner ");
    outvec.push_back ("pixsrc loop, the best source is found for each lens model.");
    outvec.push_back ("The example also demonstrates how to quickly reject bad lens models ");
    outvec.push_back ("without having to do the more costly source reconstruction.");
    outvec.push_back ("");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("\\subsection{uv-plane example}");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("");
    outvec.push_back ("This example shows how to model data taken with a radio interferometer, so that the complex visibility measurements are modeled directly.");
    outvec.push_back ("First in the analysis, an appropriate image-plane resolution is determined, which depends on the the longest baselines in data.");
    outvec.push_back ("Then, the image-plane is decomposed into radial basis functions, integrated, Fourier transformed, and evaluated at the positions of the visibility measurements.");
    outvec.push_back ("Finally, because of the large size of visibility datasets, a series of ``divide-and-conquer'' matrix algebra steps take place.");
    outvec.push_back ("The resutls of these steps are stored in binary files, which can be read in later to speed of future runs.");
    outvec.push_back ("");
    outvec.push_back ("The example shows results of the uv-plane analysis on a very small dataset, but it also performs an image-plane analysis on the image-plane data used to create the visibility measurements.");
    outvec.push_back ("In this case, as can be seen from Fig.~\ref{fig:uv}, the image-plane residuals for the uv-plane reconstruction are not as good as those of the image-plane analysis, but this is largely due to the small size of the dataset used.");
    outvec.push_back ("Fig.~\ref{fig:uvres}, on the other hand, shows that the uv-plane residuals are consistent with the model fitting the data down to the noise level.");
    outvec.push_back ("");
    outvec.push_back ("\\begin{figure}[H]");
    outvec.push_back ("\\centering");
    outvec.push_back ("\\includegraphics[width=0.95\\textwidth]{pixsrc_manual_uvplane}");
    outvec.push_back ("\\caption{Top row (image-plane analysis), left to right: data, model, reconstructed source, residuals. Bottom row (uv-plane analysis), left to right: cleaned map, model, reconstructed source, residuals. Note that the model shown for the uv-plane analysis is not the model used by the code for constraining the visibility data. A similar statement holds for the residuals shown for the uv-plane analysis, as well. These are only shown here for comparison.}");
    outvec.push_back ("\\label{fig:uv}");
    outvec.push_back ("\\end{figure}");
    outvec.push_back ("");
    outvec.push_back ("\\begin{figure}[H]");
    outvec.push_back ("\\centering");
    outvec.push_back ("\\includegraphics[width=0.5\\textwidth]{pixsrc_manual_uvres}");
    outvec.push_back ("\\caption{Standard deviation of normalized residuals ([data-model]/sigma). Each data point contains the standard deviation of 200 normalized model residuals.}");
    outvec.push_back ("\\label{fig:uvres}");
    outvec.push_back ("\\end{figure}");
    outvec.push_back ("");
    outvec.push_back ("");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("\\subsection{Lens perturbation example}");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    outvec.push_back ("");
    outvec.push_back ("This example shows how to add non-parametric components on top of a parametric lens model.");
    outvec.push_back ("One of the two lensing galaxies in the true lens model is removed, and the code is able to recover the missing galaxy, as can be seen from Fig.~\ref{fig:perturb} and \ref{fig:pertkappa}.");
    outvec.push_back ("Although this approach can be dangerous (in that it can produce unphysical features in the lens model), it can be useful to track down deficiencies in the lens model.");
    outvec.push_back ("");
    outvec.push_back ("\\begin{figure}[H]");
    outvec.push_back ("\\centering");
    outvec.push_back ("\\includegraphics[width=0.95\\textwidth]{pixsrc_manual_perturbation}");
    outvec.push_back ("\\caption{Top row (galaxy removed), left to right: data, model, reconstructed source, residuals. Bottom row (perturbations added), left to right: data, model, reconstructed source, residuals. Note that in this extreme case, the ``perturbations'' aren't really perturbations anymore; they represent a severe lacking in the lens model. However, this approach can still be used to expose something missing in the model, as long as the \"perturbations\" are physically plausible.}");
    outvec.push_back ("\\label{fig:perturb}");
    outvec.push_back ("\\end{figure}");
    outvec.push_back ("");
    outvec.push_back ("\\begin{figure}[H]");
    outvec.push_back ("\\centering");
    outvec.push_back ("\\includegraphics[width=0.5\\textwidth]{pixsrc_manual_pertkappa}");
    outvec.push_back ("\\caption{Twice the convergence map derived from the lens potential perturbations.}");
    outvec.push_back ("\\label{fig:pertkappa}");
    outvec.push_back ("\\end{figure}");
    outvec.push_back ("");
    outvec.push_back ("");
    outvec.push_back ("\\end{document}");

    // write file
    PRINTER print (0, "", "", 0, filename, outvec, 1, CONSTANT precision0, NULL);

    PRINTER print2screen("pixsrc","compiling LaTeX file.",NULL);

    // switch directories and run LaTeX
    if (chdir (CONSTANT dir_out))
        PRINTER printerror("pixsrc", "could not change directories.", NULL);

    PS_SIT usepdflatex = 1;
    cmd = "command -v pdflatex";
    if (system(cmd.c_str()))
    {
        usepdflatex = 0;
        //PRINTER printwarning("pixsrc", "pdflatex not found. Checking for latex.", NULL);
        PRINTER printwarning("pixsrc", "pdflatex not found. Only creating input files for pdflatex. Not creating pdf.", NULL);
    }

    if (usepdflatex)
    {
        cmd = "pdflatex pixsrc_manual.tex";
        for (PS_SIT mm=0; mm<2; ++mm)
        {
            if (system(cmd.c_str()))
                PRINTER printwarning("pixsrc",
                                     "problem in compiling LaTeX file; check output.", NULL);
            std::cout << std::endl << std::endl;
        }
    }
    else if (0)
    {
        cmd = "command -v latex";
        if (system(cmd.c_str()))
            PRINTER printerror("pixsrc", "created LaTeX file but latex not found", NULL);

        cmd = "latex pixsrc_manual.tex";
        for (PS_SIT mm=0; mm<2; ++mm)
        {
            if (system(cmd.c_str()))
                PRINTER printwarning("pixsrc",
                                     "problem in compiling LaTeX file; check output.", NULL);
            std::cout << std::endl << std::endl;
        }

        cmd = "command -v dvipdf";
        if (system(cmd.c_str()))
        {
            PRINTER printwarning("pixsrc", "dvi file created but dvipdf not found.", NULL);
        }
        else
        {
            cmd = "dvipdf pixsrc_manual.dvi";
            if (system(cmd.c_str()))
                PRINTER printwarning("pixsrc", "problem with running dvipdf; check output.", NULL);
        }
        std::cout << std::endl << std::endl;
    }

    PRINTER print2screen("pixsrc", "check " + string(CONSTANT dir_out) + " for manual files.", NULL);
}

vector<string> pixsrc_help::input (PS_SIT printit)
{
    string whitesp_0 ("");
    vector<string> outvec;

    outvec.push_back ("");
    outvec.push_back ("");
    outvec.push_back ("\\subsection{A few important notes}");
    outvec.push_back ("");
    outvec.push_back (whitesp_0 +            "All input files should be placed in the directory \"pixsrc\\_in\"." );
    outvec.push_back ("");
    outvec.push_back (whitesp_0 +            "When turning pixsrc on, using the" );
    outvec.push_back (whitesp_0 +            "pixsrc on" );
    outvec.push_back (whitesp_0 +            "command in lensmodel, you must specify a basename." );
    outvec.push_back (whitesp_0 +            "The basename tells pixsrc which files to read in." );
    outvec.push_back (whitesp_0 +            "Below, any occurence of " + CLR_GRE + "\\textless{}basename\\textgreater{}" + CLR_0 + " refers to this basename." );
    outvec.push_back (whitesp_0 +            "You must also supply one or more FITS files for pixsrc to read in." );
    outvec.push_back (whitesp_0 +            "Below, any occurrence of " + CLR_GRE + "\\textless{}imagename\\textgreater{}" + CLR_0 + " refers to the name of the FITS file" );
    outvec.push_back (whitesp_0 +            "without the \".fits\" extension." );
    outvec.push_back (whitesp_0 +            "For example, \\textless{}imagename\\textgreater{} corresponds to the file imagename.fits." );
    outvec.push_back ("");
    outvec.push_back (whitesp_0 +            "Any input files that end with the \".reg\" extension are DS9 region files." );
    outvec.push_back (whitesp_0 +            "In order for pixsrc to understand these files, they must be saved image coordinates" );
    outvec.push_back (whitesp_0 +            "or WCS coordinates." );
    outvec.push_back (whitesp_0 +            "If WCS, then the coordinates should be in degrees and not sexagesimal (hr:min:sec) format." );
    outvec.push_back (whitesp_0 +            "This means that there should be a line in the region file that contains the word" );
    outvec.push_back (whitesp_0 +            "\"image\" or \"fk5\"." );
    outvec.push_back (whitesp_0 +            "Furthermore, all regions should be the polygon region, as opposed to the circle or box region." );
    outvec.push_back (whitesp_0 +            "There can be multiple polygons and the polygons can self-intersect." );
    outvec.push_back (whitesp_0 +            "The one exception is the source mask described below, which only accepts circle regions." );
    outvec.push_back (whitesp_0 +            "There can be multiple circle regions." );
    outvec.push_back ("");
    outvec.push_back ("");
    outvec.push_back ("\\subsection{Input files}");
    outvec.push_back ("");
    outvec.push_back ("\\begin{description}");
    outvec.push_back ("");
    outvec.push_back (whitesp_0 + "\\item \\textbf{ " + "\\textless{}basename\\textgreater{}.include" + "} \\hfill \\\\" );
    outvec.push_back (whitesp_0 +            "This include file lists the names of the FITS file(s) to be included in the lensing analysis." );
    outvec.push_back (whitesp_0 +            "A cleaned image must be given for uv-plane modelling as well (to estimate boundaries and source/image scales)." );
    outvec.push_back (whitesp_0 +            "One file per line." );
    outvec.push_back ("");
    outvec.push_back (whitesp_0 + "\\item \\textbf{ " + "\\textless{}imagename\\textgreater{}.fits" + "} \\hfill \\\\" );
    outvec.push_back (whitesp_0 +            "This data file is a data file that contains the lensed images of the background source(s)." );
    outvec.push_back (whitesp_0 +            "There can be multiple data files, which should all be listed in the include file." );
    outvec.push_back ("");
    outvec.push_back (whitesp_0 + "\\item \\textbf{ " + "\\textless{}imagename\\textgreater{}.psf.fits" + "} \\hfill \\\\" );
    outvec.push_back (whitesp_0 +            "This PSF file contains the psf for the corresponding \\textless{}imagename\\textgreater{}.fits file." );
    outvec.push_back (whitesp_0 +            "If the psf option (see parameters) is set to \"1\", then this file must exist." );
    outvec.push_back ("");
    outvec.push_back (whitesp_0 + "\\item \\textbf{ " + "\\textless{}basename\\textgreater{}.parameters" + "} \\hfill \\\\" );
    outvec.push_back (whitesp_0 +            "This parameters file contains all the parameters for pixsrcto read in." );
    outvec.push_back (whitesp_0 +            "Please try" );
    outvec.push_back (whitesp_0 +            "pixsrc help parms" );
    outvec.push_back (whitesp_0 +            "or" );
    outvec.push_back (whitesp_0 +            "pixsrc help parms more" );
    outvec.push_back (whitesp_0 +            "for details about input parameters." );
    outvec.push_back ("");
    outvec.push_back (whitesp_0 + "\\item \\textbf{ " + "\\textless{}basename\\textgreater{}\\_\\textless{}imagename\\textgreater{}.datamask.reg" + "} \\hfill \\\\" );
    outvec.push_back (whitesp_0 +            "This data mask file masks pixels that are to be used in the lensing analysis." );
    outvec.push_back (whitesp_0 +            "Other pixels may also be included in the lensing analysis if the " );
    outvec.push_back (whitesp_0 +            "sisterpix option (see parameters) is set." );
    outvec.push_back ("");
    outvec.push_back (whitesp_0 + "\\item \\textbf{ " + "\\textless{}basename\\textgreater{}\\_\\textless{}imagename\\textgreater{}.badpixelmask.reg" + "} \\hfill \\\\" );
    outvec.push_back (whitesp_0 +            "This bad pixel mask file masks pixels that are flagged as bad and not to be used in the lensing analysis." );
    outvec.push_back (whitesp_0 +            "These pixels may, however, be used if the fillbadpix option (see parameters) is set." );
    outvec.push_back (whitesp_0 +            "They will also always be used in some penalty functions (see parameters) which do not use" );
    outvec.push_back (whitesp_0 +            "the surface brightness values." );
    outvec.push_back (whitesp_0 +            "Bad pixels may also be specified in a text file, " + "\\textbf{" + "\\textless{}basename\\textgreater{}\\_\\textless{}imagename\\textgreater{}.badpixelmask.ascii}, " );
    outvec.push_back (whitesp_0 +            "where the first column is x coordinates and the second column is y coordinates." );
    outvec.push_back (whitesp_0 +            "The pixels are numbered such that the bottom left pixel is (0,0)." );
    outvec.push_back (whitesp_0 +            "This feature is useful if you want to exclude all pixels below a certain surface brightness." );
    outvec.push_back (whitesp_0 +            "(You would have to make the pixel list yourself.)" );
    outvec.push_back ("");
    outvec.push_back (whitesp_0 + "\\item \\textbf{ " + "\\textless{}basename\\textgreater{}\\_\\textless{}imagename\\textgreater{}.magmask.reg" + "} \\hfill \\\\" );
    outvec.push_back (whitesp_0 +            "This should have N polygons, enclosing N images for which which magnifications unertainties should be calculated." );
    outvec.push_back (whitesp_0 +            "Masks should not overlap." );
    outvec.push_back (whitesp_0 +            "There will be N+1 magnification uncertainties reported." );
    outvec.push_back (whitesp_0 +            "The \"extra +1\" uncertainty calculation is for the union of all image polygons." );
    outvec.push_back (whitesp_0 +            "The first N uncertainties in the output file are ordered according to the order they appear in the input mask file." );
    outvec.push_back (whitesp_0 +            "The (N+1)th uncertainty is for the union of all images." );
    outvec.push_back ("");
    outvec.push_back (whitesp_0 + "\\item \\textbf{ " + "\\textless{}basename\\textgreater{}\\_\\textless{}imagename\\textgreater{}.chi2mask.reg" + "} \\hfill \\\\" );
    outvec.push_back (whitesp_0 +            "This $\\chi^2$ mask file masks pixels where the $\\chi^2$ term should be calculated." );
    outvec.push_back (whitesp_0 +            "No other pixels will be used to compute the $\\chi^2$ term." );
    outvec.push_back (whitesp_0 +            "However, if the Bayesian evidence is the primary statistic being computed," );
    outvec.push_back (whitesp_0 +            "then all eligible pixels will be used to constrain the source surface brightness." );
    outvec.push_back (whitesp_0 +            "Note: the value of the $\\chi^2$ term can be monitored in the output log file." );
    outvec.push_back ("");
    outvec.push_back (whitesp_0 + "\\item \\textbf{ " + "\\textless{}basename\\textgreater{}\\_\\textless{}imagename\\textgreater{}.mmimages.*.reg" + "} \\hfill \\\\" );
    outvec.push_back (whitesp_0 +            "This mmimages file masks individual lensed images of a single background source." );
    outvec.push_back (whitesp_0 +            "If there is more than one source being lensed, then multiple such files can exist." );
    outvec.push_back (whitesp_0 +            "The asterisk \"*\"in the filename can be anything and is used to distinguish background sources." );
    outvec.push_back (whitesp_0 +            "Penalty functions (see parameters) use the mmimages files to compute penalties." );
    outvec.push_back (whitesp_0 +            "Thus, if a penalty function is turned on, then at least one mmimages file must exist." );
    outvec.push_back (whitesp_0 +            "An mmimages file must contain at least two polygon regions." );
    outvec.push_back ("");
    outvec.push_back (whitesp_0 + "\\item \\textbf{ " + "\\textless{}basename\\textgreater{}\\_\\textless{}imagename\\textgreater{}.srcmask.reg" + "} \\hfill \\\\" );
    outvec.push_back (whitesp_0 +            "This source mask file masks regions of the sky within which the source grid lies." );
    outvec.push_back (whitesp_0 +            "No source pixel lies outside this region." );
    outvec.push_back (whitesp_0 +            "Unlike the other region files, the source mask only accepts circle regions; no polygons allowed." );
    outvec.push_back ("");
    outvec.push_back (whitesp_0 + "\\item \\textbf{ " + "\\textless{}analytic source file\\textgreater{}" + "} \\hfill \\\\" );
    outvec.push_back (whitesp_0 +            "This analytic source file contains analytic source models." );
    outvec.push_back (whitesp_0 +            "The filename is specified by the src option (see parameters)." );
    outvec.push_back (whitesp_0 +            "There can be multiple source models." );
    outvec.push_back (whitesp_0 +            "If there are N source models, then the first N lines specify the models." );
    outvec.push_back (whitesp_0 +            "These lines must contain nine entries each:" );
    outvec.push_back (whitesp_0 +            "model $I_0$ ra dec e PA R Z n." );
    outvec.push_back (whitesp_0 +            "model specifies the functional form of the source; options are \"none\", \"sersic\", and \"vector\\_\\textless{}vecfilename\\textgreater{}\"." );
    outvec.push_back (whitesp_0 +            "\"none\" specifies a source with zero surface brightness." );
    outvec.push_back (whitesp_0 +            "\"sersic\" specifies a S\\'{e}rsic profile." );
    outvec.push_back (whitesp_0 +            "The surface brightness, $I(r)$, along the major axis is given by:" );
    outvec.push_back (whitesp_0 +            "$I(r) = I_0 \\exp[-(r/R)^{1/n}]$, where $r$ is the distance along the major axis." );
    outvec.push_back (whitesp_0 +            "e denotes the ellipticity and PA the position angle (in degrees) of the major axis, measured east of north." );
    outvec.push_back (whitesp_0 +            "The Z parameter is unused, and positions (ra and dec) are in units of arcseconds." );
    outvec.push_back (whitesp_0 +            "\"vector\\_\\textless{}vecfilename\\textgreater{}\" allows an arbritrary surface brightness distribution to be lensed.");
    outvec.push_back (whitesp_0 +            "The file \\textless{}vecfilename\\textgreater{} contains 3 columns: right ascension, declination, surface brightness, where the positions are given in arcseconds offsets, using the coordinate system defined using \"coorsys:\".");
    outvec.push_back (whitesp_0 +            "The syntax is: vector\\_\\textless{}vecfilename\\textgreater{} $I_{0}$ ra dec m11 m12 m21 m22 Z.");
    outvec.push_back (whitesp_0 +            "$I_{0}$ is an overall scaling of the surface brightness. ra and dec are position offsets in arcseconds. Z is ignored. m?? denote elements of a linear transformation of the coordinates, so that the new coordinates are given by");
    outvec.push_back (whitesp_0 +            "\\begin{equation}");
    outvec.push_back (whitesp_0 +            "\\bigg({ra' \\atop dec'}\\bigg) = \\bigg({m11~m12 \\atop m21~m22}\\bigg) \\bigg({ra \\atop dec}\\bigg)");
    outvec.push_back (whitesp_0 +            "\\end{equation}");
    outvec.push_back (whitesp_0 +            "This allows for rotations, stretches, reflections, etc. of the coordinates to be optimized.");
    outvec.push_back (whitesp_0 +            "Then, the N+1 through 2N lines specify which parameters are to be varied." );
    outvec.push_back (whitesp_0 +            "Each of these lines must contain 8 entries." );
    outvec.push_back (whitesp_0 +            "\"0\" indicates the parameter is fixed, and \"1\" indicates the parameter is to be varied." );
    outvec.push_back ("");
    outvec.push_back (whitesp_0 + "\\item \\textbf{ " + "\\textless{}analytic source bounds file\\textgreater{}" + "} \\hfill \\\\" );
    outvec.push_back (whitesp_0 +            "This analytic source bounds file contains bounds on analytic source model parameters." );
    outvec.push_back (whitesp_0 +            "There can be multiple bounds." );
    outvec.push_back (whitesp_0 +            "Each bound is a linear combination of source model parameters." );
    outvec.push_back (whitesp_0 +            "The sum of the combination must fall within the lower and upper limits to avoid" );
    outvec.push_back (whitesp_0 +            "a very large penalty." );
    outvec.push_back (whitesp_0 +            "The file must have the following form:" );
    outvec.push_back (whitesp_0 +            "isrc iparm ai  jsrc jparm aj  ...  lo hi" );
    outvec.push_back (whitesp_0 +            "where isrc is the source number (starting from 1), iparm is the parameter" );
    outvec.push_back (whitesp_0 +            "number (starting from 1), ai is the weight of the parameter, lo is the lower bound, and" );
    outvec.push_back (whitesp_0 +            "hi is the upper bound." );
    outvec.push_back (whitesp_0 +            "Thus, each bound has the form:" );
    outvec.push_back (whitesp_0 +            "  lo $<=$ ai*pi + aj*pj + ... + aN*pN $<=$ hi" );
    outvec.push_back (whitesp_0 +            "where pj=p[jsrc][jparm] and so forth, and any source model that violates the bound" );
    outvec.push_back (whitesp_0 +            "is penalized heavily." );
    outvec.push_back (whitesp_0 +            "So, if there are N terms in a bound, then that line in the file must have 3N+2 terms." );
    outvec.push_back (whitesp_0 +            "This file is optional." );
    outvec.push_back ("");
    outvec.push_back (whitesp_0 + "\\item \\textbf{ " + "\\textless{}analytic source step sizes file\\textgreater{}" + "} \\hfill \\\\" );
    outvec.push_back (whitesp_0 +            "This analytic source step sizes file contains step sizes for setting up the initial simplex" );
    outvec.push_back (whitesp_0 +            "for source model optimization." );
    outvec.push_back (whitesp_0 +            "The number of lines in the file must equal the number of sources in the model." );
    outvec.push_back (whitesp_0 +            "There must be eight entries on each line." );
    outvec.push_back (whitesp_0 +            "Each entry specifies the step size for the corresonding parameter." );
    outvec.push_back (whitesp_0 +            "This file is optional." );
    outvec.push_back (whitesp_0 +            "Default step sizes are: 1 1 1 0.25 45 1 0 3" );
    outvec.push_back (whitesp_0 +            "Position step sizes are in units of arcseconds." );
    outvec.push_back ("");
    outvec.push_back ("");
    outvec.push_back (whitesp_0 +            "For completeness, there is an input file associated with non-parametric lens potential" );
    outvec.push_back (whitesp_0 +            "perturbations." );
    outvec.push_back (whitesp_0 +            "However, this is still under refinement and is not discussed here." );
    outvec.push_back (whitesp_0 +            "ADD SHAPEINT MASK AS WELL." );
    outvec.push_back ("");
    outvec.push_back ("\\end{description}");
    outvec.push_back ("");

    if (printit)
        for (PS_SIT i=0; i<(PS_SIT)outvec.size(); ++i)
            std::cout << outvec[i] << std::endl;

    return outvec;
}

vector<string> pixsrc_help::output (PS_SIT printit)
{
    string whitesp_0 ("");
    vector<string> outvec;

    outvec.push_back ("");
    outvec.push_back ("");
    outvec.push_back ("\\subsection{A few notes}");
    outvec.push_back ("");
    outvec.push_back (whitesp_0 +            "All output files will be placed in the directory \"pixsrc\\_out\"." );
    outvec.push_back ("");
    outvec.push_back (whitesp_0 +            "Files ending in ``.VECTOR'' contain three columns: x-pixel coordinate, y-pixel coordinate, value.");
    outvec.push_back (whitesp_0 +            "The coordinates starts from zero, with x increasing rightwards and y increasing downwards, with respect to the input data image.");
    outvec.push_back ("");
    outvec.push_back (whitesp_0 +            "Files ending in ``.fits'' are FITS files.");
    outvec.push_back ("");
    outvec.push_back (whitesp_0 +            "Files ending in any other extension are detailed below.");
    outvec.push_back ("");
    outvec.push_back (whitesp_0 +            "Additionally, the general format for output files that are NOT specific to a particular image is: ``\\textless{}basename\\textgreater{}\\_*''.");
    outvec.push_back (whitesp_0 +            "The general format for output files that are specific to a particular image is: ``\\textless{}basename\\textgreater{}\\_\\textless{}imagename\\textgreater{}\\_\\textless{}imagenumber\\textgreater{}\\_*'', where");
    outvec.push_back (whitesp_0 + "\\textless{}imagenumber\\textgreater{} is a number assigned to a particular pixsrc run (when pixsrc does multiple source reconstructions for one lens model).");
    outvec.push_back ("");
    outvec.push_back (whitesp_0 + "Lastly, in debug mode, more filetypes that those described below will be written, but these are not discussed here.");
    outvec.push_back ("");
    outvec.push_back ("\\subsection{Output files}");
    outvec.push_back ("");
    outvec.push_back ("\\begin{description}");
    outvec.push_back ("");
    outvec.push_back (whitesp_0 + "\\item \\textbf{ " + "\\textless{}basename\\textgreater{}\\_\\textless{}imagename\\textgreater{}\\_\\textless{}imagenumber\\textgreater{}\\_data.fits" + "} or \\\\ \\textbf{ " + "\\textless{}basename\\textgreater{}\\_\\textless{}imagename\\textgreater{}\\_\\textless{}imagenumber\\textgreater{}\\_data.VECTOR" + "}\\hfill \\\\" );
    outvec.push_back (whitesp_0 +            "The data pixsrc used to do the source reconstruction (except when modelling in uv-plane) after any modifications/masking." );
    outvec.push_back ("");
    outvec.push_back (whitesp_0 + "\\item \\textbf{ " + "\\textless{}basename\\textgreater{}\\_\\textless{}imagename\\textgreater{}\\_\\textless{}imagenumber\\textgreater{}\\_lensedmps.fits" + "} or \\\\ \\textbf{ " + "\\textless{}basename\\textgreater{}\\_\\textless{}imagename\\textgreater{}\\_\\textless{}imagenumber\\textgreater{}\\_lensedmps.VECTOR" + "}\\hfill \\\\" );
    outvec.push_back (whitesp_0 +            "The model that is compared to the data (except for uv-plane modelling)." );
    outvec.push_back ("");
    outvec.push_back (whitesp_0 + "\\item \\textbf{ " + "\\textless{}basename\\textgreater{}\\_\\textless{}imagename\\textgreater{}\\_\\textless{}imagenumber\\textgreater{}\\_lensedmpsnobo.fits" + "} or \\\\ \\textbf{ " + "\\textless{}basename\\textgreater{}\\_\\textless{}imagename\\textgreater{}\\_\\textless{}imagenumber\\textgreater{}\\_lensedmpsnobo.VECTOR" + "}\\hfill \\\\" );
    outvec.push_back (whitesp_0 +            "The model before convolution with the PSF." );
    outvec.push_back ("");
    outvec.push_back (whitesp_0 + "\\item \\textbf{ " + "\\textless{}basename\\textgreater{}\\_\\textless{}imagename\\textgreater{}\\_\\textless{}imagenumber\\textgreater{}\\_residuals.fits" + "} or \\\\ \\textbf{ " + "\\textless{}basename\\textgreater{}\\_\\textless{}imagename\\textgreater{}\\_\\textless{}imagenumber\\textgreater{}\\_residuals.VECTOR" + "}\\hfill \\\\" );
    outvec.push_back (whitesp_0 +            "The model subtracted from the data." );
    outvec.push_back ("");
    outvec.push_back (whitesp_0 + "\\item \\textbf{ " + "\\textless{}basename\\textgreater{}\\_\\textless{}imagename\\textgreater{}\\_\\textless{}imagenumber\\textgreater{}\\_psf.fits" + "} or \\\\ \\textbf{ " + "\\textless{}basename\\textgreater{}\\_\\textless{}imagename\\textgreater{}\\_\\textless{}imagenumber\\textgreater{}\\_psf.VECTOR" + "}\\hfill \\\\" );
    outvec.push_back (whitesp_0 +            "The point spread function used for image convolution." );
    outvec.push_back ("");
    outvec.push_back (whitesp_0 + "\\item \\textbf{ " + "\\textless{}basename\\textgreater{}\\_\\textless{}imagename\\textgreater{}\\_\\textless{}imagenumber\\textgreater{}\\_mps.fits" + "} or \\\\ \\textbf{ " + "\\textless{}basename\\textgreater{}\\_\\textless{}imagename\\textgreater{}\\_\\textless{}imagenumber\\textgreater{}\\_mps.VECTOR" + "}\\hfill \\\\" );
    outvec.push_back (whitesp_0 +            "The reconstructed source. Note that the VECTOR file contains the grid coordinates used by pixsrc. The FITS file is the source interpolated onto a regular grid." );
    outvec.push_back ("");
    outvec.push_back (whitesp_0 + "\\item \\textbf{ " + "\\textless{}basename\\textgreater{}\\_\\textless{}imagename\\textgreater{}\\_\\textless{}imagenumber\\textgreater{}\\_noise.fits" + "} or \\\\ \\textbf{ " + "\\textless{}basename\\textgreater{}\\_\\textless{}imagename\\textgreater{}\\_\\textless{}imagenumber\\textgreater{}\\_noise.VECTOR" + "}\\hfill \\\\" );
    outvec.push_back (whitesp_0 +            "The noise in the source plane. Note that the VECTOR file contains the grid coordinates used by pixsrc. The FITS file is the noise interpolated onto a regular grid." );
    outvec.push_back ("");
    outvec.push_back (whitesp_0 + "\\item \\textbf{ " + "\\textless{}basename\\textgreater{}\\_\\textless{}imagename\\textgreater{}\\_\\textless{}imagenumber\\textgreater{}\\_s2n.fits" + "} or \\\\ \\textbf{ " + "\\textless{}basename\\textgreater{}\\_\\textless{}imagename\\textgreater{}\\_\\textless{}imagenumber\\textgreater{}\\_s2n.VECTOR" + "}\\hfill \\\\" );
    outvec.push_back (whitesp_0 +            "The signal-to-noise in the source plane. Note that the VECTOR file contains the grid coordinates used by pixsrc. The FITS file is the S/N interpolated onto a regular grid." );
    outvec.push_back ("");
    outvec.push_back (whitesp_0 + "\\item \\textbf{ " + "\\textless{}basename\\textgreater{}\\_\\textless{}imagename\\textgreater{}\\_\\textless{}imagenumber\\textgreater{}\\_details.dat" + "}\\hfill \\\\" );
    outvec.push_back (whitesp_0 +            "A log file with one line per pixsrc run. This is mainly useful to check for fatal errors or to inspect the various terms that go into calculating the Bayesian evidence." );
    outvec.push_back ("");
    outvec.push_back ("\\end{description}");
    outvec.push_back ("");


if (printit)
        for (PS_SIT i=0; i<(PS_SIT)outvec.size(); ++i)
            std::cout << outvec[i] << std::endl;

    return outvec;
}

void pixsrc_help::ps_template ()
{
    string filename;
    vector<string> outvec;

    // create input directory
    // make this more portable in free time!
    string cmd = "mkdir -p " + string(CONSTANT dir_in);
    system(cmd.c_str());

    // write out first FITS file.
    filename = string(CONSTANT dir_in) + "iname.fits";
    PS_SIT imgx = 31;
    PS_SIT imgy = 31;
    PS_SIT num_pix1 = imgx*imgy;
    double arcpix = 4.0/30.0;
    double wcsinfo[10] = {16.0,16.0,
                          161.25,58.0,
                          -arcpix/3600.0 * cos(30*CONSTANT deg2rad),
                          arcpix/3600.0 * sin(30*CONSTANT deg2rad),
                          arcpix/3600.0 * sin(30*CONSTANT deg2rad),
                          arcpix/3600.0 * cos(30*CONSTANT deg2rad),
                          1, 0};
    PS_SIT *dummy1;
    double *data1;
    MEMORY ps_malloc (&dummy1,num_pix1);
    MEMORY ps_malloc (&data1,num_pix1);

    for (PS_SIT x=0; x<num_pix1; ++x)
        dummy1[x] = x;

    std::fill (data1, data1+num_pix1, 0);
    const PS_SIT num_wh_pix = 42;
    PS_SIT coorlist[num_wh_pix*2] = { 
	3,3,  4,3,  5,3,  3,4,  3,5,  3,6,  3,7,  3,8,  4,8,  5,8,	
	7,3,  7,4,  7,5,  7,6,  7,7,  7,8,  8,6,  9,6, 10,6, 10,3, 
	10,4, 10,5, 10,6, 10,7, 10,8, 12,3, 12,4, 12,5, 12,6, 12,7,
	12,8, 13,6, 14,6, 15,6, 15,3, 15,4, 15,5, 15,6, 15,7, 15,8, 
	13,3,  14,3,
    };
    for (PS_SIT dump = 0; dump<num_wh_pix; ++dump)
	data1[coorlist[dump*2+0]*imgy+coorlist[dump*2+1]] = 1;

    PRINTER printfitsimgplane (0, "", "",
                               "pixsrc", 0, filename,
                               data1, wcsinfo, imgx, imgy, 1,
                               dummy1, NULL, 1, 1, NULL);

    MEMORY ps_free (dummy1);
    MEMORY ps_free (data1);

    template_1 ();

    std::cout << std::endl <<
        "pixsrc :: A pixsrc template has been created!" << std::endl <<
        std::endl;
}

void pixsrc_help::template_1 ()
{
    string filename;
    vector<string> outvec;

    // write out parameters file.
    filename = string (CONSTANT dir_in) + "bname.parameters";
    PS_SIT size, ind;
    ps_parms_struct *pstruct;
    INIT parmscreator (NULL, &pstruct, &size);
    ind = pstruct->pindex->ps_parm_ind_coorsys;
    outvec.push_back ("coorsys: J2000 161.25 58 0 0         # " + pstruct[ind].qdescr[0]);
    for (PS_SIT i=1; i<(PS_SIT)pstruct[ind].qdescr.size(); ++i)
        outvec.push_back ("                                     # " + pstruct[ind].qdescr[i]);
    ind = pstruct->pindex->ps_parm_ind_psf;
    outvec.push_back ("psf:     0                           # " + pstruct[ind].qdescr[0]);
    for(PS_SIT i=1; i<(PS_SIT)pstruct[ind].qdescr.size(); ++i)
        outvec.push_back ("                                     # " + pstruct[ind].qdescr[i]);
    ind = pstruct->pindex->ps_parm_ind_noise;
    outvec.push_back ("noise:   0.10                        # " + pstruct[ind].qdescr[0]);
    for(PS_SIT i=1; i<(PS_SIT)pstruct[ind].qdescr.size(); ++i)
        outvec.push_back ("                                     # " + pstruct[ind].qdescr[i]);
    PRINTER print (0, "", "", 0, filename, outvec, 1, CONSTANT precision0, NULL);
    delete pstruct->pindex;
    delete [] pstruct;

    // write out include file
    filename = string (CONSTANT dir_in) + "bname.include";
    outvec.clear();
    outvec.push_back ("iname.fits");
    PRINTER print (0, "", "", 0, filename, outvec, 1, CONSTANT precision0, NULL);

    // write out lensmodel starter file
    filename = "bname.in";
    outvec.clear();
    outvec.push_back ("set chimode   = 0");
    outvec.push_back ("set gridflag  = 0");
    outvec.push_back ("");
    outvec.push_back ("pixsrc on bname");
    outvec.push_back ("");
    outvec.push_back ("setlens 1 1");
    outvec.push_back (" alpha 1e-3 0 0 0.3 90 0 0 0 0 1");
    outvec.push_back (" 0 0 0 0 0 0 0 0 0 0");
    outvec.push_back ("");
    outvec.push_back ("optimize myopt");
    outvec.push_back ("");
    outvec.push_back ("quit");
    PRINTER print (0, "", "", 0, filename, outvec, 1, CONSTANT precision0, NULL);
}

