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



#ifndef PIXSRC_PRINTER_TEMPLATES_CPP_
#define PIXSRC_PRINTER_TEMPLATES_CPP_

#include "pixsrc_printer.hpp"
#include "pixsrc_operations.hpp"
#include <stdio.h>
#include <iomanip> // for precision

template <typename T>
PS_SIT pixsrc_printer::writeoutstream( T source, std::ofstream *stream, pthread_mutex_t *lock, PS_SIT precision, vector<string>* comments )
{
    if (lock)
        pthread_mutex_lock(lock);

    if (comments)
        for (PS_SIT com=0; com<(PS_SIT)comments->size(); ++com)
            *stream << "# " << (*comments)[com] << std::endl;

    *stream << std::setprecision(precision) << source << std::endl;

    if (lock)
        pthread_mutex_unlock(lock);

    return stream->bad();
}

template <typename T>
PS_SIT pixsrc_printer::writeoutstream( T *source, PS_SIT dim1,
                                       std::ofstream *stream, pthread_mutex_t *lock, PS_SIT precision, vector<string>* comments )
{
    if (lock)
        pthread_mutex_lock(lock);

    if (comments)
        for (PS_SIT com=0; com<(PS_SIT)comments->size(); ++com)
            *stream << "# " << (*comments)[com] << std::endl;

    for(PS_SIT x=0; x<dim1; x++)
    {
        *stream << std::setprecision(precision) << source[x];
        if(x!=dim1-1)
            *stream << "\t";
    }
    *stream << std::endl;

    if (lock)
        pthread_mutex_unlock(lock);

    return stream->bad();
}

template <typename T>
PS_SIT pixsrc_printer::print(PS_SIT tracker, string basename, string imagename, bool num,
                             string filename, T source, bool nonstandardname, PS_SIT precision, vector<string>* comments)
{
    if (!nonstandardname)
    {
        imagename = imagename + "_" + OPERA tostring(tracker);
        if(num)
            filename = imagename + "_" + filename;
        filename = CONSTANT dir_out + basename + "_" + filename;
    }

    PS_SIT status = 0;

    std::ofstream file;
    file.open(filename.c_str());

    if (comments)
        for (PS_SIT com=0; com<(PS_SIT)comments->size(); ++com)
            file << "# " << (*comments)[com] << std::endl;

    file << std::setprecision(precision) << source;
    file.close();

    return file.bad();
}

template <typename T>
PS_SIT pixsrc_printer::print(PS_SIT tracker, string basename, string imagename, bool num,
                             string filename, T *source, PS_SIT dim1, bool nonstandardname, PS_SIT precision, vector<string>* comments)
{
    if (!nonstandardname)
    {
        imagename = imagename + "_" + OPERA tostring(tracker);
        if(num)
            filename = imagename + "_" + filename;
        filename = CONSTANT dir_out + basename + "_" + filename;
    }

    std::ofstream file;
    file.open(filename.c_str());

    if (comments)
        for (PS_SIT com=0; com<(PS_SIT)comments->size(); ++com)
            file << "# " << (*comments)[com] << std::endl;

    for(PS_SIT x = 0; x < dim1; x++)
    {
        file << std::setprecision(precision) << source[x];
        if(x!=dim1-1)
            file << "\n";
    }
    file.close();

    return file.bad();
}

template <typename T>
PS_SIT pixsrc_printer::print(PS_SIT tracker, string basename, string imagename, bool num,
                             string filename, vector<T>   source, bool nonstandardname, PS_SIT precision, vector<string>* comments)
{
    if (!nonstandardname)
    {
        imagename = imagename + "_" + OPERA tostring(tracker);
        if(num)
            filename = imagename + "_" + filename;
        filename = CONSTANT dir_out + basename + "_" + filename;
    }

    std::ofstream file;
    file.open(filename.c_str());

    if (comments)
        for (PS_SIT com=0; com<(PS_SIT)comments->size(); ++com)
            file << "# " << (*comments)[com] << std::endl;

    for(PS_unsignedSIT x = 0; x < source.size(); x++)
    {
        file << std::setprecision(precision) << source[x];
        if(x!=source.size()-1)
            file << "\n";
    }
    file.close();

    return file.bad();
}

template <typename T>
PS_SIT pixsrc_printer::print(PS_SIT tracker, string basename, string imagename, bool num,
                             string filename, vector< vector<T> > source, bool nonstandardname, PS_SIT precision, vector<string>* comments)
{
    if (!nonstandardname)
    {
        imagename = imagename + "_" + OPERA tostring(tracker);
        if(num)
            filename = imagename + "_" + filename;
        filename = CONSTANT dir_out + basename + "_" + filename;
    }

    std::ofstream file;
    file.open(filename.c_str());

    if (comments)
        for (PS_SIT com=0; com<(PS_SIT)comments->size(); ++com)
            file << "# " << (*comments)[com] << std::endl;

    for(PS_unsignedSIT x = 0; x < source.size(); x++)
    {
        for(PS_unsignedSIT y = 0; y < source[x].size(); y++)
        {
            file << std::setprecision(precision) << source[x][y];
            if(y!=source[x].size()-1)
                file << "\t";
        }
        if(x!=source.size()-1)
            file << "\n";
    }

    return file.bad();
}

template <typename T>
PS_SIT pixsrc_printer::print(PS_SIT tracker, string basename, string imagename, bool num,
                             string filename, T **source, PS_SIT dim1, PS_SIT dim2, bool nonstandardname, PS_SIT precision, vector<string>* comments)
{
    if (!nonstandardname)
    {
        imagename = imagename + "_" + OPERA tostring(tracker);
        if(num)
            filename = imagename + "_" + filename;
        filename = CONSTANT dir_out + basename + "_" + filename;
    }

    std::ofstream file;
    file.open(filename.c_str());

    if (comments)
        for (PS_SIT com=0; com<(PS_SIT)comments->size(); ++com)
            file << "# " << (*comments)[com] << std::endl;

    for(PS_SIT x = 0; x < dim1; x++)
    {
        for(PS_SIT y = 0; y < dim2; y++)
        {
            file << std::setprecision(precision) << *(source[x]+y);
            if(y!=dim2-1)
                file << "\t";
        }
        if(x!=dim1-1)
            file << "\n";
    }

    return file.bad();
}

template <typename T>
PS_SIT pixsrc_printer::printbinary(PS_SIT tracker, string basename, string imagename, bool num,
                                   string filename, const T *source, size_t dim1, bool nonstandardname)
{
    if (!nonstandardname)
    {
        imagename = imagename + "_" + OPERA tostring(tracker);
        if(num)
            filename = imagename + "_" + filename;
        filename = CONSTANT dir_out + basename + "_" + filename;
    }

    PS_SIT status = 0;
    FILE *file = fopen (filename.c_str(), "wb");
    if (NULL!=file)
    {
        if (dim1!=fwrite ((void*)source, sizeof(T), dim1, file))
            status = 1;
        status += fclose (file);
    }
    else
    {
        status = 1;
    }

    return status;
}

template <typename T>
PS_SIT pixsrc_printer::readbinary (string fn, T* ptr, size_t num, pthread_mutex_t *printlock)
{
    PS_SIT status = 0;
    FILE *file = fopen (fn.c_str(), "rb");
    if (NULL!=file)
    {
        if (num!=fread ((void*)ptr, sizeof(T), num, file))
            status = 1;
        status += fclose (file);
    }
    else
    {
        status = 1;
    }

    PRINTER print2screen("pixsrc",
                         fn + " loaded",
                         printlock);

    return status;
}



#endif
