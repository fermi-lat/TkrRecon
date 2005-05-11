/**
 * @file TkrException.h
 *

 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Utilities/TkrException.h,v 1.1 2003/09/27 18:53:27 burnett Exp $
 *
*/

#ifndef __TkrException_H
#define __TkrException_H

    /** @class TkrException 
         @brief hold a string
         */
    class TkrException : public std::exception{
    public: 
        TkrException(std::string error):m_what(error){}
    ~TkrException() throw() {;}
        virtual const char *what( ) const  throw() { return m_what.c_str();} 
        std::string m_what;
    };

#endif
