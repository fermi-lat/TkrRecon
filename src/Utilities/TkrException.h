/**
 * @file TkrException.h
 *

 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Utilities/TkrException.h,v 1.2 2005/05/11 04:14:35 lsrea Exp $
 *
*/

#ifndef __TkrException_H
#define __TkrException_H

    /** @class TkrException 
         @brief hold a string
         */
class TkrException : public std::exception
{
public: 
    TkrException(std::string error):m_what(error){}
    ~TkrException() throw() {;}
        
    virtual const char *what( ) const  throw() { return m_what.c_str();} 
        
    std::string m_what;
};

#endif
