/**
 * @file TkrException.h
 *

 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Utilities/TkrPoints.h,v 1.2 2002/08/30 15:44:14 atwood Exp $
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
